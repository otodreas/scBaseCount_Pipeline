from pathlib import Path

import pandas as pd
import pytest
from disease_markers.labels import (
    ControlType,
    _disease_status,
    _is_eligible,
    build_sample_label_table,
    coarse_disease_area,
    sample_labels_by_srx,
)
from study_context.models import BiologicalContext, ExperimentContext, StudyContext


def _context(
    *,
    accession: str = "SRX_TEST",
    experiment_title: str | None = None,
    tissue_type: str | None = None,
    sample_title: str | None = None,
    sample_description: str | None = None,
    attributes: dict[str, str] | None = None,
    study_title: str | None = None,
) -> ExperimentContext:
    study = None
    if study_title is not None:
        study = StudyContext(studyAccession="PRJ_TEST", studyTitle=study_title)
    return ExperimentContext(
        accession=accession,
        experimentTitle=experiment_title,
        biological=BiologicalContext(
            scientificName="Homo sapiens",
            tissueType=tissue_type,
            sampleTitle=sample_title,
            sampleDescription=sample_description,
            sampleAttributes=attributes or {},
        ),
        study=study,
    )


def test_coarse_disease_area_ipf() -> None:
    assert coarse_disease_area("idiopathic pulmonary fibrosis") == "IPF / Pulmonary Fibrosis"


def test_coarse_disease_area_normal_tissue_is_other() -> None:
    assert coarse_disease_area("normal lung tissue") == "Other"


def test_matched_adjacent_lung_retains_cancer_area() -> None:
    context = _context(
        accession="ERX11662359",
        attributes={
            "disease": "lung squamous cell carcinoma",
            "organism part": "lung",
            "sampling site": "normal tissue adjacent to tumor",
        },
        study_title="Single cell RNA-seq atlas of human NSCLC lesions and non-involved tissue",
    )

    diseased, control_type = _disease_status(context)

    assert coarse_disease_area("lung squamous cell carcinoma") == "Lung Cancer"
    assert diseased is False
    assert control_type is ControlType.MATCHED_ADJACENT
    assert _is_eligible(context, "Lung Cancer", diseased) == (True, None)


def test_protocol_text_does_not_make_tumor_a_control() -> None:
    context = _context(
        accession="ERX11876748",
        sample_description="Fresh NSCLC and unaffected autologous lung tissue were processed together.",
        attributes={
            "disease": "non-small cell lung adenocarcinoma",
            "organism part": "lung tumor central margin",
        },
        study_title="Transcriptomic profiling of the NSCLC environment",
    )

    assert _disease_status(context) == (True, None)


def test_unaffected_specimen_is_biological_control() -> None:
    context = _context(
        accession="ERX11876755",
        attributes={
            "disease": "non-small cell lung adenocarcinoma",
            "organism part": "lung unaffected",
        },
        study_title="Transcriptomic profiling of the NSCLC environment",
    )

    assert _disease_status(context) == (False, ControlType.HEALTHY)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("genotype", "WT"),
        ("treatment", "vehicle"),
        ("treatment", "no treatment"),
        ("status", "baseline"),
    ],
)
def test_experimental_comparator_alone_does_not_prove_non_disease(key: str, value: str) -> None:
    context = _context(tissue_type="lung", attributes={key: value})

    assert _disease_status(context) == (None, None)


def test_normal_lung_without_disease_cohort_is_eligible_control() -> None:
    context = _context(tissue_type="normal lung")
    diseased, control_type = _disease_status(context)

    assert diseased is False
    assert control_type is ControlType.HEALTHY
    assert coarse_disease_area("") == "Other"
    assert _is_eligible(context, "Other", diseased) == (True, None)


@pytest.mark.parametrize(
    "tissue_type",
    ["PBMC", "blood", "PBMC from lung cancer patient", "lung and lymph node", None],
)
def test_non_lung_or_unknown_tissue_is_ineligible(tissue_type: str | None) -> None:
    context = _context(
        tissue_type=tissue_type,
        attributes={"disease": "COVID-19"},
        study_title="COVID-19 study",
    )

    assert _is_eligible(context, "COVID-19 / SARS-CoV-2", True) == (False, "non_lung")


def test_pooled_donors_are_ineligible_without_per_cell_provenance() -> None:
    context = _context(
        experiment_title="Pooled BCR library for Donors 1, 2, and 3",
        tissue_type="lung",
    )

    assert _is_eligible(context, "Other", False) == (False, "mixed_sample")


def test_build_sample_label_table_emits_independent_status_fields(tmp_path: Path) -> None:
    context = _context(
        accession="ERX11662359",
        attributes={
            "disease": "lung squamous cell carcinoma",
            "organism part": "lung",
            "sampling site": "normal tissue adjacent to tumor",
        },
        study_title="Single cell RNA-seq atlas of human NSCLC lesions and non-involved tissue",
    )
    contexts_path = tmp_path / "contexts.jsonl"
    contexts_path.write_text(f"{context.model_dump_json()}\n")
    atlas_path = tmp_path / "atlas.csv"
    atlas_path.write_text("accession,status,studyAccession\nERX11662359,success,PRJ_TEST\n")

    table = build_sample_label_table(
        contexts_path,
        atlas_path,
    )

    assert table["diseased"].dtype == pd.BooleanDtype()
    assert table.to_dict(orient="records") == [
        {
            "srxAccession": "ERX11662359",
            "studyAccession": "PRJ_TEST",
            "diseaseRaw": (
                "lung squamous cell carcinoma Single cell RNA-seq atlas of human NSCLC lesions and non-involved tissue"
            ),
            "diseaseArea": "Lung Cancer",
            "diseased": False,
            "isBiologicalControl": True,
            "controlType": "matchedAdjacent",
            "eligible": True,
            "excludeReason": None,
        }
    ]
    row = sample_labels_by_srx(table)["ERX11662359"]
    assert row.diseased is False
    assert row.isBiologicalControl is True
    assert row.controlType is ControlType.MATCHED_ADJACENT
