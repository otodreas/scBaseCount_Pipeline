from pathlib import Path

import pandas as pd
import pytest
from disease_markers.labels import (
    ControlType,
    build_sample_label_table,
    coarse_disease_area,
    sample_labels_by_srx,
)
from study_context.models import BiologicalContext, ExperimentContext, StudyContext


def _context(
    *,
    accession: str,
    tissue_type: str | None = None,
    sample_title: str | None = None,
    attributes: dict[str, str] | None = None,
) -> ExperimentContext:
    return ExperimentContext(
        accession=accession,
        biological=BiologicalContext(
            scientificName="Homo sapiens",
            tissueType=tissue_type,
            sampleTitle=sample_title,
            sampleAttributes=attributes or {},
        ),
        study=StudyContext(studyAccession="PRJ_TEST", studyTitle="Test study"),
    )


def _write_label_inputs(
    tmp_path: Path,
    *,
    rows: list[dict[str, object]],
    contexts: list[ExperimentContext] | None = None,
) -> tuple[Path, Path, Path]:
    metadata = pd.DataFrame(rows)
    metadata_path = tmp_path / "sample_metadata.parquet"
    metadata.to_parquet(metadata_path, index=False)

    atlas_lines = ["accession,status,studyAccession"]
    for row in rows:
        atlas_lines.append(f"{row['srx_accession']},success,PRJ_TEST")
    atlas_path = tmp_path / "atlas.csv"
    atlas_path.write_text("\n".join(atlas_lines) + "\n")

    contexts_path = tmp_path / "contexts.jsonl"
    if contexts:
        contexts_path.write_text("".join(f"{ctx.model_dump_json()}\n" for ctx in contexts))
    else:
        contexts_path.write_text("")
    return contexts_path, atlas_path, metadata_path


def _metadata_row(
    accession: str,
    *,
    disease: str,
    tissue: str,
    cell_line: str = "none",
    organism: str = "Homo sapiens",
    perturbation: str = "none",
) -> dict[str, object]:
    return {
        "srx_accession": accession,
        "disease": disease,
        "tissue": tissue,
        "organism": organism,
        "cell_line": cell_line,
        "perturbation": perturbation,
    }


def test_coarse_disease_area_ipf() -> None:
    assert coarse_disease_area("idiopathic pulmonary fibrosis") == "IPF / Pulmonary Fibrosis"


def test_coarse_disease_area_normal_tissue_is_other() -> None:
    assert coarse_disease_area("normal lung tissue") == "Other"


def test_a549_cell_line_is_excluded_sample_type_not_non_lung(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX12285822",
                disease="SARS-CoV-2 infection",
                tissue="lung (A549 cell line, ACE2-transduced)",
                cell_line="A549 (ACE2-transduced)",
                perturbation="SARS-CoV-2 infection, MOI 0.01, 24h",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = table.iloc[0]

    assert row["diseaseArea"] == "COVID-19 / SARS-CoV-2"
    assert row["tissueRaw"].lower().startswith("lung")
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "excluded_sample_type"


def test_matched_adjacent_overrides_parquet_disease_cohort(tmp_path: Path) -> None:
    accession = "ERX11662359"
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                accession,
                disease="lung squamous cell carcinoma",
                tissue="lung",
                cell_line="unsure",
            )
        ],
        contexts=[
            _context(
                accession=accession,
                attributes={
                    "organism part": "lung",
                    "sampling site": "normal tissue adjacent to tumor",
                },
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = sample_labels_by_srx(table)[accession]

    assert row.diseaseArea == "Lung Cancer"
    assert row.diseased is False
    assert row.isBiologicalControl is True
    assert row.controlType is ControlType.MATCHED_ADJACENT
    assert row.eligible is True
    assert row.excludeReason is None


def test_non_diseased_parquet_disease_is_biological_control(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX23825222",
                disease="non-diseased",
                tissue="lung",
                cell_line="not specified",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = sample_labels_by_srx(table)["SRX23825222"]

    assert row.diseaseArea == "Other"
    assert row.diseased is False
    assert row.isBiologicalControl is True
    assert row.controlType is ControlType.HEALTHY
    assert row.eligible is True


def test_primary_lung_descriptive_cell_line_remains_eligible(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX26735834",
                disease="Healthy",
                tissue="Lung",
                cell_line="Airway epithelium",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = table.iloc[0]

    assert bool(row["diseased"]) is False
    assert bool(row["eligible"]) is True
    assert row["excludeReason"] is None


def test_pbmc_bal_mixture_is_non_lung(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX11071611",
                disease="COVID-19",
                tissue="PBMC (peripheral blood mononuclear cells), BAL (bronchoalveolar lavage)",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = table.iloc[0]

    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "non_lung"


def test_organoid_is_excluded_sample_type(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_ORG",
                disease="SARS-CoV-2 infection",
                tissue="human airway organoid",
                cell_line="none",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = table.iloc[0]

    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "excluded_sample_type"


def test_missing_parquet_coverage_fails_fast(tmp_path: Path) -> None:
    metadata_path = tmp_path / "sample_metadata.parquet"
    pd.DataFrame(
        [
            _metadata_row(
                "SRX_OTHER",
                disease="COVID-19",
                tissue="lung",
            )
        ]
    ).to_parquet(metadata_path, index=False)
    atlas_path = tmp_path / "atlas.csv"
    atlas_path.write_text("accession,status,studyAccession\nSRX_MISSING,success,PRJ_TEST\n")
    contexts_path = tmp_path / "contexts.jsonl"
    contexts_path.write_text("")

    with pytest.raises(KeyError, match="missing 1 atlas accession"):
        build_sample_label_table(contexts_path, atlas_path, metadata_path)


def test_protocol_text_in_sample_description_does_not_override_tumor(tmp_path: Path) -> None:
    accession = "ERX11876748"
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                accession,
                disease="non-small cell lung adenocarcinoma",
                tissue="lung tumor central margin",
                cell_line="fibroblast, endothelial cell",
            )
        ],
        contexts=[
            ExperimentContext(
                accession=accession,
                biological=BiologicalContext(
                    scientificName="Homo sapiens",
                    sampleDescription="Fresh NSCLC and unaffected autologous lung tissue were processed together.",
                    sampleAttributes={"organism part": "lung tumor central margin"},
                ),
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)
    row = sample_labels_by_srx(table)[accession]

    assert row.diseased is True
    assert row.isBiologicalControl is False
    assert row.controlType is None


def test_build_sample_label_table_emits_raw_source_fields(tmp_path: Path) -> None:
    contexts_path, atlas_path, metadata_path = _write_label_inputs(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_OK",
                disease="idiopathic pulmonary fibrosis",
                tissue="lung",
                cell_line="none",
            )
        ],
    )

    table = build_sample_label_table(contexts_path, atlas_path, metadata_path)

    assert table["diseased"].dtype == pd.BooleanDtype()
    assert set(table.columns) >= {
        "srxAccession",
        "studyAccession",
        "diseaseRaw",
        "tissueRaw",
        "cellLineRaw",
        "diseaseArea",
        "diseased",
        "isBiologicalControl",
        "controlType",
        "eligible",
        "excludeReason",
    }
    assert table.to_dict(orient="records") == [
        {
            "srxAccession": "SRX_OK",
            "studyAccession": "PRJ_TEST",
            "diseaseRaw": "idiopathic pulmonary fibrosis",
            "tissueRaw": "lung",
            "cellLineRaw": "none",
            "diseaseArea": "IPF / Pulmonary Fibrosis",
            "diseased": True,
            "isBiologicalControl": False,
            "controlType": None,
            "eligible": True,
            "excludeReason": None,
        }
    ]
