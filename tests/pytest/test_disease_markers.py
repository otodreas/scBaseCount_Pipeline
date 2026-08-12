import json
from pathlib import Path

import pandas as pd
import pytest
from disease_markers.labels import build_sample_label_table
from disease_markers.status import ControlType
from metadata.categorize import coarse_disease_area
from ontology_lookup import OntologyLookupConfig
from study_context.models import BiologicalContext, ExperimentContext, StudyContext

MONDO_FIXTURE: dict[str, dict[str, object]] = {
    "MONDO:0005061": {
        "label": "lung adenocarcinoma",
        "ancestors": ["MONDO:0008903", "MONDO:0005233", "MONDO:0004992"],
    },
    "MONDO:0005097": {
        "label": "squamous cell lung carcinoma",
        "ancestors": ["MONDO:0008903", "MONDO:0004992"],
    },
    "MONDO:0005233": {
        "label": "non-small cell lung carcinoma",
        "ancestors": ["MONDO:0008903", "MONDO:0004992"],
    },
    "MONDO:0008903": {"label": "lung cancer", "ancestors": ["MONDO:0004992"]},
    "MONDO:0100096": {"label": "COVID-19", "ancestors": ["MONDO:0005550"]},
    "MONDO:0800504": {
        "label": "idiopathic pulmonary fibrosis",
        "ancestors": ["MONDO:0002429", "MONDO:0015925"],
    },
    "MONDO:0002771": {"label": "pulmonary fibrosis", "ancestors": ["MONDO:0015925"]},
    "MONDO:0002429": {"label": "idiopathic interstitial pneumonia", "ancestors": ["MONDO:0015925"]},
    "MONDO:0005002": {"label": "chronic obstructive pulmonary disease", "ancestors": []},
    "MONDO:0009061": {"label": "cystic fibrosis", "ancestors": []},
    "MONDO:0015925": {"label": "interstitial lung disease", "ancestors": []},
    "MONDO:0005149": {"label": "pulmonary hypertension", "ancestors": []},
    "MONDO:0004992": {"label": "cancer", "ancestors": ["MONDO:0005070"]},
}
UBERON_FIXTURE: dict[str, dict[str, object]] = {
    "UBERON:0002048": {"label": "lung", "ancestors": ["UBERON:0001004"]},
    "UBERON:0001004": {"label": "respiratory system", "ancestors": []},
    "UBERON:0000178": {"label": "blood", "ancestors": ["UBERON:0000465"]},
    "UBERON:0039167": {"label": "bronchopulmonary lymph node", "ancestors": ["UBERON:0000029"]},
}


def _write_ontology_cache(tmp_path: Path) -> OntologyLookupConfig:
    for ontology_id, release, terms in (
        ("mondo", "2026-07-06", MONDO_FIXTURE),
        ("uberon", "2026-06-19", UBERON_FIXTURE),
    ):
        root = tmp_path / "ontologies" / ontology_id / release
        root.mkdir(parents=True)
        (root / "terms.json").write_text(json.dumps(terms, indent=2, sort_keys=True) + "\n")
        (root / "manifest.json").write_text(
            json.dumps(
                {
                    "ontologyId": ontology_id,
                    "release": release,
                    "sourceUrl": "https://example.test",
                    "generatedAt": "2026-01-01T00:00:00+00:00",
                    "termCount": len(terms),
                    "contentSha256": "test",
                },
                indent=2,
            )
            + "\n"
        )
    return OntologyLookupConfig(cacheDir=tmp_path / "ontologies")


def _context(
    *,
    accession: str,
    scientific_name: str | None = "Homo sapiens",
    tax_id: str | None = "9606",
    tissue_type: str | None = None,
    sample_title: str | None = None,
    sample_description: str | None = None,
    attributes: dict[str, str] | None = None,
    study_accession: str = "PRJ_TEST",
) -> ExperimentContext:
    return ExperimentContext(
        accession=accession,
        biological=BiologicalContext(
            scientificName=scientific_name,
            taxId=tax_id,
            tissueType=tissue_type,
            sampleTitle=sample_title,
            sampleDescription=sample_description,
            sampleAttributes=attributes or {},
        ),
        study=StudyContext(studyAccession=study_accession, studyTitle="Test study"),
    )


def _metadata_row(
    accession: str,
    *,
    disease: str,
    tissue: str,
    cell_line: str = "none",
    organism: str = "Homo sapiens",
    tissue_ontology: str | None = "UBERON:0002048",
    disease_ontology: str | None = None,
) -> dict[str, object]:
    return {
        "srx_accession": accession,
        "disease": disease,
        "disease_ontology_term_id": disease_ontology,
        "tissue": tissue,
        "tissue_ontology_term_id": tissue_ontology,
        "organism": organism,
        "cell_line": cell_line,
    }


def _write_atlas_manifest(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """Write a minimal post-cutover result manifest. rows are (accession, status, study)."""
    files = [
        {
            "accession": accession,
            "studyAccession": study,
            "r2Key": f"prefix/{accession}.h5ad",
            "status": status,
            "skipReason": None if status == "success" else "download_failed",
            "qc": None,
        }
        for accession, status, study in rows
    ]
    path.write_text(json.dumps({"files": files}) + "\n")


def _write_label_inputs(
    tmp_path: Path,
    *,
    rows: list[dict[str, object]],
    contexts: list[ExperimentContext],
) -> tuple[Path, Path, Path, OntologyLookupConfig]:
    metadata_path = tmp_path / "sample_metadata.parquet"
    pd.DataFrame(rows).to_parquet(metadata_path, index=False)

    atlas_rows = [
        (
            str(row["srx_accession"]),
            "success",
            "PRJ_TEST" if ctx.study is None else ctx.study.studyAccession,
        )
        for row, ctx in zip(rows, contexts, strict=True)
    ]
    atlas_path = tmp_path / "atlas_result.json"
    _write_atlas_manifest(atlas_path, atlas_rows)

    contexts_path = tmp_path / "contexts.jsonl"
    contexts_path.write_text("".join(f"{ctx.model_dump_json()}\n" for ctx in contexts))
    ontology_cfg = _write_ontology_cache(tmp_path)
    return contexts_path, atlas_path, metadata_path, ontology_cfg


def _build(
    tmp_path: Path,
    *,
    rows: list[dict[str, object]],
    contexts: list[ExperimentContext],
) -> pd.DataFrame:
    contexts_path, atlas_path, metadata_path, ontology_cfg = _write_label_inputs(tmp_path, rows=rows, contexts=contexts)
    return build_sample_label_table(
        contexts_path,
        atlas_path,
        metadata_path,
        ontologyConfig=ontology_cfg,
    )


def test_coarse_disease_area_ipf() -> None:
    assert coarse_disease_area("idiopathic pulmonary fibrosis") == "IPF / Pulmonary Fibrosis"


def test_coarse_disease_area_normal_tissue_is_other() -> None:
    assert coarse_disease_area("normal lung tissue") == "Other"


def test_mondo_resolves_disease_area_and_name(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_LUAD",
                disease="unsure",
                tissue="lung",
                disease_ontology="MONDO:0005061",
            )
        ],
        contexts=[_context(accession="SRX_LUAD", attributes={"tissue": "lung"})],
    )
    row = table.iloc[0]
    assert row["diseaseArea"] == "Lung Cancer"
    assert row["diseaseName"] == "lung adenocarcinoma"
    assert row["diseaseAreaSource"] == "ontology"
    assert bool(row["diseased"]) is True
    assert bool(row["eligible"]) is True


def test_generic_cancer_mondo_does_not_become_lung_cancer(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_CANCER",
                disease="cancer",
                tissue="lung",
                disease_ontology="MONDO:0004992",
            )
        ],
        contexts=[_context(accession="SRX_CANCER", attributes={"tissue": "lung"})],
    ).iloc[0]
    assert row["diseaseArea"] == "Other"
    assert row["diseaseAreaSource"] in {"unmapped_ontology", "disease_text"}


def test_text_fallback_when_ontology_missing(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[_metadata_row("SRX_TEXT", disease="lung adenocarcinoma", tissue="lung", disease_ontology=None)],
        contexts=[_context(accession="SRX_TEXT", attributes={"tissue": "lung"})],
    ).iloc[0]
    assert row["diseaseArea"] == "Lung Cancer"
    assert row["diseaseAreaSource"] == "disease_text"


def test_study_consensus_fills_missing_ontology(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_A",
                disease="lung adenocarcinoma",
                tissue="lung",
                disease_ontology="MONDO:0005061",
            ),
            _metadata_row("SRX_B", disease="none", tissue="lung", disease_ontology=None),
        ],
        contexts=[
            _context(accession="SRX_A", attributes={"tissue": "lung", "disease": "lung adenocarcinoma"}),
            _context(accession="SRX_B", attributes={"tissue": "lung", "disease": "Control"}),
        ],
    ).set_index("srxAccession")
    assert table.loc["SRX_B", "diseaseArea"] == "Lung Cancer"
    assert table.loc["SRX_B", "diseaseAreaSource"] == "study_consensus"
    assert bool(table.loc["SRX_B", "diseased"]) is False


def test_mixed_study_does_not_fill_consensus(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row("SRX_A", disease="COVID-19", tissue="lung", disease_ontology="MONDO:0100096"),
            _metadata_row(
                "SRX_B",
                disease="COPD",
                tissue="lung",
                disease_ontology="MONDO:0005002",
            ),
            _metadata_row("SRX_C", disease="none", tissue="lung", disease_ontology=None),
        ],
        contexts=[
            _context(accession="SRX_A", attributes={"tissue": "lung"}),
            _context(accession="SRX_B", attributes={"tissue": "lung"}),
            _context(accession="SRX_C", attributes={"tissue": "lung", "disease": "Control"}),
        ],
    ).set_index("srxAccession")
    assert table.loc["SRX_C", "diseaseArea"] == "Other"
    assert table.loc["SRX_C", "diseaseAreaSource"] != "study_consensus"


def test_blood_ontology_veto_even_with_lung_token(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX9058196",
                disease="COVID-19, healthy control",
                tissue="blood, bronchoalveolar lavage fluid",
                tissue_ontology="UBERON:0000178",
                disease_ontology="MONDO:0100096",
            )
        ],
        contexts=[_context(accession="SRX9058196", tissue_type="Blood", attributes={"tissue": "Blood"})],
    ).iloc[0]
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "non_lung"


def test_multi_valued_blood_plus_lung_ontology_is_non_lung(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX12708364",
                disease="severe COVID-19",
                tissue="blood (circulating leukocytes), endotracheal aspirate (ETA)",
                tissue_ontology="UBERON:0000178,UBERON:0002048",
                disease_ontology="MONDO:0100096",
            )
        ],
        contexts=[
            _context(
                accession="SRX12708364",
                tissue_type="Blood",
                attributes={"tissue": "Blood"},
            )
        ],
    ).iloc[0]
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "non_lung"


def test_context_mouse_overrides_human_parquet(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX12366723",
                disease="lung neoplasms",
                tissue="lung",
                tissue_ontology="UBERON:0002048",
                disease_ontology="MONDO:0008903",
            )
        ],
        contexts=[
            _context(
                accession="SRX12366723",
                scientific_name="Mus musculus",
                tax_id="10090",
                tissue_type="Hippocampus",
                attributes={"tissue": "Hippocampus"},
            )
        ],
    ).iloc[0]
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "non_human"


def test_context_pbmc_overrides_lung_parquet(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX13061245",
                disease="COVID-19",
                tissue="lung",
                tissue_ontology="UBERON:0002048",
                disease_ontology="MONDO:0100096",
            )
        ],
        contexts=[
            _context(
                accession="SRX13061245",
                tissue_type="PBMC",
                attributes={"tissue": "PBMC", "source_name": "CD8+ T cells"},
            )
        ],
    ).iloc[0]
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "non_lung"


def test_matched_adjacent_overrides_parquet_disease_cohort(tmp_path: Path) -> None:
    accession = "ERX11662359"
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                accession,
                disease="lung squamous cell carcinoma",
                tissue="lung",
                cell_line="unsure",
                disease_ontology="MONDO:0005097",
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
    ).iloc[0]
    assert row["diseaseArea"] == "Lung Cancer"
    assert bool(row["diseased"]) is False
    assert bool(row["isBiologicalControl"]) is True
    assert row["controlType"] == ControlType.MATCHED_ADJACENT.value
    assert bool(row["eligible"]) is True


def test_exact_status_control_in_disease_labelled_cohort(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX17412786",
                disease="idiopathic pulmonary fibrosis (IPF)",
                tissue="lung",
                disease_ontology="MONDO:0800504",
            )
        ],
        contexts=[
            _context(
                accession="SRX17412786",
                tissue_type="lung",
                attributes={"disease": "Control", "tissue": "lung"},
            )
        ],
    ).iloc[0]
    assert bool(row["diseased"]) is False
    assert bool(row["isBiologicalControl"]) is True
    assert row["controlType"] == ControlType.EXPLICIT_CONTROL.value


def test_protocol_text_in_sample_description_does_not_override_tumor(tmp_path: Path) -> None:
    accession = "ERX11876748"
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                accession,
                disease="non-small cell lung adenocarcinoma",
                tissue="lung tumor central margin",
                cell_line="fibroblast, endothelial cell",
                disease_ontology="MONDO:0005233",
            )
        ],
        contexts=[
            _context(
                accession=accession,
                sample_description="Fresh NSCLC and unaffected autologous lung tissue were processed together.",
                attributes={"organism part": "lung tumor central margin"},
            )
        ],
    ).iloc[0]
    assert bool(row["diseased"]) is True
    assert bool(row["isBiologicalControl"]) is False
    assert row["controlType"] is None


def test_paired_mock_sets_nondiseased_and_experimental_comparator(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_MOCK",
                disease="SARS-CoV-2 infection",
                tissue="lung",
                disease_ontology="MONDO:0100096",
            ),
            _metadata_row(
                "SRX_INF",
                disease="SARS-CoV-2 infection",
                tissue="lung",
                disease_ontology="MONDO:0100096",
            ),
        ],
        contexts=[
            _context(accession="SRX_MOCK", attributes={"treatment": "mock", "tissue": "lung"}),
            _context(accession="SRX_INF", attributes={"treatment": "SARS-CoV-2 infection", "tissue": "lung"}),
        ],
    ).set_index("srxAccession")
    assert bool(table.loc["SRX_MOCK", "isExperimentalComparator"]) is True
    assert bool(table.loc["SRX_MOCK", "diseased"]) is False
    assert bool(table.loc["SRX_INF", "isExperimentalComparator"]) is False
    assert bool(table.loc["SRX_INF", "diseased"]) is True


def test_paired_untreated_tumor_remains_diseased(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_UT",
                disease="lung adenocarcinoma",
                tissue="lung",
                disease_ontology="MONDO:0005061",
            ),
            _metadata_row(
                "SRX_TX",
                disease="lung adenocarcinoma",
                tissue="lung",
                disease_ontology="MONDO:0005061",
            ),
        ],
        contexts=[
            _context(accession="SRX_UT", attributes={"treatment": "untreated", "tissue": "lung"}),
            _context(
                accession="SRX_TX",
                attributes={"treatment": "kinase inhibitor treatment", "tissue": "lung"},
            ),
        ],
    ).set_index("srxAccession")
    assert bool(table.loc["SRX_UT", "isExperimentalComparator"]) is True
    assert bool(table.loc["SRX_UT", "diseased"]) is True
    assert bool(table.loc["SRX_UT", "isBiologicalControl"]) is False


def test_unpaired_comparator_is_not_experimental(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_ONLY",
                disease="COVID-19",
                tissue="lung",
                disease_ontology="MONDO:0100096",
            )
        ],
        contexts=[
            _context(
                accession="SRX_ONLY",
                attributes={"treatment": "unexposed", "tissue": "lung"},
            )
        ],
    ).iloc[0]
    assert bool(row["isExperimentalComparator"]) is False


def test_mixed_arm_value_is_rejected(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_MIX",
                disease="SARS-CoV-2 infection",
                tissue="lung",
                disease_ontology="MONDO:0100096",
            ),
            _metadata_row(
                "SRX_POS",
                disease="SARS-CoV-2 infection",
                tissue="lung",
                disease_ontology="MONDO:0100096",
            ),
        ],
        contexts=[
            _context(
                accession="SRX_MIX",
                attributes={"treatment": "4 mock + 4 infected", "tissue": "lung"},
            ),
            _context(
                accession="SRX_POS",
                attributes={"treatment": "SARS-CoV-2 infection", "tissue": "lung"},
            ),
        ],
    ).set_index("srxAccession")
    assert bool(row.loc["SRX_MIX", "isExperimentalComparator"]) is False


def test_a549_cell_line_is_excluded_sample_type(tmp_path: Path) -> None:
    row = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX12285822",
                disease="SARS-CoV-2 infection",
                tissue="lung (A549 cell line, ACE2-transduced)",
                cell_line="A549 (ACE2-transduced)",
                disease_ontology="MONDO:0100096",
            )
        ],
        contexts=[
            _context(
                accession="SRX12285822",
                attributes={"source_name": "A549-ACE2 cells", "tissue": "lung"},
            )
        ],
    ).iloc[0]
    assert bool(row["eligible"]) is False
    assert row["excludeReason"] == "excluded_sample_type"


def test_missing_parquet_coverage_fails_fast(tmp_path: Path) -> None:
    ontology_cfg = _write_ontology_cache(tmp_path)
    metadata_path = tmp_path / "sample_metadata.parquet"
    pd.DataFrame([_metadata_row("SRX_OTHER", disease="COVID-19", tissue="lung")]).to_parquet(metadata_path, index=False)
    atlas_path = tmp_path / "atlas_result.json"
    _write_atlas_manifest(atlas_path, [("SRX_MISSING", "success", "PRJ_TEST")])
    contexts_path = tmp_path / "contexts.jsonl"
    contexts_path.write_text(_context(accession="SRX_MISSING").model_dump_json() + "\n")

    with pytest.raises(KeyError, match="sample metadata missing 1 atlas accession"):
        build_sample_label_table(contexts_path, atlas_path, metadata_path, ontologyConfig=ontology_cfg)


def test_missing_context_coverage_fails_fast(tmp_path: Path) -> None:
    ontology_cfg = _write_ontology_cache(tmp_path)
    metadata_path = tmp_path / "sample_metadata.parquet"
    pd.DataFrame([_metadata_row("SRX_OK", disease="COVID-19", tissue="lung")]).to_parquet(metadata_path, index=False)
    atlas_path = tmp_path / "atlas_result.json"
    _write_atlas_manifest(atlas_path, [("SRX_OK", "success", "PRJ_TEST")])
    contexts_path = tmp_path / "contexts.jsonl"
    contexts_path.write_text("")

    with pytest.raises(KeyError, match="contexts missing 1 atlas accession"):
        build_sample_label_table(contexts_path, atlas_path, metadata_path, ontologyConfig=ontology_cfg)


def test_output_schema_and_nullable_diseased(tmp_path: Path) -> None:
    table = _build(
        tmp_path,
        rows=[
            _metadata_row(
                "SRX_OK",
                disease="idiopathic pulmonary fibrosis",
                tissue="lung",
                disease_ontology="MONDO:0800504",
            )
        ],
        contexts=[_context(accession="SRX_OK", attributes={"tissue": "lung", "disease": "IPF"})],
    )
    assert table["diseased"].dtype == pd.BooleanDtype()
    assert set(table.columns) == {
        "srxAccession",
        "studyAccession",
        "diseaseRaw",
        "diseaseOntologyTermId",
        "diseaseName",
        "tissueRaw",
        "tissueOntologyRaw",
        "cellLineRaw",
        "diseaseArea",
        "diseaseAreaSource",
        "diseased",
        "isBiologicalControl",
        "controlType",
        "isExperimentalComparator",
        "eligible",
        "excludeReason",
    }
    assert len(table) == 1
    assert table.iloc[0]["diseaseArea"] == "IPF / Pulmonary Fibrosis"
