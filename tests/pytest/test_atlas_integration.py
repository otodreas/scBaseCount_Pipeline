import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_integration import (
    AtlasIntegrationConfig,
    build_accession_study_map,
    normalize_cell_type_labels,
    study_for_accession,
)
from study_context.models import ExperimentContext, StudyContext
from study_context.utils import CONTEXTS_JSONL_PATH, load_contexts_jsonl


def test_study_for_accession_uses_study_accession() -> None:
    ctx = ExperimentContext(
        accession="SRX1",
        study=StudyContext(studyAccession="PRJNA123", studyTitle="Example study"),
    )
    assert study_for_accession("SRX1", {"SRX1": ctx}) == "PRJNA123"


def test_study_for_accession_falls_back_to_accession() -> None:
    ctx = ExperimentContext(accession="SRX2", study=None)
    assert study_for_accession("SRX2", {"SRX2": ctx}) == "SRX2"


def test_normalize_cell_type_labels_fills_missing() -> None:
    adata = sc.AnnData(X=np.zeros((3, 2)))
    adata.obs["cell_type"] = ["T cell", None, " "]
    cfg = AtlasIntegrationConfig()
    normalize_cell_type_labels(adata, cfg)
    assert adata.obs["cell_type"].tolist() == ["T cell", "unknown", "unknown"]


@pytest.mark.skipif(not CONTEXTS_JSONL_PATH.is_file(), reason="contexts.jsonl not present locally")
def test_build_accession_study_map_covers_datasets_csv() -> None:
    datasets = pd.read_csv(Path("output/metadata/datasets.csv"))
    accessions = datasets["srx_accession"].astype(str).head(10).tolist()
    mapping = build_accession_study_map(accessions, CONTEXTS_JSONL_PATH)
    assert set(mapping) == set(accessions)
    assert all(value for value in mapping.values())


@pytest.mark.skipif(not CONTEXTS_JSONL_PATH.is_file(), reason="contexts.jsonl not present locally")
def test_contexts_cover_all_datasets_when_regenerated() -> None:
    datasets = pd.read_csv(Path("output/metadata/datasets.csv"))
    contexts = load_contexts_jsonl(CONTEXTS_JSONL_PATH)
    dataset_accessions = set(datasets["srx_accession"].astype(str))
    assert dataset_accessions.issubset(set(contexts))


def test_run_metadata_roundtrip() -> None:
    payload = {
        "mergeStats": {
            "nAccessionsRequested": 2,
            "nAccessionsMerged": 2,
            "nAccessionsSkipped": 0,
            "nCellsFinal": 100,
            "nGenesFinal": 50,
            "nStudies": 1,
            "skippedAccessions": [],
        }
    }
    text = json.dumps(payload)
    loaded = json.loads(text)
    assert loaded["mergeStats"]["nCellsFinal"] == 100
