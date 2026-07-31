from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import load_and_normalize
from atlas_postprocessing.sampling import (
    SAMPLE_METHOD,
    SAMPLE_SEED,
    SAMPLE_UNS_KEY,
    _allocate_study_counts,
    sample_metadata,
    sample_study_proportional,
)


def _make_adata(study_sizes: dict[str, int]) -> sc.AnnData:
    studies: list[str] = []
    for study, size in study_sizes.items():
        studies.extend([study] * size)
    n_obs = len(studies)
    adata = sc.AnnData(
        np.ones((n_obs, 3), dtype=np.float32),
        obs=pd.DataFrame({"study_accession": studies}),
    )
    return adata


def test_allocate_study_counts_exact_and_coverage() -> None:
    sizes = {"A": 50, "B": 30, "C": 20}
    alloc = _allocate_study_counts(sizes, 10)
    assert sum(alloc.values()) == 10
    assert set(alloc) == {"A", "B", "C"}
    assert all(count >= 1 for count in alloc.values())
    assert alloc["A"] >= alloc["B"] >= alloc["C"]


def test_allocate_rejects_too_few_for_all_studies() -> None:
    with pytest.raises(ValueError, match="number of studies"):
        _allocate_study_counts({"A": 10, "B": 10, "C": 10}, 2)


def test_allocate_rejects_over_total() -> None:
    with pytest.raises(ValueError, match="exceeds available cells"):
        _allocate_study_counts({"A": 2, "B": 2}, 5)


def test_sample_study_proportional_exact_size_and_metadata() -> None:
    adata = _make_adata({"A": 40, "B": 35, "C": 25})
    sampled = sample_study_proportional(adata, n=20, stratifyKey="study_accession", seed=SAMPLE_SEED)
    assert sampled.n_obs == 20
    assert sampled.obs["study_accession"].nunique() == 3
    meta = sample_metadata(sampled)
    assert meta is not None
    assert meta["sourceCells"] == 100
    assert meta["sampleCells"] == 20
    assert meta["method"] == SAMPLE_METHOD
    assert meta["stratifyKey"] == "study_accession"
    assert meta["seed"] == SAMPLE_SEED
    assert meta["nStudies"] == 3
    assert SAMPLE_UNS_KEY in sampled.uns


def test_sample_is_deterministic() -> None:
    adata = _make_adata({"A": 40, "B": 35, "C": 25})
    first = sample_study_proportional(adata, n=20, stratifyKey="study_accession", seed=0)
    second = sample_study_proportional(adata, n=20, stratifyKey="study_accession", seed=0)
    assert list(first.obs_names) == list(second.obs_names)


def test_sample_proportional_prefers_larger_studies() -> None:
    adata = _make_adata({"A": 80, "B": 15, "C": 5})
    sampled = sample_study_proportional(adata, n=20, stratifyKey="study_accession", seed=0)
    counts = sampled.obs["study_accession"].value_counts()
    assert int(counts["A"]) > int(counts["B"])
    assert int(counts["B"]) >= int(counts["C"])
    assert int(counts["C"]) >= 1


def test_sample_rejects_missing_key() -> None:
    adata = _make_adata({"A": 5, "B": 5})
    with pytest.raises(ValueError, match="missing stratify key"):
        sample_study_proportional(adata, n=4, stratifyKey="missing")


def test_load_and_normalize_uses_preloaded_adata() -> None:
    adata = _make_adata({"A": 8, "B": 8})
    sampled = sample_study_proportional(adata, n=10, stratifyKey="study_accession", seed=0)
    cfg = AtlasPostprocessingConfig(batchKey="study_accession")
    with patch("atlas_postprocessing.core.sc.read_h5ad") as read_h5ad:
        out = load_and_normalize(cfg, adata=sampled)
    read_h5ad.assert_not_called()
    assert out.n_obs == 10
    assert out.raw is not None
    assert sample_metadata(out) is not None
