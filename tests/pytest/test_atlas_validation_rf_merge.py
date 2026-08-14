from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import capture_pre_scale_features, prepare_pca, run_postprocessing
from cluster_validation import MERGED_CLUSTER_KEY
from scipy import sparse


def _validation_counts_adata(n_cells: int = 30, n_genes: int = 40) -> sc.AnnData:
    rng = np.random.default_rng(4)
    x = np.abs(rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32))
    x[:15] += 5.0
    adata = sc.AnnData(sparse.csr_matrix(x))
    adata.obs = pd.DataFrame(
        {
            "study_accession": ["S1"] * 15 + ["S2"] * 15,
            "cell_type": ["A"] * 15 + ["B"] * 15,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_prepare_pca_can_return_pre_scale_features() -> None:
    adata = _validation_counts_adata()
    cfg = AtlasPostprocessingConfig(nTopGenes=20, nPcs=5, nPcsCompute=5, nNeighbors=5)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    prepared, features = prepare_pca(adata, cfg, returnPreScaleFeatures=True)
    assert isinstance(prepared, sc.AnnData)
    assert features.shape == (prepared.n_obs, prepared.n_vars)
    assert not np.allclose(np.asarray(prepared.X), features)


def test_capture_pre_scale_features_densifies_sparse() -> None:
    adata = sc.AnnData(sparse.csr_matrix(np.arange(12, dtype=np.float32).reshape(3, 4)))
    dense = capture_pre_scale_features(adata)
    assert isinstance(dense, np.ndarray)
    assert dense.shape == (3, 4)


def test_validation_workflow_applies_rf_merge(tmp_path: Path) -> None:
    cfg = AtlasPostprocessingConfig(
        inputH5ad=tmp_path / "unused.h5ad",
        outputH5ad=tmp_path / "out.h5ad",
        figsDir=tmp_path / "figs",
        nTopGenes=20,
        nPcs=5,
        nPcsCompute=5,
        nNeighbors=5,
        resolution=0.5,
        writePlots=False,
        nJobs=1,
        mergeThreshold=0.2,
    )
    adata = _validation_counts_adata()

    def _fake_prepare(loaded, _cfg, returnPreScaleFeatures=False):
        features = np.asarray(loaded.X.toarray() if hasattr(loaded.X, "toarray") else loaded.X).copy()
        loaded.obsm["X_pca"] = np.zeros((loaded.n_obs, _cfg.nPcsCompute), dtype=np.float32)
        if returnPreScaleFeatures:
            return loaded, features
        return loaded

    def _fake_uncorrected(loaded, _cfg):
        loaded.obs["leiden_uncorrected"] = "0"
        loaded.obsm["X_umap_uncorrected"] = np.zeros((loaded.n_obs, 2), dtype=np.float32)
        return loaded

    def _fake_harmony(loaded, _cfg, parallelUmap=False):
        del parallelUmap
        labels = np.array(["0"] * 15 + ["1"] * 15, dtype=object)
        loaded.obs["leiden_atlas"] = pd.Categorical(labels)
        loaded.obsm["X_umap"] = np.zeros((loaded.n_obs, 2), dtype=np.float32)
        return loaded

    with (
        patch("atlas_postprocessing.core.prepare_pca", side_effect=_fake_prepare),
        patch("atlas_postprocessing.core.embed_uncorrected", side_effect=_fake_uncorrected),
        patch("atlas_postprocessing.core.integrate_harmony", side_effect=_fake_harmony),
        patch("atlas_postprocessing.core.apply_thread_settings"),
        patch("atlas_postprocessing.core.save_atlas", return_value=tmp_path / "out_run.json"),
    ):
        result = run_postprocessing(cfg, adata=adata, workflow="validation")

    assert MERGED_CLUSTER_KEY in result.obs
    assert result.obs[MERGED_CLUSTER_KEY].isna().sum() == 0
    assert "rfMerge" in result.uns
    assert result.uns["rfMerge"]["featureMatrix"] == "normalized_log1p_hvg_before_scale"
    assert result.uns["rfMerge"]["nClustersPreMerge"] == 2
    assert set(result.uns["rfMerge"]["confClasses"]) == {"0", "1"}


def test_validation_workflow_rejects_undersized_clusters(tmp_path: Path) -> None:
    cfg = AtlasPostprocessingConfig(
        inputH5ad=tmp_path / "unused.h5ad",
        outputH5ad=tmp_path / "out.h5ad",
        figsDir=tmp_path / "figs",
        nTopGenes=20,
        nPcs=5,
        nPcsCompute=5,
        nNeighbors=5,
        resolution=0.5,
        writePlots=False,
        nJobs=1,
    )
    adata = _validation_counts_adata()

    def _fake_prepare(loaded, _cfg, returnPreScaleFeatures=False):
        features = np.zeros((loaded.n_obs, 10), dtype=np.float32)
        loaded.obsm["X_pca"] = np.zeros((loaded.n_obs, _cfg.nPcsCompute), dtype=np.float32)
        if returnPreScaleFeatures:
            return loaded, features
        return loaded

    def _fake_harmony(loaded, _cfg, parallelUmap=False):
        del parallelUmap
        labels = np.array(["0"] * 28 + ["1", "2"], dtype=object)
        loaded.obs["leiden_atlas"] = pd.Categorical(labels)
        return loaded

    with (
        patch("atlas_postprocessing.core.prepare_pca", side_effect=_fake_prepare),
        patch("atlas_postprocessing.core.embed_uncorrected", side_effect=lambda loaded, _cfg: loaded),
        patch("atlas_postprocessing.core.integrate_harmony", side_effect=_fake_harmony),
        patch("atlas_postprocessing.core.apply_thread_settings"),
        pytest.raises(ValueError, match="at least 3 cells"),
    ):
        run_postprocessing(cfg, adata=adata, workflow="validation")
