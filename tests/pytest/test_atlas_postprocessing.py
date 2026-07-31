from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.artifacts import reject_tuning_overrides_with_parameters_json
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import build_neighbors, integrate_harmony, validate_graph_settings
from atlas_postprocessing.plots import make_atlas_plots
from scipy.sparse import csr_matrix


def test_production_rejects_parameters_json_with_scalar_overrides() -> None:
    with pytest.raises(ValueError, match="parameters-json"):
        reject_tuning_overrides_with_parameters_json(["n_top_genes"])


def test_validate_graph_settings_rejects_impossible_neighbors() -> None:
    cfg = AtlasPostprocessingConfig(nNeighbors=100, nPcs=10, nPcsCompute=10)
    with pytest.raises(ValueError, match="nNeighbors"):
        validate_graph_settings(cfg, n_obs=50)


def test_build_neighbors_passes_harmony_graph_args() -> None:
    adata = sc.AnnData(np.zeros((20, 5), dtype=np.float32))
    adata.obsm["X_pca_harmony"] = np.random.default_rng(0).normal(size=(20, 8)).astype(np.float32)

    with patch("atlas_postprocessing.core.sc.pp.neighbors") as neighbors:
        build_neighbors(adata, nNeighbors=7, nPcs=4, useRep="X_pca_harmony")

    neighbors.assert_called_once_with(adata, n_neighbors=7, use_rep="X_pca_harmony", n_pcs=4)


def test_integrate_harmony_uses_configured_neighbors_and_pcs() -> None:
    adata = sc.AnnData(np.zeros((30, 10), dtype=np.float32))
    adata.obs = pd.DataFrame({"study_accession": ["A"] * 15 + ["B"] * 15})
    adata.obsm["X_pca"] = np.random.default_rng(1).normal(size=(30, 12)).astype(np.float32)
    cfg = AtlasPostprocessingConfig(nPcs=5, nNeighbors=8, resolution=0.4)

    fake_harmony = MagicMock()
    fake_harmony.Z_corr = np.random.default_rng(2).normal(size=(30, 5)).astype(np.float32)

    with (
        patch("atlas_postprocessing.core.harmonypy.run_harmony", return_value=fake_harmony) as run_harmony,
        patch("atlas_postprocessing.core.build_neighbors") as neighbors,
        patch("atlas_postprocessing.core.sc.tl.umap"),
        patch("atlas_postprocessing.core.run_leiden"),
    ):
        integrate_harmony(adata, cfg)

    assert run_harmony.call_args.args[0].shape == (30, 5)
    neighbors.assert_called_once_with(
        adata,
        nNeighbors=8,
        nPcs=5,
        useRep="X_pca_harmony",
    )


def test_make_atlas_plots_uses_scanpy_embedding_args(tmp_path: Path) -> None:
    n_obs = 12
    adata = sc.AnnData(csr_matrix(np.zeros((n_obs, 4), dtype=np.float32)))
    adata.obs = pd.DataFrame(
        {
            "study_accession": ["A"] * 6 + ["B"] * 6,
            "cell_type": ["t"] * n_obs,
        },
        index=[f"c{i}" for i in range(n_obs)],
    )
    coords = np.random.default_rng(0).normal(size=(n_obs, 2)).astype(np.float32)
    adata.obsm["X_umap"] = coords
    adata.obsm["X_umap_uncorrected"] = coords + 0.1
    cfg = AtlasPostprocessingConfig(figsDir=tmp_path)

    fake_fig = MagicMock()
    with patch("atlas_postprocessing.plots.sc.pl.embedding", return_value=fake_fig) as embedding:
        make_atlas_plots(adata, cfg)

    assert embedding.call_count == 4
    for call in embedding.call_args_list:
        assert call.kwargs["alpha"] == 0.25
        assert call.kwargs["size"] == 0.001
        assert call.kwargs["legend_loc"] is None
        assert call.kwargs["show"] is False
        assert call.kwargs["return_fig"] is True
    assert (tmp_path / "umap_study_accession_uncorrected.png").exists() or fake_fig.savefig.called
