import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.artifacts import reject_tuning_overrides_with_parameters_json
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import (
    build_neighbors,
    embed_uncorrected,
    integrate_harmony,
    run_postprocessing,
    scale_and_pca,
    validate_graph_settings,
)
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


@pytest.mark.parametrize("parallel_umap", [False, True])
def test_embed_uncorrected_uses_configured_umap_mode(parallel_umap: bool) -> None:
    adata = sc.AnnData(np.zeros((20, 5), dtype=np.float32))
    adata.obsm["X_pca"] = np.zeros((20, 4), dtype=np.float32)
    adata.obsm["X_umap"] = np.zeros((20, 2), dtype=np.float32)
    cfg = AtlasPostprocessingConfig(
        nPcs=4,
        nPcsCompute=4,
        nNeighbors=5,
        umapParallel=parallel_umap,
    )

    with (
        patch("atlas_postprocessing.core.build_neighbors"),
        patch("atlas_postprocessing.core.run_umap_deterministic") as deterministic,
        patch("atlas_postprocessing.core.run_umap_parallel") as parallel,
        patch("atlas_postprocessing.core.run_leiden"),
    ):
        embed_uncorrected(adata, cfg)

    if parallel_umap:
        parallel.assert_called_once_with(adata, cfg)
        deterministic.assert_not_called()
    else:
        deterministic.assert_called_once_with(adata)
        parallel.assert_not_called()
    np.testing.assert_array_equal(adata.obsm["X_umap_uncorrected"], adata.obsm["X_umap"])


def test_scale_and_pca_uses_auto_svd_solver() -> None:
    adata = sc.AnnData(np.random.default_rng(0).normal(size=(40, 25)).astype(np.float32))
    cfg = AtlasPostprocessingConfig(nPcsCompute=5, nPcs=5, nTopGenes=25)

    with (
        patch("atlas_postprocessing.core.sc.pp.scale"),
        patch("atlas_postprocessing.core.sc.pp.pca") as pca,
    ):
        scale_and_pca(adata, cfg)

    pca.assert_called_once()
    assert pca.call_args.kwargs["svd_solver"] == "auto"
    assert pca.call_args.kwargs["n_comps"] == 5


def test_integrate_harmony_uses_configured_neighbors_pcs_and_threads() -> None:
    adata = sc.AnnData(np.zeros((30, 10), dtype=np.float32))
    adata.obs = pd.DataFrame({"study_accession": ["A"] * 15 + ["B"] * 15})
    adata.obsm["X_pca"] = np.random.default_rng(1).normal(size=(30, 12)).astype(np.float32)
    cfg = AtlasPostprocessingConfig(nPcs=5, nNeighbors=8, resolution=0.4, nJobs=4, umapParallel=False)

    fake_harmony = MagicMock()
    fake_harmony.Z_corr = np.random.default_rng(2).normal(size=(30, 5)).astype(np.float32)

    with (
        patch("atlas_postprocessing.core.harmonypy.run_harmony", return_value=fake_harmony) as run_harmony,
        patch("atlas_postprocessing.core.build_neighbors") as neighbors,
        patch("atlas_postprocessing.core.run_umap_deterministic") as umap_det,
        patch("atlas_postprocessing.core.run_umap_parallel") as umap_par,
        patch("atlas_postprocessing.core.run_leiden"),
    ):
        integrate_harmony(adata, cfg)

    assert run_harmony.call_args.args[0].shape == (30, 5)
    assert run_harmony.call_args.kwargs["ncores"] == 4
    neighbors.assert_called_once_with(
        adata,
        nNeighbors=8,
        nPcs=5,
        useRep="X_pca_harmony",
    )
    umap_det.assert_called_once_with(adata)
    umap_par.assert_not_called()


def test_integrate_harmony_uses_parallel_umap_when_requested() -> None:
    adata = sc.AnnData(np.zeros((30, 10), dtype=np.float32))
    adata.obs = pd.DataFrame({"study_accession": ["A"] * 15 + ["B"] * 15})
    adata.obsm["X_pca"] = np.random.default_rng(1).normal(size=(30, 12)).astype(np.float32)
    cfg = AtlasPostprocessingConfig(nPcs=5, nNeighbors=8, resolution=0.4, nJobs=2, umapParallel=True)

    fake_harmony = MagicMock()
    fake_harmony.Z_corr = np.random.default_rng(2).normal(size=(30, 5)).astype(np.float32)

    with (
        patch("atlas_postprocessing.core.harmonypy.run_harmony", return_value=fake_harmony),
        patch("atlas_postprocessing.core.build_neighbors"),
        patch("atlas_postprocessing.core.run_umap_deterministic") as umap_det,
        patch("atlas_postprocessing.core.run_umap_parallel") as umap_par,
        patch("atlas_postprocessing.core.run_leiden"),
    ):
        integrate_harmony(adata, cfg)

    umap_par.assert_called_once_with(adata, cfg)
    umap_det.assert_not_called()


def _tiny_counts_adata(n_obs: int = 24, n_vars: int = 40) -> sc.AnnData:
    rng = np.random.default_rng(0)
    counts = csr_matrix(rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32))
    adata = sc.AnnData(
        X=counts,
        obs=pd.DataFrame(
            {
                "study_accession": (["A"] * (n_obs // 2)) + (["B"] * (n_obs - n_obs // 2)),
                "cell_type": (["t1"] * (n_obs // 2)) + (["t2"] * (n_obs - n_obs // 2)),
            },
            index=[f"c{i}" for i in range(n_obs)],
        ),
    )
    adata.var_names = [f"g{i}" for i in range(n_vars)]
    return adata


def test_production_workflow_runs_uncorrected_embedding(tmp_path: Path) -> None:
    cfg = AtlasPostprocessingConfig(
        inputH5ad=tmp_path / "unused.h5ad",
        outputH5ad=tmp_path / "out.h5ad",
        figsDir=tmp_path / "figs",
        nTopGenes=15,
        nPcs=5,
        nPcsCompute=5,
        nNeighbors=5,
        resolution=0.5,
        writePlots=False,
        nJobs=1,
    )
    adata = _tiny_counts_adata()

    def _return_loaded(loaded: sc.AnnData, _cfg: AtlasPostprocessingConfig) -> sc.AnnData:
        return loaded

    with (
        patch("atlas_postprocessing.core.embed_uncorrected", side_effect=_return_loaded) as uncorrected,
        patch("atlas_postprocessing.core.harmonypy.run_harmony") as run_harmony,
        patch("atlas_postprocessing.core.build_neighbors") as neighbors,
        patch("atlas_postprocessing.core.run_umap_parallel") as umap_par,
        patch("atlas_postprocessing.core.run_umap_deterministic") as umap_det,
        patch("atlas_postprocessing.core.run_leiden") as leiden,
        patch("atlas_postprocessing.core.prepare_pca", side_effect=lambda loaded, _cfg: loaded) as prep,
    ):
        run_harmony.return_value = MagicMock(Z_corr=np.zeros((adata.n_obs, cfg.nPcs), dtype=np.float32))
        adata.obsm["X_pca"] = np.zeros((adata.n_obs, cfg.nPcsCompute), dtype=np.float32)
        adata.obs["leiden_atlas"] = "0"
        run_postprocessing(cfg, adata=adata, workflow="production")

    prep.assert_called_once()
    uncorrected.assert_called_once_with(adata, cfg)
    neighbors.assert_called_once()
    umap_par.assert_not_called()
    umap_det.assert_called_once_with(adata)
    leiden.assert_called_once()
    summary = json.loads((tmp_path / "out_run.json").read_text())
    assert summary["workflow"] == "production"


def test_validation_workflow_keeps_both_graphs(tmp_path: Path) -> None:
    cfg = AtlasPostprocessingConfig(
        inputH5ad=tmp_path / "unused.h5ad",
        outputH5ad=tmp_path / "out.h5ad",
        figsDir=tmp_path / "figs",
        nTopGenes=15,
        nPcs=5,
        nPcsCompute=5,
        nNeighbors=5,
        resolution=0.5,
        writePlots=False,
        nJobs=1,
    )
    adata = _tiny_counts_adata()
    adata.obsm["X_pca"] = np.zeros((adata.n_obs, cfg.nPcsCompute), dtype=np.float32)

    def _fake_harmony(loaded: sc.AnnData, _cfg: AtlasPostprocessingConfig) -> sc.AnnData:
        loaded.obs["leiden_atlas"] = pd.Categorical(["0"] * loaded.n_obs)
        loaded.obs["leiden_uncorrected"] = "0"
        return loaded

    with (
        patch("atlas_postprocessing.core.prepare_pca", side_effect=lambda loaded, _cfg: loaded),
        patch("atlas_postprocessing.core.embed_uncorrected", side_effect=lambda loaded, _cfg: loaded) as uncorrected,
        patch("atlas_postprocessing.core.integrate_harmony", side_effect=_fake_harmony) as harmony,
    ):
        result = run_postprocessing(cfg, adata=adata, workflow="validation")

    uncorrected.assert_called_once()
    harmony.assert_called_once_with(adata, cfg)
    assert "leiden_merged" not in result.obs
    assert "rfMerge" not in result.uns
    summary = json.loads((tmp_path / "out_run.json").read_text())
    assert summary["workflow"] == "validation"
    assert summary["clustersMerged"] is None


def test_make_atlas_plots_production_writes_both_embeddings(tmp_path: Path) -> None:
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
        make_atlas_plots(adata, cfg, workflow="production")

    assert embedding.call_count == 4
    bases = [call.kwargs["basis"] for call in embedding.call_args_list]
    assert bases.count("X_umap_uncorrected") == 2
    assert bases.count("X_umap") == 2


def test_make_atlas_plots_validation_writes_both(tmp_path: Path) -> None:
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
        make_atlas_plots(adata, cfg, workflow="validation")

    assert embedding.call_count == 4
    bases = [call.kwargs["basis"] for call in embedding.call_args_list]
    assert bases.count("X_umap_uncorrected") == 2
    assert bases.count("X_umap") == 2
