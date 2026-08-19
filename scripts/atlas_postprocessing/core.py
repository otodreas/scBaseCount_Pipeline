import json
import logging
import time
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import harmonypy
import numba
import numpy as np
import scanpy as sc
from cluster_validation import (
    MERGED_CLUSTER_KEY,
    RF_N_ESTIMATORS,
    RF_N_SPLITS,
    RF_RANDOM_STATE,
    apply_rf_merge,
)
from shared.repo import rel_to_repo
from sklearn.utils import check_array, check_random_state
from storage.r2 import upload_to_r2, verify_upload
from umap.umap_ import find_ab_params, simplicial_set_embedding

from atlas_postprocessing.artifacts import load_approved_parameters
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.plots import make_atlas_plots, save_scree_plot
from atlas_postprocessing.sampling import sample_metadata

log = logging.getLogger(__name__)

Workflow = Literal["production", "validation"]


def timed[T](name: str, action: Callable[[], T], logger: logging.Logger | None = None) -> T:
    """Run action, logging START/DONE with elapsed seconds, and return its result."""
    active_log = logger or log
    started = time.perf_counter()
    active_log.info("START %s", name)
    result = action()
    active_log.info("DONE %s in %.1fs", name, time.perf_counter() - started)
    return result


def apply_thread_settings(cfg: AtlasPostprocessingConfig) -> None:
    """Apply configured job count to scanpy and Numba when explicitly set."""
    if cfg.nJobs > 0:
        sc.settings.n_jobs = cfg.nJobs
        numba.set_num_threads(cfg.nJobs)
        log.info("Using nJobs=%s for scanpy and Numba", cfg.nJobs)


def validate_graph_settings(cfg: AtlasPostprocessingConfig, n_obs: int) -> None:
    """Fail early for impossible PC / neighbor settings."""
    if cfg.nPcs < 1:
        raise ValueError(f"nPcs must be >= 1, got {cfg.nPcs}")
    if cfg.nPcsCompute < cfg.nPcs:
        raise ValueError(f"nPcsCompute ({cfg.nPcsCompute}) must be >= nPcs ({cfg.nPcs})")
    if cfg.nNeighbors < 2:
        raise ValueError(f"nNeighbors must be >= 2, got {cfg.nNeighbors}")
    if cfg.nNeighbors >= n_obs:
        raise ValueError(f"nNeighbors ({cfg.nNeighbors}) must be < n_obs ({n_obs})")
    if cfg.nTopGenes < 1:
        raise ValueError(f"nTopGenes must be >= 1, got {cfg.nTopGenes}")


def load_and_normalize(
    cfg: AtlasPostprocessingConfig,
    adata: sc.AnnData | None = None,
) -> sc.AnnData:
    """Read the atlas h5ad (or use a preloaded object), stash ``.raw``, and normalize + log1p."""
    if adata is None:
        adata = sc.read_h5ad(cfg.inputH5ad)
    log.info("Loaded %s cells x %s genes", f"{adata.n_obs:,}", f"{adata.n_vars:,}")
    if cfg.batchKey not in adata.obs:
        raise ValueError(f"adata.obs is missing batch key {cfg.batchKey!r}")
    log.info("Studies (%s): %d", cfg.batchKey, adata.obs[cfg.batchKey].nunique())

    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def select_hvgs(adata: sc.AnnData, cfg: AtlasPostprocessingConfig, *, nTopGenes: int | None = None) -> sc.AnnData:
    """Select batch-aware HVGs and subset to those genes."""
    n_top = cfg.nTopGenes if nTopGenes is None else nTopGenes
    if n_top > adata.n_vars:
        raise ValueError(f"nTopGenes ({n_top}) exceeds available genes ({adata.n_vars})")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, batch_key=cfg.batchKey)
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.info("Retained %s HVGs", f"{adata.n_vars:,}")
    return adata


def scale_and_pca(adata: sc.AnnData, cfg: AtlasPostprocessingConfig, *, nPcsCompute: int | None = None) -> sc.AnnData:
    """Scale expression and compute PCA."""
    n_comps = cfg.nPcsCompute if nPcsCompute is None else nPcsCompute
    if n_comps < 1:
        raise ValueError(f"nPcsCompute must be >= 1, got {n_comps}")
    if n_comps > min(adata.n_obs - 1, adata.n_vars - 1):
        raise ValueError(
            f"nPcsCompute ({n_comps}) exceeds max feasible components "
            f"({min(adata.n_obs - 1, adata.n_vars - 1)}) for shape {(adata.n_obs, adata.n_vars)}"
        )
    sc.pp.scale(adata, max_value=cfg.scaleMaxValue)
    sc.pp.pca(adata, n_comps=n_comps, svd_solver="auto")
    log.info("Computed %s PCs", n_comps)
    return adata


def capture_pre_scale_features(adata: sc.AnnData) -> np.ndarray:
    """Dense copy of the current expression matrix before scaling mutates ``.X``."""
    x = adata.X
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray())
    return np.asarray(x).copy()


def prepare_pca(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    returnPreScaleFeatures: bool = False,
) -> sc.AnnData | tuple[sc.AnnData, np.ndarray]:
    """Select HVGs, scale, and compute PCA without building a neighbor graph.

    When ``returnPreScaleFeatures`` is true, also return a dense copy of the
    normalized HVG matrix captured before ``sc.pp.scale()``.
    """
    validate_graph_settings(cfg, adata.n_obs)
    adata = select_hvgs(adata, cfg)
    features = capture_pre_scale_features(adata) if returnPreScaleFeatures else None
    adata = scale_and_pca(adata, cfg)
    if returnPreScaleFeatures:
        assert features is not None
        return adata, features
    return adata


def build_neighbors(
    adata: sc.AnnData,
    *,
    nNeighbors: int,
    nPcs: int | None = None,
    useRep: str | None = None,
) -> None:
    """Build a neighbor graph with explicit neighborhood size and optional representation."""
    kwargs: dict[str, int | str] = {"n_neighbors": nNeighbors}
    if useRep is not None:
        kwargs["use_rep"] = useRep
        if nPcs is not None:
            kwargs["n_pcs"] = nPcs
    elif nPcs is not None:
        kwargs["n_pcs"] = nPcs
    sc.pp.neighbors(adata, **kwargs)  # type: ignore[arg-type]
    log.info(
        "Computed neighbors n_neighbors=%s n_pcs=%s use_rep=%s",
        nNeighbors,
        nPcs,
        useRep or "X_pca",
    )


def run_leiden(adata: sc.AnnData, *, resolution: float, keyAdded: str) -> None:
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added=keyAdded,
    )
    log.info("%s clusters: %d", keyAdded, adata.obs[keyAdded].nunique())


def run_umap_deterministic(adata: sc.AnnData) -> None:
    """Run Scanpy UMAP with the default fixed seed (reproducible, single-threaded)."""
    sc.tl.umap(adata)
    log.info("Computed deterministic UMAP")


def run_umap_parallel(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> None:
    """Embed an existing neighbor graph with UMAP's parallel optimizer.

    Parallel UMAP is not bit-reproducible across runs. Revisit before publication and
    switch back to a seeded single-threaded embedding when coordinates must be frozen.
    """
    if "neighbors" not in adata.uns:
        raise ValueError("adata.uns['neighbors'] is required; run sc.pp.neighbors first")
    if "connectivities" not in adata.obsp:
        raise ValueError("adata.obsp['connectivities'] is required; run sc.pp.neighbors first")

    if cfg.nJobs > 0:
        numba.set_num_threads(cfg.nJobs)

    # Match Scanpy defaults used by sc.tl.umap for large graphs.
    a, b = find_ab_params(spread=1.0, min_dist=0.5)
    n_epochs = 500 if adata.n_obs <= 10_000 else 200
    # random_state=None enables parallelism; results are intentionally non-reproducible.
    random_state = check_random_state(None)
    x = np.asarray(adata.obsm.get("X_pca_harmony", adata.obsm["X_pca"]), dtype=np.float32)
    x = check_array(x, dtype="float32", accept_sparse=False)

    connectivities = adata.obsp["connectivities"]
    graph = connectivities.tocoo() if hasattr(connectivities, "tocoo") else connectivities

    x_umap, _ = simplicial_set_embedding(
        data=x,
        graph=graph,
        n_components=2,
        initial_alpha=1.0,
        a=a,
        b=b,
        gamma=1.0,
        negative_sample_rate=5,
        n_epochs=n_epochs,
        init="spectral",
        random_state=random_state,
        metric="euclidean",
        metric_kwds={},
        densmap=False,
        densmap_kwds={},
        output_dens=False,
        parallel=True,
        verbose=False,
    )
    adata.obsm["X_umap"] = x_umap
    adata.uns["umap"] = {"params": {"a": a, "b": b, "parallel": True}}
    log.info("Computed parallel UMAP (non-reproducible; revisit for publication)")


def embed_uncorrected(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> sc.AnnData:
    """Build the pre-correction neighbor graph, UMAP, and Leiden partition on existing PCA."""
    validate_graph_settings(cfg, adata.n_obs)
    if "X_pca" not in adata.obsm:
        raise ValueError("adata.obsm['X_pca'] is required before uncorrected embedding")
    timed("uncorrected neighbors", lambda: build_neighbors(adata, nNeighbors=cfg.nNeighbors, nPcs=cfg.nPcs))
    timed("uncorrected UMAP", lambda: run_umap_deterministic(adata))
    adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()
    timed(
        "uncorrected Leiden",
        lambda: run_leiden(adata, resolution=cfg.resolution, keyAdded="leiden_uncorrected"),
    )
    return adata


def run_harmony_on_pcs(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    nPcs: int,
    batchKey: str | None = None,
) -> np.ndarray:
    """Run Harmony on the first ``nPcs`` PCA dimensions and return the corrected embedding."""
    if "X_pca" not in adata.obsm:
        raise ValueError("adata.obsm['X_pca'] is required before Harmony")
    if nPcs > adata.obsm["X_pca"].shape[1]:
        raise ValueError(f"Requested nPcs ({nPcs}) exceeds computed PCs ({adata.obsm['X_pca'].shape[1]})")
    resolved_batch_key = cfg.batchKey if batchKey is None else batchKey
    if resolved_batch_key not in adata.obs:
        raise ValueError(f"adata.obs is missing batch key {resolved_batch_key!r}")
    pca_prefix = np.asarray(adata.obsm["X_pca"][:, :nPcs])
    harmony_out = harmonypy.run_harmony(
        pca_prefix,
        adata.obs,
        resolved_batch_key,
        ncores=cfg.nJobs,
    )
    corrected = np.asarray(harmony_out.Z_corr)
    log.info("Ran Harmony on %s PCs with batch key %s ncores=%s", nPcs, resolved_batch_key, cfg.nJobs)
    return corrected


def integrate_harmony(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    parallelUmap: bool = False,
) -> sc.AnnData:
    """Run Harmony on the PCA embedding and build the corrected UMAP and Leiden partition."""
    validate_graph_settings(cfg, adata.n_obs)

    def _store_harmony() -> None:
        adata.obsm["X_pca_harmony"] = run_harmony_on_pcs(adata, cfg, nPcs=cfg.nPcs)

    timed("Harmony correction", _store_harmony)
    timed(
        "Harmony neighbors",
        lambda: build_neighbors(
            adata,
            nNeighbors=cfg.nNeighbors,
            nPcs=cfg.nPcs,
            useRep="X_pca_harmony",
        ),
    )
    if parallelUmap:
        timed("Harmony UMAP", lambda: run_umap_parallel(adata, cfg))
    else:
        timed("Harmony UMAP", lambda: run_umap_deterministic(adata))
    timed(
        "Harmony Leiden",
        lambda: run_leiden(adata, resolution=cfg.resolution, keyAdded="leiden_atlas"),
    )
    return adata


def save_atlas(adata: sc.AnnData, cfg: AtlasPostprocessingConfig, *, workflow: Workflow) -> Path:
    """Write the processed atlas h5ad and a run summary JSON alongside it."""
    cfg.outputH5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(cfg.outputH5ad, compression=cfg.compression)

    calibration_summary = None
    if cfg.parametersJson is not None:
        calibration_summary = load_approved_parameters(cfg.parametersJson).calibrationSummary

    clusters_uncorrected = int(adata.obs["leiden_uncorrected"].nunique()) if "leiden_uncorrected" in adata.obs else None
    clusters_merged = int(adata.obs["leiden_merged"].nunique()) if "leiden_merged" in adata.obs else None
    summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "output": rel_to_repo(cfg.outputH5ad),
        "figsDir": rel_to_repo(cfg.figsDir),
        "workflow": workflow,
        "cells": int(adata.n_obs),
        "hvgs": int(adata.n_vars),
        "rawGenes": int(adata.raw.n_vars) if adata.raw is not None else 0,
        "studies": int(adata.obs[cfg.batchKey].nunique()),
        "clustersUncorrected": clusters_uncorrected,
        "clustersHarmony": int(adata.obs["leiden_atlas"].nunique()),
        "clustersMerged": clusters_merged,
        "resolution": cfg.resolution,
        "nPcs": cfg.nPcs,
        "nPcsCompute": cfg.nPcsCompute,
        "nTopGenes": cfg.nTopGenes,
        "nNeighbors": cfg.nNeighbors,
        "nJobs": cfg.nJobs,
        "parametersJson": rel_to_repo(cfg.parametersJson) if cfg.parametersJson else None,
        "calibrationSummary": calibration_summary,
        "sampling": sample_metadata(adata),
    }
    summary_path = cfg.outputH5ad.with_name(f"{cfg.outputH5ad.stem}_run.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    log.info("Wrote %s", rel_to_repo(cfg.outputH5ad))
    log.info("Wrote %s", rel_to_repo(summary_path))
    return summary_path


def upload_atlas(cfg: AtlasPostprocessingConfig) -> None:
    """Upload the saved atlas h5ad to R2 and confirm the object is present."""
    assert cfg.r2Key is not None
    upload_to_r2(cfg.outputH5ad, cfg.r2Key)
    if not verify_upload(cfg.r2Key):
        raise RuntimeError(f"R2 upload verification failed for {cfg.r2Key}")


def run_postprocessing(
    cfg: AtlasPostprocessingConfig,
    adata: sc.AnnData | None = None,
    *,
    workflow: Workflow = "production",
) -> sc.AnnData:
    """Run atlas postprocessing for production (Harmony-only graph) or validation (both graphs)."""
    apply_thread_settings(cfg)
    loaded = timed("load + normalize", lambda: load_and_normalize(cfg, adata=adata))
    validate_graph_settings(cfg, loaded.n_obs)

    pre_scale_features: np.ndarray | None = None
    if workflow == "validation":
        prepared = timed(
            "HVG + PCA",
            lambda: prepare_pca(loaded, cfg, returnPreScaleFeatures=True),
        )
        assert isinstance(prepared, tuple)
        loaded, pre_scale_features = prepared
        loaded = timed("uncorrected embedding", lambda: embed_uncorrected(loaded, cfg))
        loaded = timed(
            "harmony integration",
            lambda: integrate_harmony(loaded, cfg, parallelUmap=False),
        )
        if cfg.cellTypeKey not in loaded.obs:
            raise ValueError(f"adata.obs is missing weak prior key {cfg.cellTypeKey!r}")
        assert pre_scale_features is not None
        weak_prior = loaded.obs[cfg.cellTypeKey].values if cfg.rfBalanceWeakPrior else None
        loaded, merge_info = timed(
            "RF merge approved Leiden partition",
            lambda: apply_rf_merge(
                loaded,
                featureMatrix=pre_scale_features,
                clusterKey="leiden_atlas",
                mergeThreshold=cfg.mergeThreshold,
                mergedKey=MERGED_CLUSTER_KEY,
                nEstimators=RF_N_ESTIMATORS,
                nSplits=RF_N_SPLITS,
                randomState=RF_RANDOM_STATE,
                weakPriorLabels=weak_prior,
                minCellsPerCluster=3,
            ),
        )
        loaded.uns["rfMerge"] = {
            "clusterKey": "leiden_atlas",
            "mergedKey": MERGED_CLUSTER_KEY,
            "mergeThreshold": cfg.mergeThreshold,
            "rfBalanceWeakPrior": cfg.rfBalanceWeakPrior,
            "nEstimators": RF_N_ESTIMATORS,
            "nSplits": RF_N_SPLITS,
            "randomState": RF_RANDOM_STATE,
            "nClustersPreMerge": merge_info.nClustersPreMerge,
            "nClustersPostMerge": merge_info.nClustersPostMerge,
            "labelMap": merge_info.labelMap,
            "mergedGroups": merge_info.mergedGroups,
            "confClasses": [str(c) for c in merge_info.classes],
            "confMatrix": merge_info.conf.tolist(),
            "featureMatrix": "normalized_log1p_hvg_before_scale",
        }
    elif workflow == "production":
        loaded = timed("HVG + PCA", lambda: prepare_pca(loaded, cfg))
        assert isinstance(loaded, sc.AnnData)
        loaded = timed(
            "harmony integration",
            lambda: integrate_harmony(loaded, cfg, parallelUmap=True),
        )
    else:
        raise ValueError(f"Unknown workflow {workflow!r}")

    if cfg.writePlots:
        timed("scree plot", lambda: save_scree_plot(loaded, cfg))
        timed("plots", lambda: make_atlas_plots(loaded, cfg, workflow=workflow))
    timed("save atlas", lambda: save_atlas(loaded, cfg, workflow=workflow))
    if cfg.r2Key:
        timed("upload to r2", lambda: upload_atlas(cfg))
    return loaded
