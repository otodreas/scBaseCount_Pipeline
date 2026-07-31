import json
import logging
import time
from collections.abc import Callable
from pathlib import Path

import harmonypy
import numpy as np
import scanpy as sc
from shared.repo import rel_to_repo
from storage.r2 import upload_to_r2, verify_upload

from atlas_postprocessing.artifacts import load_approved_parameters
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.plots import make_atlas_plots, save_scree_plot
from atlas_postprocessing.sampling import sample_metadata

log = logging.getLogger(__name__)


def timed[T](name: str, action: Callable[[], T], logger: logging.Logger | None = None) -> T:
    """Run action, logging START/DONE with elapsed seconds, and return its result."""
    active_log = logger or log
    started = time.perf_counter()
    active_log.info("START %s", name)
    result = action()
    active_log.info("DONE %s in %.1fs", name, time.perf_counter() - started)
    return result


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
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    log.info("Computed %s PCs", n_comps)
    return adata


def build_neighbors(
    adata: sc.AnnData,
    *,
    nNeighbors: int,
    nPcs: int | None = None,
    useRep: str | None = None,
) -> None:
    """Build a neighbor graph with explicit neighborhood size and optional representation."""
    kwargs: dict[str, object] = {"n_neighbors": nNeighbors}
    if useRep is not None:
        kwargs["use_rep"] = useRep
        if nPcs is not None:
            kwargs["n_pcs"] = nPcs
    elif nPcs is not None:
        kwargs["n_pcs"] = nPcs
    sc.pp.neighbors(adata, **kwargs)
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


def embed_uncorrected(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> sc.AnnData:
    """Select HVGs, scale, run PCA, and build the pre-correction UMAP and leiden partition."""
    validate_graph_settings(cfg, adata.n_obs)
    adata = select_hvgs(adata, cfg)
    adata = scale_and_pca(adata, cfg)
    build_neighbors(adata, nNeighbors=cfg.nNeighbors, nPcs=cfg.nPcs)
    sc.tl.umap(adata)
    adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()
    log.info("Computed UMAP")
    run_leiden(adata, resolution=cfg.resolution, keyAdded="leiden_uncorrected")
    return adata


def run_harmony_on_pcs(adata: sc.AnnData, cfg: AtlasPostprocessingConfig, *, nPcs: int) -> np.ndarray:
    """Run Harmony on the first ``nPcs`` PCA dimensions and return the corrected embedding."""
    if "X_pca" not in adata.obsm:
        raise ValueError("adata.obsm['X_pca'] is required before Harmony")
    if nPcs > adata.obsm["X_pca"].shape[1]:
        raise ValueError(f"Requested nPcs ({nPcs}) exceeds computed PCs ({adata.obsm['X_pca'].shape[1]})")
    pca_prefix = np.asarray(adata.obsm["X_pca"][:, :nPcs])
    harmony_out = harmonypy.run_harmony(pca_prefix, adata.obs, cfg.batchKey)
    corrected = np.asarray(harmony_out.Z_corr)
    log.info("Ran Harmony on %s PCs with batch key %s", nPcs, cfg.batchKey)
    return corrected


def integrate_harmony(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> sc.AnnData:
    """Run Harmony on the PCA embedding and build the corrected UMAP and leiden partition."""
    validate_graph_settings(cfg, adata.n_obs)
    adata.obsm["X_pca_harmony"] = run_harmony_on_pcs(adata, cfg, nPcs=cfg.nPcs)
    build_neighbors(
        adata,
        nNeighbors=cfg.nNeighbors,
        nPcs=cfg.nPcs,
        useRep="X_pca_harmony",
    )
    sc.tl.umap(adata)
    log.info("Computed UMAP using Harmony-corrected PCA")
    run_leiden(adata, resolution=cfg.resolution, keyAdded="leiden_atlas")
    return adata


def save_atlas(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> Path:
    """Write the processed atlas h5ad and a run summary JSON alongside it."""
    cfg.outputH5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(cfg.outputH5ad, compression=cfg.compression)

    calibration_summary = None
    if cfg.parametersJson is not None:
        calibration_summary = load_approved_parameters(cfg.parametersJson).calibrationSummary

    summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "output": rel_to_repo(cfg.outputH5ad),
        "figsDir": rel_to_repo(cfg.figsDir),
        "cells": int(adata.n_obs),
        "hvgs": int(adata.n_vars),
        "rawGenes": int(adata.raw.n_vars) if adata.raw is not None else 0,
        "studies": int(adata.obs[cfg.batchKey].nunique()),
        "clustersUncorrected": int(adata.obs["leiden_uncorrected"].nunique()),
        "clustersHarmony": int(adata.obs["leiden_atlas"].nunique()),
        "resolution": cfg.resolution,
        "nPcs": cfg.nPcs,
        "nPcsCompute": cfg.nPcsCompute,
        "nTopGenes": cfg.nTopGenes,
        "nNeighbors": cfg.nNeighbors,
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
) -> sc.AnnData:
    """Run the full atlas postprocessing flow: load, embed, integrate, plot, and save."""
    loaded = timed("load + normalize", lambda: load_and_normalize(cfg, adata=adata))
    validate_graph_settings(cfg, loaded.n_obs)
    loaded = timed("HVG + PCA + pre-correction embedding", lambda: embed_uncorrected(loaded, cfg))
    loaded = timed("harmony integration", lambda: integrate_harmony(loaded, cfg))

    if cfg.writePlots:
        timed("scree plot", lambda: save_scree_plot(loaded, cfg))
        timed("plots", lambda: make_atlas_plots(loaded, cfg))
    timed("save atlas", lambda: save_atlas(loaded, cfg))
    if cfg.r2Key:
        timed("upload to r2", lambda: upload_atlas(cfg))
    return loaded
