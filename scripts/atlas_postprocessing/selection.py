import logging
import time
from typing import Any

import scanpy as sc
from cluster_validation import (
    ResolutionSelection,
    default_resolutions,
    pick_n_pcs,
    select_resolution_on_graph,
)
from shared.repo import rel_to_repo

from atlas_postprocessing.artifacts import write_json, write_metric_csv, write_parameters_template
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import (
    build_neighbors,
    load_and_normalize,
    run_harmony_on_pcs,
    scale_and_pca,
    select_hvgs,
)
from atlas_postprocessing.plots import plot_resolution_selection
from atlas_postprocessing.sampling import sample_metadata

log = logging.getLogger(__name__)

FIXED_N_TOP_GENES = 2000
FIXED_N_NEIGHBORS = 15


def _require_cell_type(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> None:
    if cfg.cellTypeKey not in adata.obs:
        raise ValueError(f"adata.obs is missing cell type key {cfg.cellTypeKey!r}")


def _resolution_rows(sel: ResolutionSelection, resolutions: list[float]) -> list[dict[str, Any]]:
    """Convert resolution selection to a list of rows for writing to a CSV file."""
    rows: list[dict[str, Any]] = []
    for idx, resolution in enumerate(resolutions):
        rows.append(
            {
                "resolution": resolution,
                "matchedJaccard": float(sel.jaccArr[idx]),
                "nClusters": int(sel.kArr[idx]),
            }
        )
    return rows


def run_calibration(
    cfg: AtlasPostprocessingConfig,
    adata: sc.AnnData | None = None,
) -> dict[str, Any]:
    """Calibrate atlas graph parameters with cluster-validation selection on a Harmony graph."""
    if cfg.nTopGenes != FIXED_N_TOP_GENES:
        raise ValueError(f"Atlas calibration fixes nTopGenes={FIXED_N_TOP_GENES}, got {cfg.nTopGenes}")
    if cfg.nNeighbors != FIXED_N_NEIGHBORS:
        raise ValueError(f"Atlas calibration fixes nNeighbors={FIXED_N_NEIGHBORS}, got {cfg.nNeighbors}")
    resolutions = list(cfg.resolutionCandidates) if cfg.resolutionCandidates else default_resolutions()
    if len(resolutions) < 2:
        raise ValueError(f"resolutionCandidates must contain at least two values, got {resolutions!r}")

    cfg.calibrationDir.mkdir(parents=True, exist_ok=True)
    metrics_dir = cfg.calibrationDir / "metrics"
    figures_dir = cfg.calibrationDir / "figures"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Allows the user to stay in the loop and approve the parameters before running the pipeline.
    log.info(
        "Weak-prior labels (%s) are not ground truth; the matched-Jaccard argmax is advisory only",
        cfg.cellTypeKey,
    )

    started_all = time.perf_counter()
    adata_norm = load_and_normalize(cfg, adata=adata)
    _require_cell_type(adata_norm, cfg)
    if cfg.nNeighbors >= adata_norm.n_obs:
        raise ValueError(f"nNeighbors ({cfg.nNeighbors}) must be < n_obs ({adata_norm.n_obs})")

    log.info("START fixed HVG selection nTopGenes=%s", cfg.nTopGenes)
    adata_hvg = select_hvgs(adata_norm, cfg, nTopGenes=cfg.nTopGenes)
    log.info("DONE fixed HVG selection")

    log.info("START scale + PCA nPcsCompute=%s", cfg.nPcsCompute)
    adata_pca = scale_and_pca(adata_hvg, cfg, nPcsCompute=cfg.nPcsCompute)
    n_pcs, cumvar = pick_n_pcs(
        adata_pca.uns["pca"]["variance_ratio"],
        nPcsMin=cfg.nPcsMin,
        nPcsCompute=cfg.nPcsCompute,
        nPcsCumvarTarget=cfg.nPcsCumvarTarget,
    )
    log.info("DONE PCA; selected nPcs=%s cumvar=%.2f%%", n_pcs, cumvar)

    log.info("START Harmony + neighbors for resolution selection")
    adata_pca.obsm["X_pca_harmony"] = run_harmony_on_pcs(adata_pca, cfg, nPcs=n_pcs)
    build_neighbors(
        adata_pca,
        nNeighbors=cfg.nNeighbors,
        nPcs=n_pcs,
        useRep="X_pca_harmony",
    )
    log.info("DONE Harmony + neighbors")

    log.info("START shared resolution selection (%s candidates)", len(resolutions))
    adata_pca, sel = select_resolution_on_graph(
        adata_pca,
        resolutions=resolutions,
        weakPriorKey=cfg.cellTypeKey,
    )
    log.info(
        "DONE resolution selection; recommended=%s matchedJaccard=%.4f nClusters=%s",
        sel.selectedResolution,
        float(sel.jaccArr[sel.bestIdx]),
        int(sel.kArr[sel.bestIdx]),
    )

    # Write resolution selection metrics to CSV.
    rows = _resolution_rows(sel, resolutions)
    write_metric_csv(
        metrics_dir / "resolution.csv",
        rows,
        ["resolution", "matchedJaccard", "nClusters"],
    )
    # Plot resolution selection metrics.
    plot_resolution_selection(
        resolutions=resolutions,
        jaccArr=sel.jaccArr.tolist(),
        kArr=sel.kArr.tolist(),
        selectedResolution=sel.selectedResolution,
        outPath=figures_dir / "resolution_matched_jaccard.png",
    )

    # Keep candidates compatible with validate_approved_against_calibration:
    # fixed graph values are singletons; any evaluated resolution may be approved.
    candidates = {
        "hvg": [cfg.nTopGenes],
        "pc": [n_pcs],
        "neighbors": [cfg.nNeighbors],
        "resolution": resolutions,
    }
    recommendation = {
        "nTopGenes": cfg.nTopGenes,
        "nPcs": n_pcs,
        "nNeighbors": cfg.nNeighbors,
        "resolution": sel.selectedResolution,
        "matchedJaccard": float(sel.jaccArr[sel.bestIdx]),
        "nClusters": int(sel.kArr[sel.bestIdx]),
        "method": "matched_jaccard_argmax",
        "note": "Advisory only; copy parameters_template.json to approved_parameters.json after review.",
    }
    summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "calibrationDir": rel_to_repo(cfg.calibrationDir),
        "cells": int(adata_norm.n_obs),
        "genes": int(adata_norm.n_vars),
        "hvgs": int(adata_hvg.n_vars),
        "sampling": sample_metadata(adata_norm),
        "batchKey": cfg.batchKey,
        "cellTypeKey": cfg.cellTypeKey,
        "graphMethod": {
            "nTopGenes": cfg.nTopGenes,
            "nNeighbors": cfg.nNeighbors,
            "nPcsCompute": cfg.nPcsCompute,
            "nPcsMin": cfg.nPcsMin,
            "nPcsCumvarTarget": cfg.nPcsCumvarTarget,
            "selectedNPcs": n_pcs,
            "selectedCumvarPercent": cumvar,
            "batchCorrection": "harmony",
            "harmonyBatchKey": cfg.batchKey,
            "useRep": "X_pca_harmony",
        },
        "baseline": {
            "nTopGenes": cfg.nTopGenes,
            "nPcs": n_pcs,
            "nNeighbors": cfg.nNeighbors,
            "resolution": sel.selectedResolution,
            "nPcsCompute": cfg.nPcsCompute,
        },
        "candidates": candidates,
        "recommendation": recommendation,
        "metrics": {
            "resolution": rel_to_repo(metrics_dir / "resolution.csv"),
        },
        "figures": {
            "resolution": rel_to_repo(figures_dir / "resolution_matched_jaccard.png"),
        },
        "timingsSeconds": round(time.perf_counter() - started_all, 3),
        "note": (
            "Cell-type labels are weak priors, not ground truth. "
            "The matched-Jaccard argmax is advisory. Copy parameters_template.json to "
            "approved_parameters.json only after review."
        ),
    }
    summary_path = cfg.calibrationDir / "calibration_summary.json"
    write_json(summary_path, summary)

    template_cfg = cfg.model_copy(
        update={
            "nTopGenes": cfg.nTopGenes,
            "nPcs": n_pcs,
            "nNeighbors": cfg.nNeighbors,
            "resolution": sel.selectedResolution,
        }
    )
    write_parameters_template(template_cfg, summary_path)
    log.info("Calibration complete. Outputs under %s", rel_to_repo(cfg.calibrationDir))
    log.info(
        "Recommended parameters: nTopGenes=%s nPcs=%s nNeighbors=%s resolution=%s",
        recommendation["nTopGenes"],
        recommendation["nPcs"],
        recommendation["nNeighbors"],
        recommendation["resolution"],
    )
    return summary
