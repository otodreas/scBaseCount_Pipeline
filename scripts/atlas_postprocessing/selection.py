import logging
import time
from typing import Any

import scanpy as sc
from cluster_validation.metrics import matched_jaccard
from shared.repo import rel_to_repo

from atlas_postprocessing.artifacts import write_json, write_metric_csv, write_parameters_template
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import (
    build_neighbors,
    load_and_normalize,
    run_harmony_on_pcs,
    run_leiden,
    scale_and_pca,
    select_hvgs,
)
from atlas_postprocessing.metrics import cross_study_macro_cell_type_neighbor_agreement, extract_plateaus
from atlas_postprocessing.plots import plot_sweep_metric

log = logging.getLogger(__name__)


def validate_candidate_lists(cfg: AtlasPostprocessingConfig) -> None:
    for name, values in (
        ("hvgCandidates", cfg.hvgCandidates),
        ("pcCandidates", cfg.pcCandidates),
        ("neighborCandidates", cfg.neighborCandidates),
        ("resolutionCandidates", cfg.resolutionCandidates),
    ):
        if len(values) < 2:
            raise ValueError(f"{name} must contain at least two values, got {values!r}")


def _require_cell_type(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> None:
    if cfg.cellTypeKey not in adata.obs:
        raise ValueError(f"adata.obs is missing cell type key {cfg.cellTypeKey!r}")


def _log_candidate(sweep: str, index: int, total: int, value: object, n_obs: int, n_vars: int) -> float:
    log.info(
        "START %s candidate %s/%s value=%s cells=%s genes=%s",
        sweep,
        index,
        total,
        value,
        f"{n_obs:,}",
        f"{n_vars:,}",
    )
    return time.perf_counter()


def _log_candidate_done(sweep: str, index: int, total: int, value: object, started: float, metric: float) -> None:
    log.info(
        "DONE %s candidate %s/%s value=%s metric=%.4f in %.1fs",
        sweep,
        index,
        total,
        value,
        metric,
        time.perf_counter() - started,
    )


def _agreement_on_harmony_graph(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    nPcs: int,
    nNeighbors: int,
) -> tuple[float, float]:
    if adata.n_obs <= nNeighbors:
        raise ValueError(f"nNeighbors ({nNeighbors}) must be < n_obs ({adata.n_obs})")
    adata.obsm["X_pca_harmony"] = run_harmony_on_pcs(adata, cfg, nPcs=nPcs)
    build_neighbors(adata, nNeighbors=nNeighbors, nPcs=nPcs, useRep="X_pca_harmony")
    return cross_study_macro_cell_type_neighbor_agreement(
        adata,
        batchKey=cfg.batchKey,
        cellTypeKey=cfg.cellTypeKey,
    )


def sweep_hvgs(adata_norm: sc.AnnData, cfg: AtlasPostprocessingConfig) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    total = len(cfg.hvgCandidates)
    for index, n_top in enumerate(cfg.hvgCandidates, start=1):
        started = _log_candidate("hvg", index, total, n_top, adata_norm.n_obs, adata_norm.n_vars)
        adata = adata_norm.copy()
        adata = select_hvgs(adata, cfg, nTopGenes=n_top)
        adata = scale_and_pca(adata, cfg, nPcsCompute=max(cfg.nPcsCompute, cfg.nPcs))
        score, coverage = _agreement_on_harmony_graph(
            adata,
            cfg,
            nPcs=cfg.nPcs,
            nNeighbors=cfg.nNeighbors,
        )
        _log_candidate_done("hvg", index, total, n_top, started, score)
        rows.append(
            {
                "nTopGenes": n_top,
                "weakPriorAgreement": score,
                "eligibleCoverage": coverage,
                "nVars": int(adata.n_vars),
            }
        )

    values = [float(r["nTopGenes"]) for r in rows]
    scores = [float(r["weakPriorAgreement"]) for r in rows]
    plateaus = extract_plateaus(values, scores, relativeThreshold=cfg.plateauRelativeThreshold)
    return {"rows": rows, "values": values, "scores": scores, "plateaus": plateaus}


def sweep_pcs(adata_norm: sc.AnnData, cfg: AtlasPostprocessingConfig) -> dict[str, Any]:
    max_pcs = max(cfg.pcCandidates)
    n_comps = max(cfg.nPcsCompute, max_pcs, cfg.nPcs)
    adata = adata_norm.copy()
    adata = select_hvgs(adata, cfg, nTopGenes=cfg.nTopGenes)
    adata = scale_and_pca(adata, cfg, nPcsCompute=n_comps)

    rows: list[dict[str, Any]] = []
    total = len(cfg.pcCandidates)
    for index, n_pcs in enumerate(cfg.pcCandidates, start=1):
        started = _log_candidate("pc", index, total, n_pcs, adata.n_obs, adata.n_vars)
        score, coverage = _agreement_on_harmony_graph(
            adata,
            cfg,
            nPcs=n_pcs,
            nNeighbors=cfg.nNeighbors,
        )
        _log_candidate_done("pc", index, total, n_pcs, started, score)
        rows.append(
            {
                "nPcs": n_pcs,
                "weakPriorAgreement": score,
                "eligibleCoverage": coverage,
            }
        )

    values = [float(r["nPcs"]) for r in rows]
    scores = [float(r["weakPriorAgreement"]) for r in rows]
    plateaus = extract_plateaus(values, scores, relativeThreshold=cfg.plateauRelativeThreshold)
    return {"rows": rows, "values": values, "scores": scores, "plateaus": plateaus}


def _baseline_harmony_adata(adata_norm: sc.AnnData, cfg: AtlasPostprocessingConfig) -> sc.AnnData:
    adata = adata_norm.copy()
    adata = select_hvgs(adata, cfg, nTopGenes=cfg.nTopGenes)
    adata = scale_and_pca(adata, cfg, nPcsCompute=max(cfg.nPcsCompute, cfg.nPcs))
    adata.obsm["X_pca_harmony"] = run_harmony_on_pcs(adata, cfg, nPcs=cfg.nPcs)
    return adata


def sweep_neighbors(adata_harmony: sc.AnnData, cfg: AtlasPostprocessingConfig) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    total = len(cfg.neighborCandidates)
    for index, n_neighbors in enumerate(cfg.neighborCandidates, start=1):
        started = _log_candidate("neighbors", index, total, n_neighbors, adata_harmony.n_obs, adata_harmony.n_vars)
        if n_neighbors >= adata_harmony.n_obs:
            raise ValueError(f"neighbor candidate {n_neighbors} must be < n_obs ({adata_harmony.n_obs})")
        build_neighbors(
            adata_harmony,
            nNeighbors=n_neighbors,
            nPcs=cfg.nPcs,
            useRep="X_pca_harmony",
        )
        score, coverage = cross_study_macro_cell_type_neighbor_agreement(
            adata_harmony,
            batchKey=cfg.batchKey,
            cellTypeKey=cfg.cellTypeKey,
        )
        _log_candidate_done("neighbors", index, total, n_neighbors, started, score)
        rows.append(
            {
                "nNeighbors": n_neighbors,
                "weakPriorAgreement": score,
                "eligibleCoverage": coverage,
            }
        )

    values = [float(r["nNeighbors"]) for r in rows]
    scores = [float(r["weakPriorAgreement"]) for r in rows]
    plateaus = extract_plateaus(values, scores, relativeThreshold=cfg.plateauRelativeThreshold)
    return {"rows": rows, "values": values, "scores": scores, "plateaus": plateaus}


def sweep_resolutions(adata_harmony: sc.AnnData, cfg: AtlasPostprocessingConfig) -> dict[str, Any]:
    build_neighbors(
        adata_harmony,
        nNeighbors=cfg.nNeighbors,
        nPcs=cfg.nPcs,
        useRep="X_pca_harmony",
    )
    ref_labels = adata_harmony.obs[cfg.cellTypeKey].values
    rows: list[dict[str, Any]] = []
    total = len(cfg.resolutionCandidates)
    for index, resolution in enumerate(cfg.resolutionCandidates, start=1):
        started = _log_candidate("resolution", index, total, resolution, adata_harmony.n_obs, adata_harmony.n_vars)
        key = f"leiden_{resolution}"
        run_leiden(adata_harmony, resolution=resolution, keyAdded=key)
        score = matched_jaccard(adata_harmony.obs[key].values, ref_labels)
        n_clusters = int(adata_harmony.obs[key].nunique())
        _log_candidate_done("resolution", index, total, resolution, started, score)
        rows.append(
            {
                "resolution": resolution,
                "matchedJaccard": score,
                "nClusters": n_clusters,
            }
        )

    values = [float(r["resolution"]) for r in rows]
    scores = [float(r["matchedJaccard"]) for r in rows]
    plateaus = extract_plateaus(values, scores, relativeThreshold=cfg.plateauRelativeThreshold)
    return {"rows": rows, "values": values, "scores": scores, "plateaus": plateaus}


def run_calibration(cfg: AtlasPostprocessingConfig) -> dict[str, Any]:
    """Sweep HVG/PC/neighbor/resolution on a representative subset and write diagnostics."""
    validate_candidate_lists(cfg)
    cfg.calibrationDir.mkdir(parents=True, exist_ok=True)
    metrics_dir = cfg.calibrationDir / "metrics"
    figures_dir = cfg.calibrationDir / "figures"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    log.info(
        "Weak-prior labels (%s) are not ground truth; optimizing them can suppress novel subdivisions",
        cfg.cellTypeKey,
    )

    started_all = time.perf_counter()
    adata_norm = load_and_normalize(cfg)
    _require_cell_type(adata_norm, cfg)
    for n_neighbors in cfg.neighborCandidates:
        if n_neighbors >= adata_norm.n_obs:
            raise ValueError(f"neighbor candidate {n_neighbors} must be < n_obs ({adata_norm.n_obs})")

    log.info("START hvg sweep (%s candidates)", len(cfg.hvgCandidates))
    hvg = sweep_hvgs(adata_norm, cfg)
    log.info("DONE hvg sweep")

    log.info("START pc sweep (%s candidates)", len(cfg.pcCandidates))
    pc = sweep_pcs(adata_norm, cfg)
    log.info("DONE pc sweep")

    log.info("START baseline Harmony embedding for neighbor/resolution sweeps")
    adata_harmony = _baseline_harmony_adata(adata_norm, cfg)
    log.info("DONE baseline Harmony embedding")

    log.info("START neighbors sweep (%s candidates)", len(cfg.neighborCandidates))
    neighbors = sweep_neighbors(adata_harmony, cfg)
    log.info("DONE neighbors sweep")

    log.info("START resolution sweep (%s candidates)", len(cfg.resolutionCandidates))
    resolution = sweep_resolutions(adata_harmony, cfg)
    log.info("DONE resolution sweep")

    write_metric_csv(
        metrics_dir / "hvg.csv",
        hvg["rows"],
        ["nTopGenes", "weakPriorAgreement", "eligibleCoverage", "nVars"],
    )
    write_metric_csv(
        metrics_dir / "pc.csv",
        pc["rows"],
        ["nPcs", "weakPriorAgreement", "eligibleCoverage"],
    )
    write_metric_csv(
        metrics_dir / "neighbors.csv",
        neighbors["rows"],
        ["nNeighbors", "weakPriorAgreement", "eligibleCoverage"],
    )
    write_metric_csv(
        metrics_dir / "resolution.csv",
        resolution["rows"],
        ["resolution", "matchedJaccard", "nClusters"],
    )

    plot_sweep_metric(
        values=hvg["values"],
        scores=hvg["scores"],
        plateaus=hvg["plateaus"],
        xlabel="nTopGenes",
        ylabel="cross-study macro cell-type neighbor agreement",
        title="HVG sweep (weak prior)",
        outPath=figures_dir / "hvg_weak_prior_agreement.png",
        baseline=float(cfg.nTopGenes),
    )
    plot_sweep_metric(
        values=pc["values"],
        scores=pc["scores"],
        plateaus=pc["plateaus"],
        xlabel="nPcs",
        ylabel="cross-study macro cell-type neighbor agreement",
        title="PC sweep (weak prior)",
        outPath=figures_dir / "pc_weak_prior_agreement.png",
        baseline=float(cfg.nPcs),
    )
    plot_sweep_metric(
        values=neighbors["values"],
        scores=neighbors["scores"],
        plateaus=neighbors["plateaus"],
        xlabel="nNeighbors",
        ylabel="cross-study macro cell-type neighbor agreement",
        title="Neighbor sweep (weak prior)",
        outPath=figures_dir / "neighbors_weak_prior_agreement.png",
        baseline=float(cfg.nNeighbors),
    )
    plot_sweep_metric(
        values=resolution["values"],
        scores=resolution["scores"],
        plateaus=resolution["plateaus"],
        xlabel="resolution",
        ylabel="matched Jaccard",
        title="Resolution sweep (weak prior)",
        outPath=figures_dir / "resolution_matched_jaccard.png",
        baseline=float(cfg.resolution),
    )

    summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "calibrationDir": rel_to_repo(cfg.calibrationDir),
        "cells": int(adata_norm.n_obs),
        "genes": int(adata_norm.n_vars),
        "batchKey": cfg.batchKey,
        "cellTypeKey": cfg.cellTypeKey,
        "baseline": {
            "nTopGenes": cfg.nTopGenes,
            "nPcs": cfg.nPcs,
            "nNeighbors": cfg.nNeighbors,
            "resolution": cfg.resolution,
            "nPcsCompute": cfg.nPcsCompute,
        },
        "candidates": {
            "hvg": cfg.hvgCandidates,
            "pc": cfg.pcCandidates,
            "neighbors": cfg.neighborCandidates,
            "resolution": cfg.resolutionCandidates,
        },
        "plateaus": {
            "hvg": hvg["plateaus"],
            "pc": pc["plateaus"],
            "neighbors": neighbors["plateaus"],
            "resolution": resolution["plateaus"],
        },
        "plateauRelativeThreshold": cfg.plateauRelativeThreshold,
        "metrics": {
            "hvg": rel_to_repo(metrics_dir / "hvg.csv"),
            "pc": rel_to_repo(metrics_dir / "pc.csv"),
            "neighbors": rel_to_repo(metrics_dir / "neighbors.csv"),
            "resolution": rel_to_repo(metrics_dir / "resolution.csv"),
        },
        "figures": {
            "hvg": rel_to_repo(figures_dir / "hvg_weak_prior_agreement.png"),
            "pc": rel_to_repo(figures_dir / "pc_weak_prior_agreement.png"),
            "neighbors": rel_to_repo(figures_dir / "neighbors_weak_prior_agreement.png"),
            "resolution": rel_to_repo(figures_dir / "resolution_matched_jaccard.png"),
        },
        "timingsSeconds": round(time.perf_counter() - started_all, 3),
        "note": (
            "Cell-type labels are weak priors, not ground truth. "
            "Inspect plateau intervals and copy parameters_template.json to approved_parameters.json "
            "only after review."
        ),
    }
    summary_path = cfg.calibrationDir / "calibration_summary.json"
    write_json(summary_path, summary)
    write_parameters_template(cfg, summary_path)
    log.info("Calibration complete. Outputs under %s", rel_to_repo(cfg.calibrationDir))
    return summary
