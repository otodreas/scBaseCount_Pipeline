from __future__ import annotations

import datetime

import scanpy as sc
from shared.logger import configure_file_logger

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.data import load_dataset
from cluster_validation.embedding import embed_dataset
from cluster_validation.merge import merge_clusters
from cluster_validation.metrics import compute_metrics
from cluster_validation.models import ClusterValidationResult
from cluster_validation.preprocess import preprocess
from cluster_validation.resolution import select_resolution_on_graph
from cluster_validation.viz import plot_all

_log = configure_file_logger("cluster_validation.log", __name__)


def _run_tag(srx: str, cfg: ClusterValidationConfig) -> str:
    return cfg.runLabel if cfg.runLabel else srx


def run_cluster_validation_on_adata(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    srx: str,
    title_suffix: str | None = None,
    *,
    write_outputs: bool = True,
    plot: bool = True,
) -> tuple[sc.AnnData, ClusterValidationResult]:
    """Run cluster validation on an in-memory AnnData object."""
    adata = adata.copy()
    adata.obs_names_make_unique()
    run_tag = _run_tag(srx, cfg)
    suffix = title_suffix if title_suffix is not None else f"{srx} ({cfg.weakPriorKey})"

    _log.info("New cluster validation run started")
    print(f"[{datetime.datetime.now().replace(microsecond=0)}] Starting cluster validation for {run_tag}")
    _log.info("starting  %s  weak_prior=%s", run_tag, cfg.weakPriorKey)

    adata, prep_stats = preprocess(adata, cfg)
    adata, n_pcs, cumvar = embed_dataset(adata, cfg)
    adata, sel = select_resolution_on_graph(
        adata,
        resolutions=cfg.resolutions,
        weakPriorKey=cfg.weakPriorKey,
    )
    adata, merge_info = merge_clusters(adata, cfg, sel)
    metric_arrays = compute_metrics(adata, cfg, sel, merge_info)

    adata_path = cfg.outputDir / f"{run_tag}_clustered.h5ad"
    if write_outputs:
        cfg.outputDir.mkdir(parents=True, exist_ok=True)
        adata.write(str(adata_path))

    result = ClusterValidationResult(
        srxAccession=srx,
        runTag=run_tag,
        weakPriorKey=cfg.weakPriorKey,
        datasetTitleSuffix=suffix,
        selectedResolution=sel.selectedResolution,
        clusterKey=sel.clusterKey,
        mergedKey=merge_info.mergedKey,
        nPcs=n_pcs,
        cumvar=cumvar,
        kPrior=prep_stats.kPrior,
        kFiltered=prep_stats.kFiltered,
        nCellsDropped=prep_stats.nDropped,
        nCellsFinal=prep_stats.nCellsFinal,
        nClustersPreMerge=merge_info.nClustersPreMerge,
        nClustersPostMerge=merge_info.nClustersPostMerge,
        adataPath=adata_path,
        labelMap=merge_info.labelMap,
        mergedGroups=merge_info.mergedGroups,
        resolutions=cfg.resolutions,
        kArr=sel.kArr.tolist(),
        jaccArr=sel.jaccArr.tolist(),
        silhouetteArr=metric_arrays.silhouetteArr,
        homogeneityArr=metric_arrays.homogeneityArr,
        completenessArr=metric_arrays.completenessArr,
        nmiArr=metric_arrays.nmiArr,
        vscoreArr=metric_arrays.vscoreArr,
        ariArr=metric_arrays.ariArr,
        confMatrix=merge_info.conf.tolist(),
        confClasses=[str(c) for c in merge_info.classes],
    )

    if plot:
        plot_all(adata, result, figs_dir=cfg.figsDir / run_tag)

    print(
        f"[{datetime.datetime.now().replace(microsecond=0)}] Done cluster validation for {run_tag} "
        f"with resolution {result.selectedResolution} and {result.nClustersPostMerge} clusters and {result.nPcs} PCs"
    )
    _log.info(
        "done      %s  resolution=%.1f  clusters=%d  pcs=%d",
        run_tag,
        result.selectedResolution,
        result.nClustersPostMerge,
        result.nPcs,
    )
    return adata, result


def run_cluster_validation(
    cfg: ClusterValidationConfig,
) -> tuple[sc.AnnData, ClusterValidationResult]:
    adata, srx, title_suffix = load_dataset(cfg)
    return run_cluster_validation_on_adata(adata, cfg, srx, title_suffix)
