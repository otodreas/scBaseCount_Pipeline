from __future__ import annotations

import datetime

import scanpy as sc
from shared.logger import configure_file_logger

from cluster_validation.clustering import sweep_leiden
from cluster_validation.config import ClusterValidationConfig
from cluster_validation.data import load_dataset
from cluster_validation.embedding import embed_dataset
from cluster_validation.merge import merge_clusters
from cluster_validation.metrics import compute_metrics
from cluster_validation.models import ClusterValidationResult
from cluster_validation.preprocess import preprocess
from cluster_validation.resolution import select_resolution
from cluster_validation.viz import plot_all

_log = configure_file_logger("cluster_validation.log", __name__)


def run_cluster_validation(
    cfg: ClusterValidationConfig,
) -> tuple[sc.AnnData, ClusterValidationResult]:
    adata, srx, title_suffix = load_dataset(cfg)
    _log.info("New cluster validation run started")
    print(f"[{datetime.datetime.now().replace(microsecond=0)}] Starting cluster validation for {srx}")
    _log.info("starting  %s", srx)
    adata, prep_stats = preprocess(adata, cfg)
    adata, n_pcs, cumvar = embed_dataset(adata, cfg)
    adata, n_clusters = sweep_leiden(adata, cfg)
    adata, sel = select_resolution(adata, cfg, n_clusters, prep_stats.kFiltered)
    adata, merge_info = merge_clusters(adata, cfg, sel)
    metric_arrays = compute_metrics(adata, cfg, sel, merge_info)

    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    adata_path = cfg.outputDir / f"{srx}_clustered.h5ad"
    adata.write(str(adata_path))

    result = ClusterValidationResult(
        srxAccession=srx,
        datasetTitleSuffix=title_suffix,
        embedding=cfg.embedding,
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

    plot_all(adata, result, figs_dir=cfg.figsDir / srx)

    if cfg.embedding == "pca":
        print(
            f"[{datetime.datetime.now().replace(microsecond=0)}] Done cluster validation for {srx} "
            f"with resolution {result.selectedResolution} and {result.nClustersPostMerge} clusters "
            f"and {result.nPcs} PCs"
        )
        _log.info(
            "done      %s  resolution=%.1f  clusters=%d  pcs=%d",
            srx,
            result.selectedResolution,
            result.nClustersPostMerge,
            result.nPcs,
        )
    else:
        print(
            f"[{datetime.datetime.now().replace(microsecond=0)}] Done cluster validation for {srx} "
            f"with resolution {result.selectedResolution} and {result.nClustersPostMerge} clusters "
            f"and {cfg.embedding} embedding"
        )
        _log.info(
            "done      %s  resolution=%.1f  clusters=%d  embedding=%s",
            srx,
            result.selectedResolution,
            result.nClustersPostMerge,
            cfg.embedding,
        )
    return adata, result
