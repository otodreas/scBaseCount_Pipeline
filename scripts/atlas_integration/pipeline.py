import json
import logging

import scanpy as sc
from umap_plots import UmapPlotConfig, plot_umap

from atlas_integration.config import AtlasIntegrationConfig
from atlas_integration.integrate import cluster_uncorrected, integrate_atlas
from atlas_integration.merge import build_merged_adata
from atlas_integration.metrics import compute_batch_mixing, compute_cluster_conservation
from atlas_integration.models import AtlasIntegrationResult, MergeStats
from atlas_integration.preprocess import preprocess_atlas

_log = logging.getLogger(__name__)


def run_atlas_integration(
    cfg: AtlasIntegrationConfig,
    adata: sc.AnnData | None = None,
    merge_stats: MergeStats | None = None,
    *,
    write_outputs: bool = True,
    plot: bool = True,
) -> tuple[sc.AnnData, AtlasIntegrationResult]:
    """Merge lung accessions, integrate with Harmony, and write atlas outputs."""
    if adata is None:
        adata, merge_stats = build_merged_adata(cfg)
    elif merge_stats is None:
        raise ValueError("merge_stats must be provided when adata is passed in")

    adata, n_pcs, cumvar = preprocess_atlas(adata, cfg)
    adata = cluster_uncorrected(adata, cfg)
    adata = integrate_atlas(adata, cfg)

    batch_mixing = compute_batch_mixing(adata, cfg)
    cluster_conservation = compute_cluster_conservation(adata, cfg)

    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    atlas_path = cfg.outputDir / "data" / cfg.atlasH5adName
    metadata_path = cfg.outputDir / "run_metadata.json"

    if write_outputs:
        atlas_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write(atlas_path)
        metadata_path.write_text(
            json.dumps(
                {
                    "atlasPath": str(atlas_path),
                    "mergeStats": merge_stats.model_dump(),
                    "batchMixing": batch_mixing.model_dump(),
                    "clusterConservation": cluster_conservation.model_dump(),
                    "nPcs": n_pcs,
                    "cumvarPct": cumvar,
                },
                indent=2,
            )
        )

    if plot:
        _plot_atlas_umaps(adata, cfg)

    result = AtlasIntegrationResult(
        atlasPath=atlas_path,
        metadataPath=metadata_path,
        mergeStats=merge_stats,
        batchMixing=batch_mixing,
        clusterConservation=cluster_conservation,
        nPcs=n_pcs,
        cumvarPct=cumvar,
    )
    _log.info(
        "Atlas integration complete: %d cells, %d studies, atlas=%s",
        adata.n_obs,
        adata.obs[cfg.batchKey].nunique(),
        atlas_path,
    )
    return adata, result


def _plot_atlas_umaps(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> None:
    cfg.figsDir.mkdir(parents=True, exist_ok=True)
    uncorrected_cfg = UmapPlotConfig(figsDir=cfg.figsDir, umapKey=cfg.umapKeyUncorrected)
    corrected_cfg = UmapPlotConfig(figsDir=cfg.figsDir, umapKey=cfg.umapKeyCorrected)

    for color_by in (cfg.batchKey, cfg.cellTypeKey):
        plot_umap(adata, color_by, uncorrected_cfg, nameSuffix="uncorrected")
        plot_umap(adata, color_by, corrected_cfg, nameSuffix="corrected")
