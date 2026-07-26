from collections.abc import Iterator

import decoupler as dc
import scanpy as sc

from disease_markers.config import DiseaseMarkersConfig


def pseudobulk_for_cluster(
    adata: sc.AnnData,
    cluster: str,
    cfg: DiseaseMarkersConfig,
) -> sc.AnnData | None:
    """Aggregate raw counts to sample-level pseudobulk for one Leiden cluster."""
    cluster_key = cfg.clusterKey
    mask = adata.obs[cluster_key].astype(str) == str(cluster)
    if not bool(mask.any()):
        return None
    sub = adata[mask].copy()
    result = dc.pp.pseudobulk(
        sub,
        sample_col=cfg.sampleKey,
        mode="sum",
        min_cells=cfg.minCellsPerProfile,
    )
    if isinstance(result, tuple):
        pdata = result[0]
    else:
        pdata = result
    if pdata.n_obs == 0:
        return None
    return pdata


def iter_cluster_pseudobulks(
    adata: sc.AnnData,
    cfg: DiseaseMarkersConfig,
) -> Iterator[tuple[str, sc.AnnData]]:
    """Yield (cluster_id, pseudobulk AnnData) for each cluster in adata."""
    cluster_key = cfg.clusterKey
    for cluster in adata.obs[cluster_key].astype(str).cat.categories:
        pdata = pseudobulk_for_cluster(adata, str(cluster), cfg)
        if pdata is not None:
            yield str(cluster), pdata
