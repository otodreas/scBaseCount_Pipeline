from __future__ import annotations

import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def sweep_leiden(
    adata: sc.AnnData, cfg: ClusterValidationConfig
) -> tuple[sc.AnnData, dict[float, int]]:
    for r in cfg.resolutions:
        sc.tl.leiden(
            adata,
            resolution=r,
            flavor="igraph",
            n_iterations=2,
            directed=False,
            key_added=f"leiden_{r}",
        )
    n_clusters = {r: adata.obs[f"leiden_{r}"].nunique() for r in cfg.resolutions}
    return adata, n_clusters
