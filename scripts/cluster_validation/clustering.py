from __future__ import annotations

import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def sweep_leiden(
    adata: sc.AnnData, cfg: ClusterValidationConfig, suffix: str = ""
) -> tuple[sc.AnnData, dict[float, int]]:
    """Run Leiden clustering at multiple resolutions and return the number of clusters at each resolution."""
    for r in cfg.resolutions:
        sc.tl.leiden(
            adata,
            resolution=r,
            flavor="igraph",
            n_iterations=2,
            directed=False,
            key_added=f"leiden_{suffix + '_' if suffix else ''}{r}",
        )
    n_clusters = {r: adata.obs[f"leiden_{suffix + '_' if suffix else ''}{r}"].nunique() for r in cfg.resolutions}
    return adata, n_clusters
