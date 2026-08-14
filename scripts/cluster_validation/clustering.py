import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def run_leiden(adata: sc.AnnData, *, resolution: float, keyAdded: str) -> None:
    """Run Leiden with the canonical cluster-validation settings."""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added=keyAdded,
    )


def sweep_leiden_resolutions(
    adata: sc.AnnData,
    resolutions: list[float],
) -> tuple[sc.AnnData, dict[float, int]]:
    """Sweep Leiden over ``resolutions`` on an existing neighbor graph."""
    for r in resolutions:
        run_leiden(adata, resolution=r, keyAdded=f"leiden_{r}")
    n_clusters = {r: int(adata.obs[f"leiden_{r}"].nunique()) for r in resolutions}
    return adata, n_clusters


def sweep_leiden(adata: sc.AnnData, cfg: ClusterValidationConfig) -> tuple[sc.AnnData, dict[float, int]]:
    return sweep_leiden_resolutions(adata, cfg.resolutions)
