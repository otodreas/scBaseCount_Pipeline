import scanpy as sc

from atlas_integration.config import AtlasIntegrationConfig


def integrate_atlas(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> sc.AnnData:
    """Run Harmony batch correction, then recompute neighbors, Leiden clusters, and UMAP."""
    sc.external.pp.harmony_integrate(
        adata,
        key=cfg.batchKey,
        basis=cfg.pcaKey,
        adjusted_basis=cfg.harmonyKey,
        max_iter_harmony=cfg.harmonyMaxIter,
    )
    sc.pp.neighbors(adata, use_rep=cfg.harmonyKey)
    sc.tl.leiden(adata, resolution=cfg.leidenResolution, key_added=cfg.leidenKeyAtlas)
    sc.tl.umap(adata)
    return adata


def cluster_uncorrected(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> sc.AnnData:
    """Cluster the uncorrected PCA embedding for conservation comparisons."""
    sc.pp.neighbors(adata, use_rep=cfg.pcaKey)
    sc.tl.leiden(adata, resolution=cfg.leidenResolution, key_added=cfg.leidenKeyUncorrected)
    return adata
