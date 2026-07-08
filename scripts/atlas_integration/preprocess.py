import numpy as np
import scanpy as sc

from atlas_integration.config import AtlasIntegrationConfig


def preprocess_atlas(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> tuple[sc.AnnData, int, float]:
    """Normalize, select batch-aware HVGs, scale, run PCA, and compute an uncorrected UMAP."""
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.nTopGenes, batch_key=cfg.batchKey)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=cfg.nPcsCompute, svd_solver="arpack")
    n_pcs, cumvar = _pick_n_pcs(adata, cfg)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    adata.obsm[cfg.umapKeyUncorrected] = adata.obsm["X_umap"].copy()
    return adata, n_pcs, cumvar


def _pick_n_pcs(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> tuple[int, float]:
    var_ratio = adata.uns["pca"]["variance_ratio"]
    min_cumvar = float(np.sum(var_ratio[: cfg.nPcsMin]))

    if min_cumvar >= cfg.nPcsCumvarTarget:
        return cfg.nPcsMin, min_cumvar * 100.0

    for i, cumvar_tail in enumerate(np.cumsum(var_ratio[cfg.nPcsMin :]), start=cfg.nPcsMin):
        cumvar = float(cumvar_tail) + min_cumvar
        if cumvar >= cfg.nPcsCumvarTarget or i == cfg.nPcsCompute - 1:
            return i + 1, cumvar * 100.0

    return cfg.nPcsCompute, float(np.sum(var_ratio)) * 100.0
