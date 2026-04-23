from __future__ import annotations

import numpy as np
import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def embed_dataset(
    adata: sc.AnnData, cfg: ClusterValidationConfig
) -> tuple[sc.AnnData, int, float]:
    sc.tl.pca(adata, n_comps=cfg.nPcsCompute, svd_solver="arpack")
    n_pcs, cumvar = _pick_n_pcs(adata, cfg)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata, n_pcs, cumvar


def _pick_n_pcs(adata: sc.AnnData, cfg: ClusterValidationConfig) -> tuple[int, float]:
    var_ratio = adata.uns["pca"]["variance_ratio"]
    min_cumvar = float(np.sum(var_ratio[: cfg.nPcsMin]))

    if min_cumvar >= cfg.nPcsCumvarTarget:
        return cfg.nPcsMin, min_cumvar * 100.0

    for i, cumvar_tail in enumerate(np.cumsum(var_ratio[cfg.nPcsMin :]), start=cfg.nPcsMin):
        cumvar = float(cumvar_tail) + min_cumvar
        if cumvar >= cfg.nPcsCumvarTarget or i == cfg.nPcsCompute - 1:
            return i + 1, cumvar * 100.0

    return cfg.nPcsCompute, float(np.sum(var_ratio)) * 100.0
