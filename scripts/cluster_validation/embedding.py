import numpy as np
import scanpy as sc
from numpy.typing import NDArray

from cluster_validation.config import ClusterValidationConfig


def pick_n_pcs(
    varianceRatio: NDArray[np.floating] | list[float],
    *,
    nPcsMin: int,
    nPcsCompute: int,
    nPcsCumvarTarget: float,
) -> tuple[int, float]:
    """Choose how many PCs to keep from explained-variance ratios.

    Returns ``(n_pcs, cumvar_percent)`` where ``cumvar_percent`` is cumulative
    explained variance of the selected PCs expressed as a percentage.
    """
    var_ratio = np.asarray(varianceRatio, dtype=np.float64)
    if nPcsMin < 1:
        raise ValueError(f"nPcsMin must be >= 1, got {nPcsMin}")
    if nPcsCompute < nPcsMin:
        raise ValueError(f"nPcsCompute ({nPcsCompute}) must be >= nPcsMin ({nPcsMin})")
    if len(var_ratio) < nPcsCompute:
        raise ValueError(f"varianceRatio length ({len(var_ratio)}) must be >= nPcsCompute ({nPcsCompute})")

    min_cumvar = float(np.sum(var_ratio[:nPcsMin]))
    if min_cumvar >= nPcsCumvarTarget:
        return nPcsMin, min_cumvar * 100.0

    for i, cumvar_tail in enumerate(np.cumsum(var_ratio[nPcsMin:]), start=nPcsMin):
        cumvar = float(cumvar_tail) + min_cumvar
        if cumvar >= nPcsCumvarTarget or i == nPcsCompute - 1:
            return i + 1, cumvar * 100.0

    return nPcsCompute, float(np.sum(var_ratio[:nPcsCompute])) * 100.0


def embed_dataset(adata: sc.AnnData, cfg: ClusterValidationConfig) -> tuple[sc.AnnData, int, float]:
    sc.tl.pca(adata, n_comps=cfg.nPcsCompute, svd_solver="arpack")
    n_pcs, cumvar = pick_n_pcs(
        adata.uns["pca"]["variance_ratio"],
        nPcsMin=cfg.nPcsMin,
        nPcsCompute=cfg.nPcsCompute,
        nPcsCumvarTarget=cfg.nPcsCumvarTarget,
    )
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata, n_pcs, cumvar
