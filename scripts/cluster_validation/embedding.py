from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def active_rep(cfg: ClusterValidationConfig) -> str:
    """Return the obsm key used for neighbors and silhouette scoring."""
    if cfg.embedding == "scgpt":
        return cfg.scgptObsmKey
    return "X_pca"


def embed_dataset(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
) -> tuple[sc.AnnData, int | None, float | None]:
    if cfg.embedding == "scgpt":
        return _embed_scgpt(adata, cfg)
    return _embed_pca(adata, cfg)


def _attach_scgpt_embedding(adata: sc.AnnData, cfg: ClusterValidationConfig) -> sc.AnnData:
    if cfg.scgptEmbedPath is None:
        raise ValueError("scgptEmbedPath is required when embedding='scgpt'")

    embed_path = Path(cfg.scgptEmbedPath)
    if not embed_path.exists():
        raise FileNotFoundError(f"scGPT embedding artifact not found: {embed_path}")

    with np.load(embed_path, allow_pickle=False) as data:
        obs_names = data["obs_names"].astype(str)
        x = data["X"].astype(np.float32)

    embed_df = pd.DataFrame(x, index=obs_names)
    aligned = embed_df.reindex(adata.obs_names.astype(str))
    missing = aligned.index[aligned.isna().all(axis=1)]
    if len(missing) > 0:
        raise ValueError(
            f"{len(missing)} cells in adata have no scGPT embedding (first missing: {missing[:3].tolist()})"
        )

    adata.obsm[cfg.scgptObsmKey] = aligned.to_numpy(dtype=np.float32)
    return adata


def _embed_scgpt(adata: sc.AnnData, cfg: ClusterValidationConfig) -> tuple[sc.AnnData, None, None]:
    adata = _attach_scgpt_embedding(adata, cfg)
    sc.pp.neighbors(adata, use_rep=cfg.scgptObsmKey)
    sc.tl.umap(adata)
    return adata, None, None


def _embed_pca(adata: sc.AnnData, cfg: ClusterValidationConfig) -> tuple[sc.AnnData, int, float]:
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
