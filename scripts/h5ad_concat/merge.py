from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad

from h5ad_concat.config import H5adConcatConfig


def concat_atlas(adatas: list[ad.AnnData], cfg: H5adConcatConfig, log: logging.Logger) -> ad.AnnData:
    """Concatenate prepared AnnData objects on the obs axis into one in-memory atlas."""
    adata = ad.concat(adatas, axis="obs", join=cfg.join)
    log.info("Concatenated %d AnnData objects into one in-memory atlas", len(adatas))
    return adata


def write_atlas(adata: ad.AnnData, cfg: H5adConcatConfig, log: logging.Logger) -> Path:
    """Write the atlas to cfg.outputPath with gzip, overwriting any existing file."""
    cfg.outputPath.parent.mkdir(parents=True, exist_ok=True)
    if cfg.outputPath.exists():
        cfg.outputPath.unlink()
    adata.write_h5ad(cfg.outputPath, compression=cfg.compression)
    log.info("Wrote atlas to %s", cfg.outputPath)
    return cfg.outputPath
