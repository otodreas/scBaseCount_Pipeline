from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad

from h5ad_concat.checkpoint import write_checkpoint
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.models import ConcatManifest


def concat_atlas(adatas: list[ad.AnnData], cfg: H5adConcatConfig) -> ad.AnnData:
    """Concatenate prepared AnnData objects on the obs axis into one in-memory atlas."""
    return ad.concat(adatas, axis="obs", join=cfg.join)


def fold_atlas(
    atlas: ad.AnnData | None,
    pending: list[ad.AnnData],
    cfg: H5adConcatConfig,
) -> ad.AnnData:
    """Collapse pending AnnData into atlas, returning a single in-memory object."""
    if atlas is None:
        return concat_atlas(pending, cfg)
    return concat_atlas([atlas, *pending], cfg)


def write_atlas(
    atlas: ad.AnnData,
    manifest: ConcatManifest,
    cfg: H5adConcatConfig,
    log: logging.Logger,
) -> Path:
    """Atomically write atlas with embedded manifest to cfg.outputPath."""
    # TODO(output): support building new atlas versions instead of overwriting outputPath
    return write_checkpoint(atlas, manifest, cfg, log)
