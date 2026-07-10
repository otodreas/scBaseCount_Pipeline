import logging
import os
from pathlib import Path

import anndata as ad

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.models import ConcatManifest

MANIFEST_KEY = "h5adConcatManifest"


def write_checkpoint(
    atlas: ad.AnnData,
    manifest: ConcatManifest,
    cfg: H5adConcatConfig,
    log: logging.Logger,
) -> Path:
    """Embed manifest in atlas.uns and atomically write to cfg.outputPath."""
    cfg.outputPath.parent.mkdir(parents=True, exist_ok=True)
    atlas.uns[MANIFEST_KEY] = manifest.model_dump_json()
    tmp_path = cfg.outputPath.with_suffix(".tmp.h5ad")
    atlas.write_h5ad(tmp_path, compression=cfg.compression)
    os.replace(tmp_path, cfg.outputPath)
    log.info(
        "Wrote checkpoint to %s (%d entries, %d obs)",
        cfg.outputPath,
        len(manifest.entries),
        atlas.n_obs,
    )
    return cfg.outputPath


def load_checkpoint(
    cfg: H5adConcatConfig,
    log: logging.Logger,
) -> tuple[ad.AnnData | None, ConcatManifest]:
    """Load atlas and manifest from cfg.outputPath when resume is enabled."""
    empty_manifest = ConcatManifest(join=cfg.join, batchKey=cfg.batchKey)
    if not cfg.resume or not cfg.outputPath.exists():
        return None, empty_manifest

    atlas = ad.read_h5ad(cfg.outputPath)
    raw_manifest = atlas.uns.get(MANIFEST_KEY)
    if raw_manifest is None:
        msg = f"Cannot resume: {cfg.outputPath} has no embedded manifest"
        raise ValueError(msg)

    manifest = ConcatManifest.model_validate_json(raw_manifest)
    if manifest.join != cfg.join or manifest.batchKey != cfg.batchKey:
        msg = (
            f"Cannot resume: manifest join/batchKey ({manifest.join}/{manifest.batchKey}) "
            f"does not match config ({cfg.join}/{cfg.batchKey})"
        )
        raise ValueError(msg)

    log.info(
        "Resumed from checkpoint: %d obs, %d processed keys",
        atlas.n_obs,
        len(manifest.processedKeys()),
    )
    return atlas, manifest
