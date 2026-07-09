from __future__ import annotations

import logging
import shutil
from pathlib import Path

from anndata.experimental import concat_on_disk
from shared.files import safe_delete

from h5ad_concat.config import H5adConcatConfig


def concat_prepared(
    prepared: list[tuple[Path, str]],
    cfg: H5adConcatConfig,
    log: logging.Logger,
) -> Path:
    """Memory-safe concat of prepared h5ads into cfg.outputPath; delete staged inputs."""
    if not prepared:
        msg = "No prepared files to concatenate"
        raise ValueError(msg)

    cfg.outputPath.parent.mkdir(parents=True, exist_ok=True)
    partial_dir = cfg.cacheDir / "partials"
    partial_dir.mkdir(parents=True, exist_ok=True)

    paths = [path for path, _ in prepared]

    accumulator: Path | None = None
    for batch_start in range(0, len(paths), cfg.mergeBatchSize):  # cfg.mergeBatchSize is step size
        batch_paths = paths[batch_start : batch_start + cfg.mergeBatchSize]
        batch_out = partial_dir / f"batch_{batch_start:05d}.h5ad"
        concat_on_disk(batch_paths, batch_out, max_loaded_elems=cfg.maxLoadedElems, join=cfg.join)
        for path in batch_paths:
            safe_delete(path, log)

        if accumulator is None:
            accumulator = batch_out
            continue

        merged = partial_dir / f"fold_{batch_start:05d}.h5ad"
        concat_on_disk([accumulator, batch_out], merged, max_loaded_elems=cfg.maxLoadedElems, join=cfg.join)
        safe_delete(accumulator, log)
        safe_delete(batch_out, log)
        accumulator = merged

    assert accumulator is not None
    if cfg.outputPath.exists():
        cfg.outputPath.unlink()
    shutil.move(str(accumulator), str(cfg.outputPath))
    return cfg.outputPath
