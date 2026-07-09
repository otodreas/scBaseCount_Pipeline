from __future__ import annotations

import logging
import shutil
from pathlib import Path

import h5py
from anndata.experimental import concat_on_disk
from shared.files import safe_delete

from h5ad_concat.config import H5adConcatConfig


def read_h5ad_shape(path: Path) -> tuple[int, int]:
    """Return (n_obs, n_vars) from an h5ad file without loading expression data."""
    with h5py.File(path, "r") as f:
        obs_group = f["obs"]
        var_group = f["var"]
        if not isinstance(obs_group, h5py.Group) or not isinstance(var_group, h5py.Group):
            msg = f"{path} is missing obs or var groups"
            raise ValueError(msg)
        obs_index = obs_group["_index"]
        var_index = var_group["_index"]
        if not isinstance(obs_index, h5py.Dataset) or not isinstance(var_index, h5py.Dataset):
            msg = f"{path} is missing obs/_index or var/_index datasets"
            raise ValueError(msg)
        return int(obs_index.shape[0]), int(var_index.shape[0])


# TODO(stream-pipeline): support incremental folding so pipeline can merge each prepared batch as it
# arrives instead of requiring the full prepared list up front. Likely fold_batch_into_accumulator()
# plus a finalize step; keep concat_prepared as a thin wrapper or replace it.


def concat_prepared(
    prepared: list[tuple[Path, str]],
    cfg: H5adConcatConfig,
    log: logging.Logger,
) -> Path:
    """
    Concatenate a list of prepared h5ad files into the atlas output file, memory-efficiently, writing to `cfg.outputPath`.
    Files are processed in memory-safe batches to limit peak memory and disk usage.

    The function works as follows:
    - Files are concatenated in batches, where the batch size is configured by `cfg.mergeBatchSize`.
    - For each batch, the input `.h5ad` files are concatenated together on disk into an intermediate "batch" file,
      using `concat_on_disk` which avoids loading all data into memory.
    - Once a batch is concatenated, its input files are deleted from disk to free up space.
    - If this is not the first batch, the output of the current batch is merged with the accumulator file
      (result of all previous batches) into a new on-disk "fold" file, and both the current batch and old accumulator are deleted.
      This "batch fold" loop continues until all input files have been consumed, effectively reducing the entire dataset
      to a single merged file in a memory- and disk-safe fashion.
    - The resulting accumulator (final merged h5ad file) is then moved to the configured `cfg.outputPath` location,
      first deleting any previous file at that location.

    All intermediate partial and fold files are stored in a dedicated partials directory, and removed as soon as possible.
    The input files and intermediates are deleted after they are no longer needed to minimize disk usage.
    """
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
    # Atlas is overwritten by default given it's likely to be quite large
    # TODO(output): support appending to existing atlas and/or building new atlas version
    if cfg.outputPath.exists():
        cfg.outputPath.unlink()
    shutil.move(str(accumulator), str(cfg.outputPath))
    return cfg.outputPath
