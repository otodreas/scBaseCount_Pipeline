"""Memory preflight and process telemetry for full-atlas aggregation."""

from __future__ import annotations

import logging
import platform
import resource
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

log = logging.getLogger(__name__)

_BYTES_PER_GIB = 1024**3


@dataclass(frozen=True)
class SparseMatrixEstimate:
    nObs: int
    nVars: int
    nnz: int
    estimatedBytes: int
    sourcePath: str


@dataclass(frozen=True)
class MemorySnapshot:
    availableBytes: int | None
    rssBytes: int
    peakRssBytes: int


def bytes_to_gib(nBytes: int | float | None) -> float | None:
    if nBytes is None:
        return None
    return float(nBytes) / _BYTES_PER_GIB


def format_bytes(nBytes: int | float | None) -> str:
    if nBytes is None:
        return "unknown"
    return f"{bytes_to_gib(nBytes):.2f} GiB"


def process_rss_bytes() -> int:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports KiB; macOS reports bytes.
    if platform.system() == "Darwin":
        return int(usage)
    return int(usage * 1024)


def peak_rss_bytes() -> int:
    return process_rss_bytes()


def available_ram_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return None
    values: dict[str, int] = {}
    for line in meminfo.read_text().splitlines():
        if ":" not in line:
            continue
        key, raw = line.split(":", 1)
        parts = raw.strip().split()
        if not parts:
            continue
        try:
            kib = int(parts[0])
        except ValueError:
            continue
        values[key] = kib * 1024
    if "MemAvailable" in values:
        return values["MemAvailable"]
    if "MemFree" in values and "Buffers" in values and "Cached" in values:
        return values["MemFree"] + values["Buffers"] + values["Cached"]
    return None


def snapshot_memory() -> MemorySnapshot:
    return MemorySnapshot(
        availableBytes=available_ram_bytes(),
        rssBytes=process_rss_bytes(),
        peakRssBytes=peak_rss_bytes(),
    )


def log_memory(stage: str, *, logger: logging.Logger | None = None) -> MemorySnapshot:
    active = logger or log
    snap = snapshot_memory()
    active.info(
        "%s memory: available=%s rss=%s peakRss=%s",
        stage,
        format_bytes(snap.availableBytes),
        format_bytes(snap.rssBytes),
        format_bytes(snap.peakRssBytes),
    )
    return snap


def _read_csr_shape_nnz(group: h5py.Group) -> tuple[int, int, int]:
    encoding = group.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    shape_attr = group.attrs.get("shape", ())
    shape = tuple(int(x) for x in list(shape_attr))
    if len(shape) != 2:
        raise ValueError(f"Expected 2D sparse matrix shape, got {shape}")
    if "indptr" not in group:
        raise ValueError("Sparse matrix group is missing indptr")
    indptr = group["indptr"]
    nnz = int(indptr[-1]) if len(indptr.shape) and indptr.shape[0] else 0
    if "data" in group:
        nnz = max(nnz, int(group["data"].shape[0]))
    return shape[0], shape[1], nnz


def estimate_raw_sparse_bytes(atlasPath: Path, *, overheadFactor: float = 2.0) -> SparseMatrixEstimate:
    """Estimate uncompressed sparse raw-count memory from H5AD metadata.

    Uses CSR layout cost: nnz * (value + index) + (nObs + 1) * indptr, then a
    configurable overhead factor for AnnData/Python object costs.
    """
    path = Path(atlasPath)
    with h5py.File(path, "r") as handle:
        if "raw" in handle and isinstance(handle["raw"], h5py.Group) and "X" in handle["raw"]:
            group = handle["raw"]["X"]
            source = "raw/X"
        elif "X" in handle:
            group = handle["X"]
            source = "X"
        else:
            raise KeyError(f"{path} has neither raw/X nor X")
        if not isinstance(group, h5py.Group):
            # Dense array encoding.
            dataset = group
            shape = tuple(int(x) for x in dataset.shape)
            if len(shape) != 2:
                raise ValueError(f"Dense matrix in {path} has unexpected shape {shape}")
            itemsize = int(np.dtype(dataset.dtype).itemsize)
            estimated = int(shape[0] * shape[1] * itemsize * overheadFactor)
            return SparseMatrixEstimate(
                nObs=int(shape[0]),
                nVars=int(shape[1]),
                nnz=int(shape[0] * shape[1]),
                estimatedBytes=estimated,
                sourcePath=source,
            )

        encoding = group.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if "array" in str(encoding):
            shape = tuple(int(x) for x in list(group.attrs.get("shape", ())))
            if len(shape) != 2:
                raise ValueError(f"Dense matrix in {path} has unexpected shape {shape}")
            estimated = int(shape[0] * shape[1] * 8 * overheadFactor)
            return SparseMatrixEstimate(
                nObs=int(shape[0]),
                nVars=int(shape[1]),
                nnz=int(shape[0] * shape[1]),
                estimatedBytes=estimated,
                sourcePath=source,
            )
        nObs, nVars, nnz = _read_csr_shape_nnz(group)
        # float64 values + int32 indices + int64 indptr is a conservative upper bound.
        valueBytes = 8
        indexBytes = 4
        indptrBytes = 8
        rawBytes = nnz * (valueBytes + indexBytes) + (nObs + 1) * indptrBytes
        estimated = int(rawBytes * overheadFactor)
        return SparseMatrixEstimate(
            nObs=nObs,
            nVars=nVars,
            nnz=nnz,
            estimatedBytes=estimated,
            sourcePath=source,
        )


def assert_memory_available(
    estimate: SparseMatrixEstimate,
    *,
    reserveBytes: int,
    logger: logging.Logger | None = None,
) -> None:
    active = logger or log
    available = available_ram_bytes()
    active.info(
        "Memory preflight: source=%s shape=%s x %s nnz=%s estimate=%s available=%s reserve=%s",
        estimate.sourcePath,
        f"{estimate.nObs:,}",
        f"{estimate.nVars:,}",
        f"{estimate.nnz:,}",
        format_bytes(estimate.estimatedBytes),
        format_bytes(available),
        format_bytes(reserveBytes),
    )
    if available is None:
        active.warning("Could not read MemAvailable; continuing without hard preflight gate")
        return
    needed = estimate.estimatedBytes + int(reserveBytes)
    if available < needed:
        raise MemoryError(
            f"Insufficient RAM for atlas load: available={format_bytes(available)}, "
            f"estimated={format_bytes(estimate.estimatedBytes)}, "
            f"reserve={format_bytes(reserveBytes)}. "
            "Lower memoryReserveBytes only if you accept less headroom."
        )
