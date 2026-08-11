"""Copy-conscious in-memory pseudobulk aggregation for atlas disease DE."""

from __future__ import annotations

import gc
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from metadata.config import MetadataConfig
from scipy import sparse

from disease_markers.candidates import annotate_obs_with_labels, same_study_contrast_support
from disease_markers.concordance import as_string_series
from disease_markers.config import AtlasDeAnalysisConfig
from disease_markers.labels import build_sample_label_table
from disease_markers.memory import (
    assert_memory_available,
    estimate_raw_sparse_bytes,
    format_bytes,
    log_memory,
    snapshot_memory,
)
from disease_markers.specificity import cluster_interpretation_table
from disease_markers.validation import filter_pseudobulk_profiles

log = logging.getLogger(__name__)


def _file_fingerprint(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "sizeBytes": int(stat.st_size),
        "mtimeNs": int(stat.st_mtime_ns),
    }


def _hash_payload(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_aggregate_fingerprint(
    cfg: AtlasDeAnalysisConfig,
    *,
    sampleMetadataPath: Path | None = None,
) -> dict[str, Any]:
    samplePath = sampleMetadataPath or MetadataConfig().sampleParquetPath
    payload = {
        "atlas": _file_fingerprint(cfg.atlasPath),
        "contexts": _file_fingerprint(cfg.contextsPath),
        "atlasCsv": _file_fingerprint(cfg.atlasCsvPath),
        "sampleMetadata": _file_fingerprint(Path(samplePath)),
        "sampleKey": cfg.sampleKey,
        "studyKey": cfg.studyKey,
        "clusterKey": cfg.clusterKey,
        "minCellsPerProfile": cfg.minCellsPerProfile,
        "geneInfo": _file_fingerprint(cfg.geneInfoPath),
    }
    return {
        "algorithm": "sparseSampleClusterSumV1",
        "sha256": _hash_payload(payload),
        "payload": payload,
    }


def reject_dense_cell_matrix(matrix: Any, *, label: str) -> None:
    if sparse.issparse(matrix):
        return
    if hasattr(matrix, "ndim") and int(getattr(matrix, "ndim", 0)) == 2:
        raise TypeError(
            f"Rejected dense cell-level matrix for {label}. Full-atlas aggregation requires a sparse raw count matrix."
        )


def sparse_raw_matrix(adata: AnnData) -> sparse.csr_matrix:
    if adata.raw is None:
        raise ValueError("Expected full-gene counts in adata.raw")
    matrix = adata.raw.X
    reject_dense_cell_matrix(matrix, label="adata.raw.X")
    return sparse.csr_matrix(matrix)


def _obs_for_groups(
    obs: pd.DataFrame,
    *,
    groupCodes: np.ndarray,
    nGroups: int,
) -> pd.DataFrame:
    """Take the first obs row for each sample x cluster group."""
    frame = obs.copy()
    frame["_groupCode"] = groupCodes
    out = frame.groupby("_groupCode", sort=True).first()
    out = out.reindex(range(nGroups))
    return out.drop(columns=["_groupCode"], errors="ignore")


def aggregate_sparse_pseudobulk(
    matrix: sparse.spmatrix | sparse.csr_matrix,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    *,
    sampleKey: str,
    clusterKey: str,
) -> AnnData:
    """Sum sparse counts by sample x cluster without densifying cell-level data."""
    reject_dense_cell_matrix(matrix, label="aggregate input")
    if matrix.shape[0] != len(obs):
        raise ValueError(f"Matrix rows ({matrix.shape[0]}) do not match obs ({len(obs)})")
    if matrix.shape[1] != len(var):
        raise ValueError(f"Matrix cols ({matrix.shape[1]}) do not match var ({len(var)})")

    sample = as_string_series(pd.Series(obs[sampleKey]))
    cluster = as_string_series(pd.Series(obs[clusterKey]))
    grouped = pd.MultiIndex.from_arrays([sample.to_numpy(), cluster.to_numpy()], names=[sampleKey, clusterKey])
    codes, uniques = pd.factorize(grouped, sort=True)
    n_groups = int(len(uniques))
    n_cells = int(matrix.shape[0])

    csr = sparse.csr_matrix(matrix)
    indicator = sparse.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (codes, np.arange(n_cells))),
        shape=(n_groups, n_cells),
    )
    counts = indicator @ csr
    detected = indicator @ (csr > 0)
    detected_array = np.asarray(detected.toarray() if sparse.issparse(detected) else detected, dtype=np.float64)
    count_array = np.asarray(counts.toarray() if sparse.issparse(counts) else counts, dtype=np.float64)

    n_cells_per_group = np.bincount(codes, minlength=n_groups).astype(np.float64)
    props = detected_array / np.maximum(n_cells_per_group[:, None], 1.0)
    total_counts = count_array.sum(axis=1)

    new_obs = _obs_for_groups(obs, groupCodes=codes, nGroups=n_groups)
    sample_values = [str(item[0]) for item in uniques]
    cluster_values = [str(item[1]) for item in uniques]
    new_obs[sampleKey] = sample_values
    new_obs[clusterKey] = cluster_values
    new_obs.index = pd.Index(
        [
            f"{sampleValue}_{clusterValue}"
            for sampleValue, clusterValue in zip(sample_values, cluster_values, strict=True)
        ],
        name="pseudobulk_id",
    )
    new_obs["psbulk_cells"] = n_cells_per_group.astype(int)
    new_obs["psbulk_counts"] = total_counts

    return AnnData(
        X=count_array,
        obs=new_obs,
        var=var.copy(),
        layers={"psbulk_props": props.astype(np.float64)},
    )


def attach_label_columns(
    pdata: AnnData,
    labelTable: pd.DataFrame,
    *,
    sampleKey: str,
    studyKey: str,
    diseaseNameKey: str,
) -> AnnData:
    labels_by_srx = labelTable.set_index("srxAccession")
    out = pdata.copy()
    sample_values = pd.Series(out.obs[sampleKey])
    out.obs["diseaseArea"] = sample_values.map(labels_by_srx["diseaseArea"]).astype("category")
    out.obs["diseased"] = sample_values.map(labels_by_srx["diseased"]).astype("boolean")
    out.obs[studyKey] = sample_values.map(labels_by_srx["studyAccession"])
    if "diseaseName" in labels_by_srx.columns:
        out.obs[diseaseNameKey] = sample_values.map(labels_by_srx["diseaseName"])
    out.obs[sampleKey] = as_string_series(sample_values)
    return out


def _sanitize_obs_for_h5ad(obs: pd.DataFrame) -> pd.DataFrame:
    """Convert mixed object columns to strings so AnnData can write them."""
    out = obs.copy()
    for col in out.columns:
        series = out[col]
        if pd.api.types.is_bool_dtype(series) or pd.api.types.is_numeric_dtype(series):
            continue
        if str(series.dtype) == "category":
            continue
        if str(series.dtype) == "boolean":
            continue
        out[col] = series.astype("string")
    return out


def write_checkpoint(
    pdata: AnnData,
    fingerprint: dict[str, Any],
    cfg: AtlasDeAnalysisConfig,
) -> None:
    cfg.checkpointsDir.mkdir(parents=True, exist_ok=True)
    pdata = pdata.copy()
    pdata.obs = _sanitize_obs_for_h5ad(pd.DataFrame(pdata.obs))
    pdata.write_h5ad(cfg.pseudobulkPath, compression=cfg.compression)
    cfg.fingerprintPath.write_text(json.dumps(fingerprint, indent=2) + "\n")
    log.info("Wrote pseudobulk checkpoint to %s", cfg.pseudobulkPath)


def load_checkpoint_if_valid(cfg: AtlasDeAnalysisConfig) -> AnnData | None:
    if not cfg.pseudobulkPath.exists() or not cfg.fingerprintPath.exists():
        return None
    expected = build_aggregate_fingerprint(cfg)
    stored = json.loads(cfg.fingerprintPath.read_text())
    if stored.get("sha256") != expected["sha256"]:
        raise ValueError(
            f"Stale pseudobulk checkpoint fingerprint. Delete {cfg.checkpointsDir} or rerun the aggregate stage."
        )
    log.info("Reusing valid pseudobulk checkpoint %s", cfg.pseudobulkPath)
    return ad.read_h5ad(cfg.pseudobulkPath)


def aggregate_atlas(cfg: AtlasDeAnalysisConfig, *, reuseCheckpoint: bool = True) -> AnnData:
    """Load the atlas once, aggregate sparse raw counts, checkpoint, and release cells."""
    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    cfg.figuresDir.mkdir(parents=True, exist_ok=True)
    cfg.checkpointsDir.mkdir(parents=True, exist_ok=True)

    fingerprint = build_aggregate_fingerprint(cfg)
    if reuseCheckpoint:
        existing = load_checkpoint_if_valid(cfg)
        if existing is not None:
            return existing

    estimate = estimate_raw_sparse_bytes(cfg.atlasPath)
    assert_memory_available(estimate, reserveBytes=cfg.memoryReserveBytes, logger=log)
    log_memory("pre-load", logger=log)

    log.info("Loading atlas from %s", cfg.atlasPath)
    adata = ad.read_h5ad(cfg.atlasPath)
    log_memory("post-load", logger=log)
    reject_dense_cell_matrix(adata.raw.X if adata.raw is not None else adata.X, label="loaded atlas")

    label_table = build_sample_label_table(cfg.contextsPath, cfg.atlasCsvPath, MetadataConfig().sampleParquetPath)
    label_table.to_csv(cfg.outputDir / "sample_labels.csv", index=False)

    annotated = annotate_obs_with_labels(
        pd.DataFrame(adata.obs),
        label_table,
        sampleKey=cfg.sampleKey,
        eligibleOnly=False,
    )
    eligible = annotated["eligible"].fillna(False).astype(bool)
    if not bool(eligible.any()):
        raise RuntimeError("No eligible cells after disease-marker labeling")

    eligible_idx = np.flatnonzero(eligible.to_numpy())
    obs_eligible = annotated.iloc[eligible_idx].copy()
    raw = sparse_raw_matrix(adata)
    # Row-slice the sparse raw matrix once; do not call raw.to_adata().
    matrix_eligible = raw[eligible_idx, :]
    var = pd.DataFrame(adata.raw.var.copy() if adata.raw is not None else adata.var.copy())
    log.info(
        "Eligible cells: %s; samples: %s; studies: %s; clusters: %s",
        f"{len(obs_eligible):,}",
        f"{obs_eligible[cfg.sampleKey].nunique():,}",
        f"{obs_eligible[cfg.studyKey].nunique():,}",
        f"{obs_eligible[cfg.clusterKey].nunique():,}",
    )
    cluster_interp = cluster_interpretation_table(
        obs_eligible,
        clusterKey=cfg.clusterKey,
        labelKey=cfg.labelKey,
        ontologyKey=cfg.ontologyKey,
        studyKey=cfg.studyKey,
        sampleKey=cfg.sampleKey,
        highPurity=cfg.highPurity,
        resolvedMinStudies=cfg.resolvedMinStudies,
    )
    cluster_interp.to_csv(cfg.outputDir / "cluster_interpretation.csv", index=False)

    contrast_support = same_study_contrast_support(
        obs_eligible,
        clusterKey=cfg.clusterKey,
        sampleKey=cfg.sampleKey,
        studyKey=cfg.studyKey,
        minCellsPerProfile=cfg.minCellsPerProfile,
    )
    contrast_support = contrast_support.merge(
        cluster_interp[["cluster", "interpretedCellType", "interpretationStatus", "topLabel", "topLabelFraction"]],
        on="cluster",
        how="left",
    )
    contrast_support.to_csv(cfg.outputDir / "same_study_contrast_support.csv", index=False)
    log_memory("post-filter", logger=log)

    # Drop the cell-level atlas before aggregation finishes writing outputs.
    del adata
    del raw
    gc.collect()
    log_memory("post-atlas-release", logger=log)

    pdata = aggregate_sparse_pseudobulk(
        matrix_eligible,
        obs_eligible,
        var,
        sampleKey=cfg.sampleKey,
        clusterKey=cfg.clusterKey,
    )
    del matrix_eligible
    del obs_eligible
    gc.collect()

    pdata = attach_label_columns(
        pdata,
        label_table,
        sampleKey=cfg.sampleKey,
        studyKey=cfg.studyKey,
        diseaseNameKey=cfg.diseaseNameKey,
    )
    pdata = filter_pseudobulk_profiles(pdata, minCellsPerProfile=cfg.minCellsPerProfile)
    pdata.obs[cfg.clusterKey] = as_string_series(pd.Series(pdata.obs[cfg.clusterKey]))
    log.info("Pseudobulk profiles after min-cell filter: %s", f"{pdata.n_obs:,}")
    pd.DataFrame(pdata.obs).to_csv(cfg.outputDir / "pseudobulk_profiles.csv", index=True)

    write_checkpoint(pdata, fingerprint, cfg)
    log_memory("post-checkpoint", logger=log)
    snap = snapshot_memory()
    log.info("Aggregation complete; peakRss=%s", format_bytes(snap.peakRssBytes))
    return pdata
