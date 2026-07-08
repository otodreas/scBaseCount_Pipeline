from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
from shared.files import safe_delete
from storage import download_from_r2
from study_context.models import ExperimentContext

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason


def accession_from_r2_key(r2_key: str) -> str:
    """Return the accession stem from an R2 object key."""
    return Path(r2_key).stem


def resolve_study_accession(accession: str, contexts: dict[str, ExperimentContext]) -> str:
    """Resolve studyAccession from contexts; raise FileRejected when missing."""
    ctx = contexts.get(accession)
    if ctx is None or ctx.study is None:
        raise FileRejected(SkipReason.missing_study)
    return ctx.study.studyAccession


def cell_type_all_missing(adata: ad.AnnData, cell_type_key: str) -> bool:
    """Return True when every cell_type value is missing or blank."""
    if cell_type_key not in adata.obs.columns:
        return True
    values = adata.obs[cell_type_key].astype("string")
    non_missing = values.notna() & (values.str.strip() != "")
    return not bool(non_missing.any())


def fill_cell_type(adata: ad.AnnData, cfg: H5adConcatConfig) -> None:
    """Fill blank or NaN cell_type entries with cfg.missingLabel."""
    values = adata.obs[cfg.cellTypeKey].astype("string")
    adata.obs[cfg.cellTypeKey] = values.fillna(cfg.missingLabel).replace("", cfg.missingLabel)


def prefix_obs_names(adata: ad.AnnData, accession: str) -> None:
    """Prefix obs_names with accession so barcodes are globally unique across studies."""
    adata.obs_names = [f"{accession}_{name}" for name in adata.obs_names]


def prepare_accession(
    r2_key: str,
    cfg: H5adConcatConfig,
    contexts: dict[str, ExperimentContext],
    log: logging.Logger,
) -> tuple[Path, str]:
    """Download, validate, enrich, and stage one h5ad; return (prepared_path, studyAccession)."""
    accession = accession_from_r2_key(r2_key)
    raw_path = cfg.cacheDir / "raw" / f"{accession}.h5ad"
    prepared_path = cfg.cacheDir / "prepared" / f"{accession}.h5ad"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    prepared_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        download_from_r2(r2_key, raw_path, verify_md5=cfg.verifyMd5)
    except ValueError as exc:
        if "MD5 mismatch" in str(exc):
            safe_delete(raw_path, log)
            raise FileRejected(SkipReason.md5_mismatch) from exc
        raise

    try:
        study_accession = resolve_study_accession(accession, contexts)
        adata = ad.read_h5ad(raw_path)

        # TODO(preprocess): when cfg.preprocess is enabled, run cluster_validation.preprocess here
        # and raise FileRejected(SkipReason.preprocess_failed) on InsufficientCellsError or other failures.

        if cell_type_all_missing(adata, cfg.cellTypeKey):
            raise FileRejected(SkipReason.cell_type_all_missing)

        adata.obs[cfg.accessionKey] = accession
        fill_cell_type(adata, cfg)
        prefix_obs_names(adata, accession)
        adata.write_h5ad(prepared_path)
        return prepared_path, study_accession
    finally:
        safe_delete(raw_path, log)
