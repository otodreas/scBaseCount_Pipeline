from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
from botocore.exceptions import BotoCoreError, ClientError
from shared.files import safe_delete
from storage import download_from_r2
from study_context.models import ExperimentContext

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason
from h5ad_concat.qc import apply_qc_gate
from h5ad_concat.reference import GeneReference, align_to_reference


def accession_from_r2_key(r2_key: str) -> str:
    """Return the accession stem from an R2 object key."""
    return Path(r2_key).stem


def resolve_batch_key(accession: str, contexts: dict[str, ExperimentContext]) -> str:
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


def validate_single_accession(adata: ad.AnnData, accession: str, cfg: H5adConcatConfig) -> None:
    """Raise FileRejected when obs accessionKey is not a single value equal to accession."""
    found = adata.obs[cfg.accessionKey].unique().tolist()
    if len(found) != 1 or found[0] != accession:
        msg = f"expected single accession {accession}, found {found}"
        raise FileRejected(SkipReason.accession_mismatch) from ValueError(msg)


def prepare_adata(
    r2_key: str,
    accession: str,
    cfg: H5adConcatConfig,
    contexts: dict[str, ExperimentContext],
    reference: GeneReference,
    log: logging.Logger,
) -> tuple[ad.AnnData, str]:
    """Download, validate, enrich one h5ad in memory; return (AnnData, studyAccession)."""
    study_accession = resolve_batch_key(accession, contexts)
    raw_path = cfg.cacheDir / "raw" / f"{accession}.h5ad"

    # Attempt to download the file and verify with MD5
    try:
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        download_from_r2(r2_key, raw_path, verify_md5=cfg.verifyMd5)
    except ValueError as exc:  # md5 mismatch is ValueError subclass
        safe_delete(raw_path, log)
        if "MD5 mismatch" in str(exc):
            raise FileRejected(SkipReason.md5_mismatch) from exc
        raise FileRejected(SkipReason.download_failed) from exc
    except (ClientError, BotoCoreError, OSError) as exc:
        safe_delete(raw_path, log)
        raise FileRejected(SkipReason.download_failed) from exc

    # Attempt to read the file and validate the AnnData object
    try:
        adata = ad.read_h5ad(raw_path)
    except Exception as exc:
        safe_delete(raw_path, log)
        raise FileRejected(SkipReason.read_failed) from exc

    # Check that cell types are present
    try:
        if cfg.preprocess:
            adata, qc_stats = apply_qc_gate(adata, cfg)
            log.info(
                "%s: QC kept %d/%d cells (%.1f%% retained)",
                accession,
                qc_stats.nCellsAfter,
                qc_stats.nCellsBefore,
                qc_stats.pctCellsAfter * 100.0,
            )

        adata, align_stats = align_to_reference(adata, reference)
        log.info(
            "%s: aligned %d genes (%d zero-filled, %d dropped)",
            accession,
            align_stats.nGenesMapped,
            align_stats.nGenesZeroFilled,
            align_stats.nGenesDropped,
        )

        validate_single_accession(adata, accession, cfg)

        if cell_type_all_missing(adata, cfg.cellTypeKey):
            raise FileRejected(SkipReason.cell_type_all_missing)

        adata.obs[cfg.batchKey] = study_accession
        fill_cell_type(adata, cfg)
        return adata, study_accession

    finally:
        safe_delete(raw_path, log)
