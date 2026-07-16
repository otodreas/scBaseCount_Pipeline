from __future__ import annotations

import pandas as pd
from shared.logger import configure_file_logger
from storage import gcs_uri_to_r2_raw_key
from study_context.utils import load_contexts_jsonl

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import H5adConcatResult, SkippedFile
from h5ad_concat.outputs import (
    append_status_row,
    finalize_outputs,
    init_status_csv,
    status_csv_path,
    write_config_manifest,
)
from h5ad_concat.prepare import accession_from_r2_key, prepare_adata
from h5ad_concat.reference import load_gene_reference

_log = configure_file_logger("h5ad_concat.log", __name__)


def resolve_r2_keys(cfg: H5adConcatConfig) -> list[str]:
    """Return explicit r2Keys or resolve them from datasets.csv file_path URIs."""
    if cfg.datasetsPath is not None:
        if not cfg.datasetsPath.is_file():
            raise FileNotFoundError(f"datasets file not found at {cfg.datasetsPath}")
        datasets = pd.read_csv(cfg.datasetsPath)
        keys = [gcs_uri_to_r2_raw_key(uri) for uri in datasets["file_path"]]
        _log.info("Resolved %d R2 key(s) from %s", len(keys), cfg.datasetsPath)
        return keys
    return cfg.r2Keys or []


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    r2_keys = resolve_r2_keys(cfg)
    _log.info("Starting h5ad_concat run: %d key(s)", len(r2_keys))

    csv_path = status_csv_path(cfg.outputPath)
    init_status_csv(csv_path)
    config_path = write_config_manifest(cfg, _log)
    contexts = load_contexts_jsonl(cfg.contextsPath)
    reference = load_gene_reference(cfg.geneInfoPath)
    skipped: list[SkippedFile] = []
    adatas = []
    studies: list[str] = []
    accessions: list[str] = []

    try:
        for r2_key in r2_keys:
            accession = accession_from_r2_key(r2_key)
            try:
                adata, study_accession = prepare_adata(r2_key, accession, cfg, contexts, reference, _log)
                adatas.append(adata)
                studies.append(study_accession)
                accessions.append(accession)
                append_status_row(csv_path, accession, r2_key, "success", "", study_accession)
            except FileRejected as exc:
                detail = f": {exc.__cause__}" if exc.__cause__ else ""
                _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
                skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))
                append_status_row(csv_path, accession, r2_key, "skip", exc.reason.value, "")

        if not adatas:
            msg = "No files passed validation; nothing to concatenate"
            raise ValueError(msg)

        try:
            atlas = concat_atlas(adatas, accessions, cfg, _log)
        except Exception:
            _log.exception("concat failed")
            raise

        try:
            output_path = write_atlas(atlas, cfg, _log)
        except Exception:
            _log.exception("write failed")
            raise

        result = H5adConcatResult(
            outputPath=output_path,
            nObs=atlas.n_obs,
            nVars=atlas.n_vars,
            nFilesConcatenated=len(adatas),
            studiesSeen=sorted(set(studies)),
            skipped=skipped,
            statusCsvPath=csv_path,
            configPath=config_path,
            atlasR2Key=cfg.atlasR2Key if cfg.uploadAtlas else None,
            conserveLayers=cfg.conserveLayers,
        )

        finalize_outputs(cfg, output_path, csv_path, result, _log)

        _log.info(
            "h5ad_concat run complete: %d concatenated, %d skipped",
            len(adatas),
            len(skipped),
        )
        return result
    except KeyboardInterrupt:
        _log.warning("h5ad_concat run interrupted (KeyboardInterrupt)")
        raise
