from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from shared.files import safe_delete
from shared.logger import configure_file_logger
from storage import gcs_uri_to_r2_raw_key, upload_to_r2, verify_upload
from storage.transfer import _MD5_METADATA_KEY, _local_md5_b64
from study_context.utils import load_contexts_jsonl

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import H5adConcatResult, SkippedFile
from h5ad_concat.prepare import accession_from_r2_key, prepare_adata

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


def _upload_atlas_and_finalize(
    cfg: H5adConcatConfig,
    output_path: Path,
    result: H5adConcatResult,
) -> None:
    """Upload the atlas to R2, write a JSON manifest, and delete the local h5ad on success."""
    if not cfg.uploadAtlas or not cfg.atlasR2Key:
        return

    r2_key = cfg.atlasR2Key
    try:
        md5 = _local_md5_b64(output_path)
        _log.info("Uploading atlas to R2 key %s", r2_key)
        upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
        if not verify_upload(r2_key):
            msg = f"Upload verification failed for {r2_key}"
            _log.error(msg)
            raise RuntimeError(msg)
        _log.info("Atlas uploaded and verified at %s", r2_key)
        manifest_path = output_path.with_suffix(".json")
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        # Local .h5ad is deleted below; point outputPath at the surviving manifest.
        result.outputPath = manifest_path
        manifest_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
        _log.info("Wrote atlas manifest to %s", manifest_path)
        safe_delete(output_path, _log)
    except Exception:
        _log.exception("Atlas upload failed for %s", r2_key)
        raise


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    r2_keys = resolve_r2_keys(cfg)
    contexts = load_contexts_jsonl(cfg.contextsPath)
    skipped: list[SkippedFile] = []
    adatas = []
    studies: list[str] = []
    accessions: list[str] = []

    for r2_key in r2_keys:
        accession = accession_from_r2_key(r2_key)
        try:
            adata, study_accession = prepare_adata(r2_key, accession, cfg, contexts, _log)
            adatas.append(adata)
            studies.append(study_accession)
            accessions.append(accession)
        except FileRejected as exc:
            detail = f": {exc.__cause__}" if exc.__cause__ else ""
            # Warn on skipped files, continue pipeline
            _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
            skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))
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
        atlasR2Key=cfg.atlasR2Key if cfg.uploadAtlas else None,
    )

    _upload_atlas_and_finalize(cfg, output_path, result)

    return result
