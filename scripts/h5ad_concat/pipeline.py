from __future__ import annotations

import csv
import json
from pathlib import Path, PurePosixPath

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
from h5ad_concat.reference import load_gene_reference

_log = configure_file_logger("h5ad_concat.log", __name__)

_STATUS_CSV_HEADER = ["accession", "r2Key", "status", "reason", "studyAccession"]


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


def _config_manifest_path(output_path: Path) -> Path:
    """Return the local config manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_config.json")


def _result_manifest_path(output_path: Path) -> Path:
    """Return the local result manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_result.json")


def _sibling_r2_key(r2_key: str, suffix: str) -> str:
    """Return an R2 key sharing the stem of r2_key with the given suffix (e.g. '.csv', '_result.json')."""
    base = PurePosixPath(r2_key)
    return str(base.with_name(f"{base.stem}{suffix}"))


def _verify_or_raise(r2_key: str) -> None:
    """Verify an R2 upload, raising RuntimeError when verification fails."""
    if not verify_upload(r2_key):
        msg = f"Upload verification failed for {r2_key}"
        _log.error(msg)
        raise RuntimeError(msg)


def _finalize_outputs(
    cfg: H5adConcatConfig,
    output_path: Path,
    csv_path: Path,
    result: H5adConcatResult,
) -> None:
    """Write config and result JSON manifests locally and optionally upload the atlas, status CSV, and manifests to R2."""
    config_path = _config_manifest_path(output_path)
    result_path = _result_manifest_path(output_path)

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(cfg.model_dump(mode="json"), indent=2))
    _log.info("Wrote config manifest to %s", config_path)

    if cfg.uploadAtlas and cfg.atlasR2Key:
        r2_key = cfg.atlasR2Key
        status_r2_key = _sibling_r2_key(r2_key, ".csv")
        config_r2_key = _sibling_r2_key(r2_key, "_config.json")
        result_r2_key = _sibling_r2_key(r2_key, "_result.json")
        try:
            md5 = _local_md5_b64(output_path)
            _log.info("Uploading atlas to R2 key %s", r2_key)
            upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
            _verify_or_raise(r2_key)
            _log.info("Atlas uploaded and verified at %s", r2_key)

            _log.info("Uploading status CSV to R2 key %s", status_r2_key)
            upload_to_r2(csv_path, status_r2_key)
            _verify_or_raise(status_r2_key)
            _log.info("Status CSV uploaded and verified at %s", status_r2_key)

            _log.info("Uploading config manifest to R2 key %s", config_r2_key)
            upload_to_r2(config_path, config_r2_key)
            _verify_or_raise(config_r2_key)
            _log.info("Config manifest uploaded and verified at %s", config_r2_key)

            result.atlasStatusR2Key = status_r2_key
            result.atlasConfigR2Key = config_r2_key
            result.atlasResultR2Key = result_r2_key
            # Local .h5ad is deleted below; point outputPath at the surviving result manifest.
            result.outputPath = result_path

            result_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
            _log.info("Wrote result manifest to %s", result_path)

            _log.info("Uploading result manifest to R2 key %s", result_r2_key)
            upload_to_r2(result_path, result_r2_key)
            _verify_or_raise(result_r2_key)
            _log.info("Result manifest uploaded and verified at %s", result_r2_key)

            safe_delete(output_path, _log)
        except Exception:
            _log.exception("Atlas upload failed for %s", r2_key)
            raise
        return

    result_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
    _log.info("Wrote result manifest to %s", result_path)


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    r2_keys = resolve_r2_keys(cfg)
    _log.info("Starting h5ad_concat run: %d key(s)", len(r2_keys))

    csv_path = cfg.outputPath.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    contexts = load_contexts_jsonl(cfg.contextsPath)
    reference = load_gene_reference(cfg.geneInfoPath)
    skipped: list[SkippedFile] = []
    adatas = []
    studies: list[str] = []
    accessions: list[str] = []

    try:
        with csv_path.open("w", newline="") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=_STATUS_CSV_HEADER)
            writer.writeheader()
            csv_file.flush()

            for r2_key in r2_keys:
                accession = accession_from_r2_key(r2_key)
                try:
                    adata, study_accession = prepare_adata(r2_key, accession, cfg, contexts, reference, _log)
                    adatas.append(adata)
                    studies.append(study_accession)
                    accessions.append(accession)
                    writer.writerow(
                        {
                            "accession": accession,
                            "r2Key": r2_key,
                            "status": "success",
                            "reason": "",
                            "studyAccession": study_accession,
                        }
                    )
                except FileRejected as exc:
                    detail = f": {exc.__cause__}" if exc.__cause__ else ""
                    _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
                    skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))
                    writer.writerow(
                        {
                            "accession": accession,
                            "r2Key": r2_key,
                            "status": "skip",
                            "reason": exc.reason.value,
                            "studyAccession": "",
                        }
                    )
                csv_file.flush()

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
            configPath=_config_manifest_path(output_path),
            atlasR2Key=cfg.atlasR2Key if cfg.uploadAtlas else None,
            conserveLayers=cfg.conserveLayers,
        )

        _finalize_outputs(cfg, output_path, csv_path, result)

        _log.info(
            "h5ad_concat run complete: %d concatenated, %d skipped",
            len(adatas),
            len(skipped),
        )
        return result
    except KeyboardInterrupt:
        _log.warning("h5ad_concat run interrupted (KeyboardInterrupt)")
        raise
