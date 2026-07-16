import csv
import json
import logging
from pathlib import Path, PurePosixPath

from shared.files import safe_delete
from storage import upload_to_r2, verify_upload
from storage.transfer import _MD5_METADATA_KEY, _local_md5_b64

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.models import H5adConcatResult

_STATUS_CSV_HEADER = ["accession", "r2Key", "status", "reason", "studyAccession"]


def status_csv_path(output_path: Path) -> Path:
    """Return the local status CSV path for the given atlas output path."""
    return output_path.with_suffix(".csv")


def config_manifest_path(output_path: Path) -> Path:
    """Return the local config manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_config.json")


def result_manifest_path(output_path: Path) -> Path:
    """Return the local result manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_result.json")


def _sibling_r2_key(r2_key: str, suffix: str) -> str:
    """Return an R2 key sharing the stem of r2_key with the given suffix (e.g. '.csv', '_result.json')."""
    base = PurePosixPath(r2_key)
    return str(base.with_name(f"{base.stem}{suffix}"))


def write_config_manifest(cfg: H5adConcatConfig, log: logging.Logger) -> Path:
    """Write the config manifest next to the atlas output at run start and return its path."""
    config_path = config_manifest_path(cfg.outputPath)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(cfg.model_dump(mode="json"), indent=2))
    log.info("Wrote config manifest to %s", config_path)
    return config_path


def init_status_csv(csv_path: Path) -> None:
    """Create or truncate the status CSV and write its header row at run start."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        csv.writer(f).writerow(_STATUS_CSV_HEADER)


def append_status_row(
    csv_path: Path, accession: str, r2_key: str, status: str, reason: str, study_accession: str
) -> None:
    """Append one per-file status row; opened in append mode so each row is flushed to disk and a partial ledger survives interruption."""
    with csv_path.open("a", newline="") as f:
        csv.DictWriter(f, fieldnames=_STATUS_CSV_HEADER).writerow(
            {
                "accession": accession,
                "r2Key": r2_key,
                "status": status,
                "reason": reason,
                "studyAccession": study_accession,
            }
        )


def _verify_or_raise(r2_key: str, log: logging.Logger) -> None:
    """Verify an R2 upload, raising RuntimeError when verification fails."""
    if not verify_upload(r2_key):
        msg = f"Upload verification failed for {r2_key}"
        log.error(msg)
        raise RuntimeError(msg)


def finalize_outputs(
    cfg: H5adConcatConfig,
    output_path: Path,
    csv_path: Path,
    result: H5adConcatResult,
    log: logging.Logger,
) -> None:
    """Write the result JSON manifest locally and, when uploadAtlas is set, upload the atlas, status CSV, and manifests to R2 then delete the local atlas."""
    result_path = result_manifest_path(output_path)

    if cfg.uploadAtlas and cfg.atlasR2Key:
        r2_key = cfg.atlasR2Key
        config_path = config_manifest_path(output_path)
        status_r2_key = _sibling_r2_key(r2_key, ".csv")
        config_r2_key = _sibling_r2_key(r2_key, "_config.json")
        result_r2_key = _sibling_r2_key(r2_key, "_result.json")
        try:
            md5 = _local_md5_b64(output_path)
            log.info("Uploading atlas to R2 key %s", r2_key)
            upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
            _verify_or_raise(r2_key, log)
            log.info("Atlas uploaded and verified at %s", r2_key)

            log.info("Uploading status CSV to R2 key %s", status_r2_key)
            upload_to_r2(csv_path, status_r2_key)
            _verify_or_raise(status_r2_key, log)
            log.info("Status CSV uploaded and verified at %s", status_r2_key)

            log.info("Uploading config manifest to R2 key %s", config_r2_key)
            upload_to_r2(config_path, config_r2_key)
            _verify_or_raise(config_r2_key, log)
            log.info("Config manifest uploaded and verified at %s", config_r2_key)

            result.atlasStatusR2Key = status_r2_key
            result.atlasConfigR2Key = config_r2_key
            result.atlasResultR2Key = result_r2_key
            # Local .h5ad is deleted below; point outputPath at the surviving result manifest.
            result.outputPath = result_path

            result_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
            log.info("Wrote result manifest to %s", result_path)

            log.info("Uploading result manifest to R2 key %s", result_r2_key)
            upload_to_r2(result_path, result_r2_key)
            _verify_or_raise(result_r2_key, log)
            log.info("Result manifest uploaded and verified at %s", result_r2_key)

            safe_delete(output_path, log)
        except Exception:
            log.exception("Atlas upload failed for %s", r2_key)
            raise
        return

    result_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
    log.info("Wrote result manifest to %s", result_path)
