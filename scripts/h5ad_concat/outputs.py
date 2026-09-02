import json
import logging
from pathlib import Path, PurePosixPath

from shared.files import safe_delete
from shared.repo import rel_to_repo
from storage import r2_key_exists, r2_object_md5, upload_to_r2, verify_upload
from storage.transfer import _MD5_METADATA_KEY, _local_md5_b64

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.models import FileRecord, H5adConcatResult


def file_log_path(output_path: Path) -> Path:
    """Return the local per-file JSONL log path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_files.jsonl")


def config_manifest_path(output_path: Path) -> Path:
    """Return the local config manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_config.json")


def result_manifest_path(output_path: Path) -> Path:
    """Return the local result manifest path for the given atlas output path."""
    return output_path.with_name(f"{output_path.stem}_result.json")


def _sibling_r2_key(r2_key: str, suffix: str) -> str:
    """Return an R2 key sharing the stem of r2_key with the given suffix (e.g. '_files.jsonl', '_result.json')."""
    base = PurePosixPath(r2_key)
    return str(base.with_name(f"{base.stem}{suffix}"))


def ensure_atlas_targets_absent(cfg: H5adConcatConfig) -> None:
    """Fail before any artifact mutation when the configured local or primary R2 atlas already exists."""
    if cfg.outputPath.exists():
        raise FileExistsError(f"Local atlas already exists: {cfg.outputPath}")

    if cfg.uploadAtlas and cfg.atlasR2Key and r2_key_exists(cfg.atlasR2Key):
        raise FileExistsError(f"R2 atlas already exists: {cfg.atlasR2Key}")


def write_config_manifest(cfg: H5adConcatConfig, log: logging.Logger) -> Path:
    """Write the config manifest next to the atlas output at run start and return its path."""
    config_path = config_manifest_path(cfg.outputPath)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    payload = cfg.model_dump(mode="json")
    for key in ("datasetsPath", "geneInfoPath", "cacheDir", "outputPath"):
        payload[key] = rel_to_repo(Path(payload[key]))
    config_path.write_text(json.dumps(payload, indent=2))
    log.info("Wrote config manifest to %s", rel_to_repo(config_path))
    return config_path


def init_file_log(log_path: Path) -> None:
    """Create or truncate the per-file JSONL log at run start. A new run overwrites any previous log at this path."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("")


def append_file_record(log_path: Path, record: FileRecord) -> None:
    """Append one complete FileRecord after the accession is processed.

    These records match H5adConcatResult.files. They are flushed before concat so completed
    files remain on disk if the run fails before the result manifest is written.

    Accessions survive if run never writes _result.json manifest.
    """
    with log_path.open("a") as handle:
        handle.write(record.model_dump_json() + "\n")


def _verify_or_raise(r2_key: str, log: logging.Logger) -> None:
    """Verify an R2 upload, raising RuntimeError when verification fails."""
    if not verify_upload(r2_key):
        msg = f"Upload verification failed for {r2_key}"
        log.error(msg)
        raise RuntimeError(msg)


def _verify_atlas_checksum(r2_key: str, local_md5: str, log: logging.Logger) -> None:
    """Confirm the R2 gcs-md5 metadata matches the pre-upload local MD5."""
    remote_md5 = r2_object_md5(r2_key)
    if remote_md5 is None:
        msg = f"Atlas upload missing gcs-md5 metadata for {r2_key}"
        log.error(msg)
        raise RuntimeError(msg)
    if remote_md5 != local_md5:
        msg = f"Atlas upload MD5 mismatch for {r2_key}: local={local_md5} remote={remote_md5}"
        log.error(msg)
        raise RuntimeError(msg)


def _result_payload(result: H5adConcatResult) -> dict:
    """Serialize the result with repository-relative local paths."""
    payload = result.model_dump(mode="json")
    for key in ("outputPath", "fileLogPath", "configPath"):
        value = payload.get(key)
        if value:
            payload[key] = rel_to_repo(Path(value))
    return payload


def finalize_outputs(
    cfg: H5adConcatConfig,
    output_path: Path,
    file_log: Path,
    result: H5adConcatResult,
    log: logging.Logger,
) -> None:
    """Write the result JSON manifest locally and, when uploadAtlas is set, upload artifacts then delete the local atlas only after checksum verification."""
    result_path = result_manifest_path(output_path)

    if cfg.uploadAtlas and cfg.atlasR2Key:
        r2_key = cfg.atlasR2Key
        config_path = config_manifest_path(output_path)
        file_log_r2_key = _sibling_r2_key(r2_key, "_files.jsonl")
        config_r2_key = _sibling_r2_key(r2_key, "_config.json")
        result_r2_key = _sibling_r2_key(r2_key, "_result.json")
        try:
            md5 = _local_md5_b64(output_path)
            log.info("Uploading atlas to R2 key %s", r2_key)
            upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
            _verify_or_raise(r2_key, log)
            _verify_atlas_checksum(r2_key, md5, log)
            log.info("Atlas uploaded and checksum-verified at %s", r2_key)

            log.info("Uploading file log to R2 key %s", file_log_r2_key)
            upload_to_r2(file_log, file_log_r2_key)
            _verify_or_raise(file_log_r2_key, log)
            log.info("File log uploaded and verified at %s", file_log_r2_key)

            log.info("Uploading config manifest to R2 key %s", config_r2_key)
            upload_to_r2(config_path, config_r2_key)
            _verify_or_raise(config_r2_key, log)
            log.info("Config manifest uploaded and verified at %s", config_r2_key)

            result.atlasFileLogR2Key = file_log_r2_key
            result.atlasConfigR2Key = config_r2_key
            result.atlasResultR2Key = result_r2_key
            # Local .h5ad is deleted below; point outputPath at the surviving result manifest.
            result.outputPath = result_path

            result_path.write_text(json.dumps(_result_payload(result), indent=2))
            log.info("Wrote result manifest to %s", rel_to_repo(result_path))

            log.info("Uploading result manifest to R2 key %s", result_r2_key)
            upload_to_r2(result_path, result_r2_key)
            _verify_or_raise(result_r2_key, log)
            log.info("Result manifest uploaded and verified at %s", result_r2_key)

            safe_delete(output_path, log)
        except Exception:
            log.exception("Atlas upload failed for %s; local atlas retained at %s", r2_key, rel_to_repo(output_path))
            raise
        return

    result_path.write_text(json.dumps(_result_payload(result), indent=2))
    log.info("Wrote result manifest to %s", rel_to_repo(result_path))
