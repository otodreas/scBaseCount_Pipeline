import json
import logging
from pathlib import Path

from shared.files import safe_delete
from storage import upload_to_r2, verify_upload
from storage.transfer import _MD5_METADATA_KEY, _local_md5_b64

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.models import H5adConcatResult


def finalize_outputs(
    cfg: H5adConcatConfig,
    output_path: Path,
    result: H5adConcatResult,
    log: logging.Logger,
) -> None:
    """Upload the atlas to R2, write a JSON manifest, and delete the local h5ad on success."""
    if not cfg.uploadAtlas or not cfg.atlasR2Key:
        return

    r2_key = cfg.atlasR2Key
    try:
        md5 = _local_md5_b64(output_path)
        log.info("Uploading atlas to R2 key %s", r2_key)
        upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
        if not verify_upload(r2_key):
            msg = f"Upload verification failed for {r2_key}"
            log.error(msg)
            raise RuntimeError(msg)
        log.info("Atlas uploaded and verified at %s", r2_key)
        manifest_path = output_path.with_suffix(".json")
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        # Local .h5ad is deleted below; point outputPath at the surviving manifest.
        result.outputPath = manifest_path
        manifest_path.write_text(json.dumps(result.model_dump(mode="json"), indent=2))
        log.info("Wrote atlas manifest to %s", manifest_path)
        safe_delete(output_path, log)
    except Exception:
        log.exception("Atlas upload failed for %s", r2_key)
        raise
