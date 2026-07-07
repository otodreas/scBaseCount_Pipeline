from __future__ import annotations

import base64
import hashlib
from pathlib import Path
from urllib.parse import urlparse

from shared.logger import configure_file_logger

from storage.r2 import r2_object_md5

_log = configure_file_logger("r2.log", __name__)

_MD5_METADATA_KEY = "gcs-md5"


def _local_md5_b64(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return base64.b64encode(h.digest()).decode()


def gcs_uri_to_r2_raw_key(gs_uri: str) -> str:
    parsed = urlparse(gs_uri)
    bucket = parsed.netloc
    blob = parsed.path.lstrip("/")
    return f"{bucket}/{blob}"


def r2_raw_matches_gcs(r2_key: str, gcs_md5: str) -> bool:
    stored = r2_object_md5(r2_key)
    if stored is None:
        return False
    matches = stored == gcs_md5
    if not matches:
        _log.warning("MD5 mismatch for r2://%s: stored=%s expected=%s", r2_key, stored, gcs_md5)
    return matches
