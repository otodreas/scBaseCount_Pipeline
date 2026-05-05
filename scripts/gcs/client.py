from __future__ import annotations

from pathlib import Path
from urllib.parse import urlparse

from google.auth.credentials import AnonymousCredentials
from google.cloud import storage as gcs

from shared.logger import configure_file_logger

_log = configure_file_logger("gcs.log", __name__)


def gcs_local_path(gs_uri: str, local_root: Path) -> Path:
    blob_name = urlparse(gs_uri).path.lstrip("/")
    return local_root / blob_name


def download_from_gcs(gs_uri: str, local_root: Path) -> Path:
    parsed = urlparse(gs_uri)
    bucket_name = parsed.netloc
    blob_name = parsed.path.lstrip("/")
    dest = local_root / blob_name
    dest.parent.mkdir(parents=True, exist_ok=True)
    client = gcs.Client(credentials=AnonymousCredentials())
    blob = client.bucket(bucket_name).blob(blob_name)
    _log.info("Downloading %s -> %s", gs_uri, dest)
    blob.download_to_filename(str(dest))
    return dest


def verify_download(gs_uri: str, local_root: Path) -> bool:
    dest = gcs_local_path(gs_uri, local_root)
    if dest.exists():
        _log.info("Download verified: %s", dest)
        return True
    _log.warning("Download verification failed: %s not found locally", dest)
    return False
