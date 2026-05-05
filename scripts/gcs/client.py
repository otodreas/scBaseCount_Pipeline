from __future__ import annotations

import os
from pathlib import Path
from urllib.parse import urlparse

from google.cloud import storage as gcs

from shared.logger import configure_file_logger

_log = configure_file_logger("gcs.log", __name__)


def download_from_gcs(gs_uri: str, dest: Path) -> None:
    parsed = urlparse(gs_uri)
    bucket_name = parsed.netloc
    blob_name = parsed.path.lstrip("/")
    dest.parent.mkdir(parents=True, exist_ok=True)
    client = gcs.Client(project=os.environ["GCP_PROJECT"])
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    _log.info("Downloading %s -> %s", gs_uri, dest)
    blob.download_to_filename(str(dest))
