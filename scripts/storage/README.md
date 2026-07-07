# storage

Cloud storage helpers for downloading from public GCS buckets and uploading to or querying an S3-compatible R2 bucket.

## GCS usage

Downloads files from public Google Cloud Storage buckets. Uses anonymous credentials; no GCP project or service account required.

```python
from storage import download_from_gcs, gcs_local_path, verify_download
from shared.repo import REPO_ROOT

GCS_LOCAL_ROOT = REPO_ROOT / "data"
gs_uri = "gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX12345678.h5ad"

local_path = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)
dest = download_from_gcs(gs_uri, GCS_LOCAL_ROOT)
verify_download(gs_uri, GCS_LOCAL_ROOT)
```

The local path mirrors the GCS blob path under `local_root`: the bucket name is stripped and the remainder of the URI becomes the relative path.

## R2 usage

Credentials are read from environment variables.

```python
from pathlib import Path

from storage import (
    download_from_r2,
    fetch_uploaded_r2_keys,
    gcs_uri_to_r2_raw_key,
    r2_key_exists,
    upload_to_r2,
    verify_upload,
)

uploaded = fetch_uploaded_r2_keys()
r2_key = gcs_uri_to_r2_raw_key("gs://bucket/path/SRX12345678.h5ad")

upload_to_r2(Path("output/cytetype/data/SRX12345678_cytetype_annotated.h5ad"), r2_key)
verify_upload(r2_key)
download_from_r2(r2_key, Path("data/r2_cache"))
```

## Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `gcs_local_path(gs_uri, local_root)` | `Path` | Derives the local destination path from the GCS URI without downloading |
| `download_from_gcs(gs_uri, local_root)` | `Path` | Downloads the file, creating parent directories as needed |
| `verify_download(gs_uri, local_root)` | `bool` | Returns `True` if the expected local file exists, `False` otherwise |
| `gcs_blob_md5(gs_uri)` | `str` | Base64 MD5 from GCS blob metadata |
| `fetch_uploaded_r2_keys(prefix=None)` | `set[str]` | Lists object keys in the bucket, optionally filtered by prefix |
| `gcs_uri_to_r2_raw_key(gs_uri)` | `str` | Maps a GCS URI to the matching R2 raw key |
| `upload_to_r2(local_path, r2_key)` | `None` | Uploads a local file to the given R2 key |
| `download_from_r2(r2_key, local_path)` | `None` | Downloads an R2 object to `local_path` |
| `verify_upload(r2_key)` | `bool` | Returns `True` if the key exists (`head_object`), `False` on 404 |
| `r2_key_exists(r2_key)` | `bool` | Same existence check as `verify_upload` |
| `r2_object_md5(r2_key)` | `str \| None` | Base64 MD5 from object metadata, if present |
| `r2_raw_matches_gcs(r2_key, gcs_md5)` | `bool` | Compares R2 object MD5 to the GCS source blob |

## Required environment variables (R2)

| Variable | Description |
|----------|-------------|
| `ENDPOINT_URL` | R2 S3-compatible endpoint URL |
| `AWS_ACCESS_KEY_ID` | R2 access key ID |
| `AWS_SECRET_ACCESS_KEY` | R2 secret access key |
| `BUCKET` | Bucket name |

## Logging

GCS steps are appended to `logs/gcs.log`. R2 steps are appended to `logs/r2.log`.
