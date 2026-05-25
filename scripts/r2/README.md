# r2

Uploads files to and queries an S3-compatible R2 bucket. Credentials are read from environment variables.

## Usage

```python
from pathlib import Path

from r2 import (
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
local_path = download_from_r2(r2_key, Path("data/r2_cache"))
```

## Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `fetch_uploaded_r2_keys()` | `set[str]` | Lists all object keys currently in the bucket |
| `gcs_uri_to_r2_raw_key(gs_uri)` | `str` | Maps a GCS URI to the matching R2 raw key |
| `upload_to_r2(local_path, r2_key)` | `None` | Uploads a local file to the given R2 key |
| `download_from_r2(r2_key, local_root)` | `Path` | Downloads an R2 object under `local_root` |
| `verify_upload(r2_key)` | `bool` | Returns `True` if the key exists (`head_object`), `False` on 404 |
| `r2_key_exists(r2_key)` | `bool` | Same existence check as `verify_upload` |
| `r2_object_md5(r2_key)` | `str \| None` | Base64 MD5 from object metadata, if present |
| `r2_raw_matches_gcs(r2_key, gs_uri)` | `bool` | Compares R2 object MD5 to the GCS source blob |

## Required environment variables

| Variable | Description |
|----------|-------------|
| `ENDPOINT_URL` | R2 S3-compatible endpoint URL |
| `AWS_ACCESS_KEY_ID` | R2 access key ID |
| `AWS_SECRET_ACCESS_KEY` | R2 secret access key |
| `BUCKET` | Bucket name |

## Logging

Steps are appended to `logs/r2.log`.
