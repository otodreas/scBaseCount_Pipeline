# gcs

Downloads files from public Google Cloud Storage buckets. Uses anonymous credentials -- no GCP project or service account required.

## Usage

```python
from gcs import download_from_gcs, gcs_local_path, verify_download
from shared.repo import REPO_ROOT

GCS_LOCAL_ROOT = REPO_ROOT / "data"
gs_uri = "gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX12345678.h5ad"

# Compute the local path without downloading
local_path = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)
# -> data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX12345678.h5ad

# Download, preserving the GCS path structure under local_root
dest = download_from_gcs(gs_uri, GCS_LOCAL_ROOT)

# Check the file landed correctly
verify_download(gs_uri, GCS_LOCAL_ROOT)
```

The local path mirrors the GCS blob path under `local_root`: the bucket name is stripped and the remainder of the URI becomes the relative path.

## Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `gcs_local_path(gs_uri, local_root)` | `Path` | Derives the local destination path from the GCS URI without downloading |
| `download_from_gcs(gs_uri, local_root)` | `Path` | Downloads the file, creating parent directories as needed |
| `verify_download(gs_uri, local_root)` | `bool` | Returns `True` if the expected local file exists, `False` otherwise |

## Logging

Steps are appended to `logs/gcs.log`.
