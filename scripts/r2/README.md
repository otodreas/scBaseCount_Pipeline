# r2

Uploads files to and queries an S3-compatible R2 bucket. Credentials are read from environment variables.

## Usage

```python
from r2 import fetch_uploaded_r2_keys, upload_to_r2, verify_upload

# List all existing keys in the bucket
uploaded = fetch_uploaded_r2_keys()

# Upload a file
upload_to_r2(Path("output/cytetype/data/SRX12345678_cytetype_annotated.h5ad"), "cytetype/SRX12345678_cytetype_annotated.h5ad")

# Verify it landed
verify_upload("cytetype/SRX12345678_cytetype_annotated.h5ad")  # returns bool
```

## Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `fetch_uploaded_r2_keys()` | `set[str]` | Lists all object keys currently in the bucket |
| `upload_to_r2(local_path, r2_key)` | `None` | Uploads a local file to the given R2 key |
| `verify_upload(r2_key)` | `bool` | Returns `True` if the key exists in the bucket (`head_object`), `False` on 404 |

## Required environment variables

| Variable | Description |
|----------|-------------|
| `ENDPOINT_URL` | R2 S3-compatible endpoint URL |
| `AWS_ACCESS_KEY_ID` | R2 access key ID |
| `AWS_SECRET_ACCESS_KEY` | R2 secret access key |
| `BUCKET` | Bucket name |

## Logging

Steps are appended to `logs/r2.log`.
