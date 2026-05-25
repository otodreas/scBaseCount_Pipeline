# scBaseCount pipeline

Master's thesis, 30 hp, Lund University

A pipeline for large-scale scRNA-seq cluster labeling assessment, built on the Arc Institute's [Virtual Cell Atlas](https://console.cloud.google.com/storage/browser/arc-institute-virtual-cell-atlas?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))).

See [`pipelines/`](pipelines/) for the runner scripts (clustering, CyteType annotation, cluster stats, and GCS-to-R2 migration), [`notebooks/README.md`](notebooks/README.md) for the analysis notebooks, and [`scripts/README.md`](scripts/README.md) for the Python packages they share.

# Requirements

- Python >= 3.12.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) (package manager)

# Data access

## Google Cloud
A Google Cloud account and project is required to download data programmatically. The [Google Cloud SDK](https://cloud.google.com/sdk/docs/install-sdk) is used via `google-cloud-storage` (locked in [`uv.lock`](uv.lock)).

## Cloudflare R2
Processed `h5ad` files are stored in Cloudflare's S3-compatible R2 storage. Credentials are required (see [`.env.example`](.env.example)).

## Optional API keys
NCBI and CyteType API keys reduce rate limiting when fetching study metadata and running annotations.

# Setup

```sh
git clone git@github.com:otodreas/scBaseCount_Pipeline.git
cd scBaseCount_Pipeline
uv sync --group dev
git config core.hooksPath .githooks
```

The `git config core.hooksPath` step enables the pre-commit hook in [`.githooks/pre-commit`](.githooks/pre-commit), which runs ruff and strips outputs from `notebooks/` and `tests/` notebooks (not `docs/`). Required once per clone.