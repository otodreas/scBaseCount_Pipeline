# scBaseCount pipeline

Master's thesis, 30 hp, Lund University

A pipeline for large-scale scRNA-seq cluster labeling assessment, built on the Arc Institute's [Virtual Cell Atlas](https://console.cloud.google.com/storage/browser/arc-institute-virtual-cell-atlas?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))).

# Requirements

- Python >= 3.12.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) (package manager)

# Repository layout
The repo splits reusable code, batch orchestration, and interactive analysis:

- **[`scripts/`](scripts/)**: Importable Python packages shared across notebooks, pipelines, and ad hoc use. See [`scripts/README.md`](scripts/README.md).
- **[`pipelines/`](pipelines/)**: Batch runners for long, unattended jobs on a server (many accessions, sustained runtime). See [`pipelines/README.md`](pipelines/README.md).
- **[`notebooks/`](notebooks/)**: Interactive workflows for one-off or short tasks, and for repeatable steps where reviewing outputs (figures, tables, spot checks) is part of the work. See [`notebooks/README.md`](notebooks/README.md).

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
uv sync --locked --group dev
git config core.hooksPath .githooks   # once per clone
```

[`.githooks/pre-commit`](.githooks/pre-commit) runs ruff and nbstripout on staged files; [`.githooks/pre-push`](.githooks/pre-push) runs the cluster validation regression test when `scripts/cluster_validation/` changed. Both are optional local help; [CI](.github/workflows/ci.yml) enforces ruff, pytest, and stripped `notebooks/` on `main`.