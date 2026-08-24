# scBaseCount pipeline

Master's project, 30 hp, Lund University

A pipeline for large-scale scRNA-seq cluster labeling assessment, built on the Arc Institute's [Virtual Cell Atlas](https://console.cloud.google.com/storage/browser/arc-institute-virtual-cell-atlas?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))).

# Requirements

- [uv](https://docs.astral.sh/uv/getting-started/installation/) (package manager)
- Python >= 3.12.12 (installed automatically by uv)

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

# On the work presented

The work presented here was done in conjunction with Nygen Analytics AB, a private, for-profit company in Lund. The work was exploratory in many regards, and therefore many tasks we embarked on did not reach the final report. For instance, these include CyteType integration efforts and differential expression analyses on the atlas. During the course of the project, reports such as [STATE vs Leiden](writeups/state_vs_leiden/README.md) were written up for internal discussions, but were ultimately deemed unnecessary or out of scope for the final report. The final report is a condensed version of the work, and therefore many details are omitted.

Furthermore, Nygen uses Cursor, an AI-powered integrated development environment. I, along with everyone else at Nygen, used Cursor to draft code and brainstorm code architecture. I did not use Cursor to write reports or analyze and interpret data. The repository presented here represents my own work.