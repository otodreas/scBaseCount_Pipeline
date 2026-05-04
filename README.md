# scBaseCount pipeline summary

Master's project, 30 hp, Lund University

This project aims to build a pipeline to run on large amounts of scRNA-seq data and assess cluster labeling performance

Uses the Arc Institute's [Virtual Cell Atlas](https://console.cloud.google.com/storage/browser/arc-institute-virtual-cell-atlas?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))

# Requirements

## Essential requirements
- Python >= 3.12.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) (package manager)

## Data management
These tools are not strictly required, but I use them to manage the large amounts of data going in and out of my programs.

### Google Cloud
A Google Cloud account and project is needed to download the data programmatically

I use the [Google Cloud SDK](https://docs.cloud.google.com/sdk/docs/install-sdk) via `google-cloud-storage` (locked in `uv.lock`) to access the dataset. Some helpful variables

```sh
export META_PATH="gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens"  # metadata path
export H5AD_PATH="gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens"      # h5ad path
```

You can use the Google Cloud CLI to get a general sense for the data in scBaseCount. Use `gcloud storage` followed by normal shell programs

```sh
gcloud storage ls "$META_PATH"                       # list files in the metadata directory
gcloud storage cp "$H5AD_PATH"/SRX123456.h5ad data/  # download an h5ad into data directory
```

### Cloudflare R2
I store processed `h5ad` files in Cloudflare's S3-compatible R2 storage. For this, you need some credentials (see `.env.example`)

### Optional API keys
I also use two API keys to get around low rate limits, namely NCBI and CyteType

# Build

To reproduce my results, paste the following commands into your terminal and run

```sh
git clone git@github.com:otodreas/scBaseCount_Pipeline.git
cd scBaseCount_Pipeline
uv sync
```