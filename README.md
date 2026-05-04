# scBaseCount Pipeline

Master's project, 30 hp, Lund University

Uses the Arc Institute's [Virtual Cell Atlas](https://console.cloud.google.com/storage/browser/arc-institute-virtual-cell-atlas?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))

## Usage

### Requirements

**Essential requirements**
- Python >= 3.12.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) (package manager)

**Optional**
- API keys
    - NCBI
    - CyteType

### Build

To reproduce my results, paste the following commands into your terminal and run

```sh
git clone git@github.com:otodreas/scBaseCount_Pipeline.git
cd scBaseCount_Pipeline
uv sync
```

## Google Cloud

I use the [Google Cloud SDK](https://docs.cloud.google.com/sdk/docs/install-sdk) to access the dataset. Some helpful variables

```sh
export META_PATH="gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens"  # metadata path
export H5AD_PATH="gs://arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens"      # h5ad path
```

Use `gcloud storage` followed by normal shell programs

```sh
gcloud storage ls "$META_PATH"                       # list files in the metadata directory
gcloud storage cp "$H5AD_PATH"/SRX123456.h5ad data/  # download an h5ad into data directory
```