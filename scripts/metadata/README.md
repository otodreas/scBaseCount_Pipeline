# metadata

Filters the scBaseCount sample metadata catalog to a lung-specific subset and exports the two CSV files that feed downstream pipeline stages. Loads sample and obs metadata from Parquet, applies regex-based disease and tissue filters, produces three summary figures, and writes `datasets.csv` (full lung intersection) and `quantiles_datasets.csv` (five representative samples used by `cluster_validation`).

## Usage

```python
from metadata import MetadataConfig, load_sample, filter_lung, export_datasets
from metadata.viz import plot_sample_breakdown, plot_disease_breakdown, plot_cell_count_distribution
from pathlib import Path

cfg = MetadataConfig(
    sampleParquetPath=Path("data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_sample_metadata.parquet"),
    obsParquetPath=Path("data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_obs_metadata.parquet"),
    outputDir=Path("output/metadata"),
)

sample = load_sample(cfg)
result = filter_lung(sample, cfg)

plot_sample_breakdown(sample, result, figs_dir=cfg.outputDir / "figs")
plot_disease_breakdown(result, figs_dir=cfg.outputDir / "figs")
plot_cell_count_distribution(result, figs_dir=cfg.outputDir / "figs")

datasets_path, quantiles_path = export_datasets(result, cfg)
```

### Accessing results

```python
result.lungIntersection          # primary analysis DataFrame (lung disease AND tissue)
result.lungUnion                 # broader set (lung disease OR tissue)
result.lungIntersectionCancer    # cancer subset of the intersection
result.sampleKnown               # all samples after dropping healthy/unknown/< minObsCount
```

### Querying obs-level data for a single sample

```python
from metadata import obs_rows_for_srx

obs = obs_rows_for_srx("SRX17412841", cfg)
obs["cell_type"].value_counts()
```

The filter is pushed into the Parquet scan, so no full-file load is needed.

## Filter logic

Filtering is applied in three steps:

1. Drop samples with fewer than `minObsCount` cells.
2. Drop samples where `disease` or `tissue` matches `NORMAL_HEALTHY_RE` (healthy controls, unknowns, unsures).
3. From the remaining `sampleKnown` set, build:
   - `lungUnion` — disease OR tissue matches a lung-related pattern
   - `lungIntersection` — disease AND tissue both match (primary analysis set; excludes samples with lung disease but non-lung tissue labels such as blood or PBMC)
   - `lungIntersectionCancer` — intersection rows where disease matches `CANCER_RE`

The intersection is used as the primary set because the union contains off-target tissue labels (blood, liver, PBMC) that are lung-disease-associated but not lung tissue.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `sampleParquetPath` | required | Path to sample-level metadata Parquet |
| `obsParquetPath` | required | Path to obs-level metadata Parquet |
| `minObsCount` | `1000` | Minimum cells per sample; samples below this are dropped before any other filtering |
| `outputDir` | `output/metadata` | Directory for `datasets.csv`, `quantiles_datasets.csv`, and figures |

## Outputs

| File | Description |
|------|-------------|
| `outputDir/datasets.csv` | Full lung intersection: `srx_accession`, `file_path`, `obs_count` |
| `outputDir/quantiles_datasets.csv` | Five rows sampled at the 25th, 33rd, 50th, 67th, and 75th percentile of `obs_count`; used by `cluster_validation` as its dataset catalog |
| `outputDir/figs/sample_breakdown.svg` | Pie chart: discarded / lung cancer / lung other / non-lung |
| `outputDir/figs/lung_disease_breakdown.svg` | Pie chart: lung intersection by disease category |
| `outputDir/figs/lung_cell_number_hist.svg` | Log-scale histogram of cells per SRX |

## Module reference

| Module | Public API |
|--------|------------|
| `config.py` | `MetadataConfig` |
| `load.py` | `load_sample(cfg)`, `obs_rows_for_srx(srx, cfg)` |
| `filter.py` | `filter_lung(sample, cfg)` → `FilterResult` |
| `export.py` | `export_datasets(result, cfg)` → `(datasets_path, quantiles_path)` |
| `viz.py` | `plot_sample_breakdown`, `plot_disease_breakdown`, `plot_cell_count_distribution` |
| `regexes.py` | `NORMAL_HEALTHY_RE`, `LUNG_DISEASE_RE`, `LUNG_TISSUE_RE`, `CANCER_RE`, `DISEASE_MAP` |
