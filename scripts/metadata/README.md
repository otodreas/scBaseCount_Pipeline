# metadata

Filters the scBaseCount sample metadata catalog to a lung-specific subset and exports the artifacts that feed downstream pipeline stages: `datasets.csv` and `quantiles_datasets.csv` for clustering, `accession_disease_categories.json` for per-accession disease labels, and (driven from the notebook) a per-cohort `datasets_subset_qc.csv` used by the cytetype evaluation runs. Also produces three summary figures.

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

### Per-accession disease categories

```python
from metadata import disease_categories_for, export_accession_disease_categories

disease_categories_for("lung adenocarcinoma")
# ['Lung Cancer', 'Non-small Cell Lung Cancer (NSCLC)', 'Lung Adenocarcinoma (LUAD)']

disease_categories_for("non-cystic fibrosis (non-CF)")
# []  (the CF label is suppressed for any string matching NON_CF_RE)

export_accession_disease_categories(result.lungIntersection, cfg)
# writes output/metadata/accession_disease_categories.json
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
2. Drop samples where `disease` or `tissue` matches `NORMAL_HEALTHY_RE`. This regex matches the usual healthy / unknown vocabulary (`normal`, `healthy`, `control`, `unknown`, `not reported`, ...) and the negation patterns `no <disease>` (e.g. `no COPD`, `no diagnosed disease`, `no donor disease`) and `non[-\s]?(disease|COPD)` (e.g. `Non-disease`). `non-cystic fibrosis` and `non-CF` are deliberately excluded from this regex; those strings represent a real disease group and are handled separately at the labelling step.
3. From the remaining `sampleKnown` set, build:
   - `lungUnion`: disease OR tissue matches a lung-related pattern
   - `lungIntersection`: disease AND tissue both match (primary analysis set; excludes samples with lung disease but non-lung tissue labels such as blood or PBMC)
   - `lungIntersectionCancer`: intersection rows where disease matches `LUNG_CANCER_RE` (same definition `plot_disease_breakdown` uses for the lung-cancer subtree)

The intersection is used as the primary set because the union contains off-target tissue labels (blood, liver, PBMC) that are lung-disease-associated but not lung tissue.

## Disease labelling

`disease_categories_for(disease)` returns the ordered list of `DISEASE_MAP` labels matched by a disease string. The label set is a mix of a nested lung-cancer subtree and a flat set of cohort labels:

```
Lung Cancer
  -> Small Cell Lung Cancer (SCLC)
  -> Non-small Cell Lung Cancer (NSCLC)
       -> Lung Adenocarcinoma (LUAD)
       -> Lung Squamous Cell Carcinoma (LUSC)
       -> Lung Large Cell Carcinoma (LCC)

IPF / Pulmonary Fibrosis
COVID-19 / SARS-CoV-2
COPD
Cystic Fibrosis
Interstitial Lung Disease
Pulmonary Hypertension
```

`DISEASE_MAP` lists parents before children, so picking `cats[-1]` yields the most-specific lung-cancer subtype for a cancer row. For the sibling cohort labels the `cats[-1]` tie-break is just the later entry in `DISEASE_MAP`; this affects only the ~5 comorbid rows across the lung intersection (a few `lung cancer, COPD` and one `COVID-19, IPF`) and is documented inside `notebooks/metadata.ipynb`.

After regex matching, `disease_categories_for` drops the `Cystic Fibrosis` label from any row whose disease string matches `NON_CF_RE` (`non-CF` or `non-cystic fibrosis`). This catches both pure non-CF controls and mixed `cystic fibrosis (CF) and non-CF` datasets; both land in `Other`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `sampleParquetPath` | required | Path to sample-level metadata Parquet |
| `obsParquetPath` | required | Path to obs-level metadata Parquet |
| `minObsCount` | `1000` | Minimum cells per sample; samples below this are dropped before any other filtering |
| `outputDir` | `output/metadata` | Directory for `datasets.csv`, `quantiles_datasets.csv`, `accession_disease_categories.json`, and figures |

All default paths are relative to the repo root.

## Outputs

| File | Description |
|------|-------------|
| `outputDir/datasets.csv` | Full lung intersection: `srx_accession`, `file_path`, `obs_count` |
| `outputDir/quantiles_datasets.csv` | Five rows sampled at the 25th, 33rd, 50th, 67th, and 75th percentile of `obs_count`; used by `cluster_validation` as its dataset catalog |
| `outputDir/accession_disease_categories.json` | `{srx_accession: {disease, categories}}` for every row in the lung intersection. `categories` is the list returned by `disease_categories_for` |
| `outputDir/datasets_subset_qc.csv` | Per-cohort QC-passing sample of up to 25 accessions across IPF / COVID-19 / COPD / Interstitial Lung Disease / Cystic Fibrosis. Built inside `notebooks/metadata.ipynb` using `compute_obs_qc` + `apply_qc` |
| `outputDir/figs/sample_breakdown.png` | Pie chart: discarded / lung cancer / lung other / non-lung |
| `outputDir/figs/lung_disease_breakdown.png` | Pie chart: lung intersection by disease category (uses `disease_categories_for` + `cats[-1]`, so categories agree with the JSON) |
| `outputDir/figs/lung_cell_number_hist.png` | Log-scale histogram of cells per SRX |

## QC helpers

`compute_obs_qc(srxAccessions, cfg)` returns per-SRX `nCellsForQc`, `medianGenesPerCell`, and `medianUmisPerCell` by streaming the obs parquet with a pushdown filter. `apply_qc(samples, qc, QcThresholds(...))` joins the metrics onto a sample frame and returns the rows that pass the configured thresholds. The cytetype QC subset CSV is written from the notebook directly because it interleaves per-cohort sampling and column projection that does not generalize.

## Module reference

| Module | Public API |
|--------|------------|
| `config.py` | `MetadataConfig` |
| `load.py` | `load_sample(cfg)`, `obs_rows_for_srx(srx, cfg)` |
| `filter.py` | `filter_lung(sample, cfg)` -> `FilterResult`, `filter_by_disease`, `available_disease_labels` |
| `categorize.py` | `disease_categories_for(disease)`, `build_accession_disease_categories(samples)`, `export_accession_disease_categories(samples, cfg)` |
| `qc.py` | `QcThresholds`, `compute_obs_qc(srxAccessions, cfg)`, `apply_qc(samples, qc, thresholds)` |
| `export.py` | `export_datasets(result, cfg)` -> `(datasets_path, quantiles_path)` |
| `viz.py` | `plot_sample_breakdown`, `plot_disease_breakdown`, `plot_cell_count_distribution` |
| `regexes.py` | `NORMAL_HEALTHY_RE`, `LUNG_DISEASE_RE`, `LUNG_TISSUE_RE`, `LUNG_CANCER_RE`, `CANCER_RE`, `NON_CF_RE`, `DISEASE_MAP` |
