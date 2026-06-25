# metadata

Filters the scBaseCount sample metadata catalog to a lung-specific subset and exports the artifacts that feed downstream pipeline stages: `datasets.csv` for clustering, `accession_disease_categories.json` for per-accession disease labels, and (driven from the notebook) a per-cohort `datasets_subset_qc.csv` used by the cytetype evaluation runs. Also produces three summary figures.

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

datasets_path = export_datasets(result, cfg)
```

### Accessing results

```python
result.lungIntersection          # primary analysis DataFrame (lung disease AND tissue)
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

### Querying metadata for a single accession

```python
from metadata import sample_row_for_srx, obs_rows_for_srx

row = sample_row_for_srx("SRX17412841", cfg)
row["disease"], row["tissue"], row["file_path"]

obs = obs_rows_for_srx("SRX17412841", cfg)
obs["cell_type"].value_counts()
```

Both lookups push the SRX filter into the Parquet scan, so no full-file load is needed.

## Filter logic

Filtering is applied in three steps:

1. Drop samples with fewer than `minObsCount` cells.
2. Drop samples where `disease` or `tissue` matches `NORMAL_HEALTHY_RE`. This regex matches the usual healthy / unknown vocabulary (`normal`, `healthy`, `control`, `unknown`, `not reported`, ...) and the negation patterns `no <disease>` (e.g. `no COPD`, `no diagnosed disease`, `no donor disease`) and `non[-\s]?(disease|COPD)` (e.g. `Non-disease`). `non-cystic fibrosis` and `non-CF` are deliberately excluded from this regex; those strings represent a real disease group and are handled separately at the labelling step.
3. From the remaining `sampleKnown` set, build:
   - `lungIntersection`: disease AND tissue both match (primary analysis set; excludes samples with lung disease but non-lung tissue labels such as blood or PBMC)
   - `lungIntersectionCancer`: intersection rows where disease matches `LUNG_CANCER_RE` (same definition `plot_disease_breakdown` uses for the lung-cancer subtree)

## Disease labelling

`disease_categories_for(disease)` returns the ordered list of `DISEASE_MAP` labels matched by a disease string. `most_specific_disease_label(disease)` is a thin wrapper that picks the last entry of that list (the most-specific match) and returns `"Other"` when nothing matched. Both the cohort tagging step in `notebooks/pipeline/metadata.ipynb` and the disjoint partition in `notebooks/analysis/clusters_to_cytetype_analysis.ipynb` go through this helper, so the per-accession bucket assignment stays consistent across notebooks and figures. The label set is a mix of a nested lung-cancer subtree and a flat set of cohort labels:

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

`DISEASE_MAP` lists parents before children, so `most_specific_disease_label` yields the most-specific lung-cancer subtype for a cancer row. For the sibling cohort labels the tie-break is just the later entry in `DISEASE_MAP`; this affects only the ~5 comorbid rows across the lung intersection (a few `lung cancer, COPD` and one `COVID-19, IPF`) and is documented inside `notebooks/pipeline/metadata.ipynb`.

After regex matching, `disease_categories_for` drops the `Cystic Fibrosis` label from any row whose disease string matches `NON_CF_RE` (`non-CF` or `non-cystic fibrosis`). This catches both pure non-CF controls and mixed `cystic fibrosis (CF) and non-CF` datasets; both land in `Other`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `sampleParquetPath` | `data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/..._sample_metadata.parquet` | Path to sample-level metadata Parquet |
| `obsParquetPath` | `data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/..._obs_metadata.parquet` | Path to obs-level metadata Parquet |
| `minObsCount` | `1000` | Minimum cells per sample; samples below this are dropped before any other filtering |
| `outputDir` | `output/metadata` | Directory for `datasets.csv`, `accession_disease_categories.json`, `datasets_subset_qc.csv`, and figures |

All default paths are relative to the repo root.

## Outputs

| File | Description |
|------|-------------|
| `outputDir/datasets.csv` | Full lung intersection: `srx_accession`, `file_path`, `obs_count`. Consumed by `clustering.ipynb` |
| `outputDir/accession_disease_categories.json` | `{srx_accession: {disease, categories}}` for every row in the lung intersection. `categories` is the list returned by `disease_categories_for` |
| `outputDir/datasets_subset_qc.csv` | Per-cohort QC-passing sample of up to 25 accessions across IPF / COVID-19 / COPD / Interstitial Lung Disease / Cystic Fibrosis. Built inside `notebooks/pipeline/metadata.ipynb` using `compute_obs_qc` + `apply_qc` |
| `outputDir/figs/sample_breakdown.png` | Pie chart: discarded / lung cancer / lung other / non-lung |
| `outputDir/figs/lung_disease_breakdown.png` | Pie chart: lung intersection by disease category (uses `most_specific_disease_label`, so categories agree with the JSON) |
| `outputDir/figs/lung_cell_number_hist.png` | Log-scale histogram of cells per SRX |

## QC helpers

`compute_obs_qc(srxAccessions, cfg)` returns per-SRX `nCellsForQc`, `medianGenesPerCell`, and `medianUmisPerCell` by streaming the obs parquet with a pushdown filter. `apply_qc(samples, qc, QcThresholds(...))` joins the metrics onto a sample frame and returns the rows that pass the configured thresholds. The cytetype QC subset CSV is written from the notebook directly because it interleaves per-cohort sampling and column projection that does not generalize.

## Module reference

| Module | Public API |
|--------|------------|
| `config.py` | `MetadataConfig` |
| `load.py` | `load_sample(cfg)`, `sample_row_for_srx(srx, cfg)`, `obs_rows_for_srx(srx, cfg)` |
| `filter.py` | `filter_lung(sample, cfg)` -> `FilterResult` |
| `categorize.py` | `disease_categories_for(disease)`, `most_specific_disease_label(disease)`, `export_accession_disease_categories(samples, cfg)` |
| `qc.py` | `QcThresholds`, `compute_obs_qc(srxAccessions, cfg)`, `apply_qc(samples, qc, thresholds)` |
| `export.py` | `export_datasets(result, cfg)` -> `datasets_path` |
| `viz.py` | `plot_sample_breakdown`, `plot_disease_breakdown`, `plot_cell_count_distribution` |
| `regexes.py` | `NORMAL_HEALTHY_RE`, `LUNG_DISEASE_RE`, `LUNG_TISSUE_RE`, `LUNG_CANCER_RE`, `NON_CF_RE`, `DISEASE_MAP` |
