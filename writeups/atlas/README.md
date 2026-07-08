# Lung atlas integration

Phase 3 work merges the lung cohort in `output/metadata/datasets.csv` (772 accessions, ~5.98M cells) into one homogeneously processed h5ad, applies Harmony batch correction on study, and compares uncorrected vs corrected embeddings.

Implementation:

- Package: [`scripts/atlas_integration/`](../../scripts/atlas_integration/)
- Server runner: [`pipelines/run_atlas_integration.py`](../../pipelines/run_atlas_integration.py)
- Analysis notebook: [`notebooks/analysis/atlas_integration.ipynb`](../../notebooks/analysis/atlas_integration.ipynb)

## Curation path

The atlas is built from accessions that already passed the metadata cascade documented in [`scripts/atlas_integration/README.md`](../../scripts/atlas_integration/README.md). In short:

1. Filter scBaseCount sample metadata to lung disease AND lung tissue (`scripts/metadata/filter.py`).
2. Re-filter with tightened `NORMAL_HEALTHY_RE` (801 to 772 accessions).
3. Mirror raw h5ad files from GCS to R2 (`pipelines/migrate_gcs_to_r2.py`).
4. Attach `studyAccession` from `output/context/contexts.jsonl`.
5. Merge, integrate, and cluster (`scripts/atlas_integration/`).

## Why the batch key is `studyAccession`

Harmony needs a batch variable that captures technical variation shared across samples from the same publication, not every SRX/ERX accession individually.

Across the current atlas cohort:

| Metric | Value |
|--------|-------|
| Accessions | 772 |
| Unique studies (`studyAccession`) | 117 |
| Median accessions per study | 2 |
| Largest study | 80 accessions |

Most studies contribute more than one experiment accession. Treating each accession as its own batch would over-split batches and leave Harmony with little signal to integrate on.

### Example: one study, two experiments

The notebook proof of concept in [`notebooks/pipeline/study_context.ipynb`](../../notebooks/pipeline/study_context.ipynb) compares `ERX11662385` and `ERX11662370`. Both map to `PRJEB68315`, but they differ at the experiment level:

| Field | ERX11662385 | ERX11662370 |
|-------|-------------|-------------|
| `study.studyAccession` | PRJEB68315 | PRJEB68315 |
| `sampleAccession` | SAMEA114591111 | SAMEA114591096 |

This pattern repeats across the cohort. Among the 68 studies with two or more accessions in `datasets.csv`:

- all 68 have distinct `sampleAccession` values across experiments (expected),
- 9 have distinct `biological.tissueType` values within the same study,
- 3 have distinct `technical.libraryStrategy` values within the same study.

Using SRX/ERX as the batch key would treat these as unrelated batches even when they come from the same publication and likely share library prep and processing context. Using `studyAccession` groups them under one batch label, which matches the Phase 3 assumption that every study is a batch.

### Fallback behavior

When an accession is missing from `contexts.jsonl` or lacks a study block, merge falls back to using the accession itself as the batch key and logs a warning. After regenerating `contexts.jsonl` from the current `datasets.csv`, all 772 atlas accessions should resolve to a study.

## Missing `cell_type` labels

Author labels in `obs["cell_type"]` are kept for post-hoc comparison against CyteType, not for atlas clustering. Missing or blank values are filled with `"unknown"` at merge time. The CyteType confusion matrix excludes `"unknown"` rows.

Across the QC subset (116 accessions), only ~2.8% of cells lack a label; 3 accessions are fully unlabeled and 105 are partially labeled. The full 772-accession atlas follows the same rule: keep all cells, fill missing labels, exclude `"unknown"` only from the reference comparison.

## Phase 3 outputs

| Phase | Deliverable |
|-------|-------------|
| 3a | Merged atlas h5ad, UMAP with and without batch correction |
| 3b | Batch-mixing and cluster-conservation metrics, CyteType on atlas clusters, `cell_type` vs CyteType confusion matrix (excluding `unknown`) |

See the analysis notebook for interactive plots and metric tables loaded from `output/atlas/run_metadata.json`.
