# scripts

Python packages for the scBaseCount pipeline. Each package is installed into the project virtualenv via `uv sync` and can be imported directly in notebooks and scripts.

## Packages

| Package | Entry point | Output |
|---------|-------------|--------|
| [`metadata/`](metadata/) | `filter_lung(sample, cfg)` | `output/metadata/datasets.csv`, `accession_disease_categories.json`, `datasets_subset_qc.csv`, figures |
| [`study_context/`](study_context/) | `pipeline_for_accession_list(accessions)` | `CONTEXTS_JSONL_PATH` → `output/context/contexts.jsonl` |
| [`cluster_validation/`](cluster_validation/) | `run_cluster_validation(cfg)` | `output/clustering/data/{srx}_clustered.h5ad`, figures |
| [`cytetype_runner/`](cytetype_runner/) | `run_cytetype(cfg, ...)` | `output/cytetype/data/{srx}_cytetype_annotated.h5ad` (job details embedded in `adata.uns["cytetype_jobDetails"]`) |
| [`cyteonto/`](cyteonto/) | `run_cyteonto(cfg)` | `output/cyteonto/runs/{run_id}.csv` |
| [`h5ad_extractor/`](h5ad_extractor/) | `extract_annotation_columns(cfg)` | `output/h5ad_extract/{stem}_{obs\|var}_columns.{parquet\|csv}` |
| [`annotation_inspector/`](annotation_inspector/) | `inspect_accession(...)`, `write_extremes_csv(...)` | Pair-level `summary.csv`, optional `extremes.csv` |
| [`umap_plots/`](umap_plots/) | `plot_umap(adata, colorBy, ...)` | `output/umap_plots/figs/umap_{colorBy}.png` |
| [`gcs/`](gcs/) | `download_from_gcs(gs_uri, local_root)` | local mirror of GCS path under `data/` |
| [`r2/`](r2/) | `upload_to_r2(local_path, r2_key)` | — (side effect: uploads to R2) |
| [`shared/`](shared/) | `REPO_ROOT`, `configure_file_logger(...)` | — (utilities only) |

Each pipeline package has its own `README.md` with usage examples, config reference, and output model. `shared` is an internal utility package; it is not called directly from notebooks.

## Pipeline order

```
metadata  →  study_context  →  cluster_validation  →  cytetype_runner  →  cyteonto
                                                              |
                                                         gcs (download)
                                                         r2  (upload)
```

`metadata` produces the dataset catalog and accession list consumed by both `study_context` and `cluster_validation`. `study_context` produces the text context strings fed into CyteType for cluster annotation. `cytetype_runner` wraps the CyteType annotation step and persists the annotated h5ad. `gcs` and `r2` handle file transfer to and from cloud storage. `cyteonto` submits CyteType-annotated h5ad label columns to the CyteOnto API and returns ontology-aware similarity scores.
