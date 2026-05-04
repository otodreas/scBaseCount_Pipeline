# scripts

Python packages for the scBaseCount pipeline. Each package is installed into the project virtualenv via `uv sync` and can be imported directly in notebooks and scripts.

## Packages

| Package | Entry point | Output |
|---------|-------------|--------|
| [`metadata/`](metadata/) | `filter_lung(sample, cfg)` | `output/metadata/datasets.csv`, `quantiles_datasets.csv`, figures |
| [`study_context/`](study_context/) | `pipeline_for_accession_list(accessions)` | `output/contexts.jsonl` |
| [`cluster_validation/`](cluster_validation/) | `run_cluster_validation(cfg)` | `output/clustering/data/{srx}_clustered.h5ad`, figures |
| [`cyteonto/`](cyteonto/) | `run_cyteonto(cfg)` | `output/cyteonto/runs/{run_id}.csv` |
| [`h5ad_extractor/`](h5ad_extractor/) | `extract_annotation_columns(cfg)` | `output/h5ad_extract/{stem}_{obs\|var}_columns.{parquet\|csv}` |
| [`shared/`](shared/) | `REPO_ROOT`, `configure_file_logger(...)` | — (utilities only) |

Each pipeline package has its own `README.md` with usage examples, config reference, and output model. `shared` is an internal utility package used by the other packages; it is not called directly from notebooks.

## Pipeline order

```
metadata  →  study_context  →  cluster_validation  →  cytetype  →  cyteonto
```

`metadata` produces the dataset catalog and accession list that both `study_context` and `cluster_validation` consume. `study_context` produces the text context strings fed into CyteType for cluster annotation. `cyteonto` submits CyteType-annotated h5ad files to the CyteOnto API and returns ontology-aware similarity scores.
