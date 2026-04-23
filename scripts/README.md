# scripts

Python packages for the scBaseCount pipeline. Each package is installed into the project virtualenv via `uv sync` and can be imported directly in notebooks and scripts.

## Packages

| Package | Entry point | Output |
|---------|-------------|--------|
| [`metadata/`](metadata/) | `filter_lung(sample, cfg)` | `output/metadata/datasets.csv`, `quantiles_datasets.csv`, figures |
| [`study_context/`](study_context/) | `pipeline_for_accession_list(accessions)` | `output/contexts.jsonl` |
| [`cluster_validation/`](cluster_validation/) | `run_cluster_validation(cfg)` | `output/clustering/data/{srx}_clustered.h5ad`, figures |

Each package has its own `README.md` with usage examples, config reference, and output model.

## Pipeline order

```
metadata  →  study_context  →  cluster_validation  →  cytetype
```

`metadata` produces the dataset catalog and accession list that both `study_context` and `cluster_validation` consume. `study_context` produces the text context strings fed into CyteType for cluster annotation.
