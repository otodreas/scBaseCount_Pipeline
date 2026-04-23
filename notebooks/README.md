# notebooks

Analysis notebooks for the scBaseCount pipeline. Run them in order — each stage consumes outputs from the previous one.

## Run order

| # | Notebook | Consumes | Produces |
|---|----------|----------|----------|
| 1 | [`metadata.ipynb`](metadata.ipynb) | scBaseCount metadata Parquet files | `output/metadata/datasets.csv`, `output/metadata/quantiles_datasets.csv`, figures |
| 2 | [`study_context.ipynb`](study_context.ipynb) | `output/metadata/datasets.csv` | `output/contexts.jsonl` |
| 3 | [`clustering.ipynb`](clustering.ipynb) | `output/metadata/quantiles_datasets.csv`, h5ad files | `output/clustering/data/{srx}_clustered.h5ad`, figures |
| 4 | [`cytetype.ipynb`](cytetype.ipynb) | `output/clustering/data/{srx}_clustered.h5ad`, `output/contexts.jsonl` | `data/other/{srx}_cytetype_annotated.h5ad` |

## Notebooks

**`metadata.ipynb`** — Loads sample metadata from Parquet, applies lung tissue and disease filters, and exports the dataset catalogs used by all downstream notebooks. Also produces three summary figures (sample breakdown, disease breakdown, cell count distribution).

**`study_context.ipynb`** — Fetches structured experiment context (study description, PubMed abstract, tissue type, library prep) from EBI ENA and NCBI for each accession in the dataset catalog. Results are cached to `output/contexts.jsonl` and can be reloaded without re-fetching. Includes field coverage, warnings, and distribution summaries.

**`clustering.ipynb`** — Runs the cluster validation pipeline across all h5ad files in the scBaseCount data directory. For each sample: preprocesses, embeds, sweeps Leiden resolutions, selects the optimal resolution via Jaccard scoring, merges over-clustered partitions using a random forest, and writes the final annotated AnnData. Produces figures per sample under `output/clustering/figs/{srx}/`.

**`cytetype.ipynb`** — Annotates cluster labels for a single sample using CyteType, driven by the study context string assembled from `output/contexts.jsonl`. Writes the annotated AnnData to `data/other/{srx}_cytetype_annotated.h5ad`.
