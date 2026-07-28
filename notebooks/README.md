# notebooks

Analysis notebooks for the scBaseCount pipeline. Run them in order — each stage consumes outputs from the previous one.

## Run order

| # | Notebook | Consumes | Produces |
|---|----------|----------|----------|
| 1 | [`metadata.ipynb`](pipeline/metadata.ipynb) | scBaseCount metadata Parquet files | `output/metadata/datasets.csv`, `output/metadata/accession_disease_categories.json`, `output/metadata/datasets_subset_qc.csv`, figures |
| 2 | [`study_context.ipynb`](pipeline/study_context.ipynb) | `output/metadata/datasets.csv` | `output/context/contexts.jsonl` |
| 3 | [`clustering.ipynb`](pipeline/clustering.ipynb) | `output/metadata/datasets.csv`, h5ad files | `output/clustering/data/{srx}_clustered.h5ad`, figures |
| 4 | [`cytetype.ipynb`](pipeline/cytetype.ipynb) | `output/clustering/data/{srx}_clustered.h5ad`, `output/context/contexts.jsonl` | `output/cytetype/data/{srx}_cytetype_annotated.h5ad` |
| 5 | [`cyteonto.ipynb`](pipeline/cyteonto.ipynb) | `output/cytetype/data/{srx}_cytetype_annotated.h5ad` | `output/cyteonto/runs/{run_id}.csv` |

## Notebooks

**`pipeline/metadata.ipynb`** — Loads sample metadata from Parquet, applies lung tissue and disease filters, and exports the artifacts used by downstream notebooks: `datasets.csv` for clustering, `accession_disease_categories.json` mapping every lung-intersection accession to its `DISEASE_MAP` labels, and `datasets_subset_qc.csv` (a per-cohort QC-passing sample of up to 25 accessions across IPF / COVID-19 / COPD / Interstitial Lung Disease / Cystic Fibrosis) used as input to cytetype evaluation runs. Also produces three summary figures (sample breakdown, disease breakdown, cell count distribution). See `scripts/metadata/README.md` for the regex catalogue and labelling rules.

**`pipeline/study_context.ipynb`** — Fetches structured experiment context (study description, PubMed abstract, tissue type, library prep) from EBI ENA and NCBI for each accession in the dataset catalog. Results are cached to `output/context/contexts.jsonl` and can be reloaded without re-fetching. Includes field coverage, warnings, and distribution summaries.

**`analysis/clusters_to_cytetype_analysis.ipynb`** — Cost-planning notebook. Loads the clustering pipeline summary (`run.csv`) and joins it against both `accession_disease_categories.json` and `datasets_subset_qc.csv` to estimate how many CyteType clusters each disease cohort would consume. Produces inclusive (parent + child overlap counted) and disjoint (most-specific label) cluster-count tables.

**`pipeline/clustering.ipynb`** — Runs the cluster validation pipeline across all h5ad files in the scBaseCount data directory. For each sample: preprocesses, embeds, sweeps Leiden resolutions, selects the optimal resolution via Jaccard scoring, merges over-clustered partitions using a random forest, and writes the final annotated AnnData. Produces figures per sample under `output/clustering/figs/{srx}/`.

**`pipeline/cytetype.ipynb`** — Annotates cluster labels for a single sample using CyteType, driven by the study context string assembled from `output/context/contexts.jsonl`. Writes the annotated AnnData to `output/cytetype/data/{srx}_cytetype_annotated.h5ad`.

**`pipeline/cyteonto.ipynb`** — Submits a CyteType-annotated h5ad file to the CyteOnto API, polls until the run completes, and returns a DataFrame of ontology-aware similarity scores (one row per `(algorithm, cell)` pair). Results are saved to `output/cyteonto/runs/{run_id}.csv`. Runs can be safely interrupted and resumed via `check_pending_runs()`.

## Utility notebooks

Not part of the main pipeline. Use these to verify connectivity or inspect h5ad files independently.

**`utility/h5ad_extractor.ipynb`** — Pulls `obs` or `var` columns from an h5ad file into Parquet or CSV using the `h5ad_extractor` package (backed read, so the full file is not loaded into memory). Useful for inspecting annotation columns without running the full pipeline.

**`utility/gcs_test.ipynb`** — Verifies that `google-cloud-storage` can authenticate and reach the scBaseCount GCS bucket. Run this first if `gcloud auth application-default login` has not been set up yet.

**`utility/r2_inspect.ipynb`** — Verifies that `boto3` can authenticate and reach the Cloudflare R2 bucket using the credentials in `.env`. Lists objects and confirms the bucket is reachable before uploading processed files.

**`analysis/cluster_stats.ipynb`** — Loads `cluster_stats.json` produced by [`pipelines/cluster_stats.py`](../pipelines/cluster_stats.py) and builds an `xarray` tensor of cell-type by cluster counts for downstream analysis.

**`analysis/annotation_inspection.ipynb`** — Loads `summary.csv` and `extremes.csv` from [`pipelines/run_annotation_inspection_pipeline.py`](../pipelines/run_annotation_inspection_pipeline.py) for pair-level groupbys and extremes review.

**`analysis/analyze_atlas_DE.ipynb`** — Loads a Harmony atlas h5ad from [`pipelines/run_atlas_harmony.py`](../pipelines/run_atlas_harmony.py), uses `.raw` for full-gene counts, joins [`disease_markers`](../scripts/disease_markers/README.md) sample labels, builds SRX-level pseudobulks with Decoupler, and runs exploratory one-vs-rest PyDESeq2 per cluster.
