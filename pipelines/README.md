# pipelines

Batch runners for server-side work across many accessions. Each script reads `output/metadata/datasets.csv` (unless noted), writes a timestamped run directory under `output/`, and logs to `logs/`.

Run from the repo root after `uv sync --group dev` and `.env` is configured (see [`.env.example`](../.env.example)).

```sh
uv run python pipelines/run_clustering_pipeline.py --help
```

Package logic lives under [`scripts/`](../scripts/README.md). Interactive exploration uses [`notebooks/`](../notebooks/README.md).

## Typical order

```
migrate_gcs_to_r2  →  run_clustering_pipeline  →  run_cytetype_pipeline  →  run_cyteonto_pipeline  →  cluster_stats
        |                      |                          |                          |
   raw h5ad in R2      clustered h5ad in R2        annotated h5ad in R2        cyteonto CSV in R2
```

`migrate_gcs_to_r2` is optional when raw files are already mirrored to R2. Notebooks can replace any batch step for a single accession.

## Shared inputs

| File | Used by |
|------|---------|
| `output/metadata/datasets.csv` | All runners (`srx_accession`, `file_path`, …) |
| `output/context/contexts.jsonl` (`study_context.CONTEXTS_JSONL_PATH`) | `run_cytetype_pipeline.py` (default; overridable with `--contexts`) |

Study context is produced by [`notebooks/pipeline/study_context.ipynb`](../notebooks/pipeline/study_context.ipynb). The cytetype runner loads summaries via `experiment_context_summary`; missing accessions proceed with empty context.

## R2 key layout

| Stage | Key pattern | Example |
|-------|-------------|---------|
| Raw mirror | `gcs_uri_to_r2_raw_key(gs_uri)` from `file_path` | Mirrors GCS path under the bucket |
| Clustering run | `{r2_prefix}/{srx}_clustered.h5ad` | `clustering_pipeline_20260511_140000/SRX…_clustered.h5ad` |
| CyteType run | `{r2_prefix}/{srx}_annotated.h5ad` | `cytetype_pipeline_20260512_090000/SRX…_annotated.h5ad` |
| CyteOnto run | `{r2_prefix}/{srx}_cyteonto.csv` | `cyteonto_pipeline_20260513_090000/SRX…_cyteonto.csv` |

`--r2-prefix` defaults to `{script_name}_{YYYYMMDD_HHMMSS}` at import time. Pass an explicit prefix when chaining runs (cytetype needs `--clustering-prefix` from the clustering run’s `metadata.json`).

---

## `migrate_gcs_to_r2.py`

Copies raw scBaseCount h5ad files from GCS to R2 for accessions present in the source CSV but absent from the baseline CSV. Skips objects already in R2 with a matching MD5. Uses a local cache under `data/` when the file is not already on disk.

**Output:** `output/migration/{timestamp}/run.csv`

| Flag | Default | Description |
|------|---------|-------------|
| `--datasets` | `output/metadata/datasets_v2.csv` | Source accession list |
| `--baseline` | `output/metadata/datasets.csv` | Accessions to exclude from migration |
| `--dry-run` | off | Log planned uploads only |

With the current metadata files, the default selection is 1,048 accessions. The process exits with status 1 if any selected row fails.

**Log:** `logs/migrate_gcs_to_r2.log`

```sh
uv run python pipelines/migrate_gcs_to_r2.py --dry-run
uv run python pipelines/migrate_gcs_to_r2.py
```

---

## `run_clustering_pipeline.py`

Runs [`cluster_validation`](../scripts/cluster_validation/README.md) per accession: download raw h5ad (local cache, else R2 raw mirror, else GCS), cluster, upload `{srx}_clustered.h5ad` to R2. Skips accessions whose output key already exists.

Per-run artifacts also land under the run directory (`figs/`, `data/`) before upload; local clustered files are deleted after a successful upload.

**Output:** `output/clustering_pipeline/{timestamp}/`

| File | Description |
|------|-------------|
| `run.csv` | Per-accession status and clustering metrics |
| `metadata.json` | Run config snapshot (`r2_prefix`, paths, optional `notes`) |
| `figs/`, `data/` | Copies of figures and clustered h5ad for this batch run |

| Flag | Default | Description |
|------|---------|-------------|
| `--datasets` | `output/metadata/datasets.csv` | Accession list |
| `--r2-prefix` | `clustering_pipeline_{timestamp}` | R2 folder for clustered outputs |
| `--metadata` | none | Free-form note stored in `metadata.json` |
| `--workers` | `1` | Parallel accession workers |

**Log:** `logs/clustering_pipeline.log`

```sh
uv run python pipelines/run_clustering_pipeline.py --workers 4
```

---

## `run_cytetype_pipeline.py`

Runs [`cytetype_runner`](../scripts/cytetype_runner/README.md) on clustered h5ad files already in R2. Requires `CYTETYPE_API_KEY`. Accessions run **serially**; optional spacing between run starts via `--min-interval`.

Downloads each input from `{clustering_prefix}/{srx}_clustered.h5ad`, annotates with study context from contexts JSONL, uploads `{r2_prefix}/{srx}_annotated.h5ad`. Per-accession metadata sent to CyteType is derived from the datasets CSV row (all columns, stringified).

**Output:** `output/cytetype_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)

| File | Description |
|------|-------------|
| `run.csv` | Per-accession status, timing, R2 keys |
| `job_details.csv` | CyteType `job_id`, `report_url`, `api_url` |
| `metadata.json` | Run config snapshot |

| Flag | Default | Description |
|------|---------|-------------|
| `--datasets` | `output/metadata/datasets.csv` | Accession list |
| `--contexts` | `output/context/contexts.jsonl` | Study context cache |
| `--clustering-prefix` | required | R2 prefix from a clustering pipeline run |
| `--r2-prefix` | `cytetype_pipeline_{timestamp}` | R2 folder for annotated outputs |
| `--metadata` | none | Run-level note in `metadata.json` only (not sent to CyteType) |
| `--min-interval` | `0` | Minimum seconds between starting consecutive accessions |
| `--dry-run` | off | Write plan CSVs without R2 or API calls |

**Log:** `logs/cytetype_pipeline.log`

```sh
uv run python pipelines/run_cytetype_pipeline.py \
  --clustering-prefix clustering_pipeline_20260511_140000 \
  --dry-run
```

---

## `run_cyteonto_pipeline.py`

Runs [`cyteonto`](../scripts/cyteonto/README.md) on annotated h5ad files already in R2. Lists every `{srx}_annotated.h5ad` under `--input-prefix` and runs accessions **serially**, blocking on each CyteOnto poll loop before starting the next.

Downloads each input from `{input_prefix}/{srx}_annotated.h5ad`, submits to the CyteOnto API, moves the result CSV to `results/{srx}_cyteonto.csv`, and uploads `{r2_prefix}/{srx}_cyteonto.csv`. Skips accessions whose output key already exists in R2.

**Output:** `output/cyteonto_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)

| File | Description |
|------|-------------|
| `run.csv` | Per-accession status, timing, R2 keys, `run_id`, local CSV path |
| `results/{srx}_cyteonto.csv` | Per-accession CyteOnto similarity CSV |
| `metadata.json` | Run config snapshot |

| Flag | Default | Description |
|------|---------|-------------|
| `--input-prefix` | required | R2 prefix containing annotated h5ads (e.g. from a cytetype pipeline run) |
| `--r2-prefix` | `cyteonto_pipeline_{timestamp}` | R2 folder for CyteOnto result CSVs |
| `--metadata` | none | Free-form note in `metadata.json` |
| `--poll-interval-s` | `10` | Seconds between CyteOnto result polls |
| `--poll-timeout-s` | `3600` | Seconds before a CyteOnto run raises `TimeoutError` |
| `--min-interval` | `0` | Minimum seconds between starting consecutive accessions |
| `--dry-run` | off | Write plan CSVs without R2 or API calls |

**Log:** `logs/cyteonto_pipeline.log`

For an in-flight accession, look up the `run_id` in `logs/cyteonto_pipeline.log` or the pending stub under `output/cyteonto/runs/`.

```sh
nohup uv run python pipelines/run_cyteonto_pipeline.py \
  --input-prefix cytetype_pipeline_20260522_175813 \
  > logs/cyteonto_pipeline.nohup.out 2>&1 &
```

---

## `cluster_stats.py`

Downloads clustered h5ad files from an R2 prefix, builds per-SRX `cell_type` × `leiden_merged` count matrices, and aggregates normalized Shannon entropy (NSE) and KL divergence (KLD) via [`cluster_validation.cell_type_metrics`](../scripts/cluster_validation/README.md).

**Output:** `output/cluster_stats/{r2_prefix}/` (default; override with `--output-dir`)

| File | Description |
|------|-------------|
| `run.csv` | Per-accession status |
| `metadata.json` | Run config snapshot |
| `cluster_stats.json` | Nested counts: `{srx: {cell_type: {cluster: count}}}` |
| `nse_matrix.csv`, `kld_matrix.csv` | Accessions × cell types |
| `cell_type_summary.csv` | Cell-type level means |
| `cell_type_metrics.png` | Summary bar chart |

| Flag | Default | Description |
|------|---------|-------------|
| `--r2-prefix` | required | Clustering run prefix (e.g. `clustering_pipeline_20260511_140000`) |
| `--output-dir` | `output/cluster_stats/{r2_prefix}` | Output root |
| `--metadata` | none | Note in `metadata.json` |
| `--workers` | `1` | Parallel accession workers |

**Log:** `logs/cluster_stats.log`

Downstream analysis: [`notebooks/analysis/cluster_stats.ipynb`](../notebooks/analysis/cluster_stats.ipynb).

```sh
uv run python pipelines/cluster_stats.py \
  --r2-prefix clustering_pipeline_20260511_140000 \
  --workers 4
```

---

## `run_annotation_inspection_pipeline.py`

Inspects CyteType-annotated h5ads from R2, joins CyteOnto cytescores, and writes a pair-level summary plus optional extremes table via [`annotation_inspector`](../scripts/annotation_inspector/README.md).

Downloads `{input_prefix}/{srx}_annotated.h5ad` and `{cyteonto_prefix}/{srx}_cyteonto.csv` (when present), inspects each accession, deletes local cache files, and streams outputs as workers finish.

**Output:** `output/annotation_inspection_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)

| File | Description |
|------|-------------|
| `summary.csv` | Pair-level rows: labels, confidence, cytescore, n_cells, report_url |
| `extremes.csv` | Top/bottom cytescore STATE types per CyteType label (when extremes enabled) |
| `run.csv` | Per-accession status and timing |
| `metadata.json` | Run config snapshot |

| Flag | Default | Description |
|------|---------|-------------|
| `--input-prefix` | required | R2 prefix containing annotated h5ads |
| `--cyteonto-prefix` | required | R2 prefix containing `{srx}_cyteonto.csv` files |
| `--workers` | `1` | Concurrent R2 fetch + inspect workers |
| `--top-n` | `10` | Top/bottom STATE cell types per CyteType label for extremes |
| `--no-extremes` | off | Skip `extremes.csv` (written by default) |
| `--from-summary` | none | Rebuild `extremes.csv` from an existing `summary.csv` (no R2 fetch) |
| `--output-dir` | summary parent | Output dir for `--from-summary` |
| `--metadata` | none | Note in `metadata.json` |
| `--dry-run` | off | Write plan CSV without R2 or inspection |

**Log:** `logs/annotation_inspection_pipeline.log`

Downstream analysis: [`notebooks/analysis/annotation_inspection.ipynb`](../notebooks/analysis/annotation_inspection.ipynb).

```sh
uv run python pipelines/run_annotation_inspection_pipeline.py \
  --input-prefix cytetype_pipeline_20260522_175813 \
  --cyteonto-prefix cyteonto_pipeline_20260526_112224 \
  --workers 4

# Regenerate extremes only
uv run python pipelines/run_annotation_inspection_pipeline.py \
  --from-summary output/annotation_inspection_pipeline/20260603_120000/summary.csv
```

---

## `run_atlas_concat.py`

Concatenates raw scBaseCount h5ads into one merged atlas via [`h5ad_concat`](../scripts/h5ad_concat/README.md): download each file from R2, validate and QC it, align to the reference gene axis, and merge the passing files. Unlike the other runners it has no CLI flags; edit the `H5adConcatConfig` in the script (for example `atlasR2Key`, and `uploadAtlas=True` to push the atlas to R2).

Reads `output/metadata/datasets.csv` by default (config `datasetsPath`).

**Output:** `output/atlas/data/` (`atlas.h5ad`, `atlas_config.json`, `atlas.csv`, `atlas_result.json`)

**Log:** `logs/h5ad_concat.log`

```sh
uv run python pipelines/run_atlas_concat.py
```

---

## `run_atlas_harmony.py`

Batch-corrects a merged atlas h5ad with Harmony on `study_accession`, then compares pre- and post-correction embeddings. Loads the atlas, normalizes and log-transforms counts, selects batch-aware HVGs, runs PCA, builds uncorrected and Harmony-corrected neighbor graphs and UMAPs, and clusters with Leiden. UMAP PNGs use [`umap_plots`](../scripts/umap_plots/README.md).

Typical input: an atlas written by [`run_atlas_concat.py`](#run_atlas_concatpy) (for example `output/atlas/data/atlas_sample20.h5ad`).

**Output:** paths below use CLI defaults; override with `--input`, `--output`, and `--figs-dir`.

| File | Description |
|------|-------------|
| `{output}.h5ad` | Processed atlas (HVGs in `.X`, full-gene counts in `.raw`) with `X_umap_uncorrected`, `X_pca_harmony`, `leiden_uncorrected`, and `leiden_atlas` |
| `{output_stem}_run.json` | Run summary (cell counts, HVGs, raw gene count, studies, cluster counts, config) |
| `{figs_dir}/umap_{batch_key}_uncorrected.png` | Pre-correction UMAP colored by batch |
| `{figs_dir}/umap_{cell_type_key}_uncorrected.png` | Pre-correction UMAP colored by cell type |
| `{figs_dir}/umap_{batch_key}_harmony.png` | Harmony-corrected UMAP colored by batch |
| `{figs_dir}/umap_{cell_type_key}_harmony.png` | Harmony-corrected UMAP colored by cell type |
| `{figs_dir}/pca_scree.png` | PCA scree plot (per-PC and cumulative variance) |

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | `output/atlas/data/atlas_sample20.h5ad` | Input atlas h5ad |
| `--output` | `output/atlas/data/atlas_sample20_harmony.h5ad` | Output atlas h5ad |
| `--figs-dir` | `output/atlas/figs` | Directory for UMAP and scree PNGs |
| `--batch-key` | `study_accession` | `obs` column for Harmony batch correction |
| `--cell-type-key` | `cell_type` | `obs` column for cell-type UMAP plots |
| `--n-top-genes` | `2000` | Number of HVGs |
| `--n-pcs` | `30` | PCs used for the neighbor graph |
| `--n-pcs-compute` | `50` | PCs computed by PCA |
| `--resolution` | `1.0` | Leiden resolution |
| `--no-plots` | off | Skip writing PNGs (still saves the h5ad and run JSON) |
| `--threads` | `0` | scanpy `n_jobs` (`0` leaves the default) |

**Log:** `logs/atlas_harmony.log`

```sh
uv run python pipelines/run_atlas_harmony.py

uv run python pipelines/run_atlas_harmony.py \
  --input output/atlas/data/atlas.h5ad \
  --output output/atlas/data/atlas_harmony.h5ad
```

Downstream atlas DE and disease-area labeling: [`notebooks/analysis/analyze_atlas_DE.ipynb`](../notebooks/analysis/analyze_atlas_DE.ipynb) loads the Harmony h5ad, uses `.raw` for full-gene pseudobulk counts, and joins labels from [`disease_markers`](../scripts/disease_markers/README.md).
