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

Study context is produced by [`notebooks/study_context.ipynb`](../notebooks/study_context.ipynb). The cytetype runner loads summaries via `experiment_context_summary`; missing accessions proceed with empty context.

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

Copies raw scBaseCount h5ad files from GCS to R2. Skips objects already in R2 with a matching MD5. Uses a local cache under `data/` when the file is not already on disk.

**Output:** `output/migration/{timestamp}/run.csv`

| Flag | Default | Description |
|------|---------|-------------|
| `--datasets` | `output/metadata/datasets.csv` | Accession list |
| `--dry-run` | off | Log planned uploads only |

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
| `--poll-interval-s` | `10` | Seconds between CyteOnto status polls |
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

Downstream analysis: [`notebooks/cluster_stats.ipynb`](../notebooks/cluster_stats.ipynb).

```sh
uv run python pipelines/cluster_stats.py \
  --r2-prefix clustering_pipeline_20260511_140000 \
  --workers 4
```
