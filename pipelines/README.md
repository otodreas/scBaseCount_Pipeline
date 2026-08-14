# pipelines

Batch runners for server-side work across many accessions. Each script reads `output/metadata/datasets.csv` (unless noted), writes a timestamped run directory under `output/`, and logs to `logs/`.

Run from the repo root after `uv sync --group dev` and `.env` is configured (see `[.env.example](../.env.example)`).

```sh
uv run python pipelines/run_clustering_pipeline.py --help
```

Package logic lives under `[scripts/](../scripts/README.md)`. Interactive exploration uses `[notebooks/](../notebooks/README.md)`.

## Typical order

```
migrate_gcs_to_r2  →  run_clustering_pipeline  →  run_cytetype_pipeline  →  run_cyteonto_pipeline  →  cluster_stats
        |                      |                          |                          |
   raw h5ad in R2      clustered h5ad in R2        annotated h5ad in R2        cyteonto CSV in R2
```

`migrate_gcs_to_r2` is optional when raw files are already mirrored to R2. Notebooks can replace any batch step for a single accession.

## Shared inputs


| File                                                                  | Used by                                                             |
| --------------------------------------------------------------------- | ------------------------------------------------------------------- |
| `output/metadata/datasets.csv`                                        | All runners (`srx_accession`, `file_path`, …)                       |
| `output/context/contexts.jsonl` (`study_context.CONTEXTS_JSONL_PATH`) | `run_cytetype_pipeline.py` (default; overridable with `--contexts`) |


Study context is produced by `[notebooks/pipeline/study_context.ipynb](../notebooks/pipeline/study_context.ipynb)`. The cytetype runner loads summaries via `experiment_context_summary`; missing accessions proceed with empty context.

## R2 key layout


| Stage          | Key pattern                                      | Example                                                   |
| -------------- | ------------------------------------------------ | --------------------------------------------------------- |
| Raw mirror     | `gcs_uri_to_r2_raw_key(gs_uri)` from `file_path` | Mirrors GCS path under the bucket                         |
| Clustering run | `{r2_prefix}/{srx}_clustered.h5ad`               | `clustering_pipeline_20260511_140000/SRX…_clustered.h5ad` |
| CyteType run   | `{r2_prefix}/{srx}_annotated.h5ad`               | `cytetype_pipeline_20260512_090000/SRX…_annotated.h5ad`   |
| CyteOnto run   | `{r2_prefix}/{srx}_cyteonto.csv`                 | `cyteonto_pipeline_20260513_090000/SRX…_cyteonto.csv`     |


`--r2-prefix` defaults to `{script_name}_{YYYYMMDD_HHMMSS}` at import time. Pass an explicit prefix when chaining runs (cytetype needs `--clustering-prefix` from the clustering run’s `metadata.json`).

---



## `migrate_gcs_to_r2.py`

Copies raw scBaseCount h5ad files from GCS to R2 for accessions present in the source CSV but absent from the baseline CSV. Skips objects already in R2 with a matching MD5. Uses a local cache under `data/` when the file is not already on disk.

**Output:** `output/migration/{timestamp}/run.csv`


| Flag         | Default                           | Description                          |
| ------------ | --------------------------------- | ------------------------------------ |
| `--datasets` | `output/metadata/datasets_v2.csv` | Source accession list                |
| `--baseline` | `output/metadata/datasets.csv`    | Accessions to exclude from migration |
| `--dry-run`  | off                               | Log planned uploads only             |


With the current metadata files, the default selection is 1,048 accessions. The process exits with status 1 if any selected row fails.

**Log:** `logs/migrate_gcs_to_r2.log`

```sh
uv run python pipelines/migrate_gcs_to_r2.py --dry-run
uv run python pipelines/migrate_gcs_to_r2.py
```

---



## `run_clustering_pipeline.py`

Runs `[cluster_validation](../scripts/cluster_validation/README.md)` per accession: download raw h5ad (local cache, else R2 raw mirror, else GCS), cluster, upload `{srx}_clustered.h5ad` to R2. Skips accessions whose output key already exists.

Per-run artifacts also land under the run directory (`figs/`, `data/`) before upload; local clustered files are deleted after a successful upload.

**Output:** `output/clustering_pipeline/{timestamp}/`


| File             | Description                                                |
| ---------------- | ---------------------------------------------------------- |
| `run.csv`        | Per-accession status and clustering metrics                |
| `metadata.json`  | Run config snapshot (`r2_prefix`, paths, optional `notes`) |
| `figs/`, `data/` | Copies of figures and clustered h5ad for this batch run    |



| Flag          | Default                           | Description                              |
| ------------- | --------------------------------- | ---------------------------------------- |
| `--datasets`  | `output/metadata/datasets.csv`    | Accession list                           |
| `--r2-prefix` | `clustering_pipeline_{timestamp}` | R2 folder for clustered outputs          |
| `--metadata`  | none                              | Free-form note stored in `metadata.json` |
| `--workers`   | `1`                               | Parallel accession workers               |


**Log:** `logs/clustering_pipeline.log`

```sh
uv run python pipelines/run_clustering_pipeline.py --workers 4
```

---



## `run_cytetype_pipeline.py`

Runs `[cytetype_runner](../scripts/cytetype_runner/README.md)` on clustered h5ad files already in R2. Requires `CYTETYPE_API_KEY`. Accessions run **serially**; optional spacing between run starts via `--min-interval`.

Downloads each input from `{clustering_prefix}/{srx}_clustered.h5ad`, annotates with study context from contexts JSONL, uploads `{r2_prefix}/{srx}_annotated.h5ad`. Per-accession metadata sent to CyteType is derived from the datasets CSV row (all columns, stringified).

**Output:** `output/cytetype_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)


| File              | Description                                |
| ----------------- | ------------------------------------------ |
| `run.csv`         | Per-accession status, timing, R2 keys      |
| `job_details.csv` | CyteType `job_id`, `report_url`, `api_url` |
| `metadata.json`   | Run config snapshot                        |



| Flag                  | Default                         | Description                                                   |
| --------------------- | ------------------------------- | ------------------------------------------------------------- |
| `--datasets`          | `output/metadata/datasets.csv`  | Accession list                                                |
| `--contexts`          | `output/context/contexts.jsonl` | Study context cache                                           |
| `--clustering-prefix` | required                        | R2 prefix from a clustering pipeline run                      |
| `--r2-prefix`         | `cytetype_pipeline_{timestamp}` | R2 folder for annotated outputs                               |
| `--metadata`          | none                            | Run-level note in `metadata.json` only (not sent to CyteType) |
| `--min-interval`      | `0`                             | Minimum seconds between starting consecutive accessions       |
| `--dry-run`           | off                             | Write plan CSVs without R2 or API calls                       |


**Log:** `logs/cytetype_pipeline.log`

```sh
uv run python pipelines/run_cytetype_pipeline.py \
  --clustering-prefix clustering_pipeline_20260511_140000 \
  --dry-run
```

---



## `run_cyteonto_pipeline.py`

Runs `[cyteonto](../scripts/cyteonto/README.md)` on annotated h5ad files already in R2. Lists every `{srx}_annotated.h5ad` under `--input-prefix` and runs accessions **serially**, blocking on each CyteOnto poll loop before starting the next.

Downloads each input from `{input_prefix}/{srx}_annotated.h5ad`, submits to the CyteOnto API, moves the result CSV to `results/{srx}_cyteonto.csv`, and uploads `{r2_prefix}/{srx}_cyteonto.csv`. Skips accessions whose output key already exists in R2.

**Output:** `output/cyteonto_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)


| File                         | Description                                                     |
| ---------------------------- | --------------------------------------------------------------- |
| `run.csv`                    | Per-accession status, timing, R2 keys, `run_id`, local CSV path |
| `results/{srx}_cyteonto.csv` | Per-accession CyteOnto similarity CSV                           |
| `metadata.json`              | Run config snapshot                                             |



| Flag                | Default                         | Description                                                              |
| ------------------- | ------------------------------- | ------------------------------------------------------------------------ |
| `--input-prefix`    | required                        | R2 prefix containing annotated h5ads (e.g. from a cytetype pipeline run) |
| `--r2-prefix`       | `cyteonto_pipeline_{timestamp}` | R2 folder for CyteOnto result CSVs                                       |
| `--metadata`        | none                            | Free-form note in `metadata.json`                                        |
| `--poll-interval-s` | `10`                            | Seconds between CyteOnto result polls                                    |
| `--poll-timeout-s`  | `3600`                          | Seconds before a CyteOnto run raises `TimeoutError`                      |
| `--min-interval`    | `0`                             | Minimum seconds between starting consecutive accessions                  |
| `--dry-run`         | off                             | Write plan CSVs without R2 or API calls                                  |


**Log:** `logs/cyteonto_pipeline.log`

For an in-flight accession, look up the `run_id` in `logs/cyteonto_pipeline.log` or the pending stub under `output/cyteonto/runs/`.

```sh
nohup uv run python pipelines/run_cyteonto_pipeline.py \
  --input-prefix cytetype_pipeline_20260522_175813 \
  > logs/cyteonto_pipeline.nohup.out 2>&1 &
```

---



## `cluster_stats.py`

Downloads clustered h5ad files from an R2 prefix, builds per-SRX `cell_type` × `leiden_merged` count matrices, and aggregates normalized Shannon entropy (NSE) and KL divergence (KLD) via `[cluster_validation.cell_type_metrics](../scripts/cluster_validation/README.md)`.

**Output:** `output/cluster_stats/{r2_prefix}/` (default; override with `--output-dir`)


| File                               | Description                                           |
| ---------------------------------- | ----------------------------------------------------- |
| `run.csv`                          | Per-accession status                                  |
| `metadata.json`                    | Run config snapshot                                   |
| `cluster_stats.json`               | Nested counts: `{srx: {cell_type: {cluster: count}}}` |
| `nse_matrix.csv`, `kld_matrix.csv` | Accessions × cell types                               |
| `cell_type_summary.csv`            | Cell-type level means                                 |
| `cell_type_metrics.png`            | Summary bar chart                                     |



| Flag           | Default                            | Description                                                        |
| -------------- | ---------------------------------- | ------------------------------------------------------------------ |
| `--r2-prefix`  | required                           | Clustering run prefix (e.g. `clustering_pipeline_20260511_140000`) |
| `--output-dir` | `output/cluster_stats/{r2_prefix}` | Output root                                                        |
| `--metadata`   | none                               | Note in `metadata.json`                                            |
| `--workers`    | `1`                                | Parallel accession workers                                         |


**Log:** `logs/cluster_stats.log`

Downstream analysis: `[notebooks/analysis/cluster_stats.ipynb](../notebooks/analysis/cluster_stats.ipynb)`.

```sh
uv run python pipelines/cluster_stats.py \
  --r2-prefix clustering_pipeline_20260511_140000 \
  --workers 4
```

---



## `run_annotation_inspection_pipeline.py`

Inspects CyteType-annotated h5ads from R2, joins CyteOnto cytescores, and writes a pair-level summary plus optional extremes table via `[annotation_inspector](../scripts/annotation_inspector/README.md)`.

Downloads `{input_prefix}/{srx}_annotated.h5ad` and `{cyteonto_prefix}/{srx}_cyteonto.csv` (when present), inspects each accession, deletes local cache files, and streams outputs as workers finish.

**Output:** `output/annotation_inspection_pipeline/{timestamp}/` (or `dry_run_{timestamp}/`)


| File            | Description                                                                 |
| --------------- | --------------------------------------------------------------------------- |
| `summary.csv`   | Pair-level rows: labels, confidence, cytescore, n_cells, report_url         |
| `extremes.csv`  | Top/bottom cytescore STATE types per CyteType label (when extremes enabled) |
| `run.csv`       | Per-accession status and timing                                             |
| `metadata.json` | Run config snapshot                                                         |



| Flag                | Default        | Description                                                         |
| ------------------- | -------------- | ------------------------------------------------------------------- |
| `--input-prefix`    | required       | R2 prefix containing annotated h5ads                                |
| `--cyteonto-prefix` | required       | R2 prefix containing `{srx}_cyteonto.csv` files                     |
| `--workers`         | `1`            | Concurrent R2 fetch + inspect workers                               |
| `--top-n`           | `10`           | Top/bottom STATE cell types per CyteType label for extremes         |
| `--no-extremes`     | off            | Skip `extremes.csv` (written by default)                            |
| `--from-summary`    | none           | Rebuild `extremes.csv` from an existing `summary.csv` (no R2 fetch) |
| `--output-dir`      | summary parent | Output dir for `--from-summary`                                     |
| `--metadata`        | none           | Note in `metadata.json`                                             |
| `--dry-run`         | off            | Write plan CSV without R2 or inspection                             |


**Log:** `logs/annotation_inspection_pipeline.log`

Downstream analysis: `[notebooks/analysis/annotation_inspection.ipynb](../notebooks/analysis/annotation_inspection.ipynb)`.

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

Concatenates raw scBaseCount h5ads into one merged atlas via `[h5ad_concat](../scripts/h5ad_concat/README.md)`: download each file from R2, validate and QC it, align to the reference gene axis, and merge the passing files. Unlike the other runners it has no CLI flags; edit the `H5adConcatConfig` in the script (for example `atlasR2Key`, and `uploadAtlas=True` to push the atlas to R2).

Reads `output/metadata/datasets.csv` by default (config `datasetsPath`).

**Output:** `output/atlas/data/` (`atlas.h5ad`, `atlas_config.json`, `atlas_files.jsonl`, `atlas_result.json`)

**Log:** `logs/h5ad_concat.log`

```sh
uv run python pipelines/run_atlas_concat.py
```

---



## Atlas postprocessing parameter selection and production

Atlas postprocessing (normalize, HVG, PCA, Harmony, neighbors, UMAP, Leiden) is split into two runners under `[scripts/atlas_postprocessing/](../scripts/atlas_postprocessing/)`:

1. `[select_atlas_parameters.py](select_atlas_parameters.py)` loads the **full atlas** (`--input`), replaces the in-memory object with a study-proportional representative sample of size `--sample-cells`, then calibrates or validates on that sample.
2. `[run_atlas_postprocessing.py](run_atlas_postprocessing.py)` runs one resolved parameter set on the full atlas (or any chosen input). It does not sample.

Sampling is deterministic (seed `0`), stratified by `--batch-key` (default `study_accession`). Every study gets at least one cell when `N` is at least the study count; remaining slots use largest-remainder proportions by study size. Requests smaller than the study count, or larger than the atlas, are rejected. The full atlas is loaded into memory once per command, so plan RAM for the whole object even though selection runs on the sample.

Harmony remains a method-specific stage (`X_pca_harmony`, Harmony UMAP plots). Overall outputs use the `post` layout under `output/atlas/<date>/post/` (legacy `v2` paths remain valid). Full-atlas production builds **only** the Harmony-corrected neighbor graph, UMAP, and Leiden partition. Subset validation still builds both the uncorrected and Harmony-corrected embeddings so scIB can compare `X_pca` vs `X_pca_harmony`, then RF-merges the approved Leiden partition.

### Recommended workflow

```sh
# 1) Calibrate on a study-proportional sample of the full atlas
uv run python pipelines/select_atlas_parameters.py calibrate \
  --input output/atlas/2026-08-12/atlas.h5ad \
  --sample-cells 100000 \
  --output-dir output/atlas/2026-08-12/post/parameter_selection/cluster_validation

# 2) Inspect metrics/figures, copy the template, optionally edit resolution
cp output/atlas/2026-08-12/post/parameter_selection/cluster_validation/parameters_template.json \
   output/atlas/2026-08-12/post/parameter_selection/cluster_validation/approved_parameters.json

# 3) Validate the approved set on a fresh sample with RF merge + scIB
uv run python pipelines/select_atlas_parameters.py validate \
  --input output/atlas/2026-08-12/atlas.h5ad \
  --sample-cells 100000 \
  --parameters-json output/atlas/2026-08-12/post/parameter_selection/cluster_validation/approved_parameters.json \
  --output-dir output/atlas/2026-08-12/post/subset_validation

# 4) After reviewing scIB and RF diagnostics, run production on the full atlas
uv run python pipelines/run_atlas_postprocessing.py \
  --input output/atlas/2026-08-12/atlas.h5ad \
  --output output/atlas/2026-08-12/post/production/atlas_post.h5ad \
  --parameters-json output/atlas/2026-08-12/post/parameter_selection/cluster_validation/approved_parameters.json
```

Calibration matches the cluster-validation control flow on a Harmony graph: fixed 2,000 batch-aware HVGs, adaptive PC count (50 computed, min 15, 50% cumvar target), fixed 15 neighbors, and a matched-Jaccard resolution sweep over `0.1..1.9` step `0.1`. The argmax is advisory only. Cell-type labels are weak priors for scoring, not ground truth. Create `approved_parameters.json` manually; the runner never writes it.

For this atlas, `100000` cells is the recommended balance between representation and calibration cost. Use `50000` for a faster preliminary run or `200000` to check stability with stronger rare-population coverage. Use the same sample size for calibration and validation.

Approved JSON shape (`AtlasPostprocessingParameters`, camelCase):

```json
{
  "nTopGenes": 2000,
  "nPcs": 18,
  "nNeighbors": 15,
  "resolution": 0.8,
  "calibrationSummary": "output/atlas/2026-08-12/post/parameter_selection/cluster_validation/calibration_summary.json"
}
```

`--parameters-json` is authoritative for those four tuning knobs. Do not combine it with `--n-top-genes`, `--n-pcs`, `--n-neighbors`, or `--resolution`. Scalar overrides remain valid when `--parameters-json` is omitted. Approved HVG/PC/neighbor values must match the singleton candidates recorded by calibration; any evaluated resolution may be approved.

### `select_atlas_parameters.py calibrate`

Loads the full atlas from `--input`, samples `--sample-cells` in memory, then runs the fixed graph + resolution selection on that sample.

**Output root default:** `output/atlas/v2/post/parameter_selection/`


| File                                     | Description                                                                                          |
| ---------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `calibration_summary.json`               | Fixed graph method, candidates, advisory recommendation, timings, artifact paths, sampling metadata |
| `parameters_template.json`               | Editable starting point seeded with resolved graph values and recommended resolution                 |
| `metrics/resolution.csv`                 | Matched Jaccard and cluster counts per resolution                                                    |
| `figures/resolution_matched_jaccard.png` | Resolution selection diagnostic with the advisory argmax marked                                      |


| Flag             | Required | Description                                              |
| ---------------- | -------- | -------------------------------------------------------- |
| `--input`        | yes      | Full atlas h5ad                                          |
| `--sample-cells` | yes      | Exact cell count for the in-memory representative sample |


Fixed method values: HVGs `2000`, neighbors `15`, PC chooser `nPcsCompute=50`, `nPcsMin=15`, `nPcsCumvarTarget=0.5`. Default resolutions follow cluster validation (`0.1` through `1.9` step `0.1`). Override the grid with `--resolution-candidates`.

`sampling` in `calibration_summary.json` records `sourceCells`, `sampleCells`, `method` (`studyProportional`), `stratifyKey`, `seed` (`0`), and `nStudies`. `input` remains the full atlas path.

**Log:** `logs/select_atlas_parameters.log`

### `select_atlas_parameters.py validate`

Loads the full atlas, draws a sample with the same `--sample-cells` policy, runs one dual-embedding validation pass with the approved JSON (uncorrected + Harmony graphs), RF-merges `leiden_atlas` into `leiden_merged` on normalized pre-scale HVGs, then scIB on `X_pca` vs `X_pca_harmony`. Review the full scIB table and RF diagnostics before production; there is no automatic pass/fail gate.

**Output root default:** `output/atlas/v2/post/subset_validation/`


| File                                             | Description                                                              |
| ------------------------------------------------ | ------------------------------------------------------------------------ |
| `atlas_pp_subset.h5ad`                           | Processed sample with `leiden_atlas` and `leiden_merged`                 |
| `atlas_pp_subset_run.json`                       | Run summary including `sampling` and `clustersMerged`                    |
| `subset_validation_summary.json`                 | Approved vs recommended resolution, RF diagnostics, sampling, scIB paths |
| `figures/`                                       | Scree + atlas-scale UMAPs                                                |
| `scib/scib_results.csv`, `scib/scib_results.svg` | scIB report                                                              |


| Flag                | Required | Description                                              |
| ------------------- | -------- | -------------------------------------------------------- |
| `--input`           | yes      | Full atlas h5ad                                          |
| `--sample-cells`    | yes      | Exact cell count for the in-memory representative sample |
| `--parameters-json` | yes      | Approved parameter JSON from calibration review          |



### `run_atlas_postprocessing.py`

Lightweight production runner. No sweeps and no scIB. Builds HVG + PCA once, then a **Harmony-only** neighbor graph, UMAP, and Leiden. Parallel production UMAP is intentionally non-reproducible for speed; revisit and freeze a seeded embedding before publication.

**Output defaults:** `output/atlas/v2/post/production/atlas_pp.h5ad` and `.../figures/`.


| File                                      | Description                                                                                                                                     |
| ----------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `{output}.h5ad`                           | Processed atlas (HVGs in `.X`, full-gene counts in `.raw`) with `X_pca_harmony`, `X_umap`, and `leiden_atlas` (no uncorrected graph artifacts) |
| `{output_stem}_run.json`                  | Run summary including `workflow`, `nNeighbors`, `nJobs`, and optional `parametersJson` (`clustersUncorrected` is `null` in production)         |
| `{figs_dir}/umap_{batch_key}_harmony.png` | Harmony UMAP colored by batch and cell type                                                                                                     |
| `{figs_dir}/pca_scree.png`                | PCA scree plot                                                                                                                                  |



| Flag                | Default                                         | Description                                                           |
| ------------------- | ----------------------------------------------- | --------------------------------------------------------------------- |
| `--input`           | `output/atlas/v2/atlas.h5ad`                    | Input atlas h5ad                                                      |
| `--output`          | `output/atlas/v2/post/production/atlas_pp.h5ad` | Output atlas h5ad                                                     |
| `--figs-dir`        | `output/atlas/v2/post/production/figures`       | Directory for UMAP and scree PNGs                                     |
| `--parameters-json` | none                                            | Approved parameter JSON (optional)                                    |
| `--batch-key`       | `study_accession`                               | `obs` column for Harmony batch correction                             |
| `--cell-type-key`   | `cell_type`                                     | `obs` column for cell-type UMAP plots                                 |
| `--n-top-genes`     | `2000`                                          | Number of HVGs (disallowed with `--parameters-json`)                  |
| `--n-pcs`           | `20`                                            | PCs used for the neighbor graph (disallowed with `--parameters-json`) |
| `--n-pcs-compute`   | `50`                                            | PCs computed by PCA                                                   |
| `--n-neighbors`     | `15`                                            | Neighbors for the graph (disallowed with `--parameters-json`)         |
| `--resolution`      | `1.0`                                           | Leiden resolution (disallowed with `--parameters-json`)               |
| `--no-plots`        | off                                             | Skip writing PNGs                                                     |
| `--threads`         | `0`                                             | Thread budget for Scanpy, Harmony, and parallel UMAP (`0` = library defaults) |


**Log:** `logs/atlas_postprocessing.log`

### `compare_atlas_batch_keys.py`

Fixed-subset comparison of Harmony batch definitions on the existing 100k validation sample. Reuses the shared `X_pca` and study `X_pca_harmony` from subset validation, joins `tech_10x` from `datasets_v2.csv`, builds study×technology and SRX Harmony embeddings, then runs scIB once per evaluation batch key across all embeddings.

**Output default:** `output/atlas/v2/post/batch_key_comparison/`


| File | Description |
| --- | --- |
| `atlas_pp_batch_comparison.h5ad` | Subset with attached tech metadata and keyed Harmony embeddings |
| `batch_key_comparison_summary.json` | Baseline checks, join audit, batch cardinalities, scIB paths, timings |
| `scib/eval_batch=*/scib_results.{csv,svg}` | One scIB report per evaluation batch key |
| `scib/scib_matrix_long.csv` | Combined long-form metric table |


| Flag | Default | Description |
| --- | --- | --- |
| `--subset-h5ad` | subset validation h5ad | Frozen sample with shared PCA |
| `--subset-run-json` | subset run JSON | Pins nTopGenes/nPcs/nNeighbors/resolution |
| `--production-run-json` | production run JSON | Confirms the same baseline |
| `--datasets` | `output/metadata/datasets_v2.csv` | Source of `tech_10x` |
| `--output-dir` | `output/atlas/v2/post/batch_key_comparison` | Comparison root |
| `--threads` | `0` | Harmony thread budget |
| `--scib-jobs` | `6` | scIB n_jobs |
| `--force-scib` | off | Re-run existing scIB evaluation dirs |
| `--skip-scib` | off | Write embeddings only |
| `--reuse-comparison-h5ad` | off | Resume from an existing comparison h5ad |


Example:

```sh
uv run python pipelines/compare_atlas_batch_keys.py \
  --subset-h5ad output/atlas/v2/post/subset_validation/atlas_pp_subset.h5ad \
  --output-dir output/atlas/v2/post/batch_key_comparison
```

**Log:** `logs/compare_atlas_batch_keys.log`

### `run_atlas_de_analysis.py`

Checkpointed disease DE and noteworthy-gene discovery on a postprocessed atlas. Aggregates sparse `.raw` counts by `SRX_accession × leiden_atlas` without building a second full-gene AnnData copy, then runs study-aware DESeq2 contrasts and adaptive shortlist ranking.

**Output default:** `output/atlas/v2/analysis/production/`

| File | Description |
| --- | --- |
| `checkpoints/pseudobulk.h5ad` | Sample x Leiden pseudobulk with `psbulk_props` |
| `noteworthy_gene_shortlist.csv` | Up to ~20 primary review candidates |
| `noteworthy_gene_extended.csv` | Extended queue capped near 60 |
| `candidate_thresholds.csv` | Per-class adaptive score cutoffs |
| `figures/` | Score distributions, heatmaps, volcanoes, evidence panels |
| `run_summary.json` | Resolved config, counts, and peak memory |

| Flag | Default | Description |
| --- | --- | --- |
| `--stage` | `all` | `aggregate`, `analyze`, or `all` |
| `--atlas` | production postprocessed atlas | Input h5ad with counts in `.raw` |
| `--output-dir` | `output/atlas/v2/analysis/production` | Output root |
| `--memory-reserve-gib` | `256` | RAM reserved beyond the estimated sparse matrix |
| `--primary-budget` | `20` | Primary shortlist size |
| `--extended-budget` | `60` | Total review queue cap |
| `--force-aggregate` | off | Ignore an existing pseudobulk checkpoint |

Example:

```bash
uv run python pipelines/run_atlas_de_analysis.py --stage all \
  --atlas output/atlas/v2/post/production/atlas_v2_post.h5ad
```

**Log:** `logs/atlas_de_analysis.log`

Downstream exploratory sample wrapper: `[notebooks/analysis/analyze_atlas_DE.py](../notebooks/analysis/analyze_atlas_DE.py)` calls the same `disease_markers` modules on the 100k sample.
