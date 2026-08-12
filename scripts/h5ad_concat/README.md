# h5ad_concat

Download scBaseCount h5ad files from R2, validate each file in-pipeline, enrich `obs`, and concatenate passing files into a local gzip-compressed atlas held in memory during merge.

## Usage

```python
from h5ad_concat import H5adConcatConfig, run_h5ad_concat

cfg = H5adConcatConfig(
    datasetsPath="output/metadata/datasets.csv",
)

result = run_h5ad_concat(cfg)
print(result.nObs, result.studiesSeen, result.skipped)
```

The required CSV columns are `file_path`, `srx_accession`, and `study_accession`; `file_path` holds `gs://` URIs mapped to R2 raw keys. To upload the completed atlas:

```python
cfg = H5adConcatConfig(
    datasetsPath="output/metadata/datasets.csv",
    uploadAtlas=True,
    atlasR2Key="atlas/2026-01-12/atlas.h5ad",
)
```

Requires R2 credentials (see `[storage/](../storage/README.md)`). Each concatenated file gains a single new obs column, `study_accession` (the `batchKey`), read from the corresponding datasets CSV row. This is the experimental batch key for downstream integration.

## Batch effect

`study_accession` is a coarse batch key. A single `studyAccession` can span multiple experiments that differ in biologically meaningful ways, including patient age, disease subtype, individual donor, and tumor grading. Grouping these under one batch key therefore folds real biological variation into a single batch, so downstream integration that treats a batch as a technical unit will treat some genuine biological signal as batch effect.

The scale of this collapse is easy to quantify. The "Batch Key Cardinality" cell in `notebooks/pipeline/study_context.ipynb` counts unique accessions and unique `studyAccession` values in `output/context/contexts.jsonl`. At the time of writing, 772 experiment accessions map to only 117 distinct studies, so using `studyAccession` as the batch key collapses those 801 experiments into 121 batches.

For a concrete case, see `notebooks/pipeline/study_context.ipynb`, which diffs experiments `ERX11662385` and `ERX11662370` from study `PRJEB68315`: they share the same `studyAccession` yet come from different patients (ages 67 and 52) with different disease subtypes (non-small cell lung cancer and lung adenocarcinoma) and different tumor grading. This is one example, not an exhaustive audit, but the pattern is inherent to how ENA study accessions aggregate samples, so we treat within-study biological heterogeneity as expected rather than exceptional.

The block below is pasted directly from the field-by-field diff output of that notebook, kept here only as a visual demonstration:

```text
accession:
  ERX11662385: ERX11662385
  ERX11662370: ERX11662370
biological.sampleAttributes.age:
  ERX11662385: 67
  ERX11662370: 52
biological.sampleAttributes.disease:
  ERX11662385: non-small cell lung cancer
  ERX11662370: lung adenocarcinoma
biological.sampleAttributes.individual:
  ERX11662385: Patient 23
  ERX11662370: Patient 24
biological.sampleAttributes.original source name:
  ERX11662385: TB21.0086
  ERX11662370: TB21.0106
biological.sampleAttributes.tumor grading:
  ERX11662385: T4N0
  ERX11662370: T3N0M0
biological.sampleTitle:
  ERX11662385: P23_B2
  ERX11662370: P24_B2
runAccessions:
  ERX11662385: ['ERR12251785', 'ERR12252043', 'ERR12251741', 'ERR12251857', 'ERR12252030', 'ERR12251724', 'ERR12251755', 'ERR12251994']
  ERX11662370: ['ERR12252007', 'ERR12251936', 'ERR12251941', 'ERR12251754']
sampleAccession:
  ERX11662385: SAMEA114591111
  ERX11662370: SAMEA114591096
```

## Validation gate

Each file is checked before it enters the concat. Failing files are recorded in `result.skipped` and excluded; the run continues.


| Check                                    | Reason                  | Active |
| ---------------------------------------- | ----------------------- | ------ |
| R2 download succeeds                       | `download_failed`       | yes    |
| MD5 matches stored `gcs-md5` metadata    | `md5_mismatch`          | yes    |
| Downloaded h5ad loads cleanly            | `read_failed`           | yes    |
| `obs` accession column is a single value matching the file accession | `accession_mismatch` | yes |
| At least one non-blank `cell_type` label | `cell_type_all_missing` | yes    |
| At least `minCellsAfterQc` cells remain after QC | `too_few_cells` | yes    |
| Cell dropout fraction at or below `minPctCellsAfterQc` | `excessive_cell_dropout` | yes |
| Gene axis maps to reference | `gene_axis_mismatch` | yes |




## Config

| Field            | Default                         | Description                                    |
| ---------------- | ------------------------------- | ---------------------------------------------- |
| `datasetsPath`   | `output/metadata/datasets.csv` | Path to datasets CSV; `file_path`, `srx_accession`, and `study_accession` define each input |
| `cellTypeKey`    | `"cell_type"`                   | Column checked and filled                      |
| `batchKey`       | `"study_accession"`             | obs column holding the batch key (ENA study accession) |
| `accessionKey`   | `"SRX_accession"`               | Existing obs column holding the per-file experiment accession |
| `missingLabel`   | `"UNKNOWN"`                     | Fill value for blank `cell_type`               |
| `geneInfoPath`   | `data/scbasecount/2026-01-12/star_references/Homo_sapiens/hg38_2020/geneInfo.tab` | STAR reference gene axis for alignment |
| `cacheDir`         | `data/h5ad_concat/cache`        | Staging for transient raw downloads             |
| `outputPath`       | `output/atlas/data/atlas.h5ad`  | Merged atlas output                             |
| `downloadBatchSize`| `8`                             | Reserved for concurrent download batching       |
| `compression`      | `"gzip"`                        | h5ad write compression for the atlas            |
| `verifyMd5`        | `True`                          | Verify download against R2 `gcs-md5` metadata   |
| `uploadAtlas`      | `False`                         | Upload merged atlas to R2 after write           |
| `atlasR2Key`       | `None`                          | R2 object key for atlas upload (required when `uploadAtlas` is true) |
| `preprocess`         | `True`                          | Run per-file QC gate before concat admission |
| `conserveLayers`     | `False`                         | Reindex every layer (e.g. STARsolo UniqueAndMult matrices) onto the reference axis instead of `X` only |
| `minGenesPerCell`    | `200`                           | Minimum genes detected per cell |
| `maxPctMito`         | `0.2`                           | Maximum mitochondrial read fraction per cell, as a fraction in (0, 1]; `1.0` keeps every cell |
| `maxPctRibo`         | `1.0`                           | Maximum ribosomal read fraction per cell, as a fraction in (0, 1]; `1.0` keeps every cell |
| `maxPctHb`           | `1.0`                           | Hemoglobin read fraction ceiling in (0, 1]; `1.0` records the metric without filtering |
| `minCellsPerGene`    | `0`                             | Minimum cells expressing a gene; set `0` to disable |
| `minCellsAfterQc`    | `100`                           | Absolute floor: minimum cells remaining after QC or file is skipped |
| `minPctCellsAfterQc` | `0.4`                           | Relative floor as a fraction in [0, 1]: reject when less than this fraction of input cells remain after QC; set `0` to disable |




## Output

`run_h5ad_concat` returns `H5adConcatResult`:

- `outputPath`: local atlas `.h5ad` path, or the result manifest path when the upload succeeds
- `nObs`, `nVars`: shape of merged object
- `nFilesConcatenated`: count of files that passed validation
- `nFilesSkipped`: count of rejected files
- `studiesSeen`: unique `studyAccession` values in the merged atlas
- `skipped`: list of rejected files with `r2Key`, `accession`, `reason`, optional `studyAccession`, and optional `qc`
- `cellFilterOrder`: sequential per-cell filter names used for drop attribution
- `qcSummary`: aggregate QC totals for concatenated files and all QC-processed files
- `files`: per-file records matching the JSONL log, including sequential `qc.nCellsDroppedByFilter`
- `fileLogPath`: local path to the append-safe per-file JSONL log written during the run
- `configPath`: local path to the config manifest written after the run
- `atlasR2Key`: R2 object key when `uploadAtlas` is enabled and upload succeeds; otherwise `None`
- `atlasFileLogR2Key`: R2 object key for the file log when `uploadAtlas` is enabled and upload succeeds; otherwise `None`
- `atlasConfigR2Key`: R2 object key for the config manifest when `uploadAtlas` is enabled and upload succeeds; otherwise `None`
- `atlasResultR2Key`: R2 object key for the result manifest when `uploadAtlas` is enabled and upload succeeds; otherwise `None`
- `conserveLayers`: whether alignment reindexed all layers onto the reference axis for this run

A run writes up to four files next to the atlas output path (`cfg.outputPath`, default `output/atlas/data/atlas.h5ad`), all locally regardless of upload. R2 upload is optional: when `uploadAtlas` is true and upload verifies, each file is uploaded to R2 under the atlas key stem and the local `.h5ad` is deleted only after the atlas `gcs-md5` metadata matches the pre-upload local MD5. The pipeline refuses to start when `outputPath` already exists, or when `uploadAtlas` is set and `atlasR2Key` already exists.

| Output | Local file | R2 key when `uploadAtlas` | Written | Contents |
| ------ | ---------- | ------------------------- | ------- | -------- |
| Config | `atlas_config.json` | `{stem}_config.json` | At run start | The `H5adConcatConfig` used for the run |
| File log | `atlas_files.jsonl` | `{stem}_files.jsonl` | One JSON object per file during the loop | Per-file record with `accession`, `studyAccession`, `r2Key`, `status`, `skipReason`, and `qc` |
| Atlas | `atlas.h5ad` (deleted after successful upload) | `atlasR2Key` (the configured key) | After concatenation | Merged gzip-compressed AnnData |
| Result | `atlas_result.json` | `{stem}_result.json` | After concatenation | The `H5adConcatResult` including `files` and `qcSummary` |

The config manifest is written up front and the file log is flushed record by record during the loop, so both survive an interrupted or failed run (including when every file is rejected). The atlas and result manifest are written only once concatenation succeeds. Together these make the run inputs and outputs readable without pulling the atlas from R2. When the atlas is uploaded, `outputPath` in the result points at the result manifest instead of the deleted local `.h5ad`.

Logs append to `logs/h5ad_concat.log`. Each run logs a start line and a completion line; `KeyboardInterrupt` (Ctrl-C) is logged as an interruption before the exception is re-raised.

## Memory and disk

Each h5ad is downloaded to a transient raw file under `cacheDir/raw`, loaded into memory with `read_h5ad`, validated and enriched, then the raw file is deleted immediately. Passing objects accumulate in RAM until all keys are processed, then `ad.concat` builds one in-memory atlas and `write_h5ad(..., compression="gzip")` writes the final output.

Peak disk usage is one raw h5ad (or a small concurrent download batch) plus the gzipped atlas. Peak RAM is roughly the sum of all loaded objects plus the concatenated result during `ad.concat`. The gzipped atlas must still fit on disk.

Per-file QC filters low-quality cells and genes while preserving raw counts in `X`. Normalization, HVG selection, and integration remain downstream on the merged atlas.

Each file is reindexed to the canonical `geneInfoPath` Ensembl-ID axis before concat, so the atlas gene space is fixed at 36,601 genes in reference order and the concat join is a no-op. Sparse-gene filtering is deferred downstream (`minCellsPerGene` defaults to `0`).

QC records `pct_counts_mt`, `pct_counts_ribo`, and `pct_counts_hb` per cell. Mitochondrial filtering is on by default (`maxPctMito`). Ribosomal and hemoglobin filtering are opt-in via `maxPctRibo` and `maxPctHb` (both default to `1.0`, which keeps every cell under the strict `<` comparison). Cell drops are attributed sequentially in the order `minGenesPerCell`, `maxPctMito`, `maxPctRibo`, `maxPctHb`.

## Cell-count gates and the dropout denominator

Per-cell QC filters low-quality cells; dataset-level gates then decide whether the filtered remainder is still trustworthy enough to merge. A file is rejected when either gate fails: fewer than `minCellsAfterQc` cells remain (absolute floor, default 100), or less than `minPctCellsAfterQc` percent cells remain  (relative floor, default 40%).

The relative gate is only meaningful when the input matrix is already cell-called. scBaseCount h5ads are produced by scRecounter with STARsolo `--soloCellFilter EmptyDrops_CR` (CellRanger-style empty-droplet removal; Youngblut et al. 2025, Methods 5.2). The denominator for `pctCellsDropped` is therefore called cells, not raw barcodes, so a high drop fraction signals systemic dataset badness rather than routine empty-droplet removal.

In the local `GeneFull/Homo_sapiens` sample, five of six files lose less than 1% of cells under default per-cell filters; the sparsest file (`SRX12708356`) loses about 49.5% (roughly half its cells fall below `minGenesPerCell=200`), which is the motivating case for the relative gate.

Reference: Youngblut et al. (2025), *scBaseCount: an AI agent-curated, uniformly processed, and continually expanding single cell data repository*, bioRxiv [10.1101/2025.02.27.640494](https://doi.org/10.1101/2025.02.27.640494). Local copy: `docs/zotero_pdfs/Youngblut et al. - 2025 - scBaseCount an AI agent-curated, uniformly processed, and continually expanding single cell data re.pdf`.

