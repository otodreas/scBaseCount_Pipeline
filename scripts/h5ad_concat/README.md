# h5ad_concat

Download scBaseCount h5ad files from R2, validate each file in-pipeline, enrich `obs`, and concatenate passing files into a local gzip-compressed atlas held in memory during merge.

## Usage

```python
from h5ad_concat import H5adConcatConfig, run_h5ad_concat

cfg = H5adConcatConfig(
    r2Keys=[
        "arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX13061245.h5ad",
        "arc-institute-virtual-cell-atlas/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX10048396.h5ad",
    ],
)

result = run_h5ad_concat(cfg)
print(result.nObs, result.studiesSeen, result.skipped)
```

Requires R2 credentials (see `[storage/](../storage/README.md)`). Each concatenated file gains a single new obs column, `study_accession` (the `batchKey`), resolved from `output/context/contexts.jsonl` via `ctx.study.studyAccession`. This is the experimental batch key for downstream integration.

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
| `studyAccession` resolves from contexts  | `missing_study`         | yes    |
| At least one non-blank `cell_type` label | `cell_type_all_missing` | yes    |
| Preprocessing passes                     | `preprocess_failed`     | future |




## Config


| Field            | Default                         | Description                                    |
| ---------------- | ------------------------------- | ---------------------------------------------- |
| `r2Keys`         | (required)                      | Explicit R2 object keys to download and concat |
| `contextsPath`   | `output/context/contexts.jsonl` | Study context lookup                           |
| `cellTypeKey`    | `"cell_type"`                   | Column checked and filled                      |
| `batchKey`       | `"study_accession"`             | obs column holding the batch key (ENA study accession) |
| `missingLabel`   | `"UNKNOWN"`                     | Fill value for blank `cell_type`               |
| `join`           | `"inner"`                       | Gene join strategy for concat                  |
| `cacheDir`         | `data/h5ad_concat/cache`        | Staging for transient raw downloads             |
| `outputPath`       | `output/atlas/data/atlas.h5ad`  | Merged atlas output                             |
| `downloadBatchSize`| `8`                             | Reserved for concurrent download batching       |
| `compression`      | `"gzip"`                        | h5ad write compression for the atlas            |
| `verifyMd5`        | `True`                          | Verify download against R2 `gcs-md5` metadata   |




## Output

`run_h5ad_concat` returns `H5adConcatResult`:

- `outputPath`: merged atlas h5ad
- `nObs`, `nVars`: shape of merged object
- `nFilesConcatenated`: count of files that passed validation
- `studiesSeen`: unique `studyAccession` values in the merged atlas
- `skipped`: list of rejected files with `r2Key`, `accession`, and `reason`

Logs append to `logs/h5ad_concat.log`.

## Memory and disk

Each h5ad is downloaded to a transient raw file under `cacheDir/raw`, loaded into memory with `read_h5ad`, validated and enriched, then the raw file is deleted immediately. Passing objects accumulate in RAM until all keys are processed, then `ad.concat` builds one in-memory atlas and `write_h5ad(..., compression="gzip")` writes the final output.

Peak disk usage is one raw h5ad (or a small concurrent download batch) plus the gzipped atlas. Peak RAM is roughly the sum of all loaded objects plus the concatenated result during `ad.concat`. The gzipped atlas must still fit on disk.

## TODO

Future work is marked in source with `# TODO(...)` comments at the call site:

- `input` (`config.py`): support `datasets.csv` as an input source (parse a column containing R2 keys) instead of requiring explicit `r2Keys`.
- `datasets-csv` (`pipeline.py`): resolve `cfg.r2Keys` from `output/metadata/datasets.csv` accessions mapped to R2 raw keys (see `pipelines/run_clustering_pipeline.py`).
- `preprocess` (`config.py`, `prepare.py`): add `preprocess: bool = True`; when enabled, run `cluster_validation.preprocess` per file and skip failures as `preprocess_failed`.
- `output` (`merge.py`): support appending to an existing atlas and/or building a new atlas version instead of always overwriting `outputPath`.
- `upload-atlas` (`pipeline.py`): upload the merged atlas to R2 via `upload_to_r2` with `_local_md5_b64` metadata.

