# h5ad_concat

Download scBaseCount h5ad files from R2, validate each file in-pipeline, enrich `obs`, and concatenate passing files into a local atlas with memory-safe on-disk merging.

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

## Validation gate

Each file is checked before it enters the concat. Failing files are recorded in `result.skipped` and excluded; the run continues.


| Check                                    | Reason                  | Active |
| ---------------------------------------- | ----------------------- | ------ |
| MD5 matches stored `gcs-md5` metadata    | `md5_mismatch`          | yes    |
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
| `cacheDir`       | `data/h5ad_concat/cache`        | Staging for downloads and partial merges       |
| `outputPath`     | `output/atlas/data/atlas.h5ad`  | Merged atlas output                            |
| `maxLoadedElems` | `100_000_000`                   | Streaming chunk size for `concat_on_disk`      |
| `mergeBatchSize` | `25`                            | Max prepared files on disk per concat batch    |
| `verifyMd5`      | `True`                          | Verify download against R2 `gcs-md5` metadata  |




## Output

`run_h5ad_concat` returns `H5adConcatResult`:

- `outputPath`: merged atlas h5ad
- `nObs`, `nVars`: shape of merged object
- `nFilesConcatenated`: count of files that passed validation
- `studiesSeen`: unique `studyAccession` values in the merged atlas
- `skipped`: list of rejected files with `r2Key`, `accession`, and `reason`

Logs append to `logs/h5ad_concat.log`.

## Memory and disk

Prepared files are concatenated with `anndata.experimental.concat_on_disk`, which streams sparse chunks instead of loading full objects into RAM. During merge, `mergeBatchSize` bounds how many prepared h5ads are folded per batch; each batch is merged and deleted before the next merge batch starts. The merged output file must still fit on disk.

Today the pipeline prepares every passing file before merge begins, so peak staging is the sum of all prepared h5ads. A later fold step can still require roughly 2x the atlas size on disk while copying the accumulator into the next fold file.

## TODO

Future work is marked in source with `# TODO(...)` comments at the call site:

- `input` (`config.py`): support `datasets.csv` as an input source (parse a column containing R2 keys) instead of requiring explicit `r2Keys`.
- `datasets-csv` (`pipeline.py`): resolve `cfg.r2Keys` from `output/metadata/datasets.csv` accessions mapped to R2 raw keys (see `pipelines/run_clustering_pipeline.py`).
- `preprocess` (`config.py`, `prepare.py`): add `preprocess: bool = True`; when enabled, run `cluster_validation.preprocess` per file and skip failures as `preprocess_failed`.
- `output` (`merge.py`): support appending to an existing atlas and/or building a new atlas version instead of always overwriting `outputPath`.
- `stream-pipeline` (`pipeline.py`, `merge.py`, `prepare.py`, `config.py`): interleave download/prepare with the batch/fold merge loop so only `mergeBatchSize` prepared files sit on disk at once, instead of staging all passing files before merge starts.
- `upload-atlas` (`pipeline.py`): upload the merged atlas to R2 via `upload_to_r2` with `_local_md5_b64` metadata.

