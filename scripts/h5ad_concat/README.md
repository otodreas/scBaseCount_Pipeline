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

Requires R2 credentials (see [`storage/`](../storage/README.md)). Study batch labels come from `output/context/contexts.jsonl` via `studyAccession`.

## Validation gate

Each file is checked before it enters the concat. Failing files are recorded in `result.skipped` and excluded; the run continues.

| Check | Reason | Active |
|-------|--------|--------|
| MD5 matches stored `gcs-md5` metadata | `md5_mismatch` | yes |
| `studyAccession` resolves from contexts | `missing_study` | yes |
| At least one non-blank `cell_type` label | `cell_type_all_missing` | yes |
| Preprocessing passes | `preprocess_failed` | future |

## Config

| Field | Default | Description |
|-------|---------|-------------|
| `r2Keys` | (required) | Explicit R2 object keys to download and concat |
| `contextsPath` | `output/context/contexts.jsonl` | Study context lookup |
| `cellTypeKey` | `"cell_type"` | Column checked and filled |
| `studyKey` | `"study"` | Batch label column written by concat |
| `accessionKey` | `"accession"` | Experiment id column written during prepare |
| `missingLabel` | `"unknown"` | Fill value for blank `cell_type` |
| `join` | `"inner"` | Gene join strategy for concat |
| `cacheDir` | `data/h5ad_concat/cache` | Staging for downloads and partial merges |
| `outputPath` | `output/atlas/data/atlas.h5ad` | Merged atlas output |
| `maxLoadedElems` | `100_000_000` | Streaming chunk size for `concat_on_disk` |
| `mergeBatchSize` | `25` | Max prepared files on disk per concat batch |
| `verifyMd5` | `True` | Verify download against R2 `gcs-md5` metadata |

## Output

`run_h5ad_concat` returns `H5adConcatResult`:

- `outputPath`: merged atlas h5ad
- `nObs`, `nVars`: shape of merged object
- `nFilesConcatenated`: count of files that passed validation
- `studiesSeen`: unique `studyAccession` values in the merged atlas
- `skipped`: list of rejected files with `r2Key`, `accession`, and `reason`

Logs append to `logs/h5ad_concat.log`.

## Memory and disk

Prepared files are concatenated with `anndata.experimental.concat_on_disk`, which streams sparse chunks instead of loading full objects into RAM. `mergeBatchSize` bounds how many prepared h5ads sit on disk at once; each batch is merged and deleted before the next batch starts. The merged output file must still fit on disk.

## TODO

Future work is marked in source with `# TODO(...)` comments at the call site:

- **`datasets-csv`** (`pipeline.py`): accept `output/metadata/datasets.csv` accessions mapped to R2 raw keys as an alternative to explicit `r2Keys`.
- **`preprocess`** (`config.py`, `prepare.py`): add `preprocess: bool = True`; when enabled, run `cluster_validation.preprocess` per file and skip failures as `preprocess_failed`.
- **`upload-atlas`** (`pipeline.py`): upload the merged atlas to R2 via `upload_to_r2` with base64 MD5 metadata.
