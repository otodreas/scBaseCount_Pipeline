# scib_benchmark

Run [scib-metrics](https://scib-metrics.readthedocs.io/) integration benchmarks on precomputed atlas embeddings and cache the results for later plotting.

Expects an h5ad with embedding keys in `obsm` (for example `X_pca`, `X_pca_harmony`) plus batch and cell-type columns in `obs`. The atlas is read in backed mode; `pre_integrated_embedding_obsm_key="X_pca"` uses uncorrected PCA as the baseline for PCR comparison and avoids recomputing PCA on `adata.X`.

## CLI

Batch runner: [`pipelines/run_scib_benchmark.py`](../../pipelines/run_scib_benchmark.py)

```sh
uv run python pipelines/run_scib_benchmark.py \
  --input output/atlas/v1/processed_1/atlas_harmony.h5ad \
  --out-dir output/atlas/v1/processed_1/scib
```

Useful flags: `--embeddings`, `--batch-key`, `--label-key`, `--n-jobs`, `--force` (re-run even when both cached outputs exist).

**Log:** `logs/scib_benchmark.log`

## Python

```python
from pathlib import Path

from scib_benchmark import run_scib_benchmark

run_scib_benchmark(
    input_h5ad=Path("output/atlas/v1/processed_1/atlas_harmony.h5ad"),
    out_dir=Path("output/atlas/v1/processed_1/scib"),
    batch_key="study_accession",
    label_key="cell_type",
    embedding_keys=["X_pca", "X_pca_harmony", "X_umap", "X_umap_uncorrected"],
    n_jobs=6,
    force=False,
)
```

## Outputs

Written under `--out-dir`:

| File | Description |
|------|-------------|
| `scib_results.csv` | Cached `Benchmarker.get_results()` table |
| `scib_results.svg` | Results table figure |

If both `scib_results.csv` and `scib_results.svg` already exist, the run exits immediately (unless `--force`). A partial cache (only one file) triggers a full re-run with a warning.

## Reload results

```python
import pandas as pd

df = pd.read_csv("output/atlas/v1/processed_1/scib/scib_results.csv", index_col=0)
df
```

The SVG is written by `Benchmarker.plot_results_table` during the benchmark run. To regenerate the figure alone, re-run with `--force` or call `plot_results_table` from a live `Benchmarker` session in a notebook.
