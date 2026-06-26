# celltypist_runner

Runs CellTypist on a raw or preprocessed h5ad and writes per-cell predicted labels onto `adata.obs`.

## Usage

```python
import scanpy as sc
from celltypist_runner import CellTypistRunnerConfig, annotate_celltypist

adata = sc.read("data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX17412841.h5ad")
cfg = CellTypistRunnerConfig(modelName="Nuclei_Lung_Airway.pkl")
adata = annotate_celltypist(adata, cfg)
```

The input `adata` is left on raw counts; normalization and log1p are applied on an internal copy only for inference. Annotation uses per-cell model predictions only (CellTypist over-clustering and majority voting are not used). Results are extracted via `predictions.to_adata()` and copied onto the input object as `obs["predicted_labels"]` (configurable via `predictedLabelKey`) and `obs["conf_score"]`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `modelName` | `Nuclei_Lung_Airway.pkl` | CellTypist model filename (downloaded when `downloadIfMissing=True`) |
| `targetSum` | `10000` | Target sum for `normalize_total` before annotation |
| `geneSymbolCol` | `gene_symbols` | `var` column used to set gene symbols for CellTypist |
| `predictedLabelKey` | `predicted_labels` | `obs` column written with CellTypist output |
| `downloadIfMissing` | `True` | Call `celltypist.models.download_models()` before loading the model |

All default paths are relative to the repo root.

## Logging

Steps are appended to `logs/celltypist_runner.log`.
