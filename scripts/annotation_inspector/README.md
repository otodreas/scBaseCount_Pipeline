# annotation_inspector

Inspect agreement between STATE labels (`obs["cell_type"]`) and CyteType labels on annotated h5ad files. Joins CyteOnto `cytescore_similarity`, maps CyteType cluster confidence, and can emit a pair-level summary CSV and an extremes table.

Batch orchestration lives in [`pipelines/run_annotation_inspection_pipeline.py`](../../pipelines/run_annotation_inspection_pipeline.py).

## Usage

```python
from pathlib import Path

from annotation_inspector import inspect_accession, top_bottom_by_cytetype, write_extremes_csv

pair_df = inspect_accession(
    "SRX17412841",
    Path("data/annotation_inspection/SRX17412841_annotated.h5ad"),
    Path("data/annotation_inspection/SRX17412841_cyteonto.csv"),
)

extremes = top_bottom_by_cytetype(pair_df, n=10)
write_extremes_csv(pair_df, n=10, output_path=Path("output/annotation_inspection_pipeline/extremes.csv"))
```

Extremes can be regenerated from an existing `summary.csv` without R2 access:

```python
import pandas as pd
from annotation_inspector import write_extremes_csv

summary = pd.read_csv("output/annotation_inspection_pipeline/20260603_120000/summary.csv")
write_extremes_csv(summary, n=10, output_path=Path("output/annotation_inspection_pipeline/20260603_120000/extremes.csv"))
```

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `inputPrefix` | required | R2 prefix for annotated h5ads |
| `cyteontoPrefix` | required | R2 prefix for `{srx}_cyteonto.csv` files |
| `topN` | `10` | Top/bottom STATE cell types per CyteType label |
| `downloadRoot` | `data/annotation_inspection` | Local cache for R2 downloads |
| `outputDir` | `output/annotation_inspection_pipeline` | Base output directory |
| `emitExtremes` | `True` | Write `extremes.csv` from accumulated summary |

## Outputs

| File | Description |
|------|-------------|
| `summary.csv` | Pair-level rows: accession, labels, confidence, cytescore, n_cells, report_url |
| `extremes.csv` | Top/bottom cytescore STATE types per CyteType label (optional; can be regenerated from summary) |
| `run.csv` | Per-accession pipeline status |

Disable extremes with `--no-extremes`, or rebuild extremes alone with `--from-summary`.
