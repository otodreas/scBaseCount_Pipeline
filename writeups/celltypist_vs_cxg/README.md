# CellTypist vs CyteType agreement with CELLxGENE labels

This writeup compares how well CellTypist and CyteType annotations agree with CELLxGENE author labels (`cell_type`, referred to as cxg) using CyteOnto cytescore.

The analysis is driven by [`agreement_inspection.ipynb`](agreement_inspection.ipynb). Intermediate outputs land under `output/celltypist_vs_cxg/`.

## Pipeline per accession

1. CellTypist inference writes `predicted_labels` to `obs`.
2. Cluster validation uses `predicted_labels` as the weak prior and produces `leiden_merged`.
3. CyteType annotates clusters and writes `cytetype_annotation_leiden_merged`.
4. CyteOnto scores both algorithms against `cell_type` in one deduplicated call.

## Comparisons

| Comparison | Algorithm key | obs column |
|------------|---------------|------------|
| CellTypist vs cxg | `celltypist` | `predicted_labels` |
| CyteType vs cxg | `cytetype` | `cytetype_annotation_leiden_merged` |

Figures in the notebook summarize per-cxg-cell-type cytescore by method and the paired delta (CellTypist minus CyteType).
