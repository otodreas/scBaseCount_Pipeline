# CellTypist vs CyteType agreement with CELLxGENE labels

This writeup compares how well cxg (`cell_type`) and CyteType annotations agree with CellTypist `predicted_labels`, using CyteOnto cytescore with CellTypist as the author reference.

## Run

Run the upstream pipeline on the server (requires `CYTETYPE_API_KEY` and CyteOnto access):

```bash
uv run python writeups/celltypist_vs_cxg/run_pipeline.py --srx SRX17412841
```

Multiple accessions:

```bash
uv run python writeups/celltypist_vs_cxg/run_pipeline.py --srx SRX17412841 SRX12345678
```

Recompute from scratch:

```bash
uv run python writeups/celltypist_vs_cxg/run_pipeline.py --srx SRX17412841 --force
```

Then open [`agreement_inspection.ipynb`](agreement_inspection.ipynb) to load cached outputs and generate figures.

## Pipeline per accession

1. CellTypist inference writes `predicted_labels` to `obs`.
2. Cluster validation uses `predicted_labels` as the weak prior and produces `leiden_merged`.
3. CyteType annotates clusters and writes `cytetype_annotation_leiden_merged`.
4. CyteOnto scores cxg and CyteType against `predicted_labels` in one deduplicated call.

## Cached outputs

All outputs live under `output/celltypist_vs_cxg/`:

| File | Description |
|------|-------------|
| `data/{srx}_celltypist_prior_clustered.h5ad` | Clustered h5ad after CellTypist weak prior |
| `data/{srx}_cytetype_annotated.h5ad` | Final annotated h5ad (notebook input) |
| `cyteonto_results/{srx}_cyteonto.csv` | CyteOnto cytescore results (notebook input) |
| `runs/{timestamp}/run.csv` | Per-run status log |
| `runs/{timestamp}/metadata.json` | Run configuration snapshot |

Reruns skip steps whose cached files already exist. Use `--force` to delete cached files and recompute.

## Comparisons

| Comparison | Algorithm key | obs column |
|------------|---------------|------------|
| cxg vs CellTypist | `cxg` | `cell_type` |
| CyteType vs CellTypist | `cytetype` | `cytetype_annotation_leiden_merged` |

Figures in the notebook summarize per-CellTypist-label cytescore by method and the paired delta (cxg minus CyteType).
