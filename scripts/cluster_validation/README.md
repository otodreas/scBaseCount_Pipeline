# cluster_validation

Standalone pipeline that takes a raw scRNA-seq `AnnData` object and returns a biologically-informed cluster partition. Resolution is selected by matching Leiden clusters to a `cell_type` weak prior using the Jaccard index and SciPy's linear sum assignment. Over-clustered partitions are reduced by merging indistinguishable clusters using a random forest out-of-fold confusion step.

## Usage

```python
from cluster_validation import ClusterValidationConfig, run_cluster_validation
from cluster_validation.viz import plot_all
from pathlib import Path

cfg = ClusterValidationConfig(
    srxAccession="SRX17412841",
    localH5adRoot=Path("data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens"),
)

adata, result = run_cluster_validation(cfg)
plot_all(adata, result, figs_dir=cfg.figsDir / cfg.srxAccession)
```

The pipeline writes the final `AnnData` to `output/clustering/data/{srx}_clustered.h5ad` and figures under `output/clustering/figs/{srx}/`. Returns a `ClusterValidationResult` with all per-resolution metrics.

## Pipeline steps

| Step | Module | Description |
|------|--------|-------------|
| Load | `data.py` | Read h5ad from local path or GCS fallback; look up dataset row from catalog CSV |
| Preprocess | `preprocess.py` | Filter rare cell types (`minCellsPerType`), QC, HVG selection, normalisation |
| Embed | `embedding.py` | PCA, select PCs by cumulative variance target, neighbors graph, UMAP |
| Sweep | `clustering.py` | Leiden clustering at each resolution in `resolutions`; one `obs` column per resolution |
| Select resolution | `resolution.py` | Jaccard matrix + SciPy linear sum assignment; pick resolution maximising matched Jaccard sum |
| Merge | `merge.py` | RF OOF confusion on HVG matrix; union-find merges pairs above `mergeThreshold`; writes `leiden_merged` |
| Metrics | `metrics.py` | Silhouette, homogeneity, completeness, NMI, V-score, ARI across all resolutions |
| Cell type metrics | `cell_type_metrics.py` | Normalized Shannon entropy and KL divergence per cell type across datasets |

## Methods

### Resolution selection

For each resolution in the sweep:

1. Build a contingency table between Leiden clusters and `cell_type` reference labels.
2. Convert each cluster-label pair to a Jaccard index (also known as intersection over union, or IoU): `J[i, j] = intersection / union`.
3. Run `scipy.optimize.linear_sum_assignment` on `-J` to find the optimal one-to-one assignment.
4. The penalised score is the sum of matched Jaccard values. The resolution that maximises this score is selected.

#### Computing the Jaccard index

For a Leiden cluster \(C\) and a weak-prior cell type \(T\):

$$
J(C,T)
= \frac{|C \cap T|}{|C \cup T|}
= \frac{|C \cap T|}{|C| + |T| - |C \cap T|}
$$

Suppose cluster 0 contains 90 cells, the macrophage label contains 85 cells, and they overlap by 80 cells. Their Jaccard index is:

$$
\frac{80}{90 + 85 - 80} = 0.84
$$

The same calculation is applied to every cluster and cell-type pair. Each matrix entry shows the overlap cell count followed by its Jaccard index in parentheses:

| Leiden cluster | Macrophage | T cell | B cell | Cluster total |
|---|---:|---:|---:|---:|
| cluster 0 | 80 (0.84) | 10 (0.06) | 0 (0.00) | 90 |
| cluster 1 | 5 (0.03) | 70 (0.67) | 5 (0.04) | 80 |
| cluster 2 | 0 (0.00) | 15 (0.10) | 60 (0.75) | 75 |
| Type total | 85 | 95 | 65 | 245 |

SciPy's [`linear_sum_assignment`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html) minimizes the total cost of a one-to-one assignment. Here the Jaccard matrix is negated before it is passed to SciPy, so minimizing the negative values is equivalent to maximizing the sum of the matched Jaccard indices. The returned row and column indices identify the optimal cluster-to-cell-type pairs, whose original positive Jaccard indices are then summed.

The selected one-to-one assignment is highlighted in the same matrix:

| Leiden cluster | Macrophage | T cell | B cell | Cluster total |
|---|---:|---:|---:|---:|
| cluster 0 | **80 (0.84)** | 10 (0.06) | 0 (0.00) | 90 |
| cluster 1 | 5 (0.03) | **70 (0.67)** | 5 (0.04) | 80 |
| cluster 2 | 0 (0.00) | 15 (0.10) | **60 (0.75)** | 75 |
| Type total | 85 | 95 | 65 | 245 |

Bold cells show the optimal one-to-one assignment. The matched-Jaccard score is the sum of those assignments:

$$
0.84 + 0.67 + 0.75 = 2.26
$$

The score can exceed 1 because it sums several Jaccard indices.

### RF cluster merging

A `RandomForestClassifier` is trained on HVG expression with stratified K-fold out-of-fold cross-validation. Pairs of clusters whose row-normalised OOF confusion exceeds `mergeThreshold` are merged. A union-find structure propagates merges transitively, so if A is confused with B and B with C, all three collapse into one cluster written as `leiden_merged`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `srxAccession` | `None` | Select dataset by SRX/ERX accession string |
| `datasetIndex` | `2` | Select dataset by row index in catalog (used when `srxAccession` is `None`) |
| `summaryPath` | `tests/quantiles_datasets.csv` | Path to catalog CSV with `srx_accession`, `file_path`, `obs_count` columns (and optional `quantile`) |
| `localH5adRoot` | `data/scbasecount/...` | Directory of local h5ad files; takes priority over `file_path` in catalog |
| `weakPriorKey` | `cell_type` | `obs` column used as the weak prior for resolution selection, type filtering, and optional RF balancing |
| `runLabel` | `None` | Optional suffix for output filenames and figure directories; defaults to the SRX accession |
| `outputDir` | `output/clustering/data` | Directory where `{run_tag}_clustered.h5ad` is written |
| `figsDir` | `output/clustering/figs` | Base directory for per-SRX figure folders (`figsDir/{srx}/`) |
| `minCellsPerType` | `20` | Minimum cells per weak-prior label; rarer types are dropped before clustering |
| `nTopGenes` | `2000` | Number of highly variable genes |
| `nPcsCompute` | `50` | Number of PCs computed |
| `nPcsMin` | `15` | Minimum PCs to use in neighbor graph regardless of variance target |
| `nPcsCumvarTarget` | `0.5` | Cumulative variance floor for PC selection |
| `resolutions` | `0.1, 0.2, ..., 1.9` | Leiden resolutions swept |
| `mergeThreshold` | `0.2` | OOF confusion threshold above which two clusters are merged |
| `rfBalanceWeakPrior` | `False` | Balance class weights in the RF by weak-prior label frequency |

Run on an in-memory `AnnData` (for example after adding CellTypist predictions to `obs`):

```python
from cluster_validation import ClusterValidationConfig, run_cluster_validation_on_adata

cfg = ClusterValidationConfig(
    weakPriorKey="predicted_labels",
    runLabel="SRX12366723_predicted_labels",
    outputDir=Path("tmp/clustering/data"),
    figsDir=Path("tmp/clustering/figs"),
)
adata, result = run_cluster_validation_on_adata(adata, cfg, srx="SRX12366723")
```

All default paths are relative to the repo root.

## Cell type metrics

`cell_type_metrics.py` exposes three functions used by the annotation pipeline after clustering is complete.

### `compute_nse_kld_row(adata, merged_key)`

Computes two scalars per cell type from a single dataset:

- **Normalized Shannon entropy (NSE)**: how fragmented the cell type is across Leiden clusters. `0` = all cells in one cluster; `1` = maximally spread across all clusters it occupies.
- **KL divergence (KLD)**: `KL(p || q)` where `p` is the cluster distribution of the cell type and `q` is the global cluster distribution of all cells. A high value means the cell type is concentrated in clusters that differ from the background, indicating coherence. Near `0` means the cell type mirrors the global distribution.

Returns `(nse_row, kld_row)`, each a `dict[str, float]` keyed by cell type name.

### `build_metric_dataframes(rows)`

Takes a list of dicts parsed from `metrics_matrix.jsonl` (one per accession, with keys `srx`, `nse`, `kld`) and returns three DataFrames:

- `nse_df`: accessions x cell types, normalized Shannon entropy values.
- `kld_df`: accessions x cell types, KL divergence values.
- `summary_df`: cell types x metrics, with columns `n_datasets`, `normalized_shannon_entropy_mean`, `kl_divergence_mean`.

### `save_metric_plot(summary_df, output_path, ...)`

Saves a two-panel horizontal bar chart to `output_path` (PNG). Left panel shows mean NSE per cell type; right panel shows mean KLD. Cell types are sorted by NSE ascending.

## Output model

```
ClusterValidationResult
├── srxAccession            str
├── runTag                  str          label used for output paths (runLabel or SRX)
├── weakPriorKey            str          obs column used as the weak prior
├── selectedResolution      float
├── clusterKey              str          obs column for the selected pre-merge partition
├── nPcs                    int
├── cumvar                  float        cumulative variance at nPcs
├── kPrior                  int          cells before QC filter
├── kFiltered               int          cells after QC filter
├── nCellsDropped           int
├── nClustersPreMerge       int
├── nClustersPostMerge      int
├── adataPath               Path
├── labelMap                dict[str, str]   merged cluster label -> representative original label
├── mergedGroups            dict[str, list[str]]
├── resolutions             list[float]
├── kArr                    list[int]    cluster count at each resolution
├── jaccArr                 list[float]  matched Jaccard score at each resolution
├── silhouetteArr           list[list]   [[resolution, value], ...]
├── homogeneityArr          list[list]
├── completenessArr         list[list]
├── nmiArr                  list[list]
├── vscoreArr               list[list]
├── ariArr                  list[list]
├── confMatrix              list[list]   RF OOF confusion matrix (pre-merge clusters)
└── confClasses             list[str]
```
