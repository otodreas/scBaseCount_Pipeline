# cluster_validation

Standalone pipeline that takes a raw scRNA-seq `AnnData` object and returns a biologically-informed cluster partition. Resolution is selected by matching Leiden clusters to a `cell_type` weak prior using Jaccard similarity and the Hungarian algorithm. Over-clustered partitions are reduced by merging indistinguishable clusters using a random forest out-of-fold confusion step.

## Usage

```python
from cluster_validation import ClusterValidationConfig, run_cluster_validation
from cluster_validation.viz import plot_all
from pathlib import Path

cfg = ClusterValidationConfig(
    srxAccession="SRX17412841",
    localH5adRoot=Path("data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens"),
    outputDir=Path("data/clustered"),
)

adata, result = run_cluster_validation(cfg)
plot_all(adata, result, figs_dir=Path("output/figs") / cfg.srxAccession)
```

The pipeline writes the final `AnnData` to `outputDir/final_adata_{srx}.h5ad` and returns a `ClusterValidationResult` with all per-resolution metrics.

## Pipeline steps

| Step | Module | Description |
|------|--------|-------------|
| Load | `data.py` | Read h5ad from local path or GCS fallback; look up dataset row from catalog CSV |
| Preprocess | `preprocess.py` | Filter rare cell types (`minCellsPerType`), QC, HVG selection, normalisation |
| Embed | `embedding.py` | PCA, select PCs by cumulative variance target, neighbors graph, UMAP |
| Sweep | `clustering.py` | Leiden clustering at each resolution in `resolutions`; one `obs` column per resolution |
| Select resolution | `resolution.py` | Jaccard matrix + Hungarian assignment; pick resolution maximising matched Jaccard sum |
| Merge | `merge.py` | RF OOF confusion on HVG matrix; union-find merges pairs above `mergeThreshold`; writes `leiden_merged` |
| Metrics | `metrics.py` | Silhouette, homogeneity, completeness, NMI, V-score, ARI across all resolutions |

## Methods

### Resolution selection

For each resolution in the sweep:

1. Build a contingency table between Leiden clusters and `cell_type` reference labels.
2. Convert each cell to a Jaccard (IoU) value: `J[i, j] = intersection / union`.
3. Run the Hungarian algorithm (`scipy.optimize.linear_sum_assignment` on `-J`) to find the optimal one-to-one assignment.
4. The penalised score is the sum of matched Jaccard values. The resolution that maximises this score is selected.

### RF cluster merging

A `RandomForestClassifier` is trained on HVG expression with stratified K-fold out-of-fold cross-validation. Pairs of clusters whose row-normalised OOF confusion exceeds `mergeThreshold` are merged. A union-find structure propagates merges transitively, so if A is confused with B and B with C, all three collapse into one cluster written as `leiden_merged`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `srxAccession` | `None` | Select dataset by SRX/ERX accession string |
| `datasetIndex` | `2` | Select dataset by row index in catalog (used when `srxAccession` is `None`) |
| `summaryPath` | `cluster_validation_sandbox/datasets_summary` | Path to catalog CSV with `srx_accession`, `file_path`, `quantile` columns |
| `localH5adRoot` | `data/scbasecount/...` | Directory of local h5ad files; takes priority over `file_path` in catalog |
| `outputDir` | `data/other` | Directory where `final_adata_{srx}.h5ad` is written |
| `minCellsPerType` | `20` | Minimum cells per `cell_type` label; rarer types are dropped before clustering |
| `nTopGenes` | `2000` | Number of highly variable genes |
| `nPcsCompute` | `50` | Number of PCs computed |
| `nPcsMin` | `15` | Minimum PCs to use in neighbor graph regardless of variance target |
| `nPcsCumvarTarget` | `0.5` | Cumulative variance floor for PC selection |
| `resolutions` | `0.1, 0.2, ..., 1.9` | Leiden resolutions swept |
| `mergeThreshold` | `0.2` | OOF confusion threshold above which two clusters are merged |
| `rfBalanceWeakPrior` | `False` | Balance class weights in the RF by `cell_type` frequency |

## Output model

```
ClusterValidationResult
├── srxAccession            str
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
├── jaccArr                 list[float]  Hungarian Jaccard score at each resolution
├── silhouetteArr           list[list]   [[resolution, value], ...]
├── homogeneityArr          list[list]
├── completenessArr         list[list]
├── nmiArr                  list[list]
├── vscoreArr               list[list]
├── ariArr                  list[list]
├── confMatrix              list[list]   RF OOF confusion matrix (pre-merge clusters)
└── confClasses             list[str]
```
