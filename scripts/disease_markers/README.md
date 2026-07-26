# disease_markers

Map Harmony Leiden clusters onto the full-gene atlas, annotate cells with coarse disease areas, and run per-cluster one-vs-rest pseudobulk differential expression.

## Usage

```python
from disease_markers import DiseaseMarkersConfig, run_disease_markers

cfg = DiseaseMarkersConfig(
    inputAtlasH5ad=Path("output/atlas/v1/atlas.h5ad"),
    harmonyAtlasH5ad=Path("output/atlas/v1/processed_1/atlas_harmony.h5ad"),
)
run_disease_markers(cfg)
```

Cluster transfer only (deterministic, no disease labels):

```python
from disease_markers.transfer import load_full_atlas_transfer_clusters

adata = load_full_atlas_transfer_clusters(cfg.inputAtlasH5ad, cfg.harmonyAtlasH5ad)
```

## Outputs

Written under `output/atlas/v1/processed_1/disease_markers/`:

| File | Description |
|------|-------------|
| `eligibility_labels.csv` | Per-SRX disease area, control flag, eligibility |
| `area_cluster_counts.csv` | Sample and study counts per cluster and disease area |
| `markers/<cluster>__<area>.csv` | Up-regulated genes for one-vs-rest tests |
| `summary.json` | Run config and high-level counts |

Optional intermediate: `output/atlas/v1/processed_1/atlas_with_clusters.h5ad` when `writeTransferredAtlas` is true (before eligibility subsetting).

## Config

See `DiseaseMarkersConfig` in `config.py`. Key fields: `clusterKey` (`leiden_atlas`), `sampleKey` (`SRX_accession`), `studyKey` (`study_accession`), `minCellsPerProfile`, `minSamplesPerArea`, `minStudiesPerArea`, `padjThreshold`, `lfcThreshold`.
