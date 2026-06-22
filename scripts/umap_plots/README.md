# umap_plots

Plot UMAPs from AnnData objects with matplotlib scatter. Expects precomputed embeddings in `adata.obsm` and colors by any `adata.obs` column, auto-detecting categorical vs continuous values.

## Usage

```python
import scanpy as sc
from umap_plots import UmapPlotConfig, plot_umap

adata = sc.read_h5ad("output/clustering/data/SRX17412841_clustered.h5ad")

fig = plot_umap(adata, "cell_type")
fig = plot_umap(adata, "leiden", nameSuffix="SRX17412841")

cfg = UmapPlotConfig(pointSize=5.0, figSize=(10.0, 8.0))
fig = plot_umap(adata, "cell_type", cfg=cfg)
```

Edit defaults in `umap_plots/config.py` or pass a custom `UmapPlotConfig`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `umapKey` | `"X_umap"` | Required key in `adata.obsm` |
| `embeddingKey` | `"X_pca"` | Required key in `adata.obsm` |
| `figsDir` | `output/umap_plots/figs` | Directory for saved PNGs |
| `figSize` | `(8.0, 6.0)` | Figure size in inches |
| `dpi` | `150` | Saved PNG resolution |
| `pointSize` | `3.0` | Scatter marker size |
| `alpha` | `0.8` | Marker opacity |
| `continuousCmap` | `"viridis"` | Colormap for continuous obs columns |
| `categoricalPalette` | `"tab20"` | Matplotlib colormap for categorical obs columns |
| `categoricalMaxUnique` | `20` | Treat numeric obs as categorical when unique count is at or below this |
| `maxCategoriesForLegend` | `30` | Skip legend when category count exceeds this |
| `legendLoc` | `"best"` | Matplotlib legend location |
| `legendFontsize` | `8.0` | Legend font size |
| `titleFontsize` | `12.0` | Title font size |

## Outputs

| File | Description |
|------|-------------|
| `output/umap_plots/figs/umap_{colorBy}.png` | UMAP colored by the given obs column |
| `output/umap_plots/figs/umap_{colorBy}_{nameSuffix}.png` | Same plot with an optional suffix |

`plot_umap` also returns the matplotlib `Figure`.
