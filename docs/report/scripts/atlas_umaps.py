import argparse
import colorsys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData, read_h5ad
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

_CELL_TYPE_KEY = "cell_type"
_UNCORRECTED_KEY = "X_umap_uncorrected"
_CORRECTED_KEY = "X_umap"
_DPI = 300
_ALPHA = 0.55
_N_LEGEND = 10
_LEGEND_TITLE = "Just the 10 most common cell types"


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the two-panel atlas UMAP composite.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="postprocessed atlas h5ad")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output PNG path")
    args = parser.parse_args()
    adata = load_atlas(args.input)
    validate_atlas(adata)
    labels = adata.obs[_CELL_TYPE_KEY]
    palette = palette_for(labels)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 6.0), facecolor="none")
    umap_custom(adata, _UNCORRECTED_KEY, "Uncorrected", axes[0], palette)
    umap_custom(adata, _CORRECTED_KEY, "Harmony", axes[1], palette)
    _label_panel(axes[0], "A")
    _label_panel(axes[1], "B")
    fig.subplots_adjust(wspace=0.04, left=0.02, right=0.98, top=0.98, bottom=0.02)
    _add_legend(fig, most_common_cell_types(labels), palette)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=_DPI, facecolor="none", transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {args.output}")


def umap_custom(
    adata: AnnData,
    basis: str,
    title: str,
    ax: Axes,
    palette: dict[str, tuple[float, float, float]],
) -> None:
    sc.pl.embedding(
        adata,
        basis=basis,
        color=_CELL_TYPE_KEY,
        palette=palette,
        ax=ax,
        show=False,
        title=title,
        legend_loc=None,
        frameon=False,
        alpha=_ALPHA,
        na_color="lightgray",
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("none")
    for spine in ax.spines.values():
        spine.set_visible(False)


def load_atlas(path: Path) -> AnnData:
    if not path.is_file():
        raise FileNotFoundError(f"atlas h5ad not found: {path}")
    return read_h5ad(path, backed="r")


def validate_atlas(adata: AnnData) -> None:
    if _CELL_TYPE_KEY not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing {_CELL_TYPE_KEY!r}")
    missing_obsm = [key for key in (_UNCORRECTED_KEY, _CORRECTED_KEY) if key not in adata.obsm]
    if missing_obsm:
        raise ValueError(f"adata.obsm is missing: {', '.join(missing_obsm)}")
    labels = adata.obs[_CELL_TYPE_KEY].astype("string")
    if labels.isna().all() or (labels.str.strip() == "").all():
        raise ValueError(f"{_CELL_TYPE_KEY} has no usable labels")


def _normalized_labels(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("missing").replace("", "missing")


def palette_for(series: pd.Series) -> dict[str, tuple[float, float, float]]:
    categories = sorted(_normalized_labels(series).unique().tolist())
    if not categories:
        raise ValueError("no cell_type categories to color")
    colors: dict[str, tuple[float, float, float]] = {}
    for i, name in enumerate(categories):
        hue = (i * 0.618033988749895) % 1.0
        sat = 0.78 if i % 2 == 0 else 0.95
        val = 0.95 if i % 3 else 0.82
        colors[name] = colorsys.hsv_to_rgb(hue, sat, val)
    return colors


def most_common_cell_types(series: pd.Series, n: int = _N_LEGEND) -> list[str]:
    ranked = _normalized_labels(series).value_counts().index.tolist()
    if not ranked:
        raise ValueError("no cell_type categories to color")
    return ranked[:n]


def _add_legend(
    fig: Figure,
    names: list[str],
    palette: dict[str, tuple[float, float, float]],
) -> None:
    missing = [name for name in names if name not in palette]
    if missing:
        raise ValueError(f"legend names missing from palette: {', '.join(missing)}")
    handles = [
        Line2D(
            [0],
            [0],
            linestyle="none",
            marker="o",
            markersize=8,
            markerfacecolor=palette[name],
            markeredgecolor="none",
            label=name,
        )
        for name in names
    ]
    fig.legend(
        handles=handles,
        title=_LEGEND_TITLE,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        borderaxespad=0.0,
        handletextpad=0.5,
        labelspacing=0.45,
    )


def _label_panel(ax: Axes, letter: str) -> None:
    ax.text(
        0.02,
        0.98,
        letter,
        transform=ax.transAxes,
        color="black",
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


if __name__ == "__main__":
    main()
