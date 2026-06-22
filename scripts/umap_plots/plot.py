from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from umap_plots.config import UmapPlotConfig

_log = logging.getLogger(__name__)


def plot_umap(
    adata: sc.AnnData,
    colorBy: str,
    cfg: UmapPlotConfig | None = None,
    nameSuffix: str | None = None,
) -> Figure:
    """Plot a UMAP colored by an obs column, save a PNG, and return the figure."""
    resolved_cfg = cfg or UmapPlotConfig()
    _validate_adata(adata, colorBy, resolved_cfg)

    coords = adata.obsm[resolved_cfg.umapKey][:, :2]
    series = adata.obs[colorBy]

    fig, ax = plt.subplots(figsize=resolved_cfg.figSize)
    if _is_categorical(series, resolved_cfg):
        _plot_categorical(ax, coords, series, colorBy, resolved_cfg)
    else:
        _plot_continuous(ax, coords, series, colorBy, resolved_cfg)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(colorBy, fontsize=resolved_cfg.titleFontsize)
    fig.tight_layout()

    output_path = _output_path(resolved_cfg.figsDir, colorBy, nameSuffix)
    resolved_cfg.figsDir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=resolved_cfg.dpi, bbox_inches="tight")
    _log.info("Saved UMAP plot to %s", output_path)
    return fig


def _validate_adata(adata: sc.AnnData, colorBy: str, cfg: UmapPlotConfig) -> None:
    missing_obsm = [key for key in (cfg.umapKey, cfg.embeddingKey) if key not in adata.obsm]
    if missing_obsm:
        raise ValueError(f"adata.obsm is missing required keys: {', '.join(missing_obsm)}")

    if colorBy not in adata.obs:
        raise ValueError(f"adata.obs is missing required column: {colorBy}")


def _is_categorical(series: pd.Series, cfg: UmapPlotConfig) -> bool:
    if pd.api.types.is_categorical_dtype(series) or pd.api.types.is_object_dtype(series):
        return True
    if pd.api.types.is_bool_dtype(series):
        return True
    if pd.api.types.is_numeric_dtype(series):
        return series.nunique(dropna=False) <= cfg.categoricalMaxUnique
    return True


def _plot_categorical(
    ax: Axes,
    coords: np.ndarray,
    series: pd.Series,
    colorBy: str,
    cfg: UmapPlotConfig,
) -> None:
    categories = series.astype(str).fillna("nan").unique().tolist()
    cmap = plt.get_cmap(cfg.categoricalPalette, max(len(categories), 1))
    colors = [cmap(i % cmap.N) for i in range(len(categories))]

    for category, color in zip(categories, colors, strict=True):
        mask = series.astype(str).fillna("nan") == category
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=cfg.pointSize,
            alpha=cfg.alpha,
            c=[color],
            label=category,
            linewidths=0,
        )

    if len(categories) <= cfg.maxCategoriesForLegend:
        ax.legend(
            loc=cfg.legendLoc,
            fontsize=cfg.legendFontsize,
            markerscale=2,
            frameon=False,
        )
    else:
        _log.warning(
            "Skipping legend for %s: %d categories exceeds maxCategoriesForLegend=%d",
            colorBy,
            len(categories),
            cfg.maxCategoriesForLegend,
        )


def _plot_continuous(
    ax: Axes,
    coords: np.ndarray,
    series: pd.Series,
    colorBy: str,
    cfg: UmapPlotConfig,
) -> None:
    values = pd.to_numeric(series, errors="coerce")
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=values,
        s=cfg.pointSize,
        alpha=cfg.alpha,
        cmap=cfg.continuousCmap,
        linewidths=0,
    )
    fig = ax.figure
    fig.colorbar(scatter, ax=ax, label=colorBy)


def _output_path(figs_dir: Path, colorBy: str, nameSuffix: str | None) -> Path:
    stem = f"umap_{colorBy}"
    if nameSuffix:
        stem = f"{stem}_{nameSuffix}"
    return figs_dir / f"{stem}.png"
