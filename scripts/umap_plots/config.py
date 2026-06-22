from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel
from shared.repo import REPO_ROOT as _REPO_ROOT


class UmapPlotConfig(BaseModel):
    umapKey: str = "X_umap"
    embeddingKey: str = "X_pca"
    figsDir: Path = _REPO_ROOT / "output" / "umap_plots" / "figs"
    figSize: tuple[float, float] = (8.0, 6.0)
    dpi: int = 150
    pointSize: float = 3.0
    alpha: float = 0.8
    continuousCmap: str = "viridis"
    categoricalPalette: str = "tab20"
    categoricalMaxUnique: int = 20
    maxCategoriesForLegend: int = 30
    legendLoc: str = "best"
    legendFontsize: float = 8.0
    titleFontsize: float = 12.0
