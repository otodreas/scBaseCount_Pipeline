import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from shared.repo import rel_to_repo

from atlas_postprocessing.config import AtlasPostprocessingConfig

log = logging.getLogger(__name__)


def make_atlas_plots(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> None:
    """Write pre- and post-correction UMAP PNGs with atlas-scale point styling."""
    cfg.figsDir.mkdir(parents=True, exist_ok=True)
    for color_by in (cfg.batchKey, cfg.cellTypeKey):
        if color_by not in adata.obs:
            log.warning("Skipping UMAP for missing obs column %s", color_by)
            continue
        _save_embedding_plot(
            adata,
            basis="X_umap_uncorrected",
            colorBy=color_by,
            outPath=cfg.figsDir / f"umap_{color_by}_uncorrected.png",
        )
        _save_embedding_plot(
            adata,
            basis="X_umap",
            colorBy=color_by,
            outPath=cfg.figsDir / f"umap_{color_by}_harmony.png",
        )


def _save_embedding_plot(
    adata: sc.AnnData,
    *,
    basis: str,
    colorBy: str,
    outPath: Path,
) -> None:
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=colorBy,
        legend_loc=None,
        show=False,
        return_fig=True,
    )
    fig.savefig(outPath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", rel_to_repo(outPath))


def save_scree_plot(adata: sc.AnnData, cfg: AtlasPostprocessingConfig) -> Path:
    """Save a PCA scree plot of per-PC and cumulative explained variance."""
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    pcs = np.arange(1, len(variance_ratio) + 1)
    cfg.figsDir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(pcs, variance_ratio * 100, color="steelblue", alpha=0.7)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)", color="steelblue")
    ax.axvline(cfg.nPcs, color="gray", linestyle="--", linewidth=1, label=f"nPcs = {cfg.nPcs}")

    ax2 = ax.twinx()
    ax2.plot(pcs, np.cumsum(variance_ratio) * 100, color="crimson", marker="o", ms=3)
    ax2.set_ylabel("Cumulative explained variance (%)", color="crimson")

    ax.set_title("PCA scree plot")
    ax.legend(loc="center right", fontsize=8)
    fig.tight_layout()

    scree_path = cfg.figsDir / "pca_scree.png"
    fig.savefig(scree_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", rel_to_repo(scree_path))
    return scree_path


def plot_sweep_metric(
    *,
    values: list[float],
    scores: list[float],
    plateaus: list[tuple[float, float]],
    xlabel: str,
    ylabel: str,
    title: str,
    outPath: Path,
    baseline: float | None = None,
) -> Path:
    """Plot one metric against candidate values and shade plateau intervals."""
    outPath.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(values, scores, marker="o", color="steelblue")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    y_min, y_max = ax.get_ylim()
    for start, end in plateaus:
        ax.axvspan(start, end, color="green", alpha=0.15, label="plateau" if start == plateaus[0][0] else None)
    if baseline is not None:
        ax.axvline(baseline, color="gray", linestyle="--", linewidth=1, label=f"baseline={baseline}")

    ax.set_ylim(y_min, y_max)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(outPath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", rel_to_repo(outPath))
    return outPath
