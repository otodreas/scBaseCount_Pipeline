import logging
from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from shared.repo import rel_to_repo

from atlas_postprocessing.config import AtlasPostprocessingConfig

log = logging.getLogger(__name__)

Workflow = Literal["production", "validation"]


def make_atlas_plots(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    workflow: Workflow = "production",
) -> None:
    """Write UMAP PNGs for the active workflow with atlas-scale point styling."""
    cfg.figsDir.mkdir(parents=True, exist_ok=True)
    for color_by in (cfg.batchKey, cfg.cellTypeKey):
        if color_by not in adata.obs:
            log.warning("Skipping UMAP for missing obs column %s", color_by)
            continue
        if workflow == "validation" and "X_umap_uncorrected" in adata.obsm:
            _save_embedding_plot(
                adata,
                basis="X_umap_uncorrected",
                colorBy=color_by,
                outPath=cfg.figsDir / f"umap_{color_by}_uncorrected.png",
            )
        if "X_umap" in adata.obsm:
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


def plot_resolution_selection(
    *,
    resolutions: list[float],
    jaccArr: list[float],
    kArr: list[int],
    selectedResolution: float,
    outPath: Path,
) -> Path:
    """Plot cluster counts and matched-Jaccard scores with the advisory argmax marked."""
    if len(resolutions) != len(jaccArr) or len(resolutions) != len(kArr):
        raise ValueError("resolutions, jaccArr, and kArr must have the same length")
    outPath.parent.mkdir(parents=True, exist_ok=True)
    best_idx = resolutions.index(selectedResolution)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(resolutions, kArr, marker="o", ms=4, color="steelblue")
    axes[0].axvline(selectedResolution, color="red", linestyle="--", label=f"argmax = {selectedResolution}")
    axes[0].set_xlabel("Resolution")
    axes[0].set_ylabel("# clusters")
    axes[0].set_title("Clusters per resolution")
    axes[0].legend(fontsize=8)

    axes[1].plot(resolutions, jaccArr, marker="o", ms=4, color="darkorange", label="matched Jaccard")
    axes[1].axvline(selectedResolution, color="red", linestyle="--", label=f"argmax = {selectedResolution}")
    axes[1].scatter([selectedResolution], [jaccArr[best_idx]], color="red", zorder=5, s=60)
    axes[1].set_xlabel("Resolution")
    axes[1].set_ylabel("Matched Jaccard")
    axes[1].set_title("Resolution selection (advisory)")
    axes[1].legend(fontsize=8)

    fig.suptitle(
        f"Atlas resolution selection\nselected = {selectedResolution}  "
        f"(k = {kArr[best_idx]}, jaccard = {jaccArr[best_idx]:.4f})"
    )
    fig.tight_layout()
    fig.savefig(outPath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", rel_to_repo(outPath))
    return outPath
