from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scanpy as sc

from cluster_validation.models import ClusterValidationResult


def _save(fig: plt.Figure, figs_dir: Path | None, name: str) -> None:
    if figs_dir is not None:
        figs_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figs_dir / name, bbox_inches="tight")


def plot_pca_cumvar(
    adata: sc.AnnData,
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    var_ratio = adata.uns["pca"]["variance_ratio"]
    n_pcs_min = result.nPcs  # floor used
    n_pcs = result.nPcs
    cumvar_target = adata.uns.get("_cv_cumvar_target", None)

    fig, ax = plt.subplots()
    ax.set_title(f"Cumulative Variance by PCs\nDataset: {result.datasetTitleSuffix}")
    ax.plot(
        range(1, len(var_ratio) + 1),
        np.cumsum(var_ratio) * 100,
    )
    ax.axvline(n_pcs, color="crimson", linestyle="--", linewidth=1, label=f"{n_pcs} PCs selected")
    ax.set_xlabel("Number of PCs")
    ax.set_ylabel("Cumulative variance (%)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, figs_dir, f"pca_cumvar_{result.srxAccession}.png")
    return fig


def plot_resolution_sweep(
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    resolutions = result.resolutions
    k_arr = result.kArr
    jacc_arr = result.jaccArr
    sel = result.selectedResolution
    best_idx = resolutions.index(sel)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    axes[0].plot(resolutions, k_arr, marker="o", ms=4, color="steelblue")
    axes[0].axhline(
        result.kFiltered, color="gray", linestyle="--", linewidth=1,
        label=f"k_filtered = {result.kFiltered}",
    )
    axes[0].axvline(sel, color="red", linestyle="--", label=f"selected = {sel}")
    axes[0].set_xlabel("Resolution")
    axes[0].set_ylabel("# clusters")
    axes[0].set_title("Clusters per resolution")
    axes[0].legend(fontsize=8)

    axes[1].plot(resolutions, jacc_arr, marker="o", ms=4, color="darkorange", label="Jaccard")
    axes[1].axvline(sel, color="red", linestyle="--", label=f"argmax = {sel}")
    axes[1].scatter([sel], [jacc_arr[best_idx]], color="red", zorder=5, s=60)
    axes[1].set_xlabel("Resolution")
    axes[1].set_ylabel("Jaccard score")
    axes[1].set_title("Jaccard score (Hungarian matched)")
    axes[1].legend(fontsize=8)

    plt.suptitle(
        f"Dataset: {result.datasetTitleSuffix}\nResolution selection\n"
        f"selected = {sel}  (k = {k_arr[best_idx]},  k_filtered = {result.kFiltered})"
    )
    plt.tight_layout()
    _save(fig, figs_dir, f"resolution_sweep_{result.srxAccession}.png")
    return fig


def plot_umap_selected(
    adata: sc.AnnData,
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> None:
    n_clusters = result.nClustersPreMerge
    sc.pl.umap(
        adata,
        color=[result.clusterKey, "cell_type"],
        ncols=2,
        legend_loc="on data",
        legend_fontsize=7,
        wspace=0.4,
        title=[
            f"Leiden {result.selectedResolution}  ({n_clusters} clusters)",
            f"cell_type filtered  ({result.kFiltered} types)",
        ],
        show=False,
    )
    if figs_dir is not None:
        figs_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(figs_dir / f"umap_selected_{result.srxAccession}.png", bbox_inches="tight")


def plot_rf_confusion(
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    conf = np.array(result.confMatrix)
    classes = result.confClasses
    fig, ax = plt.subplots(
        figsize=(max(8, len(classes) * 0.5), max(6, len(classes) * 0.45))
    )
    sns.heatmap(
        conf,
        xticklabels=classes,
        yticklabels=classes,
        cmap="Blues",
        vmin=0,
        vmax=1,
        ax=ax,
        fmt=".2f",
        linewidths=0.3,
    )
    ax.set_xlabel("Predicted cluster")
    ax.set_ylabel("True cluster")
    ax.set_title(
        f"Dataset: {result.datasetTitleSuffix}\n"
        f"RF out-of-fold confusion  ({result.clusterKey}, {len(classes)} clusters)"
    )
    plt.tight_layout()
    _save(fig, figs_dir, f"rf_confusion_{result.srxAccession}.png")
    return fig


def plot_umap_merged(
    adata: sc.AnnData,
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> None:
    merged_groups = result.mergedGroups
    multi = {m: g for m, g in merged_groups.items() if len(g) > 1}
    if multi:
        merged_str = "; ".join(
            f"[{', '.join(g)}] -> {m}" for m, g in multi.items()
        )
        umap_merged_title = (
            f"RF-merged  ({result.nClustersPostMerge} clusters)\n"
            f"Merged groups: {merged_str}"
        )
    else:
        umap_merged_title = f"RF-merged  ({result.nClustersPostMerge} clusters)\n(no merges)"

    fig = sc.pl.umap(
        adata,
        color=[result.clusterKey, "leiden_merged", "cell_type"],
        ncols=2,
        legend_loc="on data",
        legend_fontsize=7,
        wspace=0.4,
        title=[
            f"Leiden {result.selectedResolution}  ({result.nClustersPreMerge} clusters)",
            umap_merged_title,
            f"cell_type filtered  ({result.kFiltered} types)",
        ],
        show=False,
        return_fig=True,
    )
    if fig is not None:
        fig.suptitle(
            f"Dataset: {result.datasetTitleSuffix}\nRF-based cluster merging summary",
            fontsize=16,
            y=1.03,
        )
    if figs_dir is not None:
        figs_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(figs_dir / f"umap_merged_{result.srxAccession}.png", bbox_inches="tight")


def plot_composition_bars(
    adata: sc.AnnData,
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    merged_counts = (
        adata.obs["leiden_merged"].value_counts(normalize=True).sort_values(ascending=False) * 100
    )
    celltype_counts = (
        adata.obs["cell_type"].value_counts(normalize=True).sort_values(ascending=False) * 100
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    merged_counts.plot(kind="bar", ax=axes[0], color="steelblue")
    axes[0].set_title("Leiden Merged Cluster Proportions")
    axes[0].set_ylabel("Proportion of cells (%)")
    axes[0].set_xlabel("")
    axes[0].set_xticks(range(len(merged_counts)))
    axes[0].set_xticklabels([str(lbl)[:4] for lbl in merged_counts.index], rotation=90)

    celltype_counts.plot(kind="bar", ax=axes[1], color="indianred")
    axes[1].set_title("Cell Type Proportions")
    axes[1].set_xlabel("")
    axes[1].set_xticks(range(len(celltype_counts)))
    axes[1].set_xticklabels([str(lbl)[:4] for lbl in celltype_counts.index], rotation=90)

    fig.suptitle(
        f"Dataset: {result.datasetTitleSuffix}\nCluster composition", fontsize=16, y=1.03
    )
    plt.tight_layout()
    _save(fig, figs_dir, f"composition_bars_{result.srxAccession}.png")
    return fig


def plot_silhouette(
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    sil = np.array(result.silhouetteArr)
    ari = np.array(result.ariArr)
    sel = result.selectedResolution

    fig, ax = plt.subplots()
    ax.plot(sil[:, 0], sil[:, 1], marker="o", color="steelblue")
    ax.axvline(sel, color="red", linestyle="--", linewidth=1, label=f"selected = {sel}")

    valid_ari = ari[ari[:, 0] < sel] if sel > min(result.resolutions) and len(ari) > 0 else np.array([])
    if len(valid_ari) > 0:
        max_ari_res = valid_ari[valid_ari[:, 1].argmax(), 0]
        ax.axvline(
            max_ari_res,
            color="green",
            linestyle="--",
            linewidth=1,
            label=f"max ARI below selected = {max_ari_res}",
        )

    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Silhouette score")
    ax.set_title("Silhouette vs resolution (PCA embedding)")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    _save(fig, figs_dir, f"silhouette_vs_resolution_{result.srxAccession}.png")
    return fig


def plot_metrics(
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> plt.Figure:
    metrics = [
        ("Homogeneity", np.array(result.homogeneityArr)),
        ("Completeness", np.array(result.completenessArr)),
        ("V-measure", np.array(result.vscoreArr)),
        ("NMI", np.array(result.nmiArr)),
        ("ARI", np.array(result.ariArr)),
    ]
    sel = result.selectedResolution

    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    axs = axs.flatten()

    for idx, (label, arr) in enumerate(metrics):
        if arr.size > 0:
            axs[idx].plot(arr[:, 0], arr[:, 1], marker="o")
            axs[idx].axvline(
                sel, color="red", linestyle="--", linewidth=1, label=f"selected = {sel}"
            )
            axs[idx].set_title(label)
            axs[idx].set_xlabel("Leiden resolution")
            axs[idx].set_ylabel(label)
            axs[idx].legend(fontsize=8)
        else:
            axs[idx].set_visible(False)

    if len(metrics) < len(axs):
        axs[-1].set_visible(False)

    fig.suptitle("Cluster similarity metrics vs resolution (relative to merged clusters)")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    _save(fig, figs_dir, f"homogeneity_completeness_nmi_vscore_ari_vs_resolution_{result.srxAccession}.png")
    return fig


def plot_all(
    adata: sc.AnnData,
    result: ClusterValidationResult,
    figs_dir: Path | None = None,
) -> None:
    plot_pca_cumvar(adata, result, figs_dir)
    plot_resolution_sweep(result, figs_dir)
    plot_umap_selected(adata, result, figs_dir)
    plot_rf_confusion(result, figs_dir)
    plot_umap_merged(adata, result, figs_dir)
    plot_composition_bars(adata, result, figs_dir)
    plot_silhouette(result, figs_dir)
    plot_metrics(result, figs_dir)
    plt.close("all")
