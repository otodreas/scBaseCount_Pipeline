import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from cluster_validation.config import ClusterValidationConfig
from cluster_validation.embedding import embed_dataset
from cluster_validation.preprocess import preprocess
from cluster_validation.resolution import select_resolution_on_graph
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from numpy.typing import NDArray
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
    v_measure_score,
)

_EXPECTED_ACCESSION = "SRX17412841"
_CELL_TYPE_KEY = "cell_type"
_RESOLUTIONS = [i / 10 for i in range(1, 21)]
_DPI = 200


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the Leiden resolution-validation composite.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="raw SRX24313469 h5ad")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output PNG path")
    args = parser.parse_args()
    adata = load_srx(args.input)
    adata, selected, k_arr, jacc_arr, metrics = run_validation(adata)
    render(adata, selected, k_arr, jacc_arr, metrics, args.output)
    print(f"Wrote {args.output}")


def load_srx(path: Path) -> AnnData:
    if not path.is_file():
        raise FileNotFoundError(f"h5ad not found: {path}")
    if path.stem != _EXPECTED_ACCESSION:
        raise ValueError(f"expected filename stem {_EXPECTED_ACCESSION!r}, got {path.stem!r}")
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    if _CELL_TYPE_KEY not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing {_CELL_TYPE_KEY!r}")
    if "gene_symbols" not in adata.var.columns:
        raise ValueError("adata.var is missing 'gene_symbols'")
    if "SRX_accession" in adata.obs.columns:
        found = sorted(adata.obs["SRX_accession"].astype(str).unique().tolist())
        if found != [_EXPECTED_ACCESSION]:
            raise ValueError(f"SRX_accession values {found} do not match {_EXPECTED_ACCESSION}")
    return adata


def run_validation(
    adata: AnnData,
) -> tuple[AnnData, float, NDArray[np.int64], NDArray[np.float64], dict[str, list[tuple[float, float]]]]:
    cfg = ClusterValidationConfig(srxAccession=_EXPECTED_ACCESSION, resolutions=_RESOLUTIONS)
    adata, _stats = preprocess(adata, cfg)
    adata, _n_pcs, _cumvar = embed_dataset(adata, cfg)
    adata, sel = select_resolution_on_graph(
        adata,
        resolutions=cfg.resolutions,
        weakPriorKey=cfg.weakPriorKey,
    )
    metrics = _metrics_vs_cell_type(adata, cfg.resolutions, cfg.weakPriorKey)
    return adata, sel.selectedResolution, sel.kArr, sel.jaccArr, metrics


def render(
    adata: AnnData,
    selected: float,
    k_arr: NDArray[np.int64],
    jacc_arr: NDArray[np.float64],
    metrics: dict[str, list[tuple[float, float]]],
    output: Path,
) -> None:
    cluster_key = f"leiden_{selected}"
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing {cluster_key!r}")
    if "X_umap" not in adata.obsm:
        raise ValueError("adata.obsm is missing 'X_umap'")

    fig = plt.figure(figsize=(14.0, 11.0))
    outer = GridSpec(2, 2, figure=fig, wspace=0.28, hspace=0.32)

    umap_gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0], wspace=0.35)
    ax_umap_cluster = fig.add_subplot(umap_gs[0, 0])
    ax_umap_type = fig.add_subplot(umap_gs[0, 1])
    _plot_umap(
        ax_umap_cluster, adata, cluster_key, f"Leiden {selected:g} ({adata.obs[cluster_key].nunique()} clusters)"
    )
    _plot_umap(ax_umap_type, adata, _CELL_TYPE_KEY, f"cell_type ({adata.obs[_CELL_TYPE_KEY].nunique()} types)")
    _panel_label(ax_umap_cluster, "A")

    bars_gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 1], wspace=0.3)
    ax_cluster_bar = fig.add_subplot(bars_gs[0, 0])
    ax_type_bar = fig.add_subplot(bars_gs[0, 1], sharey=ax_cluster_bar)
    _plot_composition(ax_cluster_bar, adata.obs[cluster_key], "Leiden clusters", "steelblue")
    _plot_composition(ax_type_bar, adata.obs[_CELL_TYPE_KEY], "cell_type", "indianred")
    ax_cluster_bar.set_ylabel("Proportion of cells (%)")
    _panel_label(ax_cluster_bar, "B")

    sweep_gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1, 0], wspace=0.3)
    ax_k = fig.add_subplot(sweep_gs[0, 0])
    ax_jacc = fig.add_subplot(sweep_gs[0, 1])
    _plot_sweep(ax_k, _RESOLUTIONS, k_arr, selected, "# clusters", "Clusters per resolution")
    _plot_jaccard(ax_jacc, _RESOLUTIONS, jacc_arr, selected)
    _panel_label(ax_k, "C")

    metric_gs = GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[1, 1], wspace=0.35, hspace=0.45)
    metric_axes = [fig.add_subplot(metric_gs[i, j]) for i in range(2) for j in range(3)]
    names = list(metrics)
    for idx, ax in enumerate(metric_axes):
        if idx >= len(names):
            ax.set_visible(False)
            continue
        name = names[idx]
        _plot_metric(ax, metrics[name], selected, name)
    _panel_label(metric_axes[0], "D")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)


def _metrics_vs_cell_type(
    adata: AnnData,
    resolutions: list[float],
    weak_prior_key: str,
) -> dict[str, list[tuple[float, float]]]:
    ref = adata.obs[weak_prior_key]
    out: dict[str, list[tuple[float, float]]] = {
        "Homogeneity": [],
        "Completeness": [],
        "V-measure": [],
        "NMI": [],
        "ARI": [],
    }
    scorers = {
        "Homogeneity": homogeneity_score,
        "Completeness": completeness_score,
        "V-measure": v_measure_score,
        "NMI": normalized_mutual_info_score,
        "ARI": adjusted_rand_score,
    }
    for res in resolutions:
        key = f"leiden_{res}"
        labels = adata.obs[key]
        if labels.nunique() <= 1:
            continue
        for name, scorer in scorers.items():
            out[name].append((res, float(scorer(ref, labels))))
    return out


def _plot_umap(ax: Axes, adata: AnnData, color_by: str, title: str) -> None:
    coords = np.asarray(adata.obsm["X_umap"][:, :2], dtype=np.float64)
    series = adata.obs[color_by].astype(str)
    categories = series.unique().tolist()
    cmap = plt.get_cmap("tab20", max(len(categories), 1))
    for idx, category in enumerate(categories):
        mask = series.to_numpy() == category
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=4,
            alpha=0.85,
            c=[cmap(idx)],
            linewidths=0,
            rasterized=True,
        )
        centroid = coords[mask].mean(axis=0)
        ax.text(centroid[0], centroid[1], category, fontsize=6, ha="center", va="center")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal", adjustable="datalim")


def _plot_composition(ax: Axes, series: pd.Series, title: str, color: str) -> None:
    counts = series.astype(str).value_counts(normalize=True).sort_values(ascending=False) * 100
    ax.bar(range(len(counts)), counts.to_numpy(), color=color)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(list(counts.index), rotation=90, fontsize=7)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("")


def _plot_sweep(
    ax: Axes, resolutions: list[float], k_arr: NDArray[np.int64], selected: float, ylabel: str, title: str
) -> None:
    ax.plot(resolutions, k_arr, marker="o", ms=4, color="steelblue")
    ax.axvline(selected, color="red", linestyle="--", label=f"selected = {selected:g}")
    ax.set_xlabel("Resolution")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)


def _plot_jaccard(ax: Axes, resolutions: list[float], jacc_arr: NDArray[np.float64], selected: float) -> None:
    best_idx = resolutions.index(selected)
    ax.plot(resolutions, jacc_arr, marker="o", ms=4, color="darkorange", label="matched Jaccard")
    ax.axvline(selected, color="red", linestyle="--", label=f"argmax = {selected:g}")
    ax.scatter([selected], [jacc_arr[best_idx]], color="red", zorder=5, s=60)
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Matched Jaccard")
    ax.set_title("Matched Jaccard score", fontsize=9)
    ax.legend(fontsize=7)


def _plot_metric(ax: Axes, values: list[tuple[float, float]], selected: float, title: str) -> None:
    xs = [pair[0] for pair in values]
    ys = [pair[1] for pair in values]
    ax.plot(xs, ys, marker="o", ms=3, color="steelblue")
    ax.axvline(selected, color="red", linestyle="--", linewidth=1, label=f"selected = {selected:g}")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel(title)
    ax.legend(fontsize=6)


def _panel_label(ax: Axes, letter: str) -> None:
    ax.text(
        -0.18,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


if __name__ == "__main__":
    main()
