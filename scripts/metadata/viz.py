from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metadata.filter import FilterResult
from metadata.regexes import DISEASE_MAP

_TABLEAU10 = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
    "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
    "#9c755f",
]


def plot_sample_breakdown(
    sample: pd.DataFrame,
    result: FilterResult,
    figs_dir: Path,
) -> None:
    figs_dir.mkdir(parents=True, exist_ok=True)

    n_total       = len(sample)
    n_discarded   = n_total - len(result.sampleKnown)
    n_lung_cancer = len(result.lungIntersectionCancer)
    n_lung_other  = len(result.lungIntersection) - n_lung_cancer
    n_not_lung    = len(result.sampleKnown) - len(result.lungIntersection)

    labels = ["Discarded\n(unknown/healthy)", "Lung (cancer)", "Lung (other)", "Non-lung"]
    sizes  = [n_discarded, n_lung_cancer, n_lung_other, n_not_lung]
    colors = ["#d3d3d3", "#e05c5c", "#f4a261", "#5b8db8"]

    fig, ax = plt.subplots(figsize=(7, 7))
    _, _, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%\n({int(round(p * n_total / 100)):,})",
        startangle=140,
        pctdistance=0.75,
    )
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title(
        f"scBaseCount SRX breakdown by disease\n(n={n_total:,} unique SRX)",
        fontsize=13,
        pad=16,
    )
    plt.tight_layout()
    fig.savefig(figs_dir / "sample_breakdown.png", bbox_inches="tight")
    plt.show()
    plt.close(fig)


def plot_disease_breakdown(result: FilterResult, figs_dir: Path) -> None:
    figs_dir.mkdir(parents=True, exist_ok=True)

    def _categorize(d: str) -> str:
        for label, pat in DISEASE_MAP:
            if pat.search(str(d)):
                return label
        return "Other"

    disease_cats = result.lungIntersection["disease"].map(_categorize)
    cat_counts   = disease_cats.value_counts()
    total_lung   = cat_counts.sum()

    large     = cat_counts[cat_counts / total_lung >= 0.02]
    small_sum = total_lung - large.sum()
    if small_sum > 0:
        large = pd.concat([large, pd.Series({"Other": small_sum})])
    large = large.groupby(large.index).sum()

    non_other = large.drop("Other", errors="ignore").sort_values(ascending=False)
    if "Other" in large.index:
        large = pd.concat([non_other, pd.Series({"Other": large["Other"]})])
    else:
        large = non_other

    colors = [_TABLEAU10[i % len(_TABLEAU10)] for i in range(len(large) - 1)] + ["#bab0ac"]

    fig, ax = plt.subplots(figsize=(8, 8))
    _, _, autotexts = ax.pie(
        large.values,
        labels=large.index,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%\n({int(round(p * total_lung / 100)):,})",
        startangle=90,
        counterclock=False,
        pctdistance=0.75,
    )
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title(
        f"Lung disease SRX breakdown\n(n={total_lung:,} unique SRX)",
        fontsize=13,
        pad=16,
    )
    plt.tight_layout()
    fig.savefig(figs_dir / "lung_disease_breakdown.png", bbox_inches="tight")
    plt.show()
    plt.close(fig)


def plot_cell_count_distribution(result: FilterResult, figs_dir: Path) -> None:
    figs_dir.mkdir(parents=True, exist_ok=True)

    counts  = result.lungIntersection["obs_count"]
    n_srx   = len(counts)
    n_cells = int(np.sum(counts))

    plt.hist(counts / 1000, log=True)
    plt.title(
        f"Distribution of Number of Cells per Lung SRX (log scale)\n"
        f"(n={n_cells:,} cells across {n_srx:,} SRX)"
    )
    plt.xlabel("Number of Cells per SRX (thousands)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(figs_dir / "lung_cell_number_hist.png", bbox_inches="tight")
    plt.show()
    plt.close(fig)
