import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

_EMBEDDING_LABELS = {
    "X_pca": "Uncorrected",
    "X_pca_harmony": "Harmony",
}
_EMBEDDING_COLORS = {
    "X_pca": "steelblue",
    "X_pca_harmony": "darkorange",
}
_ROW_ORDER = ("X_pca", "X_pca_harmony")
_TYPE_ORDER = ("Bio conservation", "Batch correction", "Aggregate score")
_AGGREGATE_TICKS = (
    ("Bio conservation", "biological"),
    ("Batch correction", "batch correction"),
)
_INDIVIDUAL_TYPES = ("Bio conservation", "Batch correction")
_BAR_WIDTH = 0.36


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the two-panel scIB bar chart.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="scib_results.csv")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output figure path")
    args = parser.parse_args()
    scores, metric_types = load_scib(args.input)
    render(scores, metric_types, args.output)
    print(f"Wrote {args.output}")


def load_scib(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    if not path.is_file():
        raise FileNotFoundError(f"scIB CSV not found: {path}")
    raw = pd.read_csv(path, index_col=0)
    if "Metric Type" not in raw.index:
        raise ValueError(f"{path} is missing a Metric Type row")
    metric_types = raw.loc["Metric Type"].astype(str)
    scores = raw.drop(index="Metric Type").apply(pd.to_numeric, errors="raise")
    missing_rows = [name for name in _ROW_ORDER if name not in scores.index]
    if missing_rows:
        raise ValueError(f"{path} is missing embeddings: {', '.join(missing_rows)}")
    extra_rows = [name for name in scores.index if name not in _EMBEDDING_LABELS]
    if extra_rows:
        raise ValueError(f"{path} has unexpected embeddings: {', '.join(extra_rows)}")
    unknown_types = sorted(set(metric_types) - set(_TYPE_ORDER))
    if unknown_types:
        raise ValueError(f"{path} has unexpected metric types: {', '.join(unknown_types)}")
    if scores.isna().any().any():
        raise ValueError(f"{path} contains NaN scores")
    scores = scores.loc[list(_ROW_ORDER)]
    return scores, metric_types


def render(scores: pd.DataFrame, metric_types: pd.Series, output: Path) -> None:
    aggregate_cols = [column for column, _tick in _AGGREGATE_TICKS]
    missing_aggregates = [column for column in aggregate_cols if column not in scores.columns]
    if missing_aggregates:
        raise ValueError(f"missing aggregate columns: {', '.join(missing_aggregates)}")
    for column in aggregate_cols:
        if metric_types[column] != "Aggregate score":
            raise ValueError(f"{column} has metric type {metric_types[column]!r}, expected 'Aggregate score'")

    individual_cols = [col for col in scores.columns if metric_types[col] in _INDIVIDUAL_TYPES]
    if not individual_cols:
        raise ValueError("no individual scIB metrics found")
    grouped_individual: list[str] = []
    for group in _INDIVIDUAL_TYPES:
        cols = [col for col in individual_cols if metric_types[col] == group]
        if not cols:
            raise ValueError(f"no individual metrics with type {group!r}")
        grouped_individual.extend(cols)

    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(11.5, 3.8),
        gridspec_kw={"width_ratios": [1.0, 2.6]},
        sharey=True,
    )
    _grouped_bars(ax_a, scores, aggregate_cols, [tick for _column, tick in _AGGREGATE_TICKS])
    _grouped_bars(ax_b, scores, grouped_individual, grouped_individual)
    ax_a.set_ylabel("Score")
    ax_a.set_ylim(0.0, 1.0)
    ax_b.tick_params(axis="x", labelrotation=45)
    for label in ax_b.get_xticklabels():
        label.set_ha("right")
    for ax in (ax_a, ax_b):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-0.5, len(ax.get_xticks()) - 0.5)
    _panel_label(ax_a, "A")
    _panel_label(ax_b, "B")
    ax_a.legend(frameon=False, loc="upper left")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def _grouped_bars(ax: Axes, scores: pd.DataFrame, columns: list[str], tick_labels: list[str]) -> None:
    x = np.arange(len(columns))
    for i, embedding in enumerate(_ROW_ORDER):
        ax.bar(
            x + (i - 0.5) * _BAR_WIDTH,
            scores.loc[embedding, columns].to_numpy(dtype=float),
            _BAR_WIDTH,
            color=_EMBEDDING_COLORS[embedding],
            label=_EMBEDDING_LABELS[embedding],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels)


def _panel_label(ax: Axes, letter: str) -> None:
    ax.text(
        -0.08,
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
