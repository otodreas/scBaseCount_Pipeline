from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy


def compute_nse_kld_row(
    adata: Any,
    merged_key: str,
) -> tuple[dict[str, float], dict[str, float]]:
    """Compute per-cell-type normalized Shannon entropy and KL divergence for one dataset.

    For each cell type, two scalars are produced:

    - Normalized Shannon entropy (NSE): measures how fragmented the cell type is
      across Leiden clusters. A value of 0 means all cells sit in a single cluster;
      a value of 1 means cells are distributed as evenly as possible across every
      cluster they occupy.

    - KL divergence (KLD): KL(p || q) where p is the cluster distribution of the
      cell type and q is the global cluster distribution of all cells. A high value
      means the cell type is concentrated in clusters that differ strongly from the
      background, indicating high coherence. A value near 0 means the cell type
      mirrors the global distribution.

    Args:
        adata: AnnData object with obs columns "cell_type" and the column named by
            merged_key containing Leiden cluster labels.
        merged_key: Column in adata.obs that holds the merged Leiden cluster labels.

    Returns:
        A pair (nse_row, kld_row), each a dict mapping cell type name to a float.
    """
    nse_row: dict[str, float] = {}
    kld_row: dict[str, float] = {}
    global_labels = adata.obs[merged_key]
    all_clusters = global_labels.cat.categories
    q = global_labels.value_counts(normalize=True).reindex(all_clusters, fill_value=0)
    for cell_type in adata.obs["cell_type"].unique():
        labels = adata.obs[adata.obs["cell_type"] == cell_type][merged_key]
        counts = labels.value_counts(normalize=True)
        counts = counts[counts > 0]
        nse_row[cell_type] = (
            float(scipy_entropy(counts, base=2) / np.log2(len(counts)))
            if len(counts) > 1
            else 0.0
        )
        p = labels.value_counts(normalize=True).reindex(all_clusters, fill_value=0)
        kld_row[cell_type] = float(scipy_entropy(p, q, base=2))
    return nse_row, kld_row


def build_metric_dataframes(
    rows: list[dict],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build NSE matrix, KLD matrix, and cross-dataset summary from parsed JSONL rows.

    Each row is expected to have the shape::

        {"srx": "SRX123", "nse": {<cell_type>: <float>, ...}, "kld": {<cell_type>: <float>, ...}}

    Args:
        rows: List of dicts parsed from metrics_matrix.jsonl, one entry per accession.

    Returns:
        A 3-tuple (nse_df, kld_df, summary_df):

        - nse_df: DataFrame with accessions as rows and cell types as columns,
          containing normalized Shannon entropy values.
        - kld_df: DataFrame with the same shape, containing KL divergence values.
        - summary_df: DataFrame indexed by cell type with columns n_datasets,
          normalized_shannon_entropy_mean, and kl_divergence_mean.
    """
    nse_df = pd.DataFrame({r["srx"]: r["nse"] for r in rows}).T
    nse_df.index.name = "srx"

    kld_df = pd.DataFrame({r["srx"]: r["kld"] for r in rows}).T
    kld_df.index.name = "srx"

    summary_df = pd.DataFrame({
        "n_datasets": nse_df.notna().sum(),
        "normalized_shannon_entropy_mean": nse_df.mean(),
        "kl_divergence_mean": kld_df.mean(),
    })
    summary_df.index.name = "cell_type"

    return nse_df, kld_df, summary_df


def save_metric_plot(
    summary_df: pd.DataFrame,
    output_path: Path,
    sort_by: str = "normalized_shannon_entropy_mean",
    dpi: int = 150,
) -> None:
    """Save a two-panel horizontal bar chart of NSE and KLD means per cell type.

    The left panel shows normalized Shannon entropy (higher = more fragmented).
    The right panel shows mean KL divergence (higher = more coherent / distinct
    from the global cluster background). Cell types are sorted by sort_by ascending
    so the most fragmented cell type appears at the top.

    Args:
        summary_df: DataFrame produced by build_metric_dataframes, indexed by cell
            type, with columns normalized_shannon_entropy_mean and kl_divergence_mean.
        output_path: Destination path for the PNG file.
        sort_by: Column in summary_df used to sort cell types (default:
            normalized_shannon_entropy_mean).
        dpi: Resolution of the saved image (default: 150).
    """
    sorted_df = summary_df.sort_values(sort_by)
    labels = sorted_df.index.tolist()
    y = np.arange(len(labels))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(6, len(labels) * 0.5)))

    ax1.barh(y, sorted_df["normalized_shannon_entropy_mean"], color="steelblue")
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels)
    ax1.set_xlabel("normalized Shannon entropy [0, 1]")
    ax1.set_title("Normalized Shannon Entropy (↑ = more fragmented)")

    ax2.barh(y, sorted_df["kl_divergence_mean"], color="darkorange")
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("KL divergence (bits)")
    ax2.set_title("KL Divergence (↑ = more coherent)")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
