from __future__ import annotations

from pathlib import Path

import pandas as pd

CONFIDENCE_ORDER = ["High", "Moderate", "Low"]

SUMMARY_PAIR_COLUMNS = [
    "accession",
    "cell_type",
    "cytetype_annotation_leiden_merged",
    "leiden_merged",
    "cytetype_confidence",
    "cytescore_similarity",
    "n_cells",
]


def top_bottom_by_cytetype(pair_df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Return top and bottom cytescore STATE cell types per CyteType label, ordered by confidence."""
    pairs = (
        pair_df[SUMMARY_PAIR_COLUMNS]
        .drop_duplicates(["cell_type", "cytetype_annotation_leiden_merged"])
        .dropna(subset=["cytescore_similarity"])
    )

    cytetype_order = (
        pairs.groupby("cytetype_annotation_leiden_merged", observed=True)["cytetype_confidence"]
        .first()
        .reset_index()
        .assign(
            cytetype_confidence=lambda d: pd.Categorical(
                d["cytetype_confidence"], categories=CONFIDENCE_ORDER, ordered=True
            )
        )
        .sort_values(["cytetype_confidence", "cytetype_annotation_leiden_merged"])
    )

    rows: list[dict[str, object]] = []
    for _, meta in cytetype_order.iterrows():
        cytetype = meta["cytetype_annotation_leiden_merged"]
        conf = meta["cytetype_confidence"]
        sub = pairs[pairs["cytetype_annotation_leiden_merged"] == cytetype]
        n_take = min(n, len(sub))
        top = sub.nlargest(n_take, "cytescore_similarity")
        bottom = sub.nsmallest(n_take, "cytescore_similarity")
        n_rows = max(len(top), len(bottom))
        for rank in range(n_rows):
            rows.append(
                {
                    "cytetype_annotation_leiden_merged": cytetype if rank == 0 else "",
                    "cytetype_confidence": conf if rank == 0 else "",
                    "rank": rank + 1,
                    "top_cell_type": top.iloc[rank]["cell_type"] if rank < len(top) else pd.NA,
                    "top_cytescore_similarity": top.iloc[rank]["cytescore_similarity"] if rank < len(top) else pd.NA,
                    "bottom_cell_type": bottom.iloc[rank]["cell_type"] if rank < len(bottom) else pd.NA,
                    "bottom_cytescore_similarity": (
                        bottom.iloc[rank]["cytescore_similarity"] if rank < len(bottom) else pd.NA
                    ),
                }
            )

    return pd.DataFrame(rows)


def write_extremes_csv(pair_df: pd.DataFrame, n: int, output_path: Path) -> Path:
    """Write the extremes table for a pair-level summary DataFrame to CSV."""
    extremes = top_bottom_by_cytetype(pair_df, n)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    extremes.to_csv(output_path, index=False)
    return output_path
