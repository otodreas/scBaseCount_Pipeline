from __future__ import annotations

from pathlib import Path

import pandas as pd

_REQUIRED_COLUMNS = frozenset({"author_label", "algorithm_label", "cytescore_similarity"})


def dedup_table(df_path: Path, accession: str) -> pd.DataFrame | None:
    """Deduplicate CyteOnto CSVs to one row per author/algorithm label pair."""
    if not df_path.resolve().is_file():
        print("Input must be a file")
        return None

    try:
        df = pd.read_csv(df_path)
    except pd.errors.ParserError:
        print("Input path must be to a .csv file")
        return None

    missing = _REQUIRED_COLUMNS - set(df.columns)
    if missing:
        print(f"Missing required column(s) {sorted(missing)!r} in {df_path}")
        return None

    if df.empty:
        return None

    df = df.copy()
    df["pair_label"] = df["author_label"].astype(str) + df["algorithm_label"].astype(str)
    df = df.drop_duplicates(["pair_label", "cytescore_similarity"])
    df = df.dropna(subset=["pair_label", "cytescore_similarity"])

    if df.empty:
        return None

    if len(df) > df["pair_label"].nunique():
        print(f"There are multiple cytescores for at least one unique label pair for {accession}")
        df = (
            df.sort_values("cytescore_similarity", ascending=False)
            .drop_duplicates("pair_label", keep="first")
            .reset_index(drop=True)
        )

    df["accession"] = accession
    df["table_filepath"] = str(df_path)
    return df
