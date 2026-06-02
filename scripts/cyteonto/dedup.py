from pathlib import Path

import pandas as pd


def dedup_table(df_path: Path, accession: str) -> pd.DataFrame | None:
    """Deduplicate CyteOnto CSVs"""
    if df_path.resolve().is_file():
        try:
            df = pd.read_csv(df_path)
        except pd.errors.ParserError:
            print("Input path must be to a .csv file")
            return None

        df["pair_label"] = df["author_label"].astype(str) + df["algorithm_label"].astype(str)

        # Deduplicate, preserve non-duplicate cytescores, if there are any
        df = df.drop_duplicates(["pair_label", "cytescore_similarity"]).dropna()
        if len(df) > df["pair_label"].nunique():
            print(f"There are multiple cytescores for at least one unique label pair for {accession}")

        df["accession"] = accession
        df["table_filepath"] = str(df_path)
        return df

    else:
        print("Input must be a file")
        return None
