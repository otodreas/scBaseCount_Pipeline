from __future__ import annotations

import pandas as pd
import pyarrow.dataset as ds

from metadata.config import MetadataConfig


def load_sample(cfg: MetadataConfig) -> pd.DataFrame:
    return pd.read_parquet(cfg.sampleParquetPath)


def sample_row_for_srx(
    srx_accession: str,
    cfg: MetadataConfig,
    columns: list[str] | None = None,
) -> pd.Series | None:
    """Return the sample-metadata row for srx_accession, or None if not in the catalog."""
    dset = ds.dataset(str(cfg.sampleParquetPath), format="parquet")
    tbl = dset.to_table(
        columns=columns,
        filter=ds.field("srx_accession") == srx_accession,
    )
    df = tbl.to_pandas()
    if df.empty:
        return None
    if len(df) > 1:
        raise ValueError(f"Expected one row for {srx_accession!r}, found {len(df)}")
    return df.iloc[0]


def obs_rows_for_srx(
    srx_accession: str,
    cfg: MetadataConfig,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    dset = ds.dataset(str(cfg.obsParquetPath), format="parquet")
    tbl = dset.to_table(
        columns=columns,
        filter=ds.field("SRX_accession") == srx_accession,
    )
    return tbl.to_pandas()
