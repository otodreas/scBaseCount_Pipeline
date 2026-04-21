from __future__ import annotations

import pandas as pd
import pyarrow.dataset as ds

from metadata.config import MetadataConfig


def load_sample(cfg: MetadataConfig) -> pd.DataFrame:
    return pd.read_parquet(cfg.sampleParquetPath)


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
