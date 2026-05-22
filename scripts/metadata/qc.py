from __future__ import annotations

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from pydantic import BaseModel

from metadata.config import MetadataConfig

_QC_COLUMNS = ["SRX_accession", "gene_count_Unique", "umi_count_Unique"]


class QcThresholds(BaseModel):
    minMedianGenesPerCell: int | None = None
    minMedianUmisPerCell: int | None = None


def compute_obs_qc(srxAccessions: list[str], cfg: MetadataConfig) -> pd.DataFrame:
    """Compute per-SRX nCellsForQc and median genes/UMIs per cell from the obs parquet.

    Reads only rows whose SRX_accession is in srxAccessions, projects to the three
    columns needed, and returns a DataFrame with one row per SRX containing
    srx_accession, nCellsForQc, medianGenesPerCell, medianUmisPerCell.
    """
    dset = ds.dataset(str(cfg.obsParquetPath), format="parquet")
    tbl = dset.to_table(
        columns=_QC_COLUMNS,
        filter=pc.field("SRX_accession").isin(srxAccessions),
    )
    df = tbl.to_pandas()
    qc = (
        df.groupby("SRX_accession")
        .agg(
            nCellsForQc=("gene_count_Unique", "size"),
            medianGenesPerCell=("gene_count_Unique", "median"),
            medianUmisPerCell=("umi_count_Unique", "median"),
        )
        .reset_index()
        .rename(columns={"SRX_accession": "srx_accession"})
    )
    return qc


def apply_qc(samplesDf: pd.DataFrame, qcDf: pd.DataFrame, thresholds: QcThresholds) -> pd.DataFrame:
    """Left-join QC metrics onto samples and return rows meeting the configured median thresholds.

    A None threshold disables that metric. Rows missing a QC value for an enabled
    metric are filtered out (NaN >= n is False).
    """
    merged = samplesDf.merge(qcDf, on="srx_accession", how="left")
    mask = pd.Series(True, index=merged.index)
    if thresholds.minMedianGenesPerCell is not None:
        mask &= merged["medianGenesPerCell"] >= thresholds.minMedianGenesPerCell
    if thresholds.minMedianUmisPerCell is not None:
        mask &= merged["medianUmisPerCell"] >= thresholds.minMedianUmisPerCell
    return merged.loc[mask].reset_index(drop=True)
