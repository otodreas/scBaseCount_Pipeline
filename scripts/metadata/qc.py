from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from pydantic import BaseModel

from metadata.config import MetadataConfig

_QC_COLUMNS = ["SRX_accession", "gene_count_Unique", "umi_count_Unique"]
_EXPORT_COLUMNS = [
    "srx_accession",
    "file_path",
    "obs_count",
    "medianGenesPerCell",
    "medianUmisPerCell",
    "nCellsForQc",
]


class QcThresholds(BaseModel):
    minMedianGenesPerCell: int
    minMedianUmisPerCell: int


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
    """Left-join QC metrics onto samples and return rows that meet both median thresholds."""
    merged = samplesDf.merge(qcDf, on="srx_accession", how="left")
    mask = (
        (merged["medianGenesPerCell"] >= thresholds.minMedianGenesPerCell)
        & (merged["medianUmisPerCell"] >= thresholds.minMedianUmisPerCell)
    )
    return merged.loc[mask].reset_index(drop=True)


def export_datasets_qc(
    filteredDf: pd.DataFrame,
    cfg: MetadataConfig,
    filename: str = "datasets_qc.csv",
    outputPath: Path | None = None,
) -> Path:
    """Write the QC-filtered dataset CSV and return its path.

    If outputPath is provided it is used verbatim; otherwise the file is written to
    cfg.outputDir / filename. The CSV always contains srx_accession, file_path,
    obs_count, medianGenesPerCell, medianUmisPerCell, nCellsForQc.
    """
    if outputPath is None:
        cfg.outputDir.mkdir(parents=True, exist_ok=True)
        target = cfg.outputDir / filename
    else:
        outputPath.parent.mkdir(parents=True, exist_ok=True)
        target = outputPath

    filteredDf[_EXPORT_COLUMNS].to_csv(target, index=False)
    return target
