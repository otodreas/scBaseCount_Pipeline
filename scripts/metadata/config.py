from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel
from shared.repo import REPO_ROOT as _REPO_ROOT

_METADATA_DIR = _REPO_ROOT / "data" / "scbasecount" / "2026-01-12" / "metadata" / "GeneFull" / "Homo_sapiens"


class MetadataConfig(BaseModel):
    sampleParquetPath: Path = (
        _METADATA_DIR / "scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_sample_metadata.parquet"
    )
    obsParquetPath: Path = _METADATA_DIR / "scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_obs_metadata.parquet"
    minObsCount: int = 1000  # Minimum cells per sample; samples below this are dropped before any other filtering
    outputDir: Path = _REPO_ROOT / "output" / "metadata"
