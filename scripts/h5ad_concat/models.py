from __future__ import annotations

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel


class SkipReason(StrEnum):
    download_failed = "download_failed"
    read_failed = "read_failed"
    md5_mismatch = "md5_mismatch"
    missing_study = "missing_study"
    too_few_cells = "too_few_cells"
    excessive_cell_dropout = "excessive_cell_dropout"
    cell_type_all_missing = "cell_type_all_missing"
    accession_mismatch = "accession_mismatch"
    gene_axis_mismatch = "gene_axis_mismatch"


class SkippedFile(BaseModel):
    r2Key: str
    accession: str
    reason: SkipReason


class H5adConcatResult(BaseModel):
    outputPath: Path
    nObs: int
    nVars: int
    nFilesConcatenated: int
    studiesSeen: list[str]
    skipped: list[SkippedFile]
    statusCsvPath: Path
    configPath: Path
    atlasR2Key: str | None = None
    atlasStatusR2Key: str | None = None
    atlasConfigR2Key: str | None = None
    atlasResultR2Key: str | None = None
    conserveLayers: bool = False
