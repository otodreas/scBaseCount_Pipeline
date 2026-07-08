from __future__ import annotations

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel


class SkipReason(StrEnum):
    md5_mismatch = "md5_mismatch"
    missing_study = "missing_study"
    preprocess_failed = "preprocess_failed"
    cell_type_all_missing = "cell_type_all_missing"


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
