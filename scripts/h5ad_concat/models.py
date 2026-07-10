from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field


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


class ManifestEntry(BaseModel):
    r2Key: str
    accession: str
    concatenated: bool
    study: str | None = None
    reason: SkipReason | None = None


class ConcatManifest(BaseModel):
    join: Literal["inner", "outer"]
    batchKey: str
    entries: list[ManifestEntry] = Field(default_factory=list)

    def processedKeys(self) -> set[str]:
        """Return r2Keys already recorded in this manifest."""
        return {entry.r2Key for entry in self.entries}

    def skippedFiles(self) -> list[SkippedFile]:
        """Return skipped-file records derived from manifest entries."""
        return [
            SkippedFile(r2Key=entry.r2Key, accession=entry.accession, reason=entry.reason)
            for entry in self.entries
            if not entry.concatenated and entry.reason is not None
        ]

    def concatenatedCount(self) -> int:
        """Return the number of files folded into the atlas."""
        return sum(1 for entry in self.entries if entry.concatenated)

    def studiesSeen(self) -> list[str]:
        """Return sorted unique study accessions from concatenated entries."""
        studies = {entry.study for entry in self.entries if entry.concatenated and entry.study is not None}
        return sorted(studies)
