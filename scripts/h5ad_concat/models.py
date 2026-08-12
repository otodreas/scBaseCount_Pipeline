from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field, model_validator

CELL_FILTER_ORDER: tuple[str, ...] = (
    "minGenesPerCell",
    "maxPctMito",
    "maxPctRibo",
    "maxPctHb",
)


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


class PostFilterMedians(BaseModel):
    nGenesByCounts: float
    pctCountsMito: float
    pctCountsRibo: float
    pctCountsHb: float


class QcStats(BaseModel):
    nCellsBefore: int
    nCellsAfter: int
    nCellsDropped: int
    nCellsDroppedByFilter: dict[str, int]
    nGenesBefore: int
    nGenesAfter: int
    nGenesDroppedByFilter: dict[str, int]
    pctCellsAfter: float
    postFilterMedians: PostFilterMedians

    @model_validator(mode="after")
    def validate_cell_drop_invariant(self) -> "QcStats":
        expected = self.nCellsBefore - self.nCellsAfter
        filter_sum = sum(self.nCellsDroppedByFilter.values())
        if self.nCellsDropped != expected or self.nCellsDropped != filter_sum:
            msg = (
                f"QC cell-drop invariant failed: before={self.nCellsBefore} after={self.nCellsAfter} "
                f"dropped={self.nCellsDropped} filterSum={filter_sum} byFilter={self.nCellsDroppedByFilter}"
            )
            raise ValueError(msg)
        return self

    @property
    def medianGenesPerCell(self) -> float:
        return self.postFilterMedians.nGenesByCounts

    @property
    def medianPctMito(self) -> float:
        return self.postFilterMedians.pctCountsMito

    @property
    def medianPctRibo(self) -> float:
        return self.postFilterMedians.pctCountsRibo

    @property
    def medianPctHb(self) -> float:
        return self.postFilterMedians.pctCountsHb


class SkippedFile(BaseModel):
    r2Key: str
    accession: str
    reason: SkipReason
    studyAccession: str = ""
    qc: QcStats | None = None


class FileRecord(BaseModel):
    accession: str
    studyAccession: str
    r2Key: str
    status: Literal["success", "skip"]
    skipReason: SkipReason | None = None
    qc: QcStats | None = None


class QcCohortSummary(BaseModel):
    nFiles: int = 0
    nCellsBefore: int = 0
    nCellsAfter: int = 0
    nCellsDropped: int = 0
    nCellsDroppedByFilter: dict[str, int] = Field(default_factory=lambda: {name: 0 for name in CELL_FILTER_ORDER})

    def add(self, qc: QcStats) -> None:
        self.nFiles += 1
        self.nCellsBefore += qc.nCellsBefore
        self.nCellsAfter += qc.nCellsAfter
        self.nCellsDropped += qc.nCellsDropped
        for name in CELL_FILTER_ORDER:
            self.nCellsDroppedByFilter[name] = self.nCellsDroppedByFilter.get(name, 0) + qc.nCellsDroppedByFilter.get(
                name, 0
            )


class QcSummary(BaseModel):
    concatenatedFiles: QcCohortSummary = Field(default_factory=QcCohortSummary)
    allQcProcessedFiles: QcCohortSummary = Field(default_factory=QcCohortSummary)


class H5adConcatResult(BaseModel):
    outputPath: Path
    nObs: int
    nVars: int
    nFilesConcatenated: int
    nFilesSkipped: int = 0
    studiesSeen: list[str]
    skipped: list[SkippedFile]
    cellFilterOrder: list[str] = Field(default_factory=lambda: list(CELL_FILTER_ORDER))
    qcSummary: QcSummary = Field(default_factory=QcSummary)
    files: list[FileRecord] = Field(default_factory=list)
    fileLogPath: Path
    configPath: Path
    atlasR2Key: str | None = None
    atlasFileLogR2Key: str | None = None
    atlasConfigR2Key: str | None = None
    atlasResultR2Key: str | None = None
    conserveLayers: bool = False
