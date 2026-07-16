from pathlib import Path
from typing import Literal, Self

from pydantic import BaseModel, Field, model_validator
from shared.repo import REPO_ROOT
from study_context.utils import CONTEXTS_JSONL_PATH


class H5adConcatConfig(BaseModel):
    r2Keys: list[str] | None = None
    datasetsPath: Path | None = REPO_ROOT / "output" / "metadata" / "datasets.csv"
    contextsPath: Path = CONTEXTS_JSONL_PATH
    cellTypeKey: str = "cell_type"
    # Column added to obs holding the experimental batch key; value is the ENA study accession.
    batchKey: str = "study_accession"
    # Existing obs column holding the per-file experiment accession; used to make barcodes unique.
    accessionKey: str = "SRX_accession"
    missingLabel: str = "UNKNOWN"
    # Reference genes are mandatory
    geneInfoPath: Path = (
        REPO_ROOT
        / "data"
        / "scbasecount"
        / "2026-01-12"
        / "star_references"
        / "Homo_sapiens"
        / "hg38_2020"
        / "geneInfo.tab"
    )
    cacheDir: Path = REPO_ROOT / "data" / "h5ad_concat" / "cache"
    outputPath: Path = REPO_ROOT / "output" / "atlas" / "data" / "atlas.h5ad"
    # Placeholder: yet to be implemented
    downloadBatchSize: int = Field(default=8, ge=1)
    compression: Literal["gzip", "lzf"] | None = "gzip"
    verifyMd5: bool = True
    uploadAtlas: bool = True
    atlasR2Key: str | None = None
    preprocess: bool = True
    # Reindex every layer (e.g. STARsolo UniqueAndMult matrices) onto the reference axis instead of X only.
    conserveLayers: bool = False
    # Minimum genes detected per cell; 0 disables the filter.
    minGenesPerCell: int = Field(default=200, ge=0)
    # Maximum mitochondrial read fraction per cell, as a fraction in (0, 1]; 1.0 keeps every cell.
    maxPctMito: float = Field(default=0.2, gt=0.0, le=1.0)
    # Hemoglobin read fraction ceiling, as a fraction in (0, 1]; 1.0 records pct_counts_hb without filtering.
    maxPctHb: float = Field(default=1.0, gt=0.0, le=1.0)
    # Minimum cells expressing a gene; 0 disables the filter.
    minCellsPerGene: int = Field(default=0, ge=0)
    # Absolute floor on cells remaining after QC before a file is rejected.
    minCellsAfterQc: int = Field(default=100, ge=1)
    # Relative floor: reject a file when less than this fraction of input cells remain after QC; 0 disables the gate.
    minPctCellsAfterQc: float = Field(default=0.4, ge=0.0, le=1.0)

    @model_validator(mode="after")
    def validate_input_and_upload(self) -> Self:
        has_r2_keys = self.r2Keys is not None and len(self.r2Keys) > 0
        has_datasets = self.datasetsPath is not None
        if has_r2_keys and has_datasets:
            msg = "Provide r2Keys or datasetsPath, not both"
            raise ValueError(msg)
        if not has_r2_keys and not has_datasets:
            msg = "Provide either r2Keys or datasetsPath"
            raise ValueError(msg)
        if self.uploadAtlas and not self.atlasR2Key:
            msg = "atlasR2Key is required when uploadAtlas is true"
            raise ValueError(msg)
        return self
