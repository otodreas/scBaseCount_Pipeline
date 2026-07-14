from __future__ import annotations

from pathlib import Path
from typing import Literal, Self

from pydantic import BaseModel, model_validator
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
    join: Literal["inner", "outer"] = "inner"
    cacheDir: Path = REPO_ROOT / "data" / "h5ad_concat" / "cache"
    outputPath: Path = REPO_ROOT / "output" / "atlas" / "data" / "atlas.h5ad"
    downloadBatchSize: int = 8
    compression: Literal["gzip", "lzf"] | None = "gzip"
    verifyMd5: bool = True
    uploadAtlas: bool = False
    atlasR2Key: str | None = None
    # TODO(preprocess): preprocess: bool = True

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
