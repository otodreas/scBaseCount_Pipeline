from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel
from shared.repo import REPO_ROOT
from study_context.utils import CONTEXTS_JSONL_PATH


class H5adConcatConfig(BaseModel):
    r2Keys: list[str]  # TODO(input): support datasets.csv (parse file where a column contains the keys)
    contextsPath: Path = CONTEXTS_JSONL_PATH
    cellTypeKey: str = "cell_type"
    # Column added to obs holding the experimental batch key; value is the ENA study accession.
    batchKey: str = "study_accession"
    missingLabel: str = "UNKNOWN"
    join: Literal["inner", "outer"] = "inner"
    cacheDir: Path = REPO_ROOT / "data" / "h5ad_concat" / "cache"
    outputPath: Path = REPO_ROOT / "output" / "atlas" / "data" / "atlas.h5ad"
    downloadBatchSize: int = 8
    compression: Literal["gzip", "lzf"] | None = "gzip"
    verifyMd5: bool = True
    checkpointEvery: int = 25
    resume: bool = True
    # TODO(preprocess): preprocess: bool = True
