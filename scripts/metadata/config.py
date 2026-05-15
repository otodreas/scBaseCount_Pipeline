from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel
from shared.repo import REPO_ROOT as _REPO_ROOT


class MetadataConfig(BaseModel):
    sampleParquetPath: Path
    obsParquetPath: Path
    minObsCount: int = 1000
    outputDir: Path = _REPO_ROOT / "output" / "metadata"
