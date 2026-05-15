from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT as _REPO_ROOT


class H5adExtractConfig(BaseModel):
    h5adPath: Path
    columnNames: list[str] = Field(min_length=1)
    outputDir: Path = _REPO_ROOT / "output" / "h5ad_extract"
    outputPath: Path | None = None
    annotationAxis: Literal["obs", "var"] = "obs"
    outputFormat: Literal["parquet", "csv"] = "parquet"
