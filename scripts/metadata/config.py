from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel


class MetadataConfig(BaseModel):
    sampleParquetPath: Path
    obsParquetPath: Path
    minObsCount: int = 1000
    outputDir: Path = Path("output/metadata")
