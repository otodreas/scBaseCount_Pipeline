from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT as _REPO_ROOT


class CyteOntoConfig(BaseModel):
    h5adPath: Path
    authorCol: str = "cell_type"
    algorithmCols: dict[str, str] = Field(default_factory=lambda: {"algo1": "cytetype_annotation_leiden_merged"})
    payloadDir: Path = _REPO_ROOT / "output" / "cyteonto" / "payloads"
    runsDir: Path = _REPO_ROOT / "output" / "cyteonto" / "runs"
    baseUrl: str = "https://cyteonto.nygen.io"
    pollIntervalS: int = 10
    pollTimeoutS: int = 3600
