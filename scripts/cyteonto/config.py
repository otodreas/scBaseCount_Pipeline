from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel
from shared.repo import REPO_ROOT as _REPO_ROOT


class CyteOntoConfig(BaseModel):
    h5adPath: Path
    payloadDir: Path = _REPO_ROOT / "output" / "cyteonto" / "payloads"
    runsDir: Path = _REPO_ROOT / "output" / "cyteonto" / "runs"
    baseUrl: str = "https://cyteonto.nygen.io"
    pollIntervalS: int = 10
    pollTimeoutS: int = 3600
