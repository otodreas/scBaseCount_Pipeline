from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel

_REPO_ROOT = Path(__file__).parents[2]


class CyteOntoConfig(BaseModel):
    h5adPath: Path
    payloadDir: Path = _REPO_ROOT / "output" / "cyteonto" / "payloads"
    resultsDir: Path = _REPO_ROOT / "output" / "cyteonto" / "runs"
    baseUrl: str = "https://cyteonto.nygen.io"
    pollIntervalS: int = 10
    pollTimeoutS: int = 3600
