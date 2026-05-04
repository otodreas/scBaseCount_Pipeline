from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel


def _find_repo_root(start: Path) -> Path:
    for directory in [start, *start.parents]:
        if (directory / ".venv").exists() or (directory / ".git").exists():
            return directory
    raise RuntimeError(
        f"Could not locate repo root from {start}. "
        "No .venv or .git directory found in any parent."
    )


_REPO_ROOT = _find_repo_root(Path(__file__).resolve())


class CyteOntoConfig(BaseModel):
    h5adPath: Path
    payloadDir: Path = _REPO_ROOT / "output" / "cyteonto" / "payloads"
    runsDir: Path = _REPO_ROOT / "output" / "cyteonto" / "runs"
    baseUrl: str = "https://cyteonto.nygen.io"
    pollIntervalS: int = 10
    pollTimeoutS: int = 3600
