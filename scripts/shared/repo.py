from __future__ import annotations

from pathlib import Path


def _find_repo_root(start: Path) -> Path:
    for directory in [start, *start.parents]:
        if (directory / ".venv").exists() or (directory / ".git").exists():
            return directory
    raise RuntimeError(
        f"Could not locate repo root from {start}. "
        "No .venv or .git directory found in any parent."
    )


REPO_ROOT: Path = _find_repo_root(Path(__file__).resolve())
