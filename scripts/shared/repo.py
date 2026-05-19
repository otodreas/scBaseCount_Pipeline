from __future__ import annotations

from pathlib import Path


def _find_repo_root(start: Path) -> Path:
    for directory in [start, *start.parents]:
        if (directory / ".venv").exists() or (directory / ".git").exists():
            return directory
    raise RuntimeError(f"Could not locate repo root from {start}. No .venv or .git directory found in any parent.")


REPO_ROOT: Path = _find_repo_root(Path(__file__).resolve())


def rel_to_repo(p: Path) -> str:
    """Return p as a string relative to REPO_ROOT when possible, else the absolute path."""
    try:
        return str(Path(p).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(Path(p).resolve())
