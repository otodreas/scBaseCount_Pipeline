from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel

from shared.repo import REPO_ROOT as _REPO_ROOT

N_TOP_GENES = 100


class CyteTypeRunnerConfig(BaseModel):
    srxAccession: str
    outputDir: Path = _REPO_ROOT / "output" / "cytetype" / "data"
    jobDetailsDir: Path = _REPO_ROOT / "output" / "cytetype" / "job_details"
