from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel
from shared.repo import REPO_ROOT as _REPO_ROOT


class AnnotationInspectConfig(BaseModel):
    inputPrefix: str
    cyteontoPrefix: str
    topN: int = 10
    downloadRoot: Path = _REPO_ROOT / "data" / "annotation_inspection"
    outputDir: Path = _REPO_ROOT / "output" / "annotation_inspection_pipeline"
    emitUmap: bool = True
    emitExtremes: bool = True
