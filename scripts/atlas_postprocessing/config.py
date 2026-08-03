from pathlib import Path
from typing import Literal

import numpy as np
from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT


def _default_resolutions() -> list[float]:
    return [round(float(r), 1) for r in np.arange(0.2, 2.01, 0.2)]


class AtlasPostprocessingParameters(BaseModel):
    """Approved tuning knobs imported by validation and production runs."""

    nTopGenes: int
    nPcs: int
    nNeighbors: int
    resolution: float
    calibrationSummary: str | None = None


class AtlasPostprocessingConfig(BaseModel):
    inputH5ad: Path = REPO_ROOT / "output" / "atlas" / "v2" / "atlas.h5ad"
    outputH5ad: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "atlas_pp.h5ad"
    figsDir: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "figures"
    batchKey: str = "study_accession"
    cellTypeKey: str = "cell_type"
    nTopGenes: int = 2000
    nPcs: int = 20
    nPcsCompute: int = 50
    nNeighbors: int = 15
    resolution: float = 1.0
    scaleMaxValue: float = 10.0
    writePlots: bool = True
    compression: Literal["gzip", "lzf"] | None = "gzip"
    r2Key: str | None = None
    parametersJson: Path | None = None
    nJobs: int = 0
    calibrationDir: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "parameter_selection"
    validationDir: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "subset_validation"
    hvgCandidates: list[int] = Field(default_factory=lambda: [1000, 2000, 4000, 8000])
    pcCandidates: list[int] = Field(default_factory=lambda: [10, 20, 30, 50])
    neighborCandidates: list[int] = Field(default_factory=lambda: [5, 10, 15, 30, 50, 100])
    resolutionCandidates: list[float] = Field(default_factory=_default_resolutions)
    plateauRelativeThreshold: float = 0.95
