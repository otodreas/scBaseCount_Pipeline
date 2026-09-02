from pathlib import Path
from typing import Literal

from cluster_validation import default_resolutions
from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT


class AtlasPostprocessingParameters(BaseModel):
    """Approved tuning knobs imported by validation and production runs."""

    nTopGenes: int
    nPcs: int
    nNeighbors: int
    resolution: float
    calibrationSummary: str | None = None


class AtlasPostprocessingConfig(BaseModel):
    # TODO: improve default values -- input and output should not have hardcoded defaults, figs should end up in output_h5ad.parent/figs/, ...
    inputH5ad: Path = REPO_ROOT / "output" / "atlas" / "v2" / "atlas.h5ad"
    outputH5ad: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "atlas_pp.h5ad"
    figsDir: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "figures"
    batchKey: str = "study_accession"
    cellTypeKey: str = "cell_type"
    nTopGenes: int = 2000
    nPcs: int = 20
    nPcsCompute: int = 50
    nPcsMin: int = 15
    nPcsCumvarTarget: float = 0.5
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
    resolutionCandidates: list[float] = Field(default_factory=default_resolutions)
