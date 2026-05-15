from __future__ import annotations

from pathlib import Path

import numpy as np
from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT as _REPO_ROOT


class ClusterValidationConfig(BaseModel):
    srxAccession: str | None = None
    datasetIndex: int = 2
    summaryPath: Path = _REPO_ROOT / "output" / "metadata" / "quantiles_datasets.csv"
    localH5adRoot: Path = _REPO_ROOT / "data" / "scbasecount" / "2026-01-12" / "h5ad" / "GeneFull" / "Homo_sapiens"
    minCellsPerType: int = 20
    # Minimum number of cells to keep after filtering
    minCellsTotal: int = 500
    # Default to 2000 HVGs (highly variable genes). see "curse of dimensionality"
    nTopGenes: int = 2000
    nPcsCompute: int = 50
    nPcsMin: int = 15
    nPcsCumvarTarget: float = 0.5
    resolutions: list[float] = Field(default_factory=lambda: [round(r, 1) for r in np.arange(0.1, 2.0, 0.1).tolist()])
    mergeThreshold: float = 0.2
    rfBalanceWeakPrior: bool = False
    outputDir: Path = _REPO_ROOT / "output" / "clustering" / "data"
    figsDir: Path = _REPO_ROOT / "output" / "clustering" / "figs"
