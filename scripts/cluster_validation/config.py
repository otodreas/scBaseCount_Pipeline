from __future__ import annotations

from pathlib import Path

import numpy as np
from pydantic import BaseModel, Field


class ClusterValidationConfig(BaseModel):
    srxAccession: str | None = None
    datasetIndex: int = 2
    summaryPath: Path = Path("cluster_validation_sandbox/datasets_summary")
    localH5adRoot: Path = Path("data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens")
    minCellsPerType: int = 20
    nTopGenes: int = 2000
    nPcsCompute: int = 50
    nPcsMin: int = 15
    nPcsCumvarTarget: float = 0.5
    resolutions: list[float] = Field(
        default_factory=lambda: [round(r, 1) for r in np.arange(0.1, 2.0, 0.1).tolist()]
    )
    mergeThreshold: float = 0.2
    rfBalanceWeakPrior: bool = False
    outputDir: Path = Path("data/other")
