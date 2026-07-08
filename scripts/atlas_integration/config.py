from pathlib import Path

import numpy as np
from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT as _REPO_ROOT


class AtlasIntegrationConfig(BaseModel):
    datasetsCsvPath: Path = _REPO_ROOT / "output" / "metadata" / "datasets.csv"
    contextsPath: Path = _REPO_ROOT / "output" / "context" / "contexts.jsonl"
    localH5adRoot: Path = _REPO_ROOT / "data" / "scbasecount" / "2026-01-12" / "h5ad" / "GeneFull" / "Homo_sapiens"
    batchKey: str = "study"
    accessionKey: str = "accession"
    cellTypeKey: str = "cell_type"
    missingLabel: str = "unknown"
    minGenesPerCell: int = 200
    minCellsPerGene: int = 3
    maxPctCountsMt: float = 20.0
    minCellsTotal: int = 500
    nTopGenes: int = 2000
    nPcsCompute: int = 50
    nPcsMin: int = 15
    nPcsCumvarTarget: float = 0.5
    harmonyMaxIter: int = 20
    leidenResolution: float = 1.0
    leidenKeyUncorrected: str = "leiden_uncorrected"
    leidenKeyAtlas: str = "leiden_atlas"
    umapKeyUncorrected: str = "X_umap_uncorrected"
    umapKeyCorrected: str = "X_umap"
    pcaKey: str = "X_pca"
    harmonyKey: str = "X_pca_harmony"
    resolutions: list[float] = Field(default_factory=lambda: [round(r, 1) for r in np.arange(0.1, 2.0, 0.1).tolist()])
    subsampleN: int | None = None
    outputDir: Path = _REPO_ROOT / "output" / "atlas"
    figsDir: Path = _REPO_ROOT / "output" / "atlas" / "figs"
    atlasH5adName: str = "lung_atlas.h5ad"
    metricsSampleSize: int = 5000
