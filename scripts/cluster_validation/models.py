from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel


class ClusterValidationResult(BaseModel):
    srxAccession: str
    datasetTitleSuffix: str
    embedding: str = "pca"
    selectedResolution: float
    clusterKey: str
    mergedKey: str
    nPcs: int | None = None
    cumvar: float | None = None
    kPrior: int
    kFiltered: int
    nCellsDropped: int
    nCellsFinal: int
    nClustersPreMerge: int
    nClustersPostMerge: int
    adataPath: Path
    labelMap: dict[str, str]
    mergedGroups: dict[str, list[str]]
    # Per-resolution arrays stored as [[resolution, value], ...] for serialisation
    resolutions: list[float]
    kArr: list[int]
    jaccArr: list[float]
    silhouetteArr: list[list[float]]
    homogeneityArr: list[list[float]]
    completenessArr: list[list[float]]
    nmiArr: list[list[float]]
    vscoreArr: list[list[float]]
    ariArr: list[list[float]]
    confMatrix: list[list[float]]
    confClasses: list[str]
