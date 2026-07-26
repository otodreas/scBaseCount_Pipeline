from pathlib import Path
from typing import Literal

from metadata.config import MetadataConfig
from pydantic import BaseModel
from shared.repo import REPO_ROOT


class DiseaseMarkersConfig(BaseModel):
    inputAtlasH5ad: Path = REPO_ROOT / "output" / "atlas" / "v1" / "atlas.h5ad"
    harmonyAtlasH5ad: Path = REPO_ROOT / "output" / "atlas" / "v1" / "processed_1" / "atlas_harmony.h5ad"
    transferredAtlasH5ad: Path | None = (
        REPO_ROOT / "output" / "atlas" / "v1" / "processed_1" / "atlas_with_clusters.h5ad"
    )
    contextsPath: Path = REPO_ROOT / "output" / "context" / "contexts.jsonl"
    atlasCsvPath: Path = REPO_ROOT / "output" / "atlas" / "v1" / "atlas.csv"
    metadataConfig: MetadataConfig = MetadataConfig()
    outputDir: Path = REPO_ROOT / "output" / "atlas" / "v1" / "processed_1" / "disease_markers"
    clusterKey: str = "leiden_atlas"
    sampleKey: str = "SRX_accession"
    studyKey: str = "study_accession"
    minCellsPerProfile: int = 10
    minSamplesPerArea: int = 3
    minStudiesPerArea: int = 2
    padjThreshold: float = 0.05
    lfcThreshold: float = 1.0
    compression: Literal["gzip", "lzf"] | None = "gzip"
    writeTransferredAtlas: bool = True
