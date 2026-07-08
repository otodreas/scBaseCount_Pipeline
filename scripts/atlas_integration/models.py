from pathlib import Path

from pydantic import BaseModel, Field


class MergeStats(BaseModel):
    nAccessionsRequested: int
    nAccessionsMerged: int
    nAccessionsSkipped: int
    nCellsFinal: int
    nGenesFinal: int
    nStudies: int
    skippedAccessions: list[str] = Field(default_factory=list)


class BatchMixingMetrics(BaseModel):
    meanSameStudyNeighborFractionUncorrected: float
    meanSameStudyNeighborFractionCorrected: float


class ClusterConservationMetrics(BaseModel):
    ariUncorrectedVsCorrected: float
    nmiUncorrectedVsCorrected: float
    silhouetteUncorrectedEmbedding: float
    silhouetteCorrectedEmbedding: float


class AtlasIntegrationResult(BaseModel):
    atlasPath: Path
    metadataPath: Path
    mergeStats: MergeStats
    batchMixing: BatchMixingMetrics
    clusterConservation: ClusterConservationMetrics
    nPcs: int
    cumvarPct: float
