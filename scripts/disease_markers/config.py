"""Configuration for atlas disease DE and noteworthy-gene discovery."""

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT


class AtlasDeAnalysisConfig(BaseModel):
    atlasPath: Path = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "atlas_v2_post.h5ad"
    outputDir: Path = REPO_ROOT / "output" / "atlas" / "v2" / "analysis" / "production"
    contextsPath: Path = REPO_ROOT / "output" / "context" / "contexts_v2.jsonl"
    atlasManifestPath: Path = REPO_ROOT / "output" / "atlas" / "v2" / "atlas_v2_result.json"
    geneInfoPath: Path = (
        REPO_ROOT
        / "data"
        / "scbasecount"
        / "2026-01-12"
        / "star_references"
        / "Homo_sapiens"
        / "hg38_2020"
        / "geneInfo.tab"
    )

    sampleKey: str = "SRX_accession"
    studyKey: str = "study_accession"
    clusterKey: str = "leiden_atlas"
    labelKey: str = "cell_type"
    ontologyKey: str = "cell_ontology_term_id"
    diseaseNameKey: str = "diseaseName"

    minCellsPerProfile: int = Field(default=10, ge=1)
    minOverlapStudies: int = Field(default=2, ge=1)
    padj: float = Field(default=0.05, gt=0.0, le=1.0)
    lfc: float = Field(default=1.0, ge=0.0)
    minDetectionDelta: float = Field(default=0.15, ge=0.0, le=1.0)
    minTau: float = Field(default=0.8, ge=0.0, le=1.0)
    minTargetDetection: float = Field(default=0.2, ge=0.0, le=1.0)
    maxBackgroundDetection: float = Field(default=0.05, ge=0.0, le=1.0)
    minStudiesForSpecificity: int = Field(default=3, ge=1)
    minProfilesForGene: int = Field(default=5, ge=1)
    minTotalCountsForGene: int = Field(default=20, ge=0)
    highPurity: float = Field(default=0.7, ge=0.0, le=1.0)
    resolvedMinStudies: int = Field(default=3, ge=1)

    primaryBudget: int = Field(default=20, ge=1)
    extendedBudget: int = Field(default=60, ge=1)
    maxPerClassPrimary: int = Field(default=6, ge=1)
    maxPerClassExtended: int = Field(default=15, ge=1)
    maxPerGene: int = Field(default=2, ge=1)
    maxPerCluster: int = Field(default=4, ge=1)
    maxPerDiseaseArea: int = Field(default=8, ge=1)
    geneChunkSize: int = Field(default=2000, ge=100)
    maxVolcanoPlots: int = Field(default=8, ge=0)
    maxEvidencePanels: int = Field(default=20, ge=0)

    # Reserve ~256 GiB on the 2 TiB server so concurrent jobs and OS caches stay safe.
    memoryReserveBytes: int = Field(default=256 * 1024**3, ge=0)
    compression: Literal["gzip", "lzf"] | None = "gzip"

    @property
    def figuresDir(self) -> Path:
        return self.outputDir / "figures"

    @property
    def checkpointsDir(self) -> Path:
        return self.outputDir / "checkpoints"

    @property
    def deCheckpointsDir(self) -> Path:
        return self.checkpointsDir / "de_contrasts"

    @property
    def pseudobulkPath(self) -> Path:
        return self.checkpointsDir / "pseudobulk.h5ad"

    @property
    def fingerprintPath(self) -> Path:
        return self.checkpointsDir / "aggregate_fingerprint.json"
