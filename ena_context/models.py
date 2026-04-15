from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


class TechnicalContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    instrumentModel: str | None = None
    instrumentPlatform: str | None = None
    libraryStrategy: str | None = None
    librarySource: str | None = None
    librarySelection: str | None = None
    libraryLayout: str | None = None
    libraryConstructionProtocol: str | None = None


class BiologicalContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scientificName: str | None = None
    taxId: str | None = None
    strain: str | None = None
    cellType: str | None = None
    tissueType: str | None = None
    sampleTitle: str | None = None
    sampleDescription: str | None = None
    sampleAttributes: dict[str, str] = Field(default_factory=dict)


class ExperimentContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    accession: str
    studyAccession: str | None = None
    studyTitle: str | None = None
    experimentTitle: str | None = None
    sampleAccession: str | None = None
    runAccessions: list[str] = Field(default_factory=list)
    technical: TechnicalContext = Field(default_factory=TechnicalContext)
    biological: BiologicalContext = Field(default_factory=BiologicalContext)
    warnings: list[str] = Field(default_factory=list)
