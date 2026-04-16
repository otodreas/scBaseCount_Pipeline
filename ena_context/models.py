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


class PubmedArticle(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pmid: str
    title: str | None = None
    abstract: str | None = None
    journal: str | None = None
    year: int | None = None


class StudyContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    studyAccession: str
    studyTitle: str | None = None
    studyDescription: str | None = None
    geoAccession: str | None = None
    pubmedIds: list[str] = Field(default_factory=list)
    articles: list[PubmedArticle] = Field(default_factory=list)


class ExperimentContext(BaseModel):
    model_config = ConfigDict(extra="forbid")

    accession: str
    experimentTitle: str | None = None
    sampleAccession: str | None = None
    runAccessions: list[str] = Field(default_factory=list)
    technical: TechnicalContext = Field(default_factory=TechnicalContext)
    biological: BiologicalContext = Field(default_factory=BiologicalContext)
    study: StudyContext | None = None
    warnings: list[str] = Field(default_factory=list)
