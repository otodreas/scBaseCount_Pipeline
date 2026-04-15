from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


class PmidProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pmid: str
    sourceHop: str


class PubmedArticle(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pmid: str
    title: str = ""
    abstractText: str | None = None
    journal: str | None = None
    year: int | None = None


class SrxResolution(BaseModel):
    model_config = ConfigDict(extra="forbid")

    srxAccession: str
    sraUids: list[str] = Field(default_factory=list)
    pmidProvenance: list[PmidProvenance] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


class SrxAbstractBundle(BaseModel):
    model_config = ConfigDict(extra="forbid")

    srxAccession: str
    sraUids: list[str] = Field(default_factory=list)
    pmids: list[str] = Field(default_factory=list)
    pmidProvenance: list[PmidProvenance] = Field(default_factory=list)
    articles: list[PubmedArticle] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
