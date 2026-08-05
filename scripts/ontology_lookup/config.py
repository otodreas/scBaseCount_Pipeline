from pathlib import Path

from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT as _REPO_ROOT


class OntologyLookupConfig(BaseModel):
    olsBaseUrl: str = "https://www.ebi.ac.uk/ols4/api/v2"
    mondoOntologyId: str = "mondo"
    uberonOntologyId: str = "uberon"
    mondoRelease: str = "2026-07-06"
    uberonRelease: str = "2026-06-19"
    cacheDir: Path = Field(default_factory=lambda: _REPO_ROOT / "data" / "ontologies")
    requestTimeoutS: float = 60.0
