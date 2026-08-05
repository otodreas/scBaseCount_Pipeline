import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from ontology_lookup.client import assert_release, fetch_term
from ontology_lookup.config import OntologyLookupConfig


@dataclass(frozen=True)
class TermRecord:
    curie: str
    label: str
    ancestors: tuple[str, ...]

    @property
    def closure(self) -> frozenset[str]:
        return frozenset((self.curie, *self.ancestors))


class OntologyCache:
    """Offline term cache for one ontology release."""

    def __init__(
        self,
        *,
        ontologyId: str,
        release: str,
        cacheDir: Path,
        baseUrl: str,
        timeoutS: float,
        maxWorkers: int = 8,
    ) -> None:
        self.ontologyId = ontologyId
        self.release = release
        self.root = cacheDir / ontologyId / release
        self.termsPath = self.root / "terms.json"
        self.manifestPath = self.root / "manifest.json"
        self.baseUrl = baseUrl
        self.timeoutS = timeoutS
        self.maxWorkers = maxWorkers
        self._terms: dict[str, TermRecord] = {}
        if self.termsPath.exists():
            self._load()

    @classmethod
    def for_mondo(cls, cfg: OntologyLookupConfig) -> "OntologyCache":
        return cls(
            ontologyId=cfg.mondoOntologyId,
            release=cfg.mondoRelease,
            cacheDir=cfg.cacheDir,
            baseUrl=cfg.olsBaseUrl,
            timeoutS=cfg.requestTimeoutS,
        )

    @classmethod
    def for_uberon(cls, cfg: OntologyLookupConfig) -> "OntologyCache":
        return cls(
            ontologyId=cfg.uberonOntologyId,
            release=cfg.uberonRelease,
            cacheDir=cfg.cacheDir,
            baseUrl=cfg.olsBaseUrl,
            timeoutS=cfg.requestTimeoutS,
        )

    def _load(self) -> None:
        payload = json.loads(self.termsPath.read_text())
        terms: dict[str, TermRecord] = {}
        for curie, row in payload.items():
            terms[curie] = TermRecord(
                curie=curie,
                label=str(row["label"]),
                ancestors=tuple(row.get("ancestors") or ()),
            )
        self._terms = terms
        if self.manifestPath.exists():
            manifest = json.loads(self.manifestPath.read_text())
            cached_release = manifest.get("release")
            if cached_release != self.release:
                raise RuntimeError(f"cache release {cached_release!r} does not match configured {self.release!r}")

    def get(self, curie: str) -> TermRecord | None:
        return self._terms.get(curie)

    def ensure(self, curies: list[str], *, allowNetwork: bool = False) -> dict[str, TermRecord | None]:
        """Return records for curies, optionally fetching missing ones via OLS."""
        missing = [curie for curie in dict.fromkeys(curies) if curie not in self._terms]
        if missing and not allowNetwork:
            return {curie: self._terms.get(curie) for curie in curies}
        if missing:
            assert_release(
                self.ontologyId,
                self.release,
                baseUrl=self.baseUrl,
                timeoutS=self.timeoutS,
            )
            fetched = 0
            with ThreadPoolExecutor(max_workers=self.maxWorkers) as pool:
                futures = {
                    pool.submit(
                        fetch_term,
                        self.ontologyId,
                        curie,
                        baseUrl=self.baseUrl,
                        timeoutS=self.timeoutS,
                    ): curie
                    for curie in missing
                }
                for future in as_completed(futures):
                    curie = futures[future]
                    try:
                        row = future.result()
                    except KeyError:
                        continue
                    self._terms[curie] = TermRecord(
                        curie=curie,
                        label=str(row["label"]),
                        ancestors=tuple(row["ancestors"]),  # type: ignore[arg-type]
                    )
                    fetched += 1
                    if fetched % 25 == 0 or fetched == len(missing):
                        print(f"{self.ontologyId}: cached {fetched}/{len(missing)} terms", flush=True)
            self.write()
        return {curie: self._terms.get(curie) for curie in curies}

    def write(self) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        payload = {
            curie: {"label": record.label, "ancestors": list(record.ancestors)}
            for curie, record in sorted(self._terms.items())
        }
        text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        self.termsPath.write_text(text)
        digest = hashlib.sha256(text.encode()).hexdigest()
        manifest = {
            "ontologyId": self.ontologyId,
            "release": self.release,
            "sourceUrl": self.baseUrl,
            "generatedAt": datetime.now(UTC).isoformat(),
            "termCount": len(self._terms),
            "contentSha256": digest,
        }
        self.manifestPath.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
