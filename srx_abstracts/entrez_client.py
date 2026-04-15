from __future__ import annotations

import os
import threading
import time
from typing import Any

import httpx

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class EntrezClient:
    """HTTP client for NCBI E-utilities with conservative pacing and retries."""

    def __init__(
        self,
        *,
        api_key: str | None = None,
        email: str | None = None,
        tool: str | None = None,
        timeout_s: float = 60.0,
    ) -> None:
        self.api_key = api_key if api_key is not None else os.environ.get("NCBI_API_KEY")
        self.email = email if email is not None else os.environ.get("NCBI_EMAIL")
        self.tool = tool if tool is not None else os.environ.get("ENTREZ_TOOL", "srx_abstracts")
        self._min_interval = 0.11 if self.api_key else 0.35
        self._last_request_mono = 0.0
        self._lock = threading.Lock()
        self._client = httpx.Client(timeout=timeout_s)

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> EntrezClient:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def _baseParams(self) -> dict[str, str]:
        p: dict[str, str] = {"tool": self.tool}
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        return p

    def _pace(self) -> None:
        with self._lock:
            now = time.monotonic()
            wait = self._min_interval - (now - self._last_request_mono)
            if wait > 0:
                time.sleep(wait)
            self._last_request_mono = time.monotonic()

    def _get(self, endpoint: str, params: dict[str, str]) -> str:
        merged = {**self._baseParams(), **params}
        url = f"{EUTILS_BASE}/{endpoint}"
        last: httpx.Response | None = None
        for attempt in range(4):
            self._pace()
            last = self._client.get(url, params=merged)
            if last.status_code in (429, 500, 502, 503, 504):
                time.sleep(min(2**attempt, 8))
                continue
            last.raise_for_status()
            return last.text
        assert last is not None
        last.raise_for_status()
        return last.text

    def esearch(self, *, db: str, term: str, retmax: str = "20") -> str:
        return self._get(
            "esearch.fcgi",
            {"db": db, "term": term, "retmax": retmax, "retmode": "xml"},
        )

    def elink(self, *, dbfrom: str, db: str, id_param: str) -> str:
        return self._get(
            "elink.fcgi",
            {"dbfrom": dbfrom, "db": db, "id": id_param, "retmode": "xml"},
        )

    def efetch(self, *, db: str, id_param: str, rettype: str = "xml") -> str:
        return self._get(
            "efetch.fcgi",
            {"db": db, "id": id_param, "rettype": rettype, "retmode": "xml"},
        )
