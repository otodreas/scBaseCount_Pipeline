from urllib.parse import quote

import httpx

from ontology_lookup.tokens import curie_to_iri


class OntologyReleaseMismatchError(RuntimeError):
    """Raised when OLS reports a different ontology release than configured."""


def _as_label(value: object) -> str | None:
    if value is None:
        return None
    if isinstance(value, list):
        if not value:
            return None
        return _as_label(value[0])
    text = str(value).strip()
    return text or None


def fetch_ontology_release(ontologyId: str, *, baseUrl: str, timeoutS: float) -> str:
    resp = httpx.get(f"{baseUrl.rstrip('/')}/ontologies/{ontologyId}", timeout=timeoutS)
    resp.raise_for_status()
    payload = resp.json()
    version = payload.get("http://www.w3.org/2002/07/owl#versionInfo")
    if not version:
        version_iri = payload.get("http://www.w3.org/2002/07/owl#versionIRI")
        if isinstance(version_iri, str) and "/releases/" in version_iri:
            version = version_iri.split("/releases/")[1].split("/")[0]
    if not isinstance(version, str) or not version.strip():
        raise RuntimeError(f"OLS ontology {ontologyId!r} did not report a release version")
    return version.strip()


def assert_release(
    ontologyId: str,
    expectedRelease: str,
    *,
    baseUrl: str,
    timeoutS: float,
) -> str:
    observed = fetch_ontology_release(ontologyId, baseUrl=baseUrl, timeoutS=timeoutS)
    if observed != expectedRelease:
        raise OntologyReleaseMismatchError(
            f"OLS {ontologyId} release {observed!r} does not match configured {expectedRelease!r}"
        )
    return observed


def fetch_term(
    ontologyId: str,
    curie: str,
    *,
    baseUrl: str,
    timeoutS: float,
) -> dict[str, object]:
    iri = curie_to_iri(curie)
    encoded = quote(quote(iri, safe=""), safe="")
    class_url = f"{baseUrl.rstrip('/')}/ontologies/{ontologyId}/classes/{encoded}"
    resp = httpx.get(class_url, timeout=timeoutS)
    if resp.status_code == 404:
        raise KeyError(f"OLS term not found: {ontologyId}:{curie}")
    resp.raise_for_status()
    payload = resp.json()
    label = _as_label(payload.get("label"))
    if label is None:
        raise RuntimeError(f"OLS term {curie} has no preferred label")

    ancestors: list[str] = []
    page = 0
    while True:
        anc_resp = httpx.get(
            f"{class_url}/hierarchicalAncestors",
            params={"page": page, "size": 500},
            timeout=timeoutS,
        )
        anc_resp.raise_for_status()
        body = anc_resp.json()
        elements = body.get("elements") or []
        for element in elements:
            anc_curie = element.get("curie")
            if isinstance(anc_curie, str) and anc_curie.strip():
                ancestors.append(anc_curie.strip())
        total_pages = int(body.get("totalPages") or 1)
        page += 1
        if page >= total_pages:
            break

    # Stable unique order.
    seen: set[str] = set()
    unique_ancestors: list[str] = []
    for anc in ancestors:
        if anc in seen or anc == curie:
            continue
        seen.add(anc)
        unique_ancestors.append(anc)

    return {
        "curie": curie,
        "label": label,
        "ancestors": unique_ancestors,
        "iri": iri,
    }
