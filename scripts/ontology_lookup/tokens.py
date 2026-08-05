import re

CURIE_RE = re.compile(r"^(?P<prefix>[A-Za-z][A-Za-z0-9_]*):(?P<local>\d+)$")


def ontology_tokens(raw: str | None) -> list[str]:
    """Split a comma-separated ontology field into stable, unique CURIEs."""
    if raw is None:
        return []
    text = str(raw).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return []
    seen: set[str] = set()
    out: list[str] = []
    for part in text.split(","):
        token = part.strip()
        if not token:
            continue
        match = CURIE_RE.match(token)
        if match is None:
            continue
        curie = f"{match.group('prefix').upper()}:{match.group('local')}"
        if curie in seen:
            continue
        seen.add(curie)
        out.append(curie)
    return out


def curie_to_iri(curie: str) -> str:
    match = CURIE_RE.match(curie.strip())
    if match is None:
        raise ValueError(f"invalid CURIE: {curie!r}")
    prefix = match.group("prefix").upper()
    local = match.group("local")
    return f"http://purl.obolibrary.org/obo/{prefix}_{local}"
