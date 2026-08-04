import re

from ontology_lookup import OntologyCache, records_for_ids, resolve_ontology_field, respiratory_from_uberon
from ontology_lookup.tokens import ontology_tokens

from metadata.regexes import LUNG_TISSUE_RE

EXCLUDED_MODEL_RE = re.compile(
    r"\b(?:organoids?|explants?|ipscs?)\b|cell[-\s]?line|"
    r"\b(?:A549|PC-?9|WI-?38|Calu-?3|BEAS-?2B|H358|H-?23|HCC\d+|RUES2|H2228|H1975|H838)\b",
    re.IGNORECASE,
)
PRIMARY_NOT_MODEL_RE = re.compile(
    r"\bprimary\b|"
    r"(?:no|not(?:\s+an?)?)\s+(?:immortalized\s+|established\s+)?cell\s+line",
    re.IGNORECASE,
)
NON_LUNG_MIXED_RE = re.compile(
    r"\b(?:PBMCs?|peripheral\s+blood(?:\s+mononuclear\s+cells?)?|whole\s+blood|"
    r"leukocytes?\s+from\s+whole\s+blood)\b",
    re.IGNORECASE,
)
EXPLICIT_NON_RESPIRATORY_RE = re.compile(
    r"\b(?:PBMCs?|peripheral\s+blood|whole\s+blood|bone\s+marrow|lymph\s+nodes?|"
    r"hippocamp\w*|caudate|brain|cortex|cerebr\w*|olfactory|kidney|renal|liver|hepatic|"
    r"skin|colon|intestin\w*|heart|cardiac|breast|mammary|placenta|testis|testicular|"
    r"umbilical\s+cord|retina|retinal)\b",
    re.IGNORECASE,
)


def is_excluded_model(tissue: str, cellLine: str) -> bool:
    blob = f"{tissue} {cellLine}".strip()
    if not blob:
        return False
    if PRIMARY_NOT_MODEL_RE.search(blob) and not re.search(r"\b(?:organoids?|explants?|ipscs?)\b", blob, re.I):
        return False
    return bool(EXCLUDED_MODEL_RE.search(blob))


def is_respiratory_tissue(tissue: str) -> bool:
    """Text fallback when no usable UBERON term is available."""
    if not tissue or not LUNG_TISSUE_RE.search(tissue):
        return False
    return not bool(NON_LUNG_MIXED_RE.search(tissue))


def is_explicit_non_respiratory_text(text: str) -> bool:
    if not text:
        return False
    if LUNG_TISSUE_RE.search(text) and not NON_LUNG_MIXED_RE.search(text):
        return False
    return bool(EXPLICIT_NON_RESPIRATORY_RE.search(text) or NON_LUNG_MIXED_RE.search(text))


def is_respiratory_specimen(tissue: str, tissueOntology: str, uberonCache: OntologyCache) -> bool:
    """Prefer UBERON relationship closure; fall back to tissue text only when needed."""
    resolved = resolve_ontology_field(tissueOntology, uberonCache)
    records = records_for_ids(resolved.ids, uberonCache)
    respiratory = respiratory_from_uberon(records)
    if respiratory is not None:
        return respiratory
    # Partial/invalid ontology with no resolved terms: if raw tokens exist but unresolved, do not
    # invent respiratory status from missing IDs; use text fallback.
    if ontology_tokens(tissueOntology) and resolved.status == "missing":
        return is_respiratory_tissue(tissue)
    return is_respiratory_tissue(tissue)
