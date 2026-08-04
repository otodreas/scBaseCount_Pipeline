from dataclasses import dataclass

from ontology_lookup.cache import OntologyCache, TermRecord
from ontology_lookup.tokens import ontology_tokens

# Broad analysis areas keyed by stable MONDO roots. Order is priority: more specific first.
DISEASE_AREA_ROOTS: tuple[tuple[str, frozenset[str]], ...] = (
    (
        "IPF / Pulmonary Fibrosis",
        frozenset(
            {
                "MONDO:0800504",  # idiopathic pulmonary fibrosis
                "MONDO:0002771",  # pulmonary fibrosis
                "MONDO:0002429",  # idiopathic interstitial pneumonia
            }
        ),
    ),
    (
        "COVID-19 / SARS-CoV-2",
        frozenset(
            {
                "MONDO:0100096",  # COVID-19
                "MONDO:0100320",  # post-COVID-19 disorder
            }
        ),
    ),
    ("Lung Cancer", frozenset({"MONDO:0008903"})),
    ("COPD", frozenset({"MONDO:0005002"})),
    ("Cystic Fibrosis", frozenset({"MONDO:0009061"})),
    ("Pulmonary Hypertension", frozenset({"MONDO:0005149"})),
    ("Interstitial Lung Disease", frozenset({"MONDO:0015925"})),
)

RESPIRATORY_SYSTEM_ID = "UBERON:0001004"
OTHER_AREA = "Other"


@dataclass(frozen=True)
class ResolvedOntologyField:
    raw: str
    ids: tuple[str, ...]
    names: tuple[str, ...]
    status: str  # resolved | partial | missing | empty | invalid


def resolve_ontology_field(raw: str | None, cache: OntologyCache) -> ResolvedOntologyField:
    text = "" if raw is None else str(raw).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return ResolvedOntologyField(raw="", ids=(), names=(), status="empty")
    tokens = ontology_tokens(text)
    if not tokens:
        return ResolvedOntologyField(raw=text, ids=(), names=(), status="invalid")
    records = cache.ensure(tokens, allowNetwork=False)
    ids: list[str] = []
    names: list[str] = []
    missing = False
    for token in tokens:
        record = records.get(token)
        ids.append(token)
        if record is None:
            missing = True
            continue
        names.append(record.label)
    if not names:
        status = "missing"
    elif missing:
        status = "partial"
    else:
        status = "resolved"
    return ResolvedOntologyField(raw=text, ids=tuple(ids), names=tuple(names), status=status)


def area_for_record(record: TermRecord) -> str:
    closure = record.closure
    for area, roots in DISEASE_AREA_ROOTS:
        if not (closure & roots):
            continue
        # Broad neoplasm parents alone are not enough for Lung Cancer.
        if area == "Lung Cancer" and "MONDO:0008903" not in closure:
            continue
        return area
    return OTHER_AREA


def disease_area_from_records(records: list[TermRecord]) -> tuple[str, str]:
    """Return (diseaseArea, source) from resolved MONDO term records."""
    if not records:
        return OTHER_AREA, "missing_ontology"
    areas = [area_for_record(record) for record in records]
    concrete = {area for area in areas if area != OTHER_AREA}
    if not concrete:
        return OTHER_AREA, "unmapped_ontology"
    if len(concrete) > 1:
        return OTHER_AREA, "conflicting_ontology"
    return next(iter(concrete)), "ontology"


def is_respiratory_term(record: TermRecord) -> bool:
    return RESPIRATORY_SYSTEM_ID in record.closure


def respiratory_from_uberon(records: list[TermRecord]) -> bool | None:
    """Return True/False when UBERON evidence exists; None when no resolved terms."""
    if not records:
        return None
    return all(is_respiratory_term(record) for record in records)


def records_for_ids(ids: list[str] | tuple[str, ...], cache: OntologyCache) -> list[TermRecord]:
    resolved = cache.ensure(list(ids), allowNetwork=False)
    return [record for curie in ids if (record := resolved.get(curie)) is not None]
