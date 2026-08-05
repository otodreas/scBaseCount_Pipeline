import re
from enum import StrEnum

from metadata.categorize import disease_categories_for
from study_context.evidence import comparator_evidence, status_evidence, tissue_evidence
from study_context.models import ExperimentContext

MATCHED_ADJACENT_RE = re.compile(
    r"\b(?:adjacent|non[-\s]?involved|uninvolved|tumou?r[-\s]?free|para[-\s]?tumou?r|tumou?r[-\s]?distant)\b",
    re.IGNORECASE,
)
HEALTHY_RE = re.compile(
    r"\b(?:normal|healthy|unaffected|disease[-\s]?free|non[-\s]?diseased?)\b"
    r"|"
    r"\bno\s+(?:"
    r"disease|COPD|"
    r"diagnosed\s+disease|specific\s+disease|overt\s+disease|donor\s+disease|"
    r"disease\s+diagnosis|record\s+of\s+lung\s+disease"
    r")\b",
    re.IGNORECASE,
)
EXPLICIT_CONTROL_RE = re.compile(r"\bcontrol\b", re.IGNORECASE)
NONE_DISEASE_RE = re.compile(r"\bnone\b", re.IGNORECASE)
UNKNOWN_DISEASE_RE = re.compile(
    r"\b(?:unsure|unknown|not\s+specified|not\s+stated|not\s+reported|not\s+available|"
    r"none\s+specified|none\s+reported)\b",
    re.IGNORECASE,
)

WT_RE = re.compile(r"\b(?:WT|wild[-\s]?type)\b", re.IGNORECASE)
MOCK_RE = re.compile(r"\bmock\b", re.IGNORECASE)
UNINFECTED_RE = re.compile(r"\buninfected\b", re.IGNORECASE)
UNEXPOSED_RE = re.compile(r"\bunexposed\b", re.IGNORECASE)
VEHICLE_RE = re.compile(r"\bvehicle\b", re.IGNORECASE)
UNTREATED_RE = re.compile(r"\b(?:untreated|no\s+treatment)\b", re.IGNORECASE)
BASELINE_RE = re.compile(r"\bbaseline\b", re.IGNORECASE)
NAIVE_RE = re.compile(r"\bnaive\b", re.IGNORECASE)

ALTERED_GENOTYPE_RE = re.compile(
    r"\b(?:knockout|knock[- ]?in|mutant|mutation|deletion|deficient|overexpress|transgenic)\b",
    re.IGNORECASE,
)
INFECTION_RE = re.compile(
    r"\b(?:infect(?:ed|ion)?|expos(?:ed|ure)|SARS|COVID|virus|viral|MOI)\b",
    re.IGNORECASE,
)
INTERVENTION_RE = re.compile(
    r"\b(?:treat(?:ed|ment)?|drug|inhibitor|dexamethasone|chemotherapy|nivolumab|"
    r"stimulation|stimulated|dose|mg/?kg)\b",
    re.IGNORECASE,
)
FOLLOWUP_RE = re.compile(
    r"\b(?:post[-\s]?treatment|follow[-\s]?up|progression|after\s+treatment|day\s+\d+|DPI)\b",
    re.IGNORECASE,
)
ACTIVATED_RE = re.compile(
    r"\b(?:activat(?:ed|ion)|challenged|post[-\s]?infection|\d+DPI)\b",
    re.IGNORECASE,
)


class ControlType(StrEnum):
    MATCHED_ADJACENT = "matchedAdjacent"
    HEALTHY = "healthy"
    EXPLICIT_CONTROL = "explicitControl"


class ComparatorFamily(StrEnum):
    WT = "wt"
    MOCK_INFECTION = "mockInfection"
    VEHICLE_TREATMENT = "vehicleTreatment"
    BASELINE = "baseline"
    NAIVE = "naive"


def _control_type(text: str) -> ControlType | None:
    if MATCHED_ADJACENT_RE.search(text):
        return ControlType.MATCHED_ADJACENT
    if HEALTHY_RE.search(text) or NONE_DISEASE_RE.search(text):
        return ControlType.HEALTHY
    if EXPLICIT_CONTROL_RE.search(text):
        return ControlType.EXPLICIT_CONTROL
    return None


def _has_disease_category(text: str) -> bool:
    # Fallback only for unstructured disease text when ontology did not establish a disease area.
    return bool(disease_categories_for(text))


def _status_from_text(
    text: str,
    *,
    diseaseKnown: bool | None = None,
) -> tuple[bool | None, ControlType | None]:
    control_type = _control_type(text)
    has_disease = bool(diseaseKnown) if diseaseKnown is not None else _has_disease_category(text)
    if control_type is ControlType.MATCHED_ADJACENT:
        return False, control_type
    if control_type is not None and has_disease:
        return None, None
    if control_type is not None:
        return False, control_type
    if has_disease:
        return True, None
    return None, None


def _specimen_status_texts(ctx: ExperimentContext | None) -> list[str]:
    if ctx is None:
        return []
    texts = [value for _, value in tissue_evidence(ctx)]
    texts.extend(value for _, value in status_evidence(ctx))
    return texts


def infer_disease_status(
    disease: str,
    ctx: ExperimentContext | None,
    *,
    diseaseKnown: bool | None = None,
) -> tuple[bool | None, ControlType | None]:
    for text in _specimen_status_texts(ctx):
        diseased, control_type = _status_from_text(text, diseaseKnown=None)
        if diseased is not None or control_type is not None:
            return diseased, control_type

    text = disease.strip()
    if not text:
        if diseaseKnown:
            return True, None
        return None, None
    if MATCHED_ADJACENT_RE.search(text):
        return False, ControlType.MATCHED_ADJACENT
    if UNKNOWN_DISEASE_RE.search(text) and not (
        HEALTHY_RE.search(text) or EXPLICIT_CONTROL_RE.search(text) or NONE_DISEASE_RE.search(text)
    ):
        if diseaseKnown:
            return True, None
        return None, None
    return _status_from_text(text, diseaseKnown=diseaseKnown)


def _candidate_family(text: str) -> ComparatorFamily | None:
    if WT_RE.search(text):
        return ComparatorFamily.WT
    if MOCK_RE.search(text) or UNINFECTED_RE.search(text) or UNEXPOSED_RE.search(text):
        return ComparatorFamily.MOCK_INFECTION
    if VEHICLE_RE.search(text) or UNTREATED_RE.search(text):
        return ComparatorFamily.VEHICLE_TREATMENT
    if BASELINE_RE.search(text):
        return ComparatorFamily.BASELINE
    if NAIVE_RE.search(text):
        return ComparatorFamily.NAIVE
    return None


def _opposite_state(family: ComparatorFamily, text: str) -> bool:
    if family is ComparatorFamily.WT:
        remainder = WT_RE.sub(" ", text)
        return bool(ALTERED_GENOTYPE_RE.search(remainder))
    if family is ComparatorFamily.MOCK_INFECTION:
        remainder = MOCK_RE.sub(" ", text)
        remainder = UNINFECTED_RE.sub(" ", remainder)
        remainder = UNEXPOSED_RE.sub(" ", remainder)
        return bool(INFECTION_RE.search(remainder))
    if family is ComparatorFamily.VEHICLE_TREATMENT:
        remainder = VEHICLE_RE.sub(" ", UNTREATED_RE.sub(" ", text))
        return bool(INTERVENTION_RE.search(remainder))
    if family is ComparatorFamily.BASELINE:
        remainder = BASELINE_RE.sub(" ", text)
        return bool(FOLLOWUP_RE.search(remainder))
    if family is ComparatorFamily.NAIVE:
        remainder = NAIVE_RE.sub(" ", text)
        return bool(ACTIVATED_RE.search(remainder))
    return False


def _is_mixed_arm(text: str, family: ComparatorFamily) -> bool:
    return _candidate_family(text) is family and _opposite_state(family, text)


def detect_comparator_candidate(
    ctx: ExperimentContext | None,
) -> tuple[ComparatorFamily | None, str | None]:
    if ctx is None:
        return None, None
    for key, value in comparator_evidence(ctx):
        family = _candidate_family(value)
        if family is None:
            continue
        if _is_mixed_arm(value, family):
            return None, None
        return family, key
    return None, None


def has_paired_opposite(
    *,
    family: ComparatorFamily,
    fieldKey: str,
    peerContexts: list[ExperimentContext],
) -> bool:
    for peer in peerContexts:
        for key, value in comparator_evidence(peer):
            if fieldKey != "sampleTitle" and key != fieldKey:
                continue
            if fieldKey == "sampleTitle" and key != "sampleTitle":
                continue
            if _opposite_state(family, value) and not _is_mixed_arm(value, family):
                return True
    return False


def infection_comparator_sets_nondiseased(family: ComparatorFamily) -> bool:
    return family is ComparatorFamily.MOCK_INFECTION
