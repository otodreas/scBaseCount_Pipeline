import csv
import math
import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import pandas as pd
from metadata.regexes import DISEASE_MAP, LUNG_TISSUE_RE
from study_context.models import ExperimentContext
from study_context.utils import load_contexts_jsonl

COARSE_TOP_LEVEL: frozenset[str] = frozenset(
    {
        "IPF / Pulmonary Fibrosis",
        "COVID-19 / SARS-CoV-2",
        "Lung Cancer",
        "COPD",
        "Cystic Fibrosis",
        "Interstitial Lung Disease",
        "Pulmonary Hypertension",
    }
)
LUNG_CANCER_LABELS: frozenset[str] = frozenset(
    {
        "Lung Cancer",
        "Small Cell Lung Cancer (SCLC)",
        "Non-small Cell Lung Cancer (NSCLC)",
        "Lung Adenocarcinoma (LUAD)",
        "Lung Squamous Cell Carcinoma (LUSC)",
        "Lung Large Cell Carcinoma (LCC)",
    }
)
OTHER_AREA = "Other"
EXCLUDE_RE = re.compile(
    r"organoid|cell line|cell-line|ipsc|iPSC|explant|Mus musculus|\bmouse\b|embryo|olfactory|tonsil|myeloma|CRISPR|perturbation",
    re.IGNORECASE,
)
MATCHED_ADJACENT_RE = re.compile(
    r"\b(?:adjacent|non[-\s]?involved|uninvolved|tumou?r[-\s]?free|para[-\s]?tumou?r|tumou?r[-\s]?distant)\b",
    re.IGNORECASE,
)
HEALTHY_RE = re.compile(
    r"\b(?:normal|healthy|unaffected|disease[-\s]?free)\b"
    r"|"
    r"\bno\s+(?:"
    r"disease|COPD|"
    r"diagnosed\s+disease|specific\s+disease|overt\s+disease|donor\s+disease|"
    r"disease\s+diagnosis|record\s+of\s+lung\s+disease"
    r")\b"
    r"|"
    r"\bnon[-\s]?(?:disease|COPD)\b",
    re.IGNORECASE,
)
EXPLICIT_CONTROL_RE = re.compile(r"\bcontrol\b", re.IGNORECASE)
POSITIVE_SPECIMEN_RE = re.compile(
    r"\b(?:tumou?r|cancer|carcinoma|adenocarcinoma|malignan(?:t|cy)|neoplasm|fibrotic|infected|diseased)\b",
    re.IGNORECASE,
)
NON_LUNG_TISSUE_RE = re.compile(
    r"\b(?:blood|PBMCs?|peripheral blood mononuclear cells?|lymph nodes?|bone marrow)\b",
    re.IGNORECASE,
)
MIXED_SAMPLE_RE = re.compile(
    r"\bmixed sample\b|\bmultiple donors?\b|\bpooled donors?\b|\bpooled\b.*\bdonors?\b",
    re.IGNORECASE,
)
SPECIMEN_KEY_TOKENS: tuple[str, ...] = (
    "sampling site",
    "tissue",
    "source",
    "specimen",
    "sample type",
)
ANATOMY_KEY_TOKENS: tuple[str, ...] = ("organ", "anatom")
STATUS_KEY_TOKENS: tuple[str, ...] = (
    "disease",
    "diagnos",
    "condition",
    "status",
    "phenotype",
    "cohort",
    "health",
    "strain",
    "isolate",
)
MIXED_SAMPLE_ATTRIBUTE_KEYS: frozenset[str] = frozenset(
    {
        "donor",
        "donors",
        "individual",
        "pool",
        "pooled",
        "sample type",
    }
)


class ControlType(StrEnum):
    MATCHED_ADJACENT = "matchedAdjacent"
    HEALTHY = "healthy"
    EXPLICIT_CONTROL = "explicitControl"


@dataclass(frozen=True)
class SampleLabelRow:
    srxAccession: str
    studyAccession: str
    diseaseRaw: str
    diseaseArea: str
    diseased: bool | None
    isBiologicalControl: bool
    controlType: ControlType | None
    eligible: bool
    excludeReason: str | None


def _attribute_blob(ctx: ExperimentContext) -> str:
    bio = ctx.biological
    parts: list[str] = []
    if bio.sampleTitle:
        parts.append(bio.sampleTitle)
    if bio.tissueType:
        parts.append(bio.tissueType)
    if bio.sampleDescription:
        parts.append(bio.sampleDescription)
    parts.extend(str(v) for v in bio.sampleAttributes.values())
    study = ctx.study
    if study and study.studyTitle:
        parts.append(study.studyTitle)
    return " ".join(parts)


def _disease_blob(ctx: ExperimentContext) -> str:
    attrs = ctx.biological.sampleAttributes
    parts = [
        str(v) for k, v in attrs.items() if isinstance(v, str) and ("disease" in k.lower() or "diagnos" in k.lower())
    ]
    study = ctx.study
    if study and study.studyTitle:
        parts.append(study.studyTitle)
    return " ".join(parts)


def _normalized_key(key: str) -> str:
    return re.sub(r"[_-]+", " ", key).strip().lower()


def _tissue_blob(ctx: ExperimentContext) -> str:
    bio = ctx.biological
    parts: list[str] = []
    if bio.tissueType:
        parts.append(bio.tissueType)
    for key, value in bio.sampleAttributes.items():
        if not isinstance(value, str):
            continue
        normalized_key = _normalized_key(key)
        if any(token in normalized_key for token in (*SPECIMEN_KEY_TOKENS, *ANATOMY_KEY_TOKENS)):
            parts.append(value)
    return " ".join(parts)


def _control_type(text: str) -> ControlType | None:
    if MATCHED_ADJACENT_RE.search(text):
        return ControlType.MATCHED_ADJACENT
    if HEALTHY_RE.search(text):
        return ControlType.HEALTHY
    if EXPLICIT_CONTROL_RE.search(text):
        return ControlType.EXPLICIT_CONTROL
    return None


def _has_positive_disease_evidence(text: str) -> bool:
    return bool(POSITIVE_SPECIMEN_RE.search(text) or any(pattern.search(text) for _, pattern in DISEASE_MAP))


def _classify_evidence_tier(texts: list[str]) -> tuple[bool, bool | None, ControlType | None]:
    control_types: set[ControlType] = set()
    positive = False
    for text in texts:
        control_type = _control_type(text)
        if control_type is not None:
            control_types.add(control_type)
        elif _has_positive_disease_evidence(text):
            positive = True

    if not control_types and not positive:
        return False, None, None
    if control_types and positive:
        return True, None, None
    if positive:
        return True, True, None

    precedence = (
        ControlType.MATCHED_ADJACENT,
        ControlType.HEALTHY,
        ControlType.EXPLICIT_CONTROL,
    )
    control_type = next(candidate for candidate in precedence if candidate in control_types)
    return True, False, control_type


def _disease_status(ctx: ExperimentContext) -> tuple[bool | None, ControlType | None]:
    bio = ctx.biological
    specimen_texts: list[str] = []
    anatomy_texts: list[str] = []
    status_texts: list[str] = []

    if bio.tissueType:
        specimen_texts.append(bio.tissueType)
    for key, value in bio.sampleAttributes.items():
        if not isinstance(value, str):
            continue
        normalized_key = _normalized_key(key)
        if any(token in normalized_key for token in SPECIMEN_KEY_TOKENS):
            specimen_texts.append(value)
        elif any(token in normalized_key for token in ANATOMY_KEY_TOKENS):
            anatomy_texts.append(value)
        elif any(token in normalized_key for token in STATUS_KEY_TOKENS):
            status_texts.append(value)

    evidence_tiers = (specimen_texts, anatomy_texts, status_texts)
    if bio.sampleTitle:
        evidence_tiers = (*evidence_tiers, [bio.sampleTitle])

    for texts in evidence_tiers:
        has_evidence, diseased, control_type = _classify_evidence_tier(texts)
        if has_evidence:
            return diseased, control_type
    return None, None


def _is_mixed_or_pooled(ctx: ExperimentContext) -> bool:
    bio = ctx.biological
    parts = [value for value in (ctx.experimentTitle, bio.sampleTitle) if value]
    for key, value in bio.sampleAttributes.items():
        if isinstance(value, str) and _normalized_key(key) in MIXED_SAMPLE_ATTRIBUTE_KEYS:
            parts.append(value)
    return bool(MIXED_SAMPLE_RE.search(" ".join(parts)))


def coarse_disease_area(diseaseText: str, fullText: str = "") -> str:
    """Return a coarse disease area independently of specimen disease status."""
    text = diseaseText if diseaseText.strip() else fullText
    matched: list[str] = []
    for label, pattern in DISEASE_MAP:
        if pattern.search(text):
            matched.append(label)
    if not matched:
        return OTHER_AREA
    for label in matched:
        if label in COARSE_TOP_LEVEL:
            return label
    if any(label in LUNG_CANCER_LABELS for label in matched):
        return "Lung Cancer"
    return OTHER_AREA


def _is_eligible(
    ctx: ExperimentContext,
    diseaseArea: str,
    diseased: bool | None,
) -> tuple[bool, str | None]:
    if ctx.biological.scientificName and ctx.biological.scientificName.strip() != "Homo sapiens":
        return False, "non_human"
    blob = _attribute_blob(ctx)
    if EXCLUDE_RE.search(blob):
        return False, "excluded_sample_type"
    if _is_mixed_or_pooled(ctx):
        return False, "mixed_sample"
    tissue = _tissue_blob(ctx)
    if NON_LUNG_TISSUE_RE.search(tissue) or not LUNG_TISSUE_RE.search(tissue):
        return False, "non_lung"
    if diseaseArea == OTHER_AREA and diseased is not False:
        return False, "unmapped_disease"
    return True, None


def atlas_success_accessions(atlasCsvPath: Path) -> dict[str, str]:
    """Return {experiment_accession: study_accession} for rows with status success."""
    out: dict[str, str] = {}
    with atlasCsvPath.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "success":
                continue
            accession = row.get("accession") or ""
            study = row.get("studyAccession") or ""
            if accession:
                out[accession] = study
    return out


def build_sample_label_table(
    contextsPath: Path,
    atlasCsvPath: Path,
) -> pd.DataFrame:
    """Build per-SRX labels for atlas samples from contexts and coarse disease rules."""
    contexts = load_contexts_jsonl(contextsPath)
    atlas_rows = atlas_success_accessions(atlasCsvPath)
    records: list[dict[str, object]] = []
    for accession, study_accession in sorted(atlas_rows.items()):
        ctx = contexts.get(accession)
        if ctx is None:
            continue
        disease_raw = _disease_blob(ctx)
        full = _attribute_blob(ctx)
        area = coarse_disease_area(disease_raw, full)
        diseased, control_type = _disease_status(ctx)
        is_biological_control = diseased is False and control_type is not None
        eligible, exclude_reason = _is_eligible(ctx, area, diseased)
        records.append(
            {
                "srxAccession": accession,
                "studyAccession": study_accession,
                "diseaseRaw": disease_raw,
                "diseaseArea": area,
                "diseased": diseased,
                "isBiologicalControl": is_biological_control,
                "controlType": None if control_type is None else control_type.value,
                "eligible": eligible,
                "excludeReason": exclude_reason,
            }
        )
    table = pd.DataFrame.from_records(records)
    if not table.empty:
        table["diseased"] = table["diseased"].astype("boolean")
    return table


def sample_labels_by_srx(label_table: pd.DataFrame) -> dict[str, SampleLabelRow]:
    """Index label rows by SRX accession."""
    out: dict[str, SampleLabelRow] = {}
    for _, row in label_table.iterrows():
        srx = str(row["srxAccession"])
        diseased_value = row["diseased"]
        control_type_value = row["controlType"]
        exclude_reason_value = row["excludeReason"]
        out[srx] = SampleLabelRow(
            srxAccession=srx,
            studyAccession=str(row["studyAccession"]),
            diseaseRaw=str(row["diseaseRaw"]),
            diseaseArea=str(row["diseaseArea"]),
            diseased=None
            if diseased_value is None
            or diseased_value is pd.NA
            or (isinstance(diseased_value, float) and math.isnan(diseased_value))
            else bool(diseased_value),
            isBiologicalControl=bool(row["isBiologicalControl"]),
            controlType=None
            if control_type_value is None
            or control_type_value is pd.NA
            or (isinstance(control_type_value, float) and math.isnan(control_type_value))
            else ControlType(str(control_type_value)),
            eligible=bool(row["eligible"]),
            excludeReason=None
            if exclude_reason_value is None
            or exclude_reason_value is pd.NA
            or (isinstance(exclude_reason_value, float) and math.isnan(exclude_reason_value))
            else str(exclude_reason_value),
        )
    return out
