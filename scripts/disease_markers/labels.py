import csv
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from metadata.regexes import DISEASE_MAP, LUNG_DISEASE_RE, LUNG_TISSUE_RE
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
CONTROL_AREA = "Control"
OTHER_AREA = "Other"
EXCLUDE_RE = re.compile(
    r"organoid|cell line|cell-line|ipsc|iPSC|explant|Mus musculus|\bmouse\b|embryo|olfactory|tonsil|myeloma|CRISPR|perturbation",
    re.IGNORECASE,
)
CONTROL_EXTRA_RE = re.compile(
    r"adjacent|non[-\s]?involved|uninvolved|tumou?r[-\s]?free|para[-\s]?tumou?r|healthy donor",
    re.IGNORECASE,
)
# Control detection deliberately omits the missing-information tokens that live in
# NORMAL_HEALTHY_RE (none, unknown, unsure, not specified, not stated, not reported,
# not available). Those mark absent metadata, not a healthy sample, and matching them
# labelled diseased samples as Control when unrelated fields (sex, race, virus strain)
# happened to be "None" or "Unknown".
CONTROL_RE = re.compile(
    r"\b(?:"
    r"normal|healthy|control|unstimulated|naive|"
    r"uninvolved|unaffected|unexposed|vehicle|wild[-_\s]?type|"
    r"wt|no treatment|baseline"
    r")\b"
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


@dataclass(frozen=True)
class SampleLabelRow:
    srxAccession: str
    studyAccession: str
    diseaseRaw: str
    diseaseArea: str
    isControl: bool
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


def _tissue_blob(ctx: ExperimentContext) -> str:
    bio = ctx.biological
    parts: list[str] = []
    if bio.tissueType:
        parts.append(bio.tissueType)
    for key, value in bio.sampleAttributes.items():
        if not isinstance(value, str):
            continue
        if any(token in key.lower() for token in ("tissue", "organ", "source")):
            parts.append(value)
    return " ".join(parts)


def _control_blob(ctx: ExperimentContext) -> str:
    """Return sample-level text searched for control status.

    Includes disease- and diagnosis-keyed sample attributes plus the sample title and
    description. Excludes the study title, which is shared across every sample in a study
    and otherwise pulled whole mixed cohorts into Control, and excludes generic attributes
    such as sex or race that carry no disease meaning.
    """
    bio = ctx.biological
    parts: list[str] = [
        str(v)
        for k, v in bio.sampleAttributes.items()
        if isinstance(v, str) and ("disease" in k.lower() or "diagnos" in k.lower())
    ]
    if bio.sampleTitle:
        parts.append(bio.sampleTitle)
    if bio.sampleDescription:
        parts.append(bio.sampleDescription)
    return " ".join(parts)


def coarse_disease_area(diseaseText: str, fullText: str = "", controlText: str | None = None) -> str:
    """Return a coarse disease area label from disease text, optional full sample text, and optional control text.

    controlText is the sample-level text searched for control status; when None it falls back to diseaseText.
    """
    control_source = diseaseText if controlText is None else controlText
    if CONTROL_RE.search(control_source) or CONTROL_EXTRA_RE.search(control_source):
        return CONTROL_AREA

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


def _is_eligible(ctx: ExperimentContext, diseaseArea: str) -> tuple[bool, str | None]:
    if ctx.biological.scientificName and ctx.biological.scientificName.strip() != "Homo sapiens":
        return False, "non_human"
    blob = _attribute_blob(ctx)
    if EXCLUDE_RE.search(blob):
        return False, "excluded_sample_type"
    tissue = _tissue_blob(ctx)
    disease = _disease_blob(ctx)
    lung = bool(LUNG_TISSUE_RE.search(tissue) or LUNG_DISEASE_RE.search(disease) or LUNG_TISSUE_RE.search(blob))
    if not lung:
        return False, "non_lung"
    if diseaseArea == OTHER_AREA:
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
        control_text = _control_blob(ctx)
        area = coarse_disease_area(disease_raw, full, control_text)
        is_control = area == CONTROL_AREA
        eligible, exclude_reason = _is_eligible(ctx, area)
        records.append(
            {
                "srxAccession": accession,
                "studyAccession": study_accession,
                "diseaseRaw": disease_raw,
                "diseaseArea": area,
                "isControl": is_control,
                "eligible": eligible,
                "excludeReason": exclude_reason,
            }
        )
    return pd.DataFrame.from_records(records)


def sample_labels_by_srx(label_table: pd.DataFrame) -> dict[str, SampleLabelRow]:
    """Index label rows by SRX accession."""
    out: dict[str, SampleLabelRow] = {}
    for _, row in label_table.iterrows():
        srx = str(row["srxAccession"])
        out[srx] = SampleLabelRow(
            srxAccession=srx,
            studyAccession=str(row["studyAccession"]),
            diseaseRaw=str(row["diseaseRaw"]),
            diseaseArea=str(row["diseaseArea"]),
            isControl=bool(row["isControl"]),
            eligible=bool(row["eligible"]),
            excludeReason=None if pd.isna(row["excludeReason"]) else str(row["excludeReason"]),
        )
    return out
