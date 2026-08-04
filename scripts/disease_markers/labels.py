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
REQUIRED_METADATA_COLUMNS: frozenset[str] = frozenset(
    {
        "srx_accession",
        "disease",
        "tissue",
        "organism",
        "cell_line",
        "perturbation",
    }
)
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
ENA_SPECIMEN_KEY_TOKENS: tuple[str, ...] = (
    "sampling site",
    "tissue",
    "organism part",
    "source name",
    "source",
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
    tissueRaw: str
    cellLineRaw: str
    diseaseArea: str
    diseased: bool | None
    isBiologicalControl: bool
    controlType: ControlType | None
    eligible: bool
    excludeReason: str | None


def _normalized_key(key: str) -> str:
    return re.sub(r"[_-]+", " ", key).strip().lower()


def _as_text(value: object) -> str:
    if value is None or value is pd.NA:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _control_type(text: str) -> ControlType | None:
    if MATCHED_ADJACENT_RE.search(text):
        return ControlType.MATCHED_ADJACENT
    if HEALTHY_RE.search(text) or NONE_DISEASE_RE.search(text):
        return ControlType.HEALTHY
    if EXPLICIT_CONTROL_RE.search(text):
        return ControlType.EXPLICIT_CONTROL
    return None


def _has_disease_category(text: str) -> bool:
    return any(pattern.search(text) for _, pattern in DISEASE_MAP)


def _ena_specimen_texts(ctx: ExperimentContext | None) -> list[str]:
    if ctx is None:
        return []
    bio = ctx.biological
    parts: list[str] = []
    if bio.tissueType:
        parts.append(bio.tissueType)
    if bio.sampleTitle:
        parts.append(bio.sampleTitle)
    for key, value in bio.sampleAttributes.items():
        if not isinstance(value, str):
            continue
        normalized_key = _normalized_key(key)
        if any(token in normalized_key for token in ENA_SPECIMEN_KEY_TOKENS):
            parts.append(value)
    return parts


def _disease_status_from_text(text: str) -> tuple[bool | None, ControlType | None]:
    control_type = _control_type(text)
    has_disease = _has_disease_category(text)
    if control_type is ControlType.MATCHED_ADJACENT:
        return False, control_type
    if control_type is not None and has_disease:
        return None, None
    if control_type is not None:
        return False, control_type
    if has_disease:
        return True, None
    return None, None


def _disease_status_from_parquet(disease: str) -> tuple[bool | None, ControlType | None]:
    text = disease.strip()
    if not text:
        return None, None
    if MATCHED_ADJACENT_RE.search(text):
        return False, ControlType.MATCHED_ADJACENT
    if UNKNOWN_DISEASE_RE.search(text) and not (
        HEALTHY_RE.search(text) or EXPLICIT_CONTROL_RE.search(text) or NONE_DISEASE_RE.search(text)
    ):
        return None, None
    return _disease_status_from_text(text)


def _disease_status(
    disease: str,
    ctx: ExperimentContext | None,
) -> tuple[bool | None, ControlType | None]:
    for text in _ena_specimen_texts(ctx):
        diseased, control_type = _disease_status_from_text(text)
        if diseased is not None or control_type is not None:
            return diseased, control_type
    return _disease_status_from_parquet(disease)


def _is_excluded_model(tissue: str, cell_line: str) -> bool:
    blob = f"{tissue} {cell_line}".strip()
    if not blob:
        return False
    if PRIMARY_NOT_MODEL_RE.search(blob) and not re.search(r"\b(?:organoids?|explants?|ipscs?)\b", blob, re.I):
        return False
    return bool(EXCLUDED_MODEL_RE.search(blob))


def _is_lung_tissue(tissue: str) -> bool:
    if not tissue or not LUNG_TISSUE_RE.search(tissue):
        return False
    return not bool(NON_LUNG_MIXED_RE.search(tissue))


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
    *,
    organism: str,
    tissue: str,
    cell_line: str,
    diseaseArea: str,
    diseased: bool | None,
) -> tuple[bool, str | None]:
    if organism and organism.strip() != "Homo sapiens":
        return False, "non_human"
    if _is_excluded_model(tissue, cell_line):
        return False, "excluded_sample_type"
    if not _is_lung_tissue(tissue):
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


def _load_sample_metadata(sampleMetadataPath: Path) -> pd.DataFrame:
    metadata = pd.read_parquet(sampleMetadataPath)
    missing = REQUIRED_METADATA_COLUMNS - set(metadata.columns)
    if missing:
        raise KeyError(f"sample metadata missing required columns: {sorted(missing)}")
    if not metadata["srx_accession"].is_unique:
        duplicated = metadata.loc[metadata["srx_accession"].duplicated(), "srx_accession"].tolist()
        raise ValueError(f"sample metadata has duplicate srx_accession values: {duplicated[:5]}")
    return metadata.set_index("srx_accession", drop=False)


def build_sample_label_table(
    contextsPath: Path,
    atlasCsvPath: Path,
    sampleMetadataPath: Path,
) -> pd.DataFrame:
    """Build per-SRX labels from scBaseCount sample metadata and atlas success rows."""
    contexts = load_contexts_jsonl(contextsPath)
    atlas_rows = atlas_success_accessions(atlasCsvPath)
    metadata = _load_sample_metadata(sampleMetadataPath)

    missing = sorted(accession for accession in atlas_rows if accession not in metadata.index)
    if missing:
        raise KeyError(f"sample metadata missing {len(missing)} atlas accession(s); examples: {missing[:5]}")

    records: list[dict[str, object]] = []
    for accession, study_accession in sorted(atlas_rows.items()):
        row = metadata.loc[accession]
        disease_raw = _as_text(row["disease"])
        tissue_raw = _as_text(row["tissue"])
        cell_line_raw = _as_text(row["cell_line"])
        organism = _as_text(row["organism"])
        area = coarse_disease_area(disease_raw)
        diseased, control_type = _disease_status(disease_raw, contexts.get(accession))
        is_biological_control = diseased is False and control_type is not None
        eligible, exclude_reason = _is_eligible(
            organism=organism,
            tissue=tissue_raw,
            cell_line=cell_line_raw,
            diseaseArea=area,
            diseased=diseased,
        )
        records.append(
            {
                "srxAccession": accession,
                "studyAccession": study_accession,
                "diseaseRaw": disease_raw,
                "tissueRaw": tissue_raw,
                "cellLineRaw": cell_line_raw,
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
            tissueRaw=str(row["tissueRaw"]),
            cellLineRaw=str(row["cellLineRaw"]),
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
