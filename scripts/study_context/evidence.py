import re

from metadata.categorize import disease_categories_for
from metadata.specimen import is_explicit_non_respiratory_text

from study_context.models import ExperimentContext

TISSUE_ATTRIBUTE_KEYS: frozenset[str] = frozenset(
    {
        "tissue",
        "organism part",
        "sampling site",
    }
)
STATUS_ATTRIBUTE_KEYS: frozenset[str] = frozenset(
    {
        "disease",
        "disease state",
        "disease state/diagnosis",
        "diagnosis",
        "condition",
        "clinic status",
        "subject status",
        "health state",
        "infection status",
        "infection",
    }
)
COMPARATOR_ATTRIBUTE_KEYS: frozenset[str] = frozenset(
    {
        "genotype",
        "phenotype",
        "treatment",
        "kinase inhibitor treatment",
        "agent",
        "status",
        "time",
        "infection status",
        "cohort",
        "group",
    }
)
HUMAN_MARKERS: frozenset[str] = frozenset({"homo sapiens", "9606", "human"})
COMPATIBLE_DISEASE_AREAS: frozenset[frozenset[str]] = frozenset(
    {
        frozenset({"IPF / Pulmonary Fibrosis", "Interstitial Lung Disease"}),
    }
)
COARSE_AREAS: frozenset[str] = frozenset(
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


def normalized_key(key: str) -> str:
    return re.sub(r"[_-]+", " ", key).strip().lower()


def attribute_values(ctx: ExperimentContext, keys: frozenset[str]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for key, value in ctx.biological.sampleAttributes.items():
        if not isinstance(value, str) or not value.strip():
            continue
        normalized = normalized_key(key)
        if normalized in keys:
            out.append((normalized, value.strip()))
    return out


def tissue_evidence(ctx: ExperimentContext) -> list[tuple[str, str]]:
    parts: list[tuple[str, str]] = []
    if ctx.biological.tissueType and ctx.biological.tissueType.strip():
        parts.append(("tissueType", ctx.biological.tissueType.strip()))
    parts.extend(attribute_values(ctx, TISSUE_ATTRIBUTE_KEYS))
    return parts


def status_evidence(ctx: ExperimentContext) -> list[tuple[str, str]]:
    return attribute_values(ctx, STATUS_ATTRIBUTE_KEYS)


def comparator_evidence(ctx: ExperimentContext) -> list[tuple[str, str]]:
    parts: list[tuple[str, str]] = []
    if ctx.biological.sampleTitle and ctx.biological.sampleTitle.strip():
        parts.append(("sampleTitle", ctx.biological.sampleTitle.strip()))
    parts.extend(attribute_values(ctx, COMPARATOR_ATTRIBUTE_KEYS))
    return parts


def context_is_nonhuman(ctx: ExperimentContext) -> bool:
    markers = {
        value.strip().lower()
        for value in (ctx.biological.scientificName, ctx.biological.taxId)
        if value and value.strip()
    }
    return bool(markers) and markers.isdisjoint(HUMAN_MARKERS)


def coarse_areas_from_text(text: str) -> set[str]:
    labels = set(disease_categories_for(text))
    out = {label for label in labels if label in COARSE_AREAS}
    if labels & LUNG_CANCER_LABELS:
        out.add("Lung Cancer")
    return out


def disease_areas_compatible(left: set[str], right: set[str]) -> bool:
    if not left or not right:
        return True
    if left & right:
        return True
    for pair in COMPATIBLE_DISEASE_AREAS:
        if left <= pair and right <= pair:
            return True
    return False


def context_exclude_reason(
    *,
    parquetOrganism: str,
    parquetTissue: str,
    parquetDisease: str,
    parquetDiseaseArea: str,
    ctx: ExperimentContext,
) -> str | None:
    del parquetOrganism, parquetTissue
    if context_is_nonhuman(ctx):
        return "non_human"

    for _, value in tissue_evidence(ctx):
        if is_explicit_non_respiratory_text(value):
            return "non_lung"

    if parquetDiseaseArea and parquetDiseaseArea != "Other":
        parquet_areas = {parquetDiseaseArea}
    else:
        parquet_areas = coarse_areas_from_text(parquetDisease)
    context_areas: set[str] = set()
    for _, value in status_evidence(ctx):
        context_areas |= coarse_areas_from_text(value)
    if parquet_areas and context_areas and not disease_areas_compatible(parquet_areas, context_areas):
        return "disease_mismatch"
    return None
