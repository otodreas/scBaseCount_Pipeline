import json
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from metadata.categorize import coarse_disease_area
from metadata.specimen import is_excluded_model, is_respiratory_specimen
from ontology_lookup import (
    OntologyCache,
    OntologyLookupConfig,
    disease_area_from_records,
    records_for_ids,
    resolve_ontology_field,
)
from study_context.evidence import context_exclude_reason
from study_context.models import ExperimentContext
from study_context.utils import load_contexts_jsonl

from disease_markers.status import (
    ComparatorFamily,
    ControlType,
    detect_comparator_candidate,
    has_paired_opposite,
    infection_comparator_sets_nondiseased,
    infer_disease_status,
)

OTHER_AREA = "Other"
REQUIRED_METADATA_COLUMNS: frozenset[str] = frozenset(
    {
        "srx_accession",
        "disease",
        "disease_ontology_term_id",
        "tissue",
        "tissue_ontology_term_id",
        "organism",
        "cell_line",
    }
)


@dataclass
class _DraftLabel:
    srxAccession: str
    studyAccession: str
    diseaseRaw: str
    diseaseOntologyTermId: str
    diseaseName: str
    tissueRaw: str
    tissueOntologyRaw: str
    cellLineRaw: str
    diseaseArea: str
    diseaseAreaSource: str
    diseased: bool | None
    controlType: str | None
    excludeReasonDraft: str | None
    comparatorFamily: ComparatorFamily | None
    comparatorField: str | None
    context: ExperimentContext


def _as_text(value: object) -> str:
    if value is None or value is pd.NA:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def atlas_success_accessions(atlasManifestPath: Path) -> dict[str, str]:
    """Return {experiment_accession: study_accession} for successful files in the result manifest."""
    payload = json.loads(atlasManifestPath.read_text())
    if "files" not in payload:
        raise ValueError(f"{atlasManifestPath} missing files[]; expected a post-cutover atlas result manifest")
    out: dict[str, str] = {}
    for row in payload.get("files", []):
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


def _specimen_exclude_reason(
    *,
    organism: str,
    tissue: str,
    tissueOntology: str,
    cellLine: str,
    disease: str,
    diseaseArea: str,
    ctx: ExperimentContext,
    uberonCache: OntologyCache,
) -> str | None:
    if organism and organism.strip() != "Homo sapiens":
        return "non_human"
    context_reason = context_exclude_reason(
        parquetOrganism=organism,
        parquetTissue=tissue,
        parquetDisease=disease,
        parquetDiseaseArea=diseaseArea,
        ctx=ctx,
    )
    if context_reason is not None:
        return context_reason
    if is_excluded_model(tissue, cellLine):
        return "excluded_sample_type"
    if not is_respiratory_specimen(tissue, tissueOntology, uberonCache):
        return "non_lung"
    return None


def _eligibility(
    *,
    excludeReason: str | None,
    diseaseArea: str,
    diseased: bool | None,
) -> tuple[bool, str | None]:
    if excludeReason is not None:
        return False, excludeReason
    if diseaseArea == OTHER_AREA and diseased is not False:
        return False, "unmapped_disease"
    return True, None


def _derive_disease_area(
    *,
    diseaseOntologyRaw: str,
    diseaseRaw: str,
    mondoCache: OntologyCache,
) -> tuple[str, str, str, tuple[str, ...]]:
    resolved = resolve_ontology_field(diseaseOntologyRaw, mondoCache)
    disease_name = "; ".join(resolved.names)
    records = records_for_ids(resolved.ids, mondoCache)
    if records:
        area, source = disease_area_from_records(records)
        if area != OTHER_AREA:
            return area, source, disease_name, resolved.ids
        if source == "conflicting_ontology":
            return OTHER_AREA, source, disease_name, resolved.ids
    if diseaseRaw:
        text_area = coarse_disease_area(diseaseRaw)
        if text_area != OTHER_AREA:
            return text_area, "disease_text", disease_name, resolved.ids
    if resolved.status in {"empty", "invalid", "missing"}:
        return OTHER_AREA, "missing_ontology", disease_name, resolved.ids
    return OTHER_AREA, "unmapped_ontology", disease_name, resolved.ids


def _apply_study_consensus(draft: list[_DraftLabel]) -> None:
    by_study: dict[str, list[_DraftLabel]] = {}
    for record in draft:
        by_study.setdefault(record.studyAccession, []).append(record)

    for study_records in by_study.values():
        peer_areas = {
            record.diseaseArea
            for record in study_records
            if record.diseaseArea != OTHER_AREA and record.diseaseAreaSource in {"ontology", "disease_text"}
        }
        if len(peer_areas) != 1:
            continue
        consensus = next(iter(peer_areas))
        for record in study_records:
            if record.diseaseArea == OTHER_AREA and record.diseaseAreaSource in {
                "missing_ontology",
                "unmapped_ontology",
            }:
                record.diseaseArea = consensus
                record.diseaseAreaSource = "study_consensus"


def build_sample_label_table(
    contextsPath: Path,
    atlasManifestPath: Path,
    sampleMetadataPath: Path,
    *,
    ontologyConfig: OntologyLookupConfig | None = None,
) -> pd.DataFrame:
    """Build per-SRX labels from ontology metadata, with narrow text fallbacks."""
    contexts = load_contexts_jsonl(contextsPath)
    atlas_rows = atlas_success_accessions(atlasManifestPath)
    metadata = _load_sample_metadata(sampleMetadataPath)
    cfg = ontologyConfig or OntologyLookupConfig()
    mondo_cache = OntologyCache.for_mondo(cfg)
    uberon_cache = OntologyCache.for_uberon(cfg)

    missing_metadata = sorted(accession for accession in atlas_rows if accession not in metadata.index)
    if missing_metadata:
        raise KeyError(
            f"sample metadata missing {len(missing_metadata)} atlas accession(s); examples: {missing_metadata[:5]}"
        )
    missing_contexts = sorted(accession for accession in atlas_rows if accession not in contexts)
    if missing_contexts:
        raise KeyError(f"contexts missing {len(missing_contexts)} atlas accession(s); examples: {missing_contexts[:5]}")

    draft: list[_DraftLabel] = []
    for accession, study_accession in sorted(atlas_rows.items()):
        ctx = contexts[accession]
        context_study = None if ctx.study is None else ctx.study.studyAccession
        if context_study is not None and context_study != study_accession:
            raise ValueError(
                f"study accession mismatch for {accession}: atlas={study_accession!r} context={context_study!r}"
            )

        row = metadata.loc[accession]
        disease_raw = _as_text(row["disease"])
        disease_ontology_raw = _as_text(row["disease_ontology_term_id"])
        tissue_raw = _as_text(row["tissue"])
        tissue_ontology_raw = _as_text(row["tissue_ontology_term_id"])
        cell_line_raw = _as_text(row["cell_line"])
        organism = _as_text(row["organism"])

        area, area_source, disease_name, _mondo_ids = _derive_disease_area(
            diseaseOntologyRaw=disease_ontology_raw,
            diseaseRaw=disease_raw,
            mondoCache=mondo_cache,
        )
        exclude_reason = _specimen_exclude_reason(
            organism=organism,
            tissue=tissue_raw,
            tissueOntology=tissue_ontology_raw,
            cellLine=cell_line_raw,
            disease=disease_raw,
            diseaseArea=area,
            ctx=ctx,
            uberonCache=uberon_cache,
        )
        diseased, control_type = infer_disease_status(
            disease_raw,
            ctx,
            diseaseKnown=area != OTHER_AREA,
        )
        if exclude_reason == "disease_mismatch":
            diseased, control_type = None, None
        family, field_key = detect_comparator_candidate(ctx)
        draft.append(
            _DraftLabel(
                srxAccession=accession,
                studyAccession=study_accession,
                diseaseRaw=disease_raw,
                diseaseOntologyTermId=disease_ontology_raw,
                diseaseName=disease_name,
                tissueRaw=tissue_raw,
                tissueOntologyRaw=tissue_ontology_raw,
                cellLineRaw=cell_line_raw,
                diseaseArea=area,
                diseaseAreaSource=area_source,
                diseased=diseased,
                controlType=None if control_type is None else control_type.value,
                excludeReasonDraft=exclude_reason,
                comparatorFamily=family,
                comparatorField=field_key,
                context=ctx,
            )
        )

    _apply_study_consensus(draft)

    by_study: dict[str, list[_DraftLabel]] = {}
    for record in draft:
        by_study.setdefault(record.studyAccession, []).append(record)

    records: list[dict[str, object]] = []
    for study_records in by_study.values():
        for record in study_records:
            diseased = record.diseased
            control_type = record.controlType
            is_experimental_comparator = False
            if record.comparatorFamily is not None and record.comparatorField is not None:
                peers = [peer.context for peer in study_records if peer.srxAccession != record.srxAccession]
                if has_paired_opposite(
                    family=record.comparatorFamily,
                    fieldKey=record.comparatorField,
                    peerContexts=peers,
                ):
                    is_experimental_comparator = True
                    if infection_comparator_sets_nondiseased(record.comparatorFamily) and diseased is not False:
                        diseased = False
                        control_type = None

            eligible, exclude_reason = _eligibility(
                excludeReason=record.excludeReasonDraft,
                diseaseArea=record.diseaseArea,
                diseased=diseased,
            )
            records.append(
                {
                    "srxAccession": record.srxAccession,
                    "studyAccession": record.studyAccession,
                    "diseaseRaw": record.diseaseRaw,
                    "diseaseOntologyTermId": record.diseaseOntologyTermId or None,
                    "diseaseName": record.diseaseName or None,
                    "tissueRaw": record.tissueRaw,
                    "tissueOntologyRaw": record.tissueOntologyRaw,
                    "cellLineRaw": record.cellLineRaw,
                    "diseaseArea": record.diseaseArea,
                    "diseaseAreaSource": record.diseaseAreaSource,
                    "diseased": diseased,
                    "isBiologicalControl": diseased is False and control_type is not None,
                    "controlType": control_type,
                    "isExperimentalComparator": is_experimental_comparator,
                    "eligible": eligible,
                    "excludeReason": exclude_reason,
                }
            )

    table = pd.DataFrame.from_records(records)
    if not table.empty:
        table = table.sort_values("srxAccession").reset_index(drop=True)
        table["diseased"] = table["diseased"].astype("boolean")
    return table


__all__ = [
    "ControlType",
    "OTHER_AREA",
    "atlas_success_accessions",
    "build_sample_label_table",
]
