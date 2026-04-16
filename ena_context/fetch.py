from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET

from ena_context.ena_client import curl_get
from ena_context.models import BiologicalContext, ExperimentContext, StudyContext, TechnicalContext

PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
BROWSER_BASE = "https://www.ebi.ac.uk/ena/browser/api"


def _str(val: str | None) -> str | None:
    if not val or not val.strip():
        return None
    return val.strip()


def _parse_pubmed_ids(tag: str) -> list[str]:
    return re.findall(r"xref:PubMed:(\d+)", tag)


def _parse_sample_attributes(xml_text: str) -> dict[str, str]:
    root = ET.fromstring(xml_text)
    attrs: dict[str, str] = {}
    for attr in root.iter("SAMPLE_ATTRIBUTE"):
        tag_el = attr.find("TAG")
        val_el = attr.find("VALUE")
        if tag_el is not None and tag_el.text:
            attrs[tag_el.text] = val_el.text.strip() if val_el is not None and val_el.text else ""
    return attrs


def _fetch_study_context(study_accession: str, warnings: list[str]) -> StudyContext | None:
    url = (
        f"{PORTAL_BASE}/filereport"
        f"?accession={study_accession}&result=study&fields=all&format=json"
    )
    try:
        raw = curl_get(url)
        records: list[dict[str, str]] = json.loads(raw)
    except Exception as exc:
        warnings.append(f"study_api_failed:{exc}")
        return None

    if not records:
        warnings.append(f"study_api_empty_response:{study_accession}")
        return None

    r = records[0]
    return StudyContext(
        studyAccession=study_accession,
        studyTitle=_str(r.get("study_title")),
        studyDescription=_str(r.get("study_description")),
        geoAccession=_str(r.get("geo_accession")),
        pubmedIds=_parse_pubmed_ids(r.get("tag", "")),
    )


def fetch_experiment_context(accession: str) -> ExperimentContext:
    print(f"[DEBUG] Fetching experiment context for accession: {accession}")
    warnings: list[str] = []

    url = (
        f"{PORTAL_BASE}/filereport"
        f"?accession={accession}&result=read_experiment&fields=all&format=json"
    )
    try:
        raw = curl_get(url)
        records: list[dict[str, str]] = json.loads(raw)
    except Exception as exc:
        warnings.append(f"portal_api_failed:{exc}")
        return ExperimentContext(accession=accession, warnings=warnings)

    if not records:
        warnings.append("portal_api_empty_response")
        return ExperimentContext(accession=accession, warnings=warnings)

    first = records[0]
    run_accessions = [r["run_accession"] for r in records if r.get("run_accession")]

    technical = TechnicalContext(
        instrumentModel=_str(first.get("instrument_model")),
        instrumentPlatform=_str(first.get("instrument_platform")),
        libraryStrategy=_str(first.get("library_strategy")),
        librarySource=_str(first.get("library_source")),
        librarySelection=_str(first.get("library_selection")),
        libraryLayout=_str(first.get("library_layout")),
        libraryConstructionProtocol=_str(first.get("library_construction_protocol")),
    )

    sample_accession = _str(first.get("sample_accession"))

    biological = BiologicalContext(
        scientificName=_str(first.get("scientific_name")),
        taxId=_str(first.get("tax_id")),
        strain=_str(first.get("strain")),
        cellType=_str(first.get("cell_type")),
        tissueType=_str(first.get("tissue_type")),
        sampleTitle=_str(first.get("sample_title")),
        sampleDescription=_str(first.get("sample_description")),
    )

    if sample_accession:
        try:
            xml_text = curl_get(f"{BROWSER_BASE}/xml/{sample_accession}")
            biological = biological.model_copy(
                update={"sampleAttributes": _parse_sample_attributes(xml_text)}
            )
        except Exception as exc:
            warnings.append(f"sample_xml_failed:{exc}")

    study_accession = _str(first.get("study_accession"))
    study = _fetch_study_context(study_accession, warnings) if study_accession else None

    return ExperimentContext(
        accession=accession,
        experimentTitle=_str(first.get("experiment_title")),
        sampleAccession=sample_accession,
        runAccessions=run_accessions,
        technical=technical,
        biological=biological,
        study=study,
        warnings=warnings,
    )
