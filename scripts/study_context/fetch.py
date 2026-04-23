from __future__ import annotations

import datetime
import json
import logging
import os
import re
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import httpx

from study_context.models import BiologicalContext, ExperimentContext, StudyContext, TechnicalContext

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
BROWSER_BASE = "https://www.ebi.ac.uk/ena/browser/api"

_LOG_PATH = Path(__file__).parents[2] / "logs" / "study_context.log"
_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(_LOG_PATH, mode="a", encoding="utf-8"),
    ],
)
_log = logging.getLogger(__name__)

_http = httpx.Client(timeout=30.0, follow_redirects=True)



def _http_get(url: str, *, retries: int = 3) -> str:
    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            r = _http.get(url)
            r.raise_for_status()
            return r.text
        except Exception as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"HTTP GET failed for {url!r}: {last_exc}")


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


def _fetch_pubmed_abstract(pmids: list[str], warnings: list[str]) -> str | None:
    if not pmids:
        return None

    print(f"[{datetime.datetime.now().replace(microsecond=0)}] Fetching PubMed abstracts for PMIDs: {pmids}")
    _log.info("Fetching PubMed abstracts for PMIDs: %s", pmids)
    api_key = os.environ.get("NCBI_API_KEY", "")
    key_param = f"&api_key={api_key}" if api_key else ""
    url = (
        f"{NCBI_EUTILS_BASE}/efetch.fcgi"
        f"?db=pubmed&id={','.join(pmids)}&rettype=xml&retmode=xml{key_param}"
    )
    try:
        xml_text = _http_get(url)
        root = ET.fromstring(xml_text)
    except Exception as exc:
        warnings.append(f"pubmed_fetch_failed:{exc}")
        return None

    abstract_el = root.find(".//AbstractText")
    if abstract_el is None:
        return None
    return "".join(abstract_el.itertext()).strip() or None


def _fetch_study_context(study_accession: str, warnings: list[str]) -> StudyContext | None:
    url = (
        f"{PORTAL_BASE}/filereport"
        f"?accession={study_accession}&result=study&fields=all&format=json"
    )
    try:
        raw = _http_get(url)
        records: list[dict[str, str]] = json.loads(raw)
    except Exception as exc:
        warnings.append(f"study_api_failed:{exc}")
        return None

    if not records:
        warnings.append(f"study_api_empty_response:{study_accession}")
        return None

    r = records[0]
    pubmed_ids = _parse_pubmed_ids(r.get("tag", ""))
    pubmed_abstract = _fetch_pubmed_abstract(pubmed_ids, warnings)

    return StudyContext(
        studyAccession=study_accession,
        studyTitle=_str(r.get("study_title")),
        studyDescription=_str(r.get("study_description")),
        geoAccession=_str(r.get("geo_accession")),
        pubmedIds=pubmed_ids,
        pubmedAbstract=pubmed_abstract,
    )


def fetch_experiment_context(accession: str) -> ExperimentContext:
    print(f"[{datetime.datetime.now().replace(microsecond=0)}] Fetching experiment context for accession: {accession}")
    _log.info("Fetching experiment context for accession: %s", accession)
    warnings: list[str] = []

    url = (
        f"{PORTAL_BASE}/filereport"
        f"?accession={accession}&result=read_experiment&fields=all&format=json"
    )
    try:
        raw = _http_get(url)
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

    study_accession = _str(first.get("study_accession"))

    with ThreadPoolExecutor(max_workers=2) as pool:
        xml_fut = (
            pool.submit(_http_get, f"{BROWSER_BASE}/xml/{sample_accession}")
            if sample_accession else None
        )
        study_fut = (
            pool.submit(_fetch_study_context, study_accession, warnings)
            if study_accession else None
        )

        if xml_fut is not None:
            try:
                xml_text = xml_fut.result()
                biological = biological.model_copy(
                    update={"sampleAttributes": _parse_sample_attributes(xml_text)}
                )
            except Exception as exc:
                warnings.append(f"sample_xml_failed:{exc}")

        study = study_fut.result() if study_fut is not None else None

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
