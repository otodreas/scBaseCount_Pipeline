from __future__ import annotations

import xml.etree.ElementTree as ET

from srx_abstracts.entrez_client import EntrezClient
from srx_abstracts.models import PmidProvenance, SrxResolution


def _stripNs(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", maxsplit=1)[-1]
    return tag


def _parseEsearchIdsRobust(xmlText: str) -> list[str]:
    root = ET.fromstring(xmlText)
    for err in root.iter():
        if _stripNs(err.tag) == "ERROR" and (err.text or "").strip():
            raise ValueError(f"esearch ERROR: {err.text.strip()}")
    out: list[str] = []
    for el in root.iter():
        if _stripNs(el.tag) != "IdList":
            continue
        for child in el:
            if _stripNs(child.tag) == "Id" and child.text:
                out.append(child.text.strip())
    return out


def _linkSetDbTargets(xmlText: str, *, dbTo: str) -> set[str]:
    root = ET.fromstring(xmlText)
    for err in root.iter():
        if _stripNs(err.tag) == "ERROR" and (err.text or "").strip():
            raise ValueError(f"elink ERROR: {err.text.strip()}")
    want = dbTo.strip().lower()
    out: set[str] = set()
    for lsd in root.iter():
        if _stripNs(lsd.tag) != "LinkSetDb":
            continue
        dbToEl: str | None = None
        for child in lsd:
            t = _stripNs(child.tag)
            if t == "DbTo" and child.text:
                dbToEl = child.text.strip().lower()
                break
        if dbToEl != want:
            continue
        for el in lsd.iter():
            if el is lsd:
                continue
            if _stripNs(el.tag) == "Id" and el.text:
                out.add(el.text.strip())
    return out


def _parseElinkPubmedIds(xmlText: str) -> set[str]:
    return _linkSetDbTargets(xmlText, dbTo="pubmed")


def _parseElinkLinkedIds(xmlText: str, *, dbTo: str) -> set[str]:
    return _linkSetDbTargets(xmlText, dbTo=dbTo)


def _dedupePreserveOrder(pmids: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for p in pmids:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def resolve_pmids(client: EntrezClient, srx_accession: str) -> SrxResolution:
    acc = srx_accession.strip()
    if not acc:
        return SrxResolution(srxAccession=srx_accession, warnings=["empty_accession"])

    warnings: list[str] = []
    provenance: list[PmidProvenance] = []

    term = f"{acc}[accn]"
    try:
        searchXml = client.esearch(db="sra", term=term, retmax="50")
        sraUids = _parseEsearchIdsRobust(searchXml)
    except (ET.ParseError, ValueError) as e:
        return SrxResolution(
            srxAccession=acc,
            warnings=[f"esearch_failed:{e}"],
        )

    if not sraUids:
        warnings.append("no_sra_uid_for_accession")
        return SrxResolution(
            srxAccession=acc,
            sraUids=[],
            pmidProvenance=[],
            warnings=warnings,
        )

    idParam = ",".join(sraUids)

    try:
        xml = client.elink(dbfrom="sra", db="pubmed", id_param=idParam)
        for pm in _parseElinkPubmedIds(xml):
            provenance.append(PmidProvenance(pmid=pm, sourceHop="sra_pubmed"))
    except (ET.ParseError, ValueError) as e:
        warnings.append(f"sra_pubmed_elink_failed:{e}")

    bioprojectIds: set[str] = set()
    try:
        xmlBp = client.elink(dbfrom="sra", db="bioproject", id_param=idParam)
        bioprojectIds = _parseElinkLinkedIds(xmlBp, dbTo="bioproject")
    except (ET.ParseError, ValueError) as e:
        warnings.append(f"sra_bioproject_elink_failed:{e}")

    if bioprojectIds:
        bpParam = ",".join(sorted(bioprojectIds))
        try:
            xmlPm = client.elink(dbfrom="bioproject", db="pubmed", id_param=bpParam)
            for pm in _parseElinkPubmedIds(xmlPm):
                provenance.append(PmidProvenance(pmid=pm, sourceHop="bioproject_pubmed"))
        except (ET.ParseError, ValueError) as e:
            warnings.append(f"bioproject_pubmed_elink_failed:{e}")

    gdsIds: set[str] = set()
    try:
        xmlGds = client.elink(dbfrom="sra", db="gds", id_param=idParam)
        gdsIds = _parseElinkLinkedIds(xmlGds, dbTo="gds")
    except (ET.ParseError, ValueError) as e:
        warnings.append(f"sra_gds_elink_failed:{e}")

    if gdsIds:
        gdsParam = ",".join(sorted(gdsIds))
        try:
            xmlGdsPm = client.elink(dbfrom="gds", db="pubmed", id_param=gdsParam)
            for pm in _parseElinkPubmedIds(xmlGdsPm):
                provenance.append(PmidProvenance(pmid=pm, sourceHop="gds_pubmed"))
        except (ET.ParseError, ValueError) as e:
            warnings.append(f"gds_pubmed_elink_failed:{e}")

    if not provenance:
        warnings.append("no_pmids_from_elink")

    return SrxResolution(
        srxAccession=acc,
        sraUids=sraUids,
        pmidProvenance=provenance,
        warnings=warnings,
    )


def pmids_from_resolution(resolution: SrxResolution) -> list[str]:
    return _dedupePreserveOrder([p.pmid for p in resolution.pmidProvenance])
