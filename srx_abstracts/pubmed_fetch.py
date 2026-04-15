from __future__ import annotations

import xml.etree.ElementTree as ET

from srx_abstracts.entrez_client import EntrezClient
from srx_abstracts.models import PubmedArticle


def _stripNs(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", maxsplit=1)[-1]
    return tag


def _textJoin(parts: list[str]) -> str:
    return "\n".join(p for p in parts if p).strip()


def _findFirstDesc(parent: ET.Element, name: str) -> ET.Element | None:
    for el in parent.iter():
        if el is parent:
            continue
        if _stripNs(el.tag) == name:
            return el
    return None


def _parsePubmedArticle(node: ET.Element) -> PubmedArticle | None:
    mc = None
    for ch in node:
        if _stripNs(ch.tag) == "MedlineCitation":
            mc = ch
            break
    if mc is None:
        return None

    pmid: str | None = None
    for ch in mc:
        if _stripNs(ch.tag) == "PMID" and ch.text:
            pmid = ch.text.strip()
            break
    if not pmid:
        return None

    artEl = None
    for ch in mc:
        if _stripNs(ch.tag) == "Article":
            artEl = ch
            break

    title = ""
    abstractPieces: list[str] = []
    journalTitle: str | None = None
    year: int | None = None

    if artEl is not None:
        at = _findFirstDesc(artEl, "ArticleTitle")
        if at is not None and at.text:
            title = at.text.strip()

        for ch in artEl:
            if _stripNs(ch.tag) != "Abstract":
                continue
            for ab in ch:
                if _stripNs(ab.tag) != "AbstractText":
                    continue
                label = ab.attrib.get("Label")
                chunk = (ab.text or "").strip()
                if label and chunk:
                    abstractPieces.append(f"{label}: {chunk}")
                elif chunk:
                    abstractPieces.append(chunk)

        for ch in artEl:
            if _stripNs(ch.tag) != "Journal":
                continue
            jt = _findFirstDesc(ch, "Title")
            if jt is not None and jt.text:
                journalTitle = jt.text.strip()
            ji = None
            for jch in ch:
                if _stripNs(jch.tag) == "JournalIssue":
                    ji = jch
                    break
            if ji is not None:
                pd = None
                for jich in ji:
                    if _stripNs(jich.tag) == "PubDate":
                        pd = jich
                        break
                if pd is not None:
                    for z in pd:
                        if _stripNs(z.tag) == "Year" and z.text and z.text.strip().isdigit():
                            year = int(z.text.strip())
                            break

    abstractText = _textJoin(abstractPieces) or None
    return PubmedArticle(
        pmid=pmid,
        title=title,
        abstractText=abstractText,
        journal=journalTitle,
        year=year,
    )


def _parsePubmedBatch(xmlText: str) -> list[PubmedArticle]:
    root = ET.fromstring(xmlText)
    for err in root.iter():
        if _stripNs(err.tag) == "ERROR" and (err.text or "").strip():
            raise ValueError(f"efetch ERROR: {err.text.strip()}")

    articles: list[PubmedArticle] = []
    for el in root:
        if _stripNs(el.tag) != "PubmedArticle":
            continue
        parsed = _parsePubmedArticle(el)
        if parsed is not None:
            articles.append(parsed)
    return articles


def fetch_pubmed_records(
    client: EntrezClient,
    pmids: list[str],
    *,
    chunk_size: int = 100,
) -> list[PubmedArticle]:
    out: list[PubmedArticle] = []
    clean = [p.strip() for p in pmids if p.strip()]
    for i in range(0, len(clean), chunk_size):
        chunk = clean[i : i + chunk_size]
        id_param = ",".join(chunk)
        xmlText = client.efetch(db="pubmed", id_param=id_param, rettype="xml")
        out.extend(_parsePubmedBatch(xmlText))
    byPmid = {a.pmid: a for a in out}
    return [byPmid[p] for p in clean if p in byPmid]
