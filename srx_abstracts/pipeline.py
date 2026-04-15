from __future__ import annotations

from srx_abstracts.entrez_client import EntrezClient
from srx_abstracts.models import SrxAbstractBundle
from srx_abstracts.pubmed_fetch import fetch_pubmed_records
from srx_abstracts.resolve import pmids_from_resolution, resolve_pmids


def pipeline_for_srx(client: EntrezClient, srx_accession: str) -> SrxAbstractBundle:
    res = resolve_pmids(client, srx_accession)
    pmids = pmids_from_resolution(res)
    articles = fetch_pubmed_records(client, pmids) if pmids else []
    fetched = {a.pmid for a in articles}
    warnings = list(res.warnings)
    missing = [p for p in pmids if p not in fetched]
    if missing:
        warnings.append(f"efetch_missing_pmids:{','.join(missing)}")
    return SrxAbstractBundle(
        srxAccession=res.srxAccession,
        sraUids=res.sraUids,
        pmids=pmids,
        pmidProvenance=res.pmidProvenance,
        articles=articles,
        warnings=warnings,
    )


def pipeline_for_srx_list(client: EntrezClient, srx_accessions: list[str]) -> list[SrxAbstractBundle]:
    return [pipeline_for_srx(client, s) for s in srx_accessions]
