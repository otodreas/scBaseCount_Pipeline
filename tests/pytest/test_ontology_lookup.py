import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from ontology_lookup import (
    OntologyCache,
    OntologyLookupConfig,
    OntologyReleaseMismatchError,
    area_for_record,
    assert_release,
    disease_area_from_records,
    ontology_tokens,
    resolve_ontology_field,
    respiratory_from_uberon,
)
from ontology_lookup.cache import TermRecord
from ontology_lookup.client import fetch_term


def test_ontology_tokens_single_and_multi() -> None:
    assert ontology_tokens("MONDO:0005061") == ["MONDO:0005061"]
    assert ontology_tokens("MONDO:0008903,MONDO:0005002") == ["MONDO:0008903", "MONDO:0005002"]
    assert ontology_tokens("mondo:0005061, MONDO:0005061") == ["MONDO:0005061"]
    assert ontology_tokens("") == []
    assert ontology_tokens("not-a-curie") == []


def test_area_for_record_lung_cancer_and_generic() -> None:
    luad = TermRecord(
        curie="MONDO:0005061",
        label="lung adenocarcinoma",
        ancestors=("MONDO:0008903", "MONDO:0005233", "MONDO:0004992"),
    )
    generic = TermRecord(
        curie="MONDO:0004992",
        label="cancer",
        ancestors=("MONDO:0000001",),
    )
    assert area_for_record(luad) == "Lung Cancer"
    assert area_for_record(generic) == "Other"


def test_disease_area_conflict_and_agreement() -> None:
    covid = TermRecord(curie="MONDO:0100096", label="COVID-19", ancestors=())
    copd = TermRecord(curie="MONDO:0005002", label="COPD", ancestors=())
    assert disease_area_from_records([covid, covid]) == ("COVID-19 / SARS-CoV-2", "ontology")
    assert disease_area_from_records([covid, copd])[0] == "Other"
    assert disease_area_from_records([covid, copd])[1] == "conflicting_ontology"


def test_respiratory_from_uberon() -> None:
    lung = TermRecord(curie="UBERON:0002048", label="lung", ancestors=("UBERON:0001004",))
    blood = TermRecord(curie="UBERON:0000178", label="blood", ancestors=("UBERON:0000465",))
    assert respiratory_from_uberon([lung]) is True
    assert respiratory_from_uberon([lung, blood]) is False
    assert respiratory_from_uberon([]) is None


def _write_cache(tmp_path: Path, ontology_id: str, release: str, terms: dict[str, dict]) -> OntologyCache:
    root = tmp_path / ontology_id / release
    root.mkdir(parents=True)
    text = json.dumps(terms, indent=2, sort_keys=True) + "\n"
    (root / "terms.json").write_text(text)
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "ontologyId": ontology_id,
                "release": release,
                "sourceUrl": "https://example.test",
                "generatedAt": "2026-01-01T00:00:00+00:00",
                "termCount": len(terms),
                "contentSha256": "abc",
            },
            indent=2,
        )
        + "\n"
    )
    return OntologyCache(
        ontologyId=ontology_id,
        release=release,
        cacheDir=tmp_path,
        baseUrl="https://example.test",
        timeoutS=1.0,
    )


def test_resolve_ontology_field_statuses(tmp_path: Path) -> None:
    cache = _write_cache(
        tmp_path,
        "mondo",
        "2026-07-06",
        {
            "MONDO:0005061": {"label": "lung adenocarcinoma", "ancestors": ["MONDO:0008903"]},
        },
    )
    resolved = resolve_ontology_field("MONDO:0005061", cache)
    assert resolved.status == "resolved"
    assert resolved.names == ("lung adenocarcinoma",)
    partial = resolve_ontology_field("MONDO:0005061,MONDO:9999999", cache)
    assert partial.status == "partial"
    missing = resolve_ontology_field("MONDO:9999999", cache)
    assert missing.status == "missing"
    empty = resolve_ontology_field("", cache)
    assert empty.status == "empty"


def test_assert_release_mismatch() -> None:
    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"http://www.w3.org/2002/07/owl#versionInfo": "1999-01-01"}
    with patch("ontology_lookup.client.httpx.get", return_value=mock_resp):
        with pytest.raises(OntologyReleaseMismatchError):
            assert_release("mondo", "2026-07-06", baseUrl="https://example.test", timeoutS=1.0)


def test_cache_reuse_without_network(tmp_path: Path) -> None:
    cache = _write_cache(
        tmp_path,
        "mondo",
        "2026-07-06",
        {"MONDO:0005061": {"label": "lung adenocarcinoma", "ancestors": ["MONDO:0008903"]}},
    )
    with patch("ontology_lookup.cache.fetch_term") as fetch:
        out = cache.ensure(["MONDO:0005061"], allowNetwork=False)
        fetch.assert_not_called()
    assert out["MONDO:0005061"] is not None
    assert out["MONDO:0005061"].label == "lung adenocarcinoma"


def test_fetch_term_parses_label_and_ancestors() -> None:
    class_resp = MagicMock()
    class_resp.status_code = 200
    class_resp.raise_for_status = MagicMock()
    class_resp.json.return_value = {"label": ["lung adenocarcinoma"], "curie": "MONDO:0005061"}
    anc_resp = MagicMock()
    anc_resp.raise_for_status = MagicMock()
    anc_resp.json.return_value = {
        "elements": [{"curie": "MONDO:0008903"}, {"curie": "MONDO:0005061"}],
        "totalPages": 1,
    }

    def _get(url: str, **kwargs):
        del kwargs
        if url.endswith("/hierarchicalAncestors") or "hierarchicalAncestors" in url:
            return anc_resp
        return class_resp

    with patch("ontology_lookup.client.httpx.get", side_effect=_get):
        row = fetch_term("mondo", "MONDO:0005061", baseUrl="https://example.test", timeoutS=1.0)
    assert row["label"] == "lung adenocarcinoma"
    assert row["ancestors"] == ["MONDO:0008903"]


def test_config_defaults() -> None:
    cfg = OntologyLookupConfig()
    assert cfg.mondoRelease == "2026-07-06"
    assert cfg.uberonRelease == "2026-06-19"
