from unittest.mock import patch

from study_context.fetch import (
    _fetch_ena_sample_attributes,
    _fetch_ncbi_biosample_attributes,
    _fetch_sample_attributes,
    _is_ncbi_biosample,
    _parse_biosample_attributes,
    fetch_experiment_context,
)

BIOSAMPLE_XML = """\
<BioSampleSet>
  <BioSample accession="SAMN40627613">
    <Attributes>
      <Attribute attribute_name="tissue" harmonized_name="tissue">Lung parenchyma</Attribute>
      <Attribute attribute_name="disease" harmonized_name="disease">COPD</Attribute>
      <Attribute harmonized_name="cell_type">lung cells</Attribute>
    </Attributes>
  </BioSample>
</BioSampleSet>
"""

ENA_SAMPLE_XML = """\
<SAMPLE_SET>
  <SAMPLE accession="SAMEA114591064">
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE>
        <TAG>disease</TAG>
        <VALUE>lung adenocarcinoma</VALUE>
      </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</SAMPLE_SET>
"""

PORTAL_RECORD = {
    "run_accession": "SRR123",
    "sample_accession": "SAMN40627613",
    "study_accession": "PRJNA1092571",
    "scientific_name": "Homo sapiens",
}


def test_is_ncbi_biosample() -> None:
    assert _is_ncbi_biosample("SAMN40627613")
    assert _is_ncbi_biosample("samn40627613")
    assert not _is_ncbi_biosample("SAMEA114591064")


def test_parse_biosample_attributes() -> None:
    attrs = _parse_biosample_attributes(BIOSAMPLE_XML)
    assert attrs["tissue"] == "Lung parenchyma"
    assert attrs["disease"] == "COPD"
    assert attrs["cell_type"] == "lung cells"


def test_fetch_sample_attributes_routes_samn_to_ncbi() -> None:
    warnings: list[str] = []
    with patch("study_context.fetch._fetch_ncbi_biosample_attributes", return_value={"tissue": "lung"}) as ncbi:
        with patch("study_context.fetch._fetch_ena_sample_attributes") as ena:
            attrs = _fetch_sample_attributes("SAMN40627613", warnings)
    ncbi.assert_called_once_with("SAMN40627613", warnings)
    ena.assert_not_called()
    assert attrs == {"tissue": "lung"}


def test_fetch_sample_attributes_routes_samea_to_ena() -> None:
    warnings: list[str] = []
    with patch("study_context.fetch._fetch_ncbi_biosample_attributes") as ncbi:
        with patch("study_context.fetch._fetch_ena_sample_attributes", return_value={"disease": "COPD"}) as ena:
            attrs = _fetch_sample_attributes("SAMEA114591064", warnings)
    ena.assert_called_once_with("SAMEA114591064", warnings)
    ncbi.assert_not_called()
    assert attrs == {"disease": "COPD"}


def test_fetch_ncbi_biosample_attributes_success() -> None:
    warnings: list[str] = []
    with patch("study_context.fetch._http_get_ncbi", return_value=BIOSAMPLE_XML):
        attrs = _fetch_ncbi_biosample_attributes("SAMN40627613", warnings)
    assert attrs["disease"] == "COPD"
    assert warnings == []


def test_fetch_ncbi_biosample_attributes_failure() -> None:
    warnings: list[str] = []
    with patch("study_context.fetch._http_get_ncbi", side_effect=RuntimeError("network down")):
        attrs = _fetch_ncbi_biosample_attributes("SAMN40627613", warnings)
    assert attrs == {}
    assert len(warnings) == 1
    assert warnings[0].startswith("biosample_fetch_failed:")


def test_fetch_ena_sample_attributes_success() -> None:
    warnings: list[str] = []
    with patch("study_context.fetch._http_get", return_value=ENA_SAMPLE_XML):
        attrs = _fetch_ena_sample_attributes("SAMEA114591064", warnings)
    assert attrs["disease"] == "lung adenocarcinoma"
    assert warnings == []


def test_fetch_experiment_context_uses_ncbi_for_samn() -> None:
    def fake_http_get_ncbi(url: str, *, retries: int = 3) -> str:
        if "db=biosample" in url:
            return BIOSAMPLE_XML
        raise AssertionError(f"unexpected URL for SAMN sample: {url}")

    with patch("study_context.fetch.fetch_read_experiment_records", return_value=[PORTAL_RECORD]):
        with patch("study_context.fetch._fetch_study_context", return_value=None):
            with patch("study_context.fetch._http_get_ncbi", side_effect=fake_http_get_ncbi):
                ctx = fetch_experiment_context("SRX24073315")

    assert ctx.biological.sampleAttributes["disease"] == "COPD"
    assert not any(w.startswith("sample_xml_failed") for w in ctx.warnings)


def test_fetch_experiment_context_uses_ena_for_samea() -> None:
    record = {**PORTAL_RECORD, "sample_accession": "SAMEA114591064"}

    def fake_http_get(url: str, *, retries: int = 3) -> str:
        if "/xml/SAMEA114591064" in url:
            return ENA_SAMPLE_XML
        raise AssertionError(f"unexpected URL for SAMEA sample: {url}")

    with patch("study_context.fetch.fetch_read_experiment_records", return_value=[record]):
        with patch("study_context.fetch._fetch_study_context", return_value=None):
            with patch("study_context.fetch._http_get", side_effect=fake_http_get):
                ctx = fetch_experiment_context("ERX11662338")

    assert ctx.biological.sampleAttributes["disease"] == "lung adenocarcinoma"
    assert not any("biosample" in w for w in ctx.warnings)


def test_ncbi_rate_limiter_spaces_requests(monkeypatch) -> None:
    import study_context.fetch as fetch

    fetch._ncbi_last_request_at = 0.0
    clock = [0.0]

    monkeypatch.setattr(fetch.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(fetch.time, "sleep", lambda seconds: clock.__setitem__(0, clock[0] + seconds))
    monkeypatch.delenv("NCBI_API_KEY", raising=False)

    with patch("study_context.fetch._http_get", return_value="ok") as http_get:
        fetch._http_get_ncbi("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=SAMN1")
        fetch._http_get_ncbi("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=SAMN2")

    assert http_get.call_count == 2
    assert clock[0] >= fetch._ncbi_min_interval()
