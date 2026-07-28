import pytest
from disease_markers.labels import build_sample_label_table, coarse_disease_area
from shared.repo import REPO_ROOT


def test_coarse_disease_area_ipf() -> None:
    assert coarse_disease_area("idiopathic pulmonary fibrosis") == "IPF / Pulmonary Fibrosis"


def test_coarse_disease_area_control() -> None:
    assert coarse_disease_area("normal lung tissue") == "Control"


@pytest.mark.skipif(
    not (REPO_ROOT / "output/context/contexts.jsonl").is_file(),
    reason="contexts.jsonl not present",
)
def test_build_sample_label_table_non_empty() -> None:
    table = build_sample_label_table(
        REPO_ROOT / "output/context/contexts.jsonl",
        REPO_ROOT / "output/atlas/v1/atlas.csv",
    )
    assert len(table) > 0
    assert "eligible" in table.columns
