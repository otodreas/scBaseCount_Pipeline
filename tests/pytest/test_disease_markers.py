import anndata as ad
import numpy as np
import pandas as pd
import pytest
from disease_markers.labels import build_sample_label_table, coarse_disease_area
from disease_markers.transfer import transfer_leiden_clusters
from shared.repo import REPO_ROOT


def test_coarse_disease_area_ipf() -> None:
    assert coarse_disease_area("idiopathic pulmonary fibrosis") == "IPF / Pulmonary Fibrosis"


def test_coarse_disease_area_control() -> None:
    assert coarse_disease_area("normal lung tissue") == "Control"


def test_transfer_leiden_clusters_by_obs_names() -> None:
    n = 5
    obs_names = [f"cell_{i}" for i in range(n)]
    full = ad.AnnData(X=np.ones((n, 3)), obs=pd.DataFrame(index=obs_names))
    harmony = ad.AnnData(
        X=np.ones((n, 2)),
        obs=pd.DataFrame({"leiden_atlas": ["0", "0", "1", "1", "1"]}, index=obs_names),
    )
    out = transfer_leiden_clusters(full, harmony, clusterKey="leiden_atlas")
    assert list(out.obs["leiden_atlas"].astype(str)) == ["0", "0", "1", "1", "1"]


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
