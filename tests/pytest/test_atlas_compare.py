"""Tests for bounded-memory atlas candidate comparison."""

import anndata as ad
import numpy as np
import pandas as pd
from h5ad_concat.compare import AtlasCompareConfig, compare_atlases

ad.settings.allow_write_nullable_strings = True


def test_compare_atlases_reports_baseline_subset(tmp_path) -> None:
    genes = ["GAPDH", "RPS18"]
    counts = np.array([[80, 20], [40, 60]], dtype=np.float32)
    obs = pd.DataFrame(
        {
            "SRX_accession": ["SRX1", "SRX1"],
            "cell_type": ["T cell", "T cell"],
            "pct_counts_ribo": [20.0, 60.0],
        },
        index=["cell0_SRX1", "cell1_SRX1"],
    )
    var = pd.DataFrame(index=genes)
    baseline = ad.AnnData(X=counts, obs=obs, var=var)
    candidate = baseline[:1].copy()
    baseline_path = tmp_path / "baseline.h5ad"
    candidate_path = tmp_path / "candidate.h5ad"
    baseline.write_h5ad(baseline_path)
    candidate.write_h5ad(candidate_path)

    report = compare_atlases(
        AtlasCompareConfig(
            baselinePath=baseline_path,
            candidatePath=candidate_path,
            reportPath=tmp_path / "report.json",
            chunkSize=1,
        )
    )

    assert report.byteIdentical is False
    assert report.fullLogicalIdentical is False
    assert report.nCandidateOnlyCells == 0
    assert report.nBaselineOnlyCells == 1
    assert report.commonMatrixEqual is True
    assert report.baselineOnlyDirectRiboRemovals == 1
    assert (tmp_path / "report.json").exists()
