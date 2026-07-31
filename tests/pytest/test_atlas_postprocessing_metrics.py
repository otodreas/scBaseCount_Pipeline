import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.metrics import (
    cross_study_macro_cell_type_neighbor_agreement,
    extract_plateaus,
)
from cluster_validation.metrics import matched_jaccard
from scipy import sparse


def test_extract_plateaus_relative_to_max() -> None:
    values = [5.0, 10.0, 15.0, 30.0, 50.0, 100.0]
    scores = [0.50, 0.92, 0.95, 0.94, 0.70, 0.60]
    plateaus = extract_plateaus(values, scores, relativeThreshold=0.95)
    assert plateaus == [(10.0, 30.0)]


def test_extract_plateaus_requires_aligned_inputs() -> None:
    try:
        extract_plateaus([1.0, 2.0], [0.5])
    except ValueError as exc:
        assert "same length" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_matched_jaccard_perfect_alignment() -> None:
    clusters = np.array(["a", "a", "b", "b"])
    labels = np.array(["x", "x", "y", "y"])
    assert matched_jaccard(clusters, labels) == pytest.approx(2.0)


def test_matched_jaccard_empty() -> None:
    assert matched_jaccard([], []) == 0.0


def _adata_with_distances() -> sc.AnnData:
    # 6 cells: two batches, two cell types. Each cell neighbors the other batch same-type cell.
    adata = sc.AnnData(np.zeros((6, 3), dtype=np.float32))
    adata.obs = pd.DataFrame(
        {
            "study_accession": ["S1", "S1", "S1", "S2", "S2", "S2"],
            "cell_type": ["A", "A", None, "A", "B", "B"],
        },
        index=[f"c{i}" for i in range(6)],
    )
    # distances graph (undirected pairs): 0-3 (A cross), 1-3 (A cross), 4-5 same batch ignored for eligibility alone
    rows = [0, 3, 1, 3, 4, 5]
    cols = [3, 0, 3, 1, 5, 4]
    data = [1.0] * len(rows)
    adata.obsp["distances"] = sparse.csr_matrix((data, (rows, cols)), shape=(6, 6))
    return adata


def test_cross_study_macro_agreement_and_coverage() -> None:
    adata = _adata_with_distances()
    score, coverage = cross_study_macro_cell_type_neighbor_agreement(
        adata,
        batchKey="study_accession",
        cellTypeKey="cell_type",
    )
    # eligible cells: 0,1,3 (type A). Type B cells 4,5 only have same-batch neighbors -> ineligible.
    assert coverage == 3 / 6
    assert score == 1.0
