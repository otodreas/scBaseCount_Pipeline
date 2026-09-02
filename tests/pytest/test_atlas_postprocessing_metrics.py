import numpy as np
import pytest
from cluster_validation.metrics import matched_jaccard


def test_matched_jaccard_perfect_alignment() -> None:
    clusters = np.array(["a", "a", "b", "b"])
    labels = np.array(["x", "x", "y", "y"])
    assert matched_jaccard(clusters, labels) == pytest.approx(2.0)


def test_matched_jaccard_empty() -> None:
    assert matched_jaccard([], []) == 0.0
