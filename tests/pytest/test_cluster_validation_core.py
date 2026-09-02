import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from cluster_validation import (
    apply_rf_merge,
    default_resolutions,
    merge_by_confusion,
    pick_n_pcs,
    rf_pairwise_confusion,
    select_resolution_on_graph,
)
from cluster_validation.clustering import sweep_leiden_resolutions
from cluster_validation.resolution import score_resolutions


def test_default_resolutions_match_canonical_grid() -> None:
    assert default_resolutions() == [round(r, 1) for r in np.arange(0.1, 2.0, 0.1).tolist()]


def test_pick_n_pcs_returns_minimum_when_target_already_met() -> None:
    # First 15 PCs already explain >= 50%.
    ratios = np.full(50, 0.04, dtype=np.float64)
    n_pcs, cumvar = pick_n_pcs(ratios, nPcsMin=15, nPcsCompute=50, nPcsCumvarTarget=0.5)
    assert n_pcs == 15
    assert cumvar == pytest.approx(60.0)


def test_pick_n_pcs_grows_until_cumvar_target() -> None:
    ratios = np.full(50, 0.01, dtype=np.float64)
    ratios[:10] = 0.03
    n_pcs, cumvar = pick_n_pcs(ratios, nPcsMin=15, nPcsCompute=50, nPcsCumvarTarget=0.5)
    assert n_pcs > 15
    assert cumvar >= 50.0
    assert cumvar == pytest.approx(float(np.sum(ratios[:n_pcs])) * 100.0)


def test_pick_n_pcs_rejects_invalid_bounds() -> None:
    with pytest.raises(ValueError, match="nPcsMin"):
        pick_n_pcs([0.1] * 10, nPcsMin=0, nPcsCompute=10, nPcsCumvarTarget=0.5)
    with pytest.raises(ValueError, match="nPcsCompute"):
        pick_n_pcs([0.1] * 10, nPcsMin=15, nPcsCompute=10, nPcsCumvarTarget=0.5)


def _toy_graph_adata(n_cells: int = 60) -> sc.AnnData:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(n_cells, 8)).astype(np.float32)
    adata = sc.AnnData(x)
    labels = np.array(["A"] * (n_cells // 2) + ["B"] * (n_cells - n_cells // 2))
    adata.obs["cell_type"] = labels
    # Two well-separated blobs so Leiden finds stable partitions.
    x[: n_cells // 2] += 5.0
    adata.X = x
    sc.pp.neighbors(adata, n_neighbors=5, n_pcs=None, use_rep="X")
    return adata


def test_select_resolution_on_graph_argmax_and_tie_behavior() -> None:
    adata = _toy_graph_adata()
    resolutions = [0.2, 0.4, 0.6]
    adata, sel = select_resolution_on_graph(
        adata,
        resolutions=resolutions,
        weakPriorKey="cell_type",
    )
    assert sel.bestIdx == int(np.argmax(sel.jaccArr))
    assert sel.selectedResolution == resolutions[sel.bestIdx]
    assert sel.clusterKey == f"leiden_{sel.selectedResolution}"
    assert sel.kArr.shape == (3,)
    assert sel.jaccArr.shape == (3,)
    for r in resolutions:
        assert f"leiden_{r}" in adata.obs


def test_score_resolutions_matches_manual_sweep() -> None:
    adata = _toy_graph_adata()
    resolutions = [0.3, 0.5]
    adata, n_clusters = sweep_leiden_resolutions(adata, resolutions)
    sel = score_resolutions(
        adata,
        resolutions=resolutions,
        weakPriorKey="cell_type",
        nClusters=n_clusters,
    )
    assert list(sel.kArr) == [n_clusters[r] for r in resolutions]
    assert sel.bestIdx == int(np.argmax(sel.jaccArr))


def test_merge_by_confusion_directional_and_transitive() -> None:
    labels = np.array(["0", "0", "1", "1", "2", "2"])
    classes = np.array(["0", "1", "2"])
    # 0 <-> 1 above threshold; 1 <-> 2 above threshold => all merge transitively.
    conf = np.array(
        [
            [0.7, 0.3, 0.0],
            [0.25, 0.5, 0.25],
            [0.0, 0.3, 0.7],
        ],
        dtype=np.float64,
    )
    merged, label_map = merge_by_confusion(labels, conf, classes, threshold=0.2)
    assert len(set(label_map.values())) == 1
    assert len(set(merged.tolist())) == 1


def test_merge_by_confusion_requires_strictly_greater_than_threshold() -> None:
    labels = np.array(["0", "0", "1", "1"])
    classes = np.array(["0", "1"])
    conf = np.array(
        [
            [0.8, 0.2],
            [0.2, 0.8],
        ],
        dtype=np.float64,
    )
    merged, label_map = merge_by_confusion(labels, conf, classes, threshold=0.2)
    assert label_map == {"0": "0", "1": "1"}
    assert set(merged.tolist()) == {"0", "1"}


def test_rf_pairwise_confusion_is_deterministic_and_row_normalized() -> None:
    rng = np.random.default_rng(1)
    x = rng.normal(size=(40, 5)).astype(np.float32)
    labels = np.array(["0"] * 20 + ["1"] * 20)
    conf_a, classes_a = rf_pairwise_confusion(x, labels)
    conf_b, classes_b = rf_pairwise_confusion(x, labels)
    assert np.array_equal(classes_a, classes_b)
    assert np.allclose(conf_a, conf_b)
    assert np.allclose(conf_a.sum(axis=1), np.ones(len(classes_a)))


def test_apply_rf_merge_rejects_undersized_clusters() -> None:
    rng = np.random.default_rng(2)
    x = rng.normal(size=(10, 4)).astype(np.float32)
    adata = sc.AnnData(x)
    adata.obs["leiden_atlas"] = pd.Categorical(["0"] * 7 + ["1"] * 2 + ["2"])
    with pytest.raises(ValueError, match="at least 3 cells"):
        apply_rf_merge(
            adata,
            featureMatrix=x,
            clusterKey="leiden_atlas",
            minCellsPerCluster=3,
        )
