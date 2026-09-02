from dataclasses import dataclass

import numpy as np
import pandas as pd
import scanpy as sc
from numpy.typing import NDArray
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.resolution import ResolutionSelection

MERGED_CLUSTER_KEY = "leiden_merged"
RF_N_ESTIMATORS = 300
RF_N_SPLITS = 3
RF_RANDOM_STATE = 42


@dataclass
class MergeInfo:
    conf: NDArray[np.float64]
    classes: NDArray
    labelMap: dict[str, str]
    mergedGroups: dict[str, list[str]]
    nClustersPreMerge: int
    nClustersPostMerge: int
    mergedKey: str = MERGED_CLUSTER_KEY


def rf_pairwise_confusion(
    X: NDArray[np.floating],
    clusterLabels: NDArray | list,
    *,
    nEstimators: int = RF_N_ESTIMATORS,
    nSplits: int = RF_N_SPLITS,
    randomState: int = RF_RANDOM_STATE,
    weakPriorLabels: NDArray | None = None,
) -> tuple[NDArray[np.float64], NDArray]:
    """Out-of-fold Random Forest confusion among cluster labels.

    Rows are true clusters and columns are predicted clusters. Each row is
    normalized to sum to one (or left as zeros when a class has no cells).
    """
    X_arr = np.asarray(X)
    le = LabelEncoder()
    y = le.fit_transform(np.asarray(clusterLabels))
    n_classes = len(le.classes_)
    min_class_size = int(np.bincount(y).min())
    if min_class_size < 2:
        raise ValueError(
            f"Each cluster must have at least 2 cells for stratified folds; min class size is {min_class_size}"
        )
    n_splits = min(nSplits, min_class_size)

    w_full: NDArray[np.float64] | None = None
    if weakPriorLabels is not None:
        s = pd.Series(np.asarray(weakPriorLabels))
        vc = s.value_counts(dropna=False)
        cnt = s.map(vc).replace(0, 1).astype(float)
        w_full = (1.0 / cnt).to_numpy(dtype=np.float64)
        w_full *= len(w_full) / w_full.sum()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=randomState)
    oof_preds = np.zeros(len(y), dtype=int)

    for train_idx, test_idx in skf.split(X_arr, y):
        rf = RandomForestClassifier(n_estimators=nEstimators, n_jobs=-1, random_state=randomState)
        if w_full is None:
            rf.fit(X_arr[train_idx], y[train_idx])
        else:
            rf.fit(X_arr[train_idx], y[train_idx], sample_weight=w_full[train_idx])
        oof_preds[test_idx] = rf.predict(X_arr[test_idx])

    conf = np.zeros((n_classes, n_classes), dtype=np.float64)
    for true, pred in zip(y, oof_preds, strict=True):
        conf[true, pred] += 1
    row_sums = conf.sum(axis=1, keepdims=True)
    conf = conf / np.where(row_sums == 0, 1, row_sums)

    return conf, le.classes_


def merge_by_confusion(
    clusterLabels: NDArray | list,
    conf: NDArray[np.float64],
    classes: NDArray,
    threshold: float,
) -> tuple[NDArray, dict]:
    """Transitively merge clusters when either directional confusion exceeds ``threshold``."""
    n = len(classes)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        parent[find(x)] = find(y)

    for i in range(n):
        for j in range(i + 1, n):
            if conf[i, j] > threshold or conf[j, i] > threshold:
                union(i, j)

    root_map: dict[int, str] = {}
    counter = 0
    for idx in range(n):
        root = find(idx)
        if root not in root_map:
            root_map[root] = str(counter)
            counter += 1

    label_to_merged = {classes[i]: root_map[find(i)] for i in range(n)}
    merged = np.vectorize(label_to_merged.get)(np.asarray(clusterLabels))
    return merged, label_to_merged


def apply_rf_merge(
    adata: sc.AnnData,
    *,
    featureMatrix: NDArray[np.floating],
    clusterKey: str,
    mergeThreshold: float = 0.2,
    mergedKey: str = MERGED_CLUSTER_KEY,
    nEstimators: int = RF_N_ESTIMATORS,
    nSplits: int = RF_N_SPLITS,
    randomState: int = RF_RANDOM_STATE,
    weakPriorLabels: NDArray | None = None,
    minCellsPerCluster: int = 3,
) -> tuple[sc.AnnData, MergeInfo]:
    """Merge clusters with RF out-of-fold confusion on an explicit feature matrix."""
    if clusterKey not in adata.obs:
        raise ValueError(f"adata.obs is missing cluster key {clusterKey!r}")
    if featureMatrix.shape[0] != adata.n_obs:
        raise ValueError(f"featureMatrix rows ({featureMatrix.shape[0]}) must equal adata.n_obs ({adata.n_obs})")

    cluster_labels = adata.obs[clusterKey].values
    counts = pd.Series(cluster_labels).value_counts()
    undersized = counts[counts < minCellsPerCluster]
    if not undersized.empty:
        detail = ", ".join(f"{label}={int(n)}" for label, n in undersized.items())
        raise ValueError(
            f"Each cluster must have at least {minCellsPerCluster} cells for RF merge; undersized: {detail}"
        )

    conf, classes = rf_pairwise_confusion(
        featureMatrix,
        cluster_labels,
        nEstimators=nEstimators,
        nSplits=nSplits,
        randomState=randomState,
        weakPriorLabels=weakPriorLabels,
    )
    merged_labels, label_map = merge_by_confusion(cluster_labels, conf, classes, mergeThreshold)
    adata.obs[mergedKey] = pd.Categorical(merged_labels)

    merged_groups: dict[str, list[str]] = {}
    for original, merged in label_map.items():
        merged_groups.setdefault(str(merged), []).append(str(original))

    return adata, MergeInfo(
        conf=conf,
        classes=classes,
        labelMap={str(k): str(v) for k, v in label_map.items()},
        mergedGroups=merged_groups,
        nClustersPreMerge=len(classes),
        nClustersPostMerge=int(adata.obs[mergedKey].nunique()),
        mergedKey=mergedKey,
    )


def merge_clusters(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    sel: ResolutionSelection,
) -> tuple[sc.AnnData, MergeInfo]:
    X_hvg = adata.X[:, adata.var.highly_variable.values]
    if hasattr(X_hvg, "toarray"):
        X_hvg = X_hvg.toarray()

    weak_prior = adata.obs[cfg.weakPriorKey].values if cfg.rfBalanceWeakPrior else None
    return apply_rf_merge(
        adata,
        featureMatrix=np.asarray(X_hvg),
        clusterKey=sel.clusterKey,
        mergeThreshold=cfg.mergeThreshold,
        weakPriorLabels=weak_prior,
        # Historical per-dataset path allowed StratifiedKFold to shrink to min class size.
        minCellsPerCluster=2,
    )
