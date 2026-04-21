from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
import scanpy as sc

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.resolution import ResolutionSelection


@dataclass
class MergeInfo:
    conf: NDArray[np.float64]
    classes: NDArray
    labelMap: dict[str, str]
    mergedGroups: dict[str, list[str]]
    nClustersPreMerge: int
    nClustersPostMerge: int


def merge_clusters(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    sel: ResolutionSelection,
) -> tuple[sc.AnnData, MergeInfo]:
    X_hvg = adata.X[:, adata.var.highly_variable.values]
    if hasattr(X_hvg, "toarray"):
        X_hvg = X_hvg.toarray()

    weak_prior = adata.obs["cell_type"].values if cfg.rfBalanceWeakPrior else None
    conf, classes = _rf_pairwise_confusion(
        X_hvg,
        adata.obs[sel.clusterKey].values,
        weak_prior_labels=weak_prior,
    )

    merged_labels, label_map = _merge_by_confusion(
        adata.obs[sel.clusterKey].values, conf, classes, cfg.mergeThreshold
    )
    adata.obs["leiden_merged"] = pd.Categorical(merged_labels)

    merged_groups: dict[str, list[str]] = {}
    for original, merged in label_map.items():
        merged_groups.setdefault(str(merged), []).append(str(original))

    return adata, MergeInfo(
        conf=conf,
        classes=classes,
        labelMap={str(k): str(v) for k, v in label_map.items()},
        mergedGroups=merged_groups,
        nClustersPreMerge=len(classes),
        nClustersPostMerge=int(adata.obs["leiden_merged"].nunique()),
    )


def _rf_pairwise_confusion(
    X: NDArray[np.float32],
    cluster_labels: NDArray[np.str_],
    n_estimators: int = 300,
    n_splits: int = 3,
    random_state: int = 42,
    weak_prior_labels: NDArray | None = None,
) -> tuple[NDArray[np.float64], NDArray]:
    le = LabelEncoder()
    y = le.fit_transform(cluster_labels)
    n_classes = len(le.classes_)
    min_class_size = int(np.bincount(y).min())
    n_splits = min(n_splits, min_class_size)

    w_full: NDArray[np.float64] | None = None
    if weak_prior_labels is not None:
        s = pd.Series(np.asarray(weak_prior_labels))
        vc = s.value_counts(dropna=False)
        cnt = s.map(vc).replace(0, 1).astype(float)
        w_full = (1.0 / cnt).to_numpy(dtype=np.float64)
        w_full *= len(w_full) / w_full.sum()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    oof_preds = np.zeros(len(y), dtype=int)

    for train_idx, test_idx in skf.split(X, y):
        rf = RandomForestClassifier(
            n_estimators=n_estimators, n_jobs=-1, random_state=random_state
        )
        if w_full is None:
            rf.fit(X[train_idx], y[train_idx])
        else:
            rf.fit(X[train_idx], y[train_idx], sample_weight=w_full[train_idx])
        oof_preds[test_idx] = rf.predict(X[test_idx])

    conf = np.zeros((n_classes, n_classes))
    for true, pred in zip(y, oof_preds):
        conf[true, pred] += 1
    row_sums = conf.sum(axis=1, keepdims=True)
    conf = conf / np.where(row_sums == 0, 1, row_sums)

    return conf, le.classes_


def _merge_by_confusion(
    cluster_labels: NDArray[np.str_],
    conf: NDArray[np.float64],
    classes: NDArray,
    threshold: float,
) -> tuple[NDArray[np.str_], dict]:
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
    merged = np.vectorize(label_to_merged.get)(np.asarray(cluster_labels))
    return merged, label_to_merged
