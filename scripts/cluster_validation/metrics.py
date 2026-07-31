from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
    silhouette_score,
    v_measure_score,
)

from cluster_validation.config import ClusterValidationConfig


def matched_jaccard(
    clusterLabels: NDArray[Any] | list[Any],
    refLabels: NDArray[Any] | list[Any],
) -> float:
    """Hungarian-matched sum of pairwise Jaccard scores between clusters and reference labels."""
    cluster_labels = np.asarray(clusterLabels)
    ref_labels = np.asarray(refLabels)
    if cluster_labels.size == 0 or ref_labels.size == 0:
        return 0.0

    clusters = np.unique(cluster_labels)
    celltypes = np.unique(ref_labels)
    k, m = len(clusters), len(celltypes)
    if k == 0 or m == 0:
        return 0.0

    cl_idx = {c: i for i, c in enumerate(clusters)}
    ct_idx = {t: j for j, t in enumerate(celltypes)}

    contingency = np.zeros((k, m), dtype=np.float64)
    for cl, ct in zip(cluster_labels, ref_labels, strict=True):
        contingency[cl_idx[cl], ct_idx[ct]] += 1

    cl_sizes = contingency.sum(axis=1)
    ct_sizes = contingency.sum(axis=0)
    jaccard = contingency / (cl_sizes[:, None] + ct_sizes[None, :] - contingency + 1e-10)

    row_ind, col_ind = linear_sum_assignment(-jaccard)
    return float(jaccard[row_ind, col_ind].sum())


@dataclass
class MetricArrays:
    silhouetteArr: list[list[float]]
    homogeneityArr: list[list[float]]
    completenessArr: list[list[float]]
    nmiArr: list[list[float]]
    vscoreArr: list[list[float]]
    ariArr: list[list[float]]


def compute_metrics(
    adata: Any,
    cfg: ClusterValidationConfig,
    sel: Any,
    merge_info: Any,
) -> MetricArrays:
    del sel, merge_info
    silhouette: list[list[float]] = []
    homogeneity: list[list[float]] = []
    completeness: list[list[float]] = []
    nmi: list[list[float]] = []
    vscore: list[list[float]] = []
    ari: list[list[float]] = []

    for res in cfg.resolutions:
        key = f"leiden_{res}"
        labels = adata.obs[key]
        if labels.nunique() <= 1:
            continue
        sil = float(silhouette_score(adata.obsm["X_pca"], labels))
        silhouette.append([res, sil])
        homogeneity.append([res, float(homogeneity_score(adata.obs[key], adata.obs["leiden_merged"]))])
        completeness.append([res, float(completeness_score(adata.obs[key], adata.obs["leiden_merged"]))])
        nmi.append([res, float(normalized_mutual_info_score(adata.obs[key], adata.obs["leiden_merged"]))])
        vscore.append([res, float(v_measure_score(adata.obs[key], adata.obs["leiden_merged"]))])
        ari.append([res, float(adjusted_rand_score(adata.obs[key], adata.obs["leiden_merged"]))])

    return MetricArrays(
        silhouetteArr=silhouette,
        homogeneityArr=homogeneity,
        completenessArr=completeness,
        nmiArr=nmi,
        vscoreArr=vscore,
        ariArr=ari,
    )
