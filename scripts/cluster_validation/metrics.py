from __future__ import annotations

from dataclasses import dataclass

import scanpy as sc
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
    silhouette_score,
    v_measure_score,
)

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.merge import MergeInfo
from cluster_validation.resolution import ResolutionSelection


@dataclass
class MetricArrays:
    silhouetteArr: list[list[float]]
    homogeneityArr: list[list[float]]
    completenessArr: list[list[float]]
    nmiArr: list[list[float]]
    vscoreArr: list[list[float]]
    ariArr: list[list[float]]


def compute_metrics(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    sel: ResolutionSelection,
    merge_info: MergeInfo,
) -> MetricArrays:
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
