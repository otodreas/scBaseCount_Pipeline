from dataclasses import dataclass

import numpy as np
import scanpy as sc
from numpy.typing import NDArray

from cluster_validation.clustering import sweep_leiden_resolutions
from cluster_validation.config import ClusterValidationConfig
from cluster_validation.metrics import matched_jaccard


@dataclass
class ResolutionSelection:
    selectedResolution: float
    clusterKey: str
    jaccArr: NDArray[np.float64]
    kArr: NDArray[np.int64]
    bestIdx: int


def score_resolutions(
    adata: sc.AnnData,
    *,
    resolutions: list[float],
    weakPriorKey: str,
    nClusters: dict[float, int],
) -> ResolutionSelection:
    """Score precomputed ``leiden_{r}`` columns and pick the matched-Jaccard argmax."""
    if weakPriorKey not in adata.obs:
        raise ValueError(f"adata.obs is missing weak prior key {weakPriorKey!r}")
    if not resolutions:
        raise ValueError("resolutions must be non-empty")

    ref_labels = adata.obs[weakPriorKey].values
    k_arr = np.array([nClusters[r] for r in resolutions], dtype=np.int64)
    jacc_arr = np.zeros(len(resolutions), dtype=np.float64)

    for idx, r in enumerate(resolutions):
        key = f"leiden_{r}"
        if key not in adata.obs:
            raise ValueError(f"adata.obs is missing Leiden column {key!r}")
        jacc_arr[idx] = matched_jaccard(adata.obs[key].values, ref_labels)

    best_idx = int(np.argmax(jacc_arr))
    selected = resolutions[best_idx]
    return ResolutionSelection(
        selectedResolution=selected,
        clusterKey=f"leiden_{selected}",
        jaccArr=jacc_arr,
        kArr=k_arr,
        bestIdx=best_idx,
    )


def select_resolution_on_graph(
    adata: sc.AnnData,
    *,
    resolutions: list[float],
    weakPriorKey: str,
) -> tuple[sc.AnnData, ResolutionSelection]:
    """Sweep Leiden on the active neighbor graph, score, and return the argmax selection."""
    adata, n_clusters = sweep_leiden_resolutions(adata, resolutions)
    sel = score_resolutions(
        adata,
        resolutions=resolutions,
        weakPriorKey=weakPriorKey,
        nClusters=n_clusters,
    )
    return adata, sel


def select_resolution(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    n_clusters: dict[float, int],
    k_filtered: int,
) -> tuple[sc.AnnData, ResolutionSelection]:
    del k_filtered  # retained for call-site compatibility
    sel = score_resolutions(
        adata,
        resolutions=cfg.resolutions,
        weakPriorKey=cfg.weakPriorKey,
        nClusters=n_clusters,
    )
    return adata, sel
