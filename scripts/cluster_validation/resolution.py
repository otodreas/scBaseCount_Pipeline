from dataclasses import dataclass

import numpy as np
import scanpy as sc
from numpy.typing import NDArray

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.metrics import matched_jaccard


@dataclass
class ResolutionSelection:
    selectedResolution: float
    clusterKey: str
    jaccArr: NDArray[np.float64]
    kArr: NDArray[np.int64]
    bestIdx: int


def select_resolution(
    adata: sc.AnnData,
    cfg: ClusterValidationConfig,
    n_clusters: dict[float, int],
    k_filtered: int,
) -> tuple[sc.AnnData, ResolutionSelection]:
    del k_filtered  # retained for call-site compatibility
    resolutions = cfg.resolutions
    ref_labels = adata.obs[cfg.weakPriorKey].values
    k_arr = np.array([n_clusters[r] for r in resolutions], dtype=np.int64)
    jacc_arr = np.zeros(len(resolutions))

    for idx, r in enumerate(resolutions):
        jacc_arr[idx] = matched_jaccard(adata.obs[f"leiden_{r}"].values, ref_labels)

    best_idx = int(np.argmax(jacc_arr))
    selected = resolutions[best_idx]
    cluster_key = f"leiden_{selected}"

    return adata, ResolutionSelection(
        selectedResolution=selected,
        clusterKey=cluster_key,
        jaccArr=jacc_arr,
        kArr=k_arr,
        bestIdx=best_idx,
    )
