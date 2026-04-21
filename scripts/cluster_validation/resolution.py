from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linear_sum_assignment
import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


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
    resolutions = cfg.resolutions
    ref_labels = adata.obs["cell_type"].values
    celltypes = np.unique(ref_labels)
    k_arr = np.array([n_clusters[r] for r in resolutions], dtype=np.int64)
    jacc_arr = np.zeros(len(resolutions))

    for idx, r in enumerate(resolutions):
        leiden_labels = adata.obs[f"leiden_{r}"].values
        clusters = np.unique(leiden_labels)
        k, m = len(clusters), len(celltypes)

        cl_idx = {c: i for i, c in enumerate(clusters)}
        ct_idx = {t: j for j, t in enumerate(celltypes)}

        C = np.zeros((k, m), dtype=np.float64)
        for cl, ct in zip(leiden_labels, ref_labels):
            C[cl_idx[cl], ct_idx[ct]] += 1

        cl_sizes = C.sum(axis=1)
        ct_sizes = C.sum(axis=0)
        J = C / (cl_sizes[:, None] + ct_sizes[None, :] - C + 1e-10)

        row_ind, col_ind = linear_sum_assignment(-J)
        jacc_arr[idx] = J[row_ind, col_ind].sum()

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
