import numpy as np
import scanpy as sc
from numpy.typing import NDArray
from scipy import sparse


def extract_plateaus(
    values: list[float] | NDArray[np.floating],
    scores: list[float] | NDArray[np.floating],
    *,
    relativeThreshold: float = 0.95,
) -> list[tuple[float, float]]:
    """Return contiguous value ranges whose scores stay within relativeThreshold of the max."""
    if len(values) == 0 or len(scores) == 0:
        return []
    if len(values) != len(scores):
        raise ValueError("values and scores must have the same length")
    if not (0.0 < relativeThreshold <= 1.0):
        raise ValueError("relativeThreshold must be in (0, 1]")

    value_arr = np.asarray(values, dtype=np.float64)
    score_arr = np.asarray(scores, dtype=np.float64)
    max_score = float(np.max(score_arr))
    if not np.isfinite(max_score):
        return []

    threshold = max_score * relativeThreshold
    in_plateau = score_arr >= threshold

    plateaus: list[tuple[float, float]] = []
    start: int | None = None
    for idx, flag in enumerate(in_plateau):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            plateaus.append((float(value_arr[start]), float(value_arr[idx - 1])))
            start = None
    if start is not None:
        plateaus.append((float(value_arr[start]), float(value_arr[-1])))
    return plateaus


def _neighbor_indices(distances: sparse.spmatrix) -> list[NDArray[np.int64]]:
    """Return neighbor index arrays per cell from a distances sparse matrix."""
    mat = distances.tocsr()
    out: list[NDArray[np.int64]] = []
    for row in range(mat.shape[0]):
        start, end = mat.indptr[row], mat.indptr[row + 1]
        out.append(np.asarray(mat.indices[start:end], dtype=np.int64))
    return out


def cross_study_macro_cell_type_neighbor_agreement(
    adata: sc.AnnData,
    *,
    batchKey: str,
    cellTypeKey: str,
) -> tuple[float, float]:
    """Macro-average cross-study same-label neighbor agreement and eligible-cell coverage.

    For each cell with a non-null cell type, consider neighbors from other batches that also
    have a non-null cell type. The cell score is the fraction of those neighbors sharing its
    label. Scores are averaged within each cell type, then macro-averaged across cell types.
    Coverage is the fraction of cells that contribute at least one such neighbor.
    """
    if "distances" not in adata.obsp:
        raise ValueError("adata.obsp['distances'] is required; run sc.pp.neighbors first")
    if batchKey not in adata.obs:
        raise ValueError(f"adata.obs is missing batch key {batchKey!r}")
    if cellTypeKey not in adata.obs:
        raise ValueError(f"adata.obs is missing cell type key {cellTypeKey!r}")

    batches = adata.obs[batchKey].astype(str).to_numpy()
    labels = adata.obs[cellTypeKey]
    label_vals = labels.astype(object).to_numpy()
    label_is_null = labels.isna().to_numpy() if hasattr(labels, "isna") else np.array([x is None for x in label_vals])

    neighbor_lists = _neighbor_indices(adata.obsp["distances"])
    per_type_scores: dict[str, list[float]] = {}
    eligible = 0

    for cell_idx, neighbors in enumerate(neighbor_lists):
        if label_is_null[cell_idx]:
            continue
        cell_label = label_vals[cell_idx]
        cell_batch = batches[cell_idx]

        same = 0
        cross = 0
        for neighbor_idx in neighbors:
            if label_is_null[neighbor_idx]:
                continue
            if batches[neighbor_idx] == cell_batch:
                continue
            cross += 1
            if label_vals[neighbor_idx] == cell_label:
                same += 1

        if cross == 0:
            continue

        eligible += 1
        key = str(cell_label)
        per_type_scores.setdefault(key, []).append(same / cross)

    coverage = eligible / adata.n_obs if adata.n_obs else 0.0
    if not per_type_scores:
        return 0.0, coverage

    type_means = [float(np.mean(scores)) for scores in per_type_scores.values()]
    return float(np.mean(type_means)), float(coverage)
