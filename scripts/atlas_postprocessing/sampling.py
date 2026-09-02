import logging
from collections.abc import Mapping
from typing import Any

import numpy as np
import scanpy as sc

log = logging.getLogger(__name__)

# SAMPLE_SEED = 0
SAMPLE_METHOD = "studyProportional"
SAMPLE_UNS_KEY = "parameterSelectionSample"


def _allocate_study_counts(sizes: Mapping[str, int], n: int) -> dict[str, int]:
    """Allocate exactly ``n`` cells across studies with one guaranteed slot each.

    Remaining slots use largest-remainder proportional allocation by study size,
    capped by available cells per study.
    """
    studies = sorted(sizes)
    n_studies = len(studies)
    total = sum(sizes[study] for study in studies)
    if n < 1:
        raise ValueError(f"sampleCells must be >= 1, got {n}")
    if n_studies < 1:
        raise ValueError("Cannot sample: no studies present in stratify key")
    if n < n_studies:
        raise ValueError(f"sampleCells ({n}) must be >= number of studies ({n_studies}) to include every study")
    if n > total:
        raise ValueError(f"sampleCells ({n}) exceeds available cells ({total})")

    remaining = n - n_studies
    capacity = {study: sizes[study] - 1 for study in studies}
    if remaining == 0:
        return {study: 1 for study in studies}

    raw = {study: remaining * sizes[study] / total for study in studies}
    floors = {study: int(raw[study]) for study in studies}
    leftover = remaining - sum(floors.values())
    frac_order = sorted(studies, key=lambda study: (-(raw[study] - floors[study]), study))
    for study in frac_order:
        if leftover == 0:
            break
        floors[study] += 1
        leftover -= 1

    overflow = 0
    for study in studies:
        if floors[study] > capacity[study]:
            overflow += floors[study] - capacity[study]
            floors[study] = capacity[study]

    while overflow > 0:
        receivers = sorted(
            [study for study in studies if floors[study] < capacity[study]],
            key=lambda study: (-(capacity[study] - floors[study]), study),
        )
        if not receivers:
            raise ValueError("Unable to allocate sampleCells under study size caps")
        for study in receivers:
            if overflow == 0:
                break
            floors[study] += 1
            overflow -= 1

    return {study: 1 + floors[study] for study in studies}


def sample_metadata(adata: sc.AnnData) -> dict[str, Any] | None:
    """Return sampling audit fields stored on ``adata.uns``, if present."""
    payload = adata.uns.get(SAMPLE_UNS_KEY)
    if payload is None:
        return None
    return dict(payload)


def _study_positions(adata: sc.AnnData, stratifyKey: str) -> dict[str, np.ndarray]:
    labels = [str(value) for value in adata.obs[stratifyKey].tolist()]
    buckets: dict[str, list[int]] = {}
    for index, label in enumerate(labels):
        buckets.setdefault(label, []).append(index)
    return {label: np.asarray(positions, dtype=np.int64) for label, positions in buckets.items()}


def sample_study_proportional(
    adata: sc.AnnData,
    *,
    n: int,
    stratifyKey: str,
    seed: int,
) -> sc.AnnData:  # mentioned in methods
    """Return an exact-size study-proportional sample without replacement.

    Every study receives at least one cell when ``n`` is at least the study count.
    Remaining cells are allocated by deterministic largest-remainder proportions.
    """
    if stratifyKey not in adata.obs:
        raise ValueError(f"adata.obs is missing stratify key {stratifyKey!r}")
    if adata.n_obs < 1:
        raise ValueError("Cannot sample an empty AnnData")

    groups = _study_positions(adata, stratifyKey)
    sizes = {study: int(positions.size) for study, positions in groups.items()}
    alloc = _allocate_study_counts(sizes, n)

    rng = np.random.default_rng(seed)
    chosen: list[np.ndarray] = []
    for study in sorted(alloc):
        positions = groups[study]
        take = alloc[study]
        if take == positions.size:
            selected = positions
        else:
            selected = rng.choice(positions, size=take, replace=False)
        chosen.append(np.sort(selected))

    indices = np.sort(np.concatenate(chosen))
    sampled = adata[indices].copy()
    sampled.uns[SAMPLE_UNS_KEY] = {
        "sourceCells": int(adata.n_obs),
        "sampleCells": int(sampled.n_obs),
        "method": SAMPLE_METHOD,
        "stratifyKey": stratifyKey,
        "seed": int(seed),
        "nStudies": int(len(sizes)),
    }
    log.info(
        "Sampled %s/%s cells across %s studies (%s, seed=%s)",
        f"{sampled.n_obs:,}",
        f"{adata.n_obs:,}",
        len(sizes),
        SAMPLE_METHOD,
        seed,
    )
    return sampled
