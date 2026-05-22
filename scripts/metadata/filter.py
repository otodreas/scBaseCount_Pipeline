from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from metadata.config import MetadataConfig
from metadata.regexes import LUNG_CANCER_RE, LUNG_DISEASE_RE, LUNG_TISSUE_RE, NORMAL_HEALTHY_RE


@dataclass
class FilterResult:
    sampleKnown: pd.DataFrame
    lungIntersection: pd.DataFrame
    lungIntersectionCancer: pd.DataFrame


def filter_lung(sample: pd.DataFrame, cfg: MetadataConfig) -> FilterResult:
    """Apply the lung-intersection cascade and return sampleKnown, lungIntersection, lungIntersectionCancer."""
    sample = sample[sample["obs_count"] >= cfg.minObsCount].copy()

    unknown_mask = sample["disease"].str.contains(NORMAL_HEALTHY_RE, na=True) | sample["tissue"].str.contains(
        NORMAL_HEALTHY_RE, na=True
    )
    sample_known = sample.loc[~unknown_mask]

    lung_disease_mask = sample_known["disease"].str.contains(LUNG_DISEASE_RE, regex=True, na=False)
    lung_tissue_mask = sample_known["tissue"].str.contains(LUNG_TISSUE_RE, regex=True, na=False)

    lung_intersection = sample_known.loc[lung_disease_mask & lung_tissue_mask].reset_index(drop=True)

    cancer_mask = lung_intersection["disease"].str.contains(LUNG_CANCER_RE, regex=True, na=False)
    lung_intersection_cancer = lung_intersection.loc[cancer_mask].reset_index(drop=True)

    return FilterResult(
        sampleKnown=sample_known.reset_index(drop=True),
        lungIntersection=lung_intersection,
        lungIntersectionCancer=lung_intersection_cancer,
    )
