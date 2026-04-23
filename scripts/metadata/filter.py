from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from metadata.config import MetadataConfig
from metadata.regexes import CANCER_RE, LUNG_DISEASE_RE, LUNG_TISSUE_RE, NORMAL_HEALTHY_RE


@dataclass
class FilterResult:
    sampleKnown: pd.DataFrame
    lungUnion: pd.DataFrame
    lungIntersection: pd.DataFrame
    lungIntersectionCancer: pd.DataFrame


def filter_lung(sample: pd.DataFrame, cfg: MetadataConfig) -> FilterResult:
    sample = sample[sample["obs_count"] >= cfg.minObsCount].copy()

    unknown_mask = (
        sample[["disease", "tissue"]]
        .apply(lambda col: col.str.contains(NORMAL_HEALTHY_RE, na=True))
        .any(axis=1)
    )
    sample_known = sample.loc[~unknown_mask]

    lung_disease_mask = sample_known["disease"].str.contains(LUNG_DISEASE_RE, regex=True, na=False)
    lung_tissue_mask  = sample_known["tissue"].str.contains(LUNG_TISSUE_RE,   regex=True, na=False)

    lung_union        = sample_known.loc[lung_disease_mask | lung_tissue_mask].reset_index(drop=True)
    lung_intersection = sample_known.loc[lung_disease_mask & lung_tissue_mask].reset_index(drop=True)

    cancer_mask            = lung_intersection["disease"].str.contains(CANCER_RE, regex=True, na=False)
    lung_intersection_cancer = lung_intersection.loc[cancer_mask].reset_index(drop=True)

    return FilterResult(
        sampleKnown=sample_known.reset_index(drop=True),
        lungUnion=lung_union,
        lungIntersection=lung_intersection,
        lungIntersectionCancer=lung_intersection_cancer,
    )
