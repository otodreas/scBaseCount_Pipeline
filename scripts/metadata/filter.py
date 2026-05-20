from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

from metadata.config import MetadataConfig
from metadata.regexes import CANCER_RE, DISEASE_MAP, LUNG_DISEASE_RE, LUNG_TISSUE_RE, NORMAL_HEALTHY_RE


@dataclass
class FilterResult:
    sampleKnown: pd.DataFrame
    lungUnion: pd.DataFrame
    lungIntersection: pd.DataFrame
    lungIntersectionCancer: pd.DataFrame


def filter_lung(sample: pd.DataFrame, cfg: MetadataConfig) -> FilterResult:
    sample = sample[sample["obs_count"] >= cfg.minObsCount].copy()

    unknown_mask = (
        sample[["disease", "tissue"]].apply(lambda col: col.str.contains(NORMAL_HEALTHY_RE, na=True)).any(axis=1)
    )
    sample_known = sample.loc[~unknown_mask]

    lung_disease_mask = sample_known["disease"].str.contains(LUNG_DISEASE_RE, regex=True, na=False)
    lung_tissue_mask = sample_known["tissue"].str.contains(LUNG_TISSUE_RE, regex=True, na=False)

    lung_union = sample_known.loc[lung_disease_mask | lung_tissue_mask].reset_index(drop=True)
    lung_intersection = sample_known.loc[lung_disease_mask & lung_tissue_mask].reset_index(drop=True)

    cancer_mask = lung_intersection["disease"].str.contains(CANCER_RE, regex=True, na=False)
    lung_intersection_cancer = lung_intersection.loc[cancer_mask].reset_index(drop=True)

    return FilterResult(
        sampleKnown=sample_known.reset_index(drop=True),
        lungUnion=lung_union,
        lungIntersection=lung_intersection,
        lungIntersectionCancer=lung_intersection_cancer,
    )


def available_disease_labels() -> list[str]:
    """Return the named disease labels recognised by filter_by_disease (from DISEASE_MAP)."""
    return [label for label, _ in DISEASE_MAP]


def filter_by_disease(
    samplesDf: pd.DataFrame,
    diseaseLabel: str | None = None,
    diseaseRegex: str | re.Pattern[str] | None = None,
) -> pd.DataFrame:
    """Return rows of samplesDf whose disease column matches the chosen filter.

    Pass diseaseLabel for a named pattern from DISEASE_MAP, or diseaseRegex for a
    custom pattern (string is compiled case-insensitive). When both are None the
    frame is returned unchanged. Passing both raises ValueError. An unknown
    diseaseLabel raises KeyError listing the valid labels.
    """
    if diseaseLabel is not None and diseaseRegex is not None:
        raise ValueError("Pass either diseaseLabel or diseaseRegex, not both.")

    if diseaseLabel is None and diseaseRegex is None:
        return samplesDf

    if diseaseLabel is not None:
        labels = available_disease_labels()
        if diseaseLabel not in labels:
            raise KeyError(f"Unknown diseaseLabel {diseaseLabel!r}. Valid labels: {labels}")
        pattern: re.Pattern[str] = next(p for label, p in DISEASE_MAP if label == diseaseLabel)
    else:
        pattern = re.compile(diseaseRegex, re.IGNORECASE) if isinstance(diseaseRegex, str) else diseaseRegex

    mask = samplesDf["disease"].str.contains(pattern, regex=True, na=False)
    return samplesDf.loc[mask].reset_index(drop=True)
