"""Compare integrated Leiden clusters with preserved source cell-type labels."""

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
    v_measure_score,
)

MISSING_LABELS: frozenset[str] = frozenset({"", "UNKNOWN", "unknown", "NA", "nan", "None"})


def as_string_series(values: pd.Series) -> pd.Series:
    return values.astype("string").fillna("").str.strip()


def is_usable_label(values: pd.Series) -> pd.Series:
    text = as_string_series(values)
    return text.ne("") & ~text.isin(MISSING_LABELS)


def source_label_contingency(
    obs: pd.DataFrame,
    *,
    clusterKey: str,
    labelKey: str,
) -> pd.DataFrame:
    """Return cluster-by-label cell counts for usable source labels."""
    for key in (clusterKey, labelKey):
        if key not in obs.columns:
            raise KeyError(f"obs is missing required column {key!r}")

    usable = is_usable_label(obs[labelKey])
    subset = obs.loc[usable, [clusterKey, labelKey]].copy()
    subset[clusterKey] = as_string_series(subset[clusterKey])
    subset[labelKey] = as_string_series(subset[labelKey])
    return pd.crosstab(subset[clusterKey], subset[labelKey])


def _metric_row(clusters: pd.Series, labels: pd.Series) -> dict[str, float]:
    cluster_vals = as_string_series(clusters).to_numpy()
    label_vals = as_string_series(labels).to_numpy()
    if len(np.unique(cluster_vals)) < 2 or len(np.unique(label_vals)) < 2:
        return {
            "nmi": float("nan"),
            "ari": float("nan"),
            "homogeneity": float("nan"),
            "completeness": float("nan"),
            "vMeasure": float("nan"),
        }
    return {
        "nmi": float(normalized_mutual_info_score(label_vals, cluster_vals)),
        "ari": float(adjusted_rand_score(label_vals, cluster_vals)),
        "homogeneity": float(homogeneity_score(label_vals, cluster_vals)),
        "completeness": float(completeness_score(label_vals, cluster_vals)),
        "vMeasure": float(v_measure_score(label_vals, cluster_vals)),
    }


def concordance_summary(
    obs: pd.DataFrame,
    *,
    clusterKey: str,
    labelKey: str,
    studyKey: str,
) -> pd.DataFrame:
    """Global and study-macro concordance between clusters and source labels.

    Source labels are treated as weak references from individual datasets, not
    ground truth. Study-macro means average metric values across studies with
    usable labels so large studies cannot dominate the summary.
    """
    for key in (clusterKey, labelKey, studyKey):
        if key not in obs.columns:
            raise KeyError(f"obs is missing required column {key!r}")

    usable = is_usable_label(obs[labelKey])
    subset = obs.loc[usable, [clusterKey, labelKey, studyKey]].copy()
    if subset.empty:
        return pd.DataFrame(
            columns=[
                "scope",
                "studyAccession",
                "nCells",
                "nClusters",
                "nLabels",
                "nmi",
                "ari",
                "homogeneity",
                "completeness",
                "vMeasure",
            ]
        )

    subset[clusterKey] = as_string_series(subset[clusterKey])
    subset[labelKey] = as_string_series(subset[labelKey])
    subset[studyKey] = as_string_series(subset[studyKey])

    rows: list[dict[str, object]] = []
    global_metrics = _metric_row(subset[clusterKey], subset[labelKey])
    rows.append(
        {
            "scope": "global",
            "studyAccession": None,
            "nCells": int(len(subset)),
            "nClusters": int(subset[clusterKey].nunique()),
            "nLabels": int(subset[labelKey].nunique()),
            **global_metrics,
        }
    )

    study_rows: list[dict[str, object]] = []
    for study, study_obs in subset.groupby(studyKey, observed=True):
        metrics = _metric_row(study_obs[clusterKey], study_obs[labelKey])
        study_rows.append(
            {
                "scope": "study",
                "studyAccession": str(study),
                "nCells": int(len(study_obs)),
                "nClusters": int(study_obs[clusterKey].nunique()),
                "nLabels": int(study_obs[labelKey].nunique()),
                **metrics,
            }
        )
    rows.extend(study_rows)

    if study_rows:
        study_frame = pd.DataFrame(study_rows)
        macro = {
            "scope": "studyMacro",
            "studyAccession": None,
            "nCells": int(study_frame["nCells"].sum()),
            "nClusters": int(subset[clusterKey].nunique()),
            "nLabels": int(subset[labelKey].nunique()),
        }
        for metric in ("nmi", "ari", "homogeneity", "completeness", "vMeasure"):
            macro[metric] = float(study_frame[metric].mean(skipna=True))
        rows.append(macro)

    return pd.DataFrame(rows)


def cluster_label_purity(
    obs: pd.DataFrame,
    *,
    clusterKey: str,
    labelKey: str,
    studyKey: str | None = None,
) -> pd.DataFrame:
    """Per-cluster purity against source labels, plus support counts."""
    for key in (clusterKey, labelKey):
        if key not in obs.columns:
            raise KeyError(f"obs is missing required column {key!r}")
    if studyKey is not None and studyKey not in obs.columns:
        raise KeyError(f"obs is missing required column {studyKey!r}")

    usable = is_usable_label(obs[labelKey])
    columns = [clusterKey, labelKey] if studyKey is None else [clusterKey, labelKey, studyKey]
    subset = obs.loc[usable, columns].copy()
    if subset.empty:
        return pd.DataFrame(
            columns=[
                "cluster",
                "nCellsLabeled",
                "nLabels",
                "nStudies",
                "topLabel",
                "topLabelFraction",
                "labelEntropy",
            ]
        )

    subset[clusterKey] = as_string_series(subset[clusterKey])
    subset[labelKey] = as_string_series(subset[labelKey])
    if studyKey is not None:
        subset[studyKey] = as_string_series(subset[studyKey])

    rows: list[dict[str, object]] = []
    for cluster, group in subset.groupby(clusterKey, observed=True):
        counts = group[labelKey].value_counts(dropna=False)
        total = int(counts.sum())
        probs = counts.to_numpy(dtype=np.float64) / total
        positive = probs[probs > 0]
        entropy = float(-(positive * np.log2(positive)).sum()) if positive.size else 0.0
        rows.append(
            {
                "cluster": str(cluster),
                "nCellsLabeled": total,
                "nLabels": int(counts.size),
                "nStudies": int(group[studyKey].nunique()) if studyKey is not None else None,
                "topLabel": str(counts.index[0]),
                "topLabelFraction": float(counts.iloc[0] / total),
                "labelEntropy": entropy,
            }
        )
    return pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)
