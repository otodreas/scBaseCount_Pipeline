import numpy as np
import pandas as pd
import scanpy as sc

from disease_markers.config import DiseaseMarkersConfig


def _sample_means(pdata: sc.AnnData, sample_ids: pd.Index, gene: str) -> float:
    if gene not in pdata.var_names:
        return float("nan")
    sub = pdata[list(sample_ids), gene]
    matrix = sub.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    values = np.asarray(matrix).ravel()
    if values.size == 0:
        return float("nan")
    return float(np.log1p(values).mean())


def flag_study_consistency(
    pdata: sc.AnnData,
    area: str,
    de_hits: pd.DataFrame,
    cfg: DiseaseMarkersConfig,
) -> pd.DataFrame:
    """Add studyConsistency column via per-study area-vs-rest mean log1p sign agreement."""
    if de_hits.empty:
        return de_hits

    meta = pdata.obs
    area_rows = meta[meta["diseaseArea"].astype(str) == area]
    studies = sorted(area_rows[cfg.studyKey].astype(str).unique())
    rest_samples = meta.index[meta["diseaseArea"].astype(str) != area]

    flags: list[bool] = []
    for gene in de_hits["gene"].astype(str):
        study_signs: list[int] = []
        for study in studies:
            area_samples = area_rows.index[area_rows[cfg.studyKey].astype(str) == study]
            if len(area_samples) == 0:
                continue
            area_mean = _sample_means(pdata, area_samples, gene)
            rest_mean = _sample_means(pdata, rest_samples, gene)
            if np.isnan(area_mean) or np.isnan(rest_mean):
                continue
            diff = area_mean - rest_mean
            if diff == 0:
                study_signs.append(0)
            else:
                study_signs.append(1 if diff > 0 else -1)
        non_zero = [s for s in study_signs if s != 0]
        if len(non_zero) <= 1:
            flags.append(True)
        else:
            flags.append(len(set(non_zero)) == 1 and non_zero[0] > 0)

    out = de_hits.copy()
    out["studyConsistent"] = flags
    return out
