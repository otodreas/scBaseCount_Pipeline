import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from disease_markers.config import DiseaseMarkersConfig


def _counts_dataframe(pdata: sc.AnnData) -> pd.DataFrame:
    """Return integer sample-by-gene counts from pseudobulk AnnData."""
    matrix = pdata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    counts = np.asarray(matrix, dtype=np.float64)
    counts = np.rint(counts).astype(int)
    return pd.DataFrame(counts, index=pdata.obs_names, columns=pdata.var_names)


def _area_meets_depth(pdata: sc.AnnData, area: str, cfg: DiseaseMarkersConfig) -> bool:
    area_samples = pdata.obs[pdata.obs["diseaseArea"].astype(str) == area]
    if area_samples.empty:
        return False
    if area_samples[cfg.sampleKey].nunique() < cfg.minSamplesPerArea:
        return False
    return area_samples[cfg.studyKey].nunique() >= cfg.minStudiesPerArea


def one_vs_rest_de(
    pdata: sc.AnnData,
    area: str,
    cfg: DiseaseMarkersConfig,
) -> pd.DataFrame:
    """Run pydeseq2 one-vs-rest DE for a disease area on pseudobulk profiles."""
    if not _area_meets_depth(pdata, area, cfg):
        return pd.DataFrame()

    metadata = pdata.obs.copy()
    metadata["group"] = np.where(metadata["diseaseArea"].astype(str) == area, area, "rest")
    if metadata["group"].nunique() < 2:
        return pd.DataFrame()
    if (metadata["group"] == area).sum() < 2:
        return pd.DataFrame()

    counts = _counts_dataframe(pdata)
    if counts.shape[0] < 4:
        return pd.DataFrame()

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="group",
        ref_level=["group", "rest"],
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", area, "rest"], quiet=True)
    stats.summary()
    results = stats.results_df.copy()
    results["gene"] = results.index.astype(str)
    up = results[
        (results["padj"] <= cfg.padjThreshold)
        & (results["log2FoldChange"] >= cfg.lfcThreshold)
        & results["padj"].notna()
    ]
    return up.reset_index(drop=True)


def de_areas_in_cluster(
    pdata: sc.AnnData,
    cfg: DiseaseMarkersConfig,
) -> dict[str, pd.DataFrame]:
    """Run one-vs-rest DE for each disease area present in a cluster pseudobulk."""
    areas = sorted(pdata.obs["diseaseArea"].astype(str).unique())
    out: dict[str, pd.DataFrame] = {}
    for area in areas:
        if area == "Control":
            continue
        hits = one_vs_rest_de(pdata, area, cfg)
        if not hits.empty:
            out[area] = hits
    return out
