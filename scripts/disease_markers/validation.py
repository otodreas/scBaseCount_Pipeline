"""Study-aware abundance and pseudobulk DE helpers for atlas cluster candidates."""

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from disease_markers.candidates import OTHER_AREA
from disease_markers.concordance import as_string_series


def _inverse_variance_mean(effects: list[float], variances: list[float]) -> tuple[float, float]:
    weights = np.asarray([1.0 / max(variance, 1e-12) for variance in variances], dtype=np.float64)
    effect_arr = np.asarray(effects, dtype=np.float64)
    combined = float(np.sum(weights * effect_arr) / np.sum(weights))
    combined_se = float(np.sqrt(1.0 / np.sum(weights)))
    return combined, combined_se


def counts_dataframe(adata) -> pd.DataFrame:
    """Return integer count matrix as a DataFrame indexed like ``adata``."""
    matrix = adata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    counts = np.rint(np.asarray(matrix, dtype=np.float64)).astype(int)
    return pd.DataFrame(counts, index=adata.obs_names, columns=adata.var_names)


def filter_pseudobulk_profiles(
    pdata,
    *,
    minCellsPerProfile: int = 10,
    cellsKey: str = "psbulk_cells",
):
    """Keep pseudobulk profiles with enough contributing cells."""
    if cellsKey not in pdata.obs.columns:
        raise KeyError(f"pseudobulk obs is missing {cellsKey!r}")
    keep = pdata.obs[cellsKey].astype(float) >= float(minCellsPerProfile)
    return pdata[keep.to_numpy()].copy()


def same_study_case_control_profiles(
    pdata,
    *,
    area: str,
    cluster: str | None = None,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
    diseaseAreaKey: str = "diseaseArea",
    diseasedKey: str = "diseased",
) -> pd.DataFrame:
    """Return a boolean mask for same-study case/control profiles for one area."""
    obs = pdata.obs.copy()
    for key in (studyKey, diseaseAreaKey, diseasedKey):
        if key not in obs.columns:
            raise KeyError(f"pseudobulk obs is missing {key!r}")
    if cluster is not None:
        if clusterKey not in obs.columns:
            raise KeyError(f"pseudobulk obs is missing {clusterKey!r}")
        obs = obs[as_string_series(obs[clusterKey]) == str(cluster)]

    diseased = obs[diseasedKey].astype("boolean")
    area_mask = as_string_series(obs[diseaseAreaKey]) == str(area)
    case_mask = diseased.eq(True).fillna(False) & area_mask
    control_mask = diseased.eq(False).fillna(False)
    case_studies = set(as_string_series(obs.loc[case_mask, studyKey]))
    control_studies = set(as_string_series(obs.loc[control_mask, studyKey]))
    overlap = case_studies & control_studies
    selected = (case_mask | control_mask) & as_string_series(obs[studyKey]).isin(overlap)
    return selected.reindex(pdata.obs_names, fill_value=False)


def filter_two_sided_de(
    results: pd.DataFrame,
    *,
    padj: float = 0.05,
    lfc: float = 1.0,
) -> pd.DataFrame:
    """Keep significant up- and down-regulated genes."""
    if results.empty:
        return results.copy()
    required = {"padj", "log2FoldChange"}
    missing = required - set(results.columns)
    if missing:
        raise KeyError(f"DE results missing required columns: {sorted(missing)}")
    keep = results["padj"].notna() & (results["padj"] <= padj) & (results["log2FoldChange"].abs() >= lfc)
    return results.loc[keep].copy().reset_index(drop=True)


def disease_vs_control_deseq2(
    pdata,
    *,
    area: str,
    cluster: str,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
    diseaseAreaKey: str = "diseaseArea",
    diseasedKey: str = "diseased",
    minProfilesPerGroup: int = 2,
    includeStudy: bool = True,
) -> pd.DataFrame:
    """Run DESeq2 for one cluster and disease area using same-study controls.

    The design includes study when more than one overlapping study is present.
    Distinct disease areas are never pooled into one case class.
    """
    if str(area) == OTHER_AREA:
        return pd.DataFrame()

    selected = same_study_case_control_profiles(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=clusterKey,
        studyKey=studyKey,
        diseaseAreaKey=diseaseAreaKey,
        diseasedKey=diseasedKey,
    )
    if not bool(selected.any()):
        return pd.DataFrame()

    analysis = pdata[selected.to_numpy()].copy()
    meta = analysis.obs.copy()
    diseased = meta[diseasedKey].astype("boolean").eq(True).fillna(False)
    meta["group"] = np.where(diseased.to_numpy(dtype=bool), str(area), "nonDiseased")
    meta["group"] = meta["group"].astype(str)
    meta["study"] = as_string_series(meta[studyKey]).astype(str)

    # A study can lose one arm after min-cell filtering; keep only balanced studies.
    balanced_studies: list[str] = []
    for study, study_meta in meta.groupby("study", observed=True):
        groups = set(study_meta["group"].tolist())
        if str(area) in groups and "nonDiseased" in groups:
            balanced_studies.append(str(study))
    if not balanced_studies:
        return pd.DataFrame()
    keep = meta["study"].isin(balanced_studies)
    analysis = analysis[keep.to_numpy()].copy()
    meta = meta.loc[keep].copy()

    group_counts = meta["group"].value_counts()
    if group_counts.get(str(area), 0) < minProfilesPerGroup:
        return pd.DataFrame()
    if group_counts.get("nonDiseased", 0) < minProfilesPerGroup:
        return pd.DataFrame()

    n_studies = int(meta["study"].nunique())
    design_factors = ["group", "study"] if includeStudy and n_studies >= 2 else ["group"]
    design_meta = meta.loc[:, design_factors].copy()
    for col in design_factors:
        design_meta[col] = design_meta[col].astype(str)
    counts = counts_dataframe(analysis)
    dds = DeseqDataSet(
        counts=counts,
        metadata=design_meta,
        design_factors=design_factors,
        ref_level=["group", "nonDiseased"],
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", str(area), "nonDiseased"], quiet=True)
    stats.summary()
    results = stats.results_df.copy()
    results["gene"] = results.index.astype(str)
    results.insert(0, "cluster", str(cluster))
    results.insert(1, "diseaseArea", str(area))
    results["nCaseProfiles"] = int(group_counts.get(str(area), 0))
    results["nControlProfiles"] = int(group_counts.get("nonDiseased", 0))
    results["nStudies"] = n_studies
    return results.reset_index(drop=True)


def de_supported_candidates(
    pdata,
    candidates: pd.DataFrame,
    contrastSupport: pd.DataFrame,
    *,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
    padj: float = 0.05,
    lfc: float = 1.0,
    minOverlapStudies: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run study-aware DE for disease-associated cluster candidates only."""
    if candidates.empty or contrastSupport.empty:
        empty = pd.DataFrame()
        return empty, empty, empty

    disease_clusters = set(
        candidates.loc[candidates["isDiseaseAssociatedCandidate"].astype(bool), "cluster"].astype(str)
    )
    eligible = contrastSupport[
        contrastSupport["cluster"].astype(str).isin(disease_clusters)
        & (contrastSupport["nOverlapStudies"] >= minOverlapStudies)
        & contrastSupport["eligibleForContrast"].astype(bool)
    ].copy()

    summary_rows: list[dict[str, object]] = []
    result_frames: list[pd.DataFrame] = []
    hit_frames: list[pd.DataFrame] = []

    for row in eligible.itertuples(index=False):
        full = disease_vs_control_deseq2(
            pdata,
            area=str(row.diseaseArea),
            cluster=str(row.cluster),
            clusterKey=clusterKey,
            studyKey=studyKey,
        )
        hits = filter_two_sided_de(full, padj=padj, lfc=lfc) if not full.empty else full
        summary_rows.append(
            {
                "cluster": str(row.cluster),
                "diseaseArea": str(row.diseaseArea),
                "nOverlapStudies": int(row.nOverlapStudies),
                "nCaseProfiles": int(row.nCaseProfiles),
                "nControlProfiles": int(row.nControlProfiles),
                "nHits": int(len(hits)),
                "ran": not full.empty,
            }
        )
        if full.empty:
            continue
        result_frames.append(full)
        if not hits.empty:
            hit_frames.append(hits)

    summary = pd.DataFrame(summary_rows)
    results = pd.concat(result_frames, ignore_index=True) if result_frames else pd.DataFrame()
    hits = pd.concat(hit_frames, ignore_index=True) if hit_frames else pd.DataFrame()
    return summary, hits, results


def sample_cluster_proportions(
    obs: pd.DataFrame,
    *,
    clusterKey: str = "leiden_atlas",
    sampleKey: str = "SRX_accession",
    studyKey: str = "study_accession",
    diseaseAreaKey: str = "diseaseArea",
    diseasedKey: str = "diseased",
) -> pd.DataFrame:
    """Compute per-sample cluster proportions for abundance testing."""
    required = {clusterKey, sampleKey, studyKey, diseaseAreaKey, diseasedKey}
    missing = required - set(obs.columns)
    if missing:
        raise KeyError(f"obs is missing required columns: {sorted(missing)}")

    frame = obs.copy()
    frame[clusterKey] = as_string_series(frame[clusterKey])
    frame[sampleKey] = as_string_series(frame[sampleKey])
    frame[studyKey] = as_string_series(frame[studyKey])
    frame[diseaseAreaKey] = as_string_series(frame[diseaseAreaKey])
    frame[diseasedKey] = frame[diseasedKey].astype("boolean")

    sample_meta = (
        frame.groupby(sampleKey, observed=True)
        .agg(
            studyAccession=(studyKey, "first"),
            diseaseArea=(diseaseAreaKey, "first"),
            diseased=(diseasedKey, "first"),
            nCells=(clusterKey, "size"),
        )
        .reset_index()
        .rename(columns={sampleKey: "srxAccession"})
    )
    counts = (
        frame.groupby([sampleKey, clusterKey], observed=True)
        .size()
        .rename("nClusterCells")
        .reset_index()
        .rename(columns={sampleKey: "srxAccession", clusterKey: "cluster"})
    )
    merged = counts.merge(sample_meta, on="srxAccession", how="left")
    merged["proportion"] = merged["nClusterCells"] / merged["nCells"].astype(float)
    return merged


def differential_abundance_by_study(
    proportions: pd.DataFrame,
    *,
    area: str,
    cluster: str,
    minSamplesPerArm: int = 2,
) -> pd.DataFrame:
    """Study-stratified mean proportion difference for one cluster and disease area.

    Cells are never treated as independent replicates. Each study contributes one
    case-vs-control mean-proportion difference. When multiple studies are present,
    a fixed-effect combination is reported when statsmodels is available.
    """
    if proportions.empty:
        return pd.DataFrame()

    required = {
        "cluster",
        "srxAccession",
        "studyAccession",
        "diseaseArea",
        "diseased",
        "proportion",
    }
    missing = required - set(proportions.columns)
    if missing:
        raise KeyError(f"proportions missing required columns: {sorted(missing)}")

    frame = proportions.copy()
    frame["cluster"] = as_string_series(frame["cluster"])
    frame["diseaseArea"] = as_string_series(frame["diseaseArea"])
    frame["studyAccession"] = as_string_series(frame["studyAccession"])
    frame["diseased"] = frame["diseased"].astype("boolean")
    subset = frame[frame["cluster"] == str(cluster)].copy()
    if subset.empty:
        return pd.DataFrame()

    case_mask = subset["diseased"].eq(True).fillna(False) & (subset["diseaseArea"] == str(area))
    control_mask = subset["diseased"].eq(False).fillna(False)
    cases = subset.loc[case_mask]
    controls = subset.loc[control_mask]
    overlap = sorted(set(cases["studyAccession"]) & set(controls["studyAccession"]))
    if not overlap:
        return pd.DataFrame()

    study_rows: list[dict[str, object]] = []
    effects: list[float] = []
    variances: list[float] = []
    for study in overlap:
        case_vals = cases.loc[cases["studyAccession"] == study, "proportion"].astype(float)
        control_vals = controls.loc[controls["studyAccession"] == study, "proportion"].astype(float)
        if len(case_vals) < minSamplesPerArm or len(control_vals) < minSamplesPerArm:
            continue
        effect = float(case_vals.mean() - control_vals.mean())
        var = float(case_vals.var(ddof=1) / len(case_vals) + control_vals.var(ddof=1) / len(control_vals))
        if not np.isfinite(var) or var <= 0:
            var = 1e-8
        study_rows.append(
            {
                "cluster": str(cluster),
                "diseaseArea": str(area),
                "studyAccession": study,
                "nCaseSamples": int(len(case_vals)),
                "nControlSamples": int(len(control_vals)),
                "meanCaseProportion": float(case_vals.mean()),
                "meanControlProportion": float(control_vals.mean()),
                "effect": effect,
                "variance": var,
            }
        )
        effects.append(effect)
        variances.append(var)

    if not study_rows:
        return pd.DataFrame()

    out = pd.DataFrame(study_rows)
    combined_effect, combined_se = _inverse_variance_mean(effects, variances)
    out["combinedEffect"] = combined_effect
    out["combinedSe"] = combined_se
    out["nStudiesCombined"] = len(effects)
    return out


def shared_direction_genes(
    deHits: pd.DataFrame,
    *,
    minDiseaseAreas: int = 2,
) -> pd.DataFrame:
    """Genes with the same sign of effect across multiple disease areas."""
    if deHits.empty:
        return pd.DataFrame(columns=["gene", "nDiseaseAreas", "nClusters", "meanLog2FoldChange", "direction"])
    required = {"gene", "diseaseArea", "cluster", "log2FoldChange"}
    missing = required - set(deHits.columns)
    if missing:
        raise KeyError(f"DE hits missing required columns: {sorted(missing)}")

    frame = deHits.copy()
    frame["gene"] = frame["gene"].astype(str)
    frame["diseaseArea"] = frame["diseaseArea"].astype(str)
    frame["direction"] = np.sign(frame["log2FoldChange"].astype(float)).astype(int)
    frame = frame[frame["direction"] != 0]

    rows: list[dict[str, object]] = []
    for gene, group in frame.groupby("gene", observed=True):
        # Require a consistent direction within each disease area first.
        area_dirs = group.groupby("diseaseArea", observed=True)["direction"].agg(
            lambda values: int(values.mode().iloc[0]) if not values.mode().empty else 0
        )
        area_dirs = area_dirs[area_dirs != 0]
        if area_dirs.nunique() != 1:
            continue
        if int(area_dirs.size) < minDiseaseAreas:
            continue
        direction = int(area_dirs.iloc[0])
        rows.append(
            {
                "gene": str(gene),
                "nDiseaseAreas": int(area_dirs.size),
                "nClusters": int(group["cluster"].nunique()),
                "meanLog2FoldChange": float(group["log2FoldChange"].mean()),
                "direction": "up" if direction > 0 else "down",
            }
        )
    if not rows:
        return pd.DataFrame(columns=["gene", "nDiseaseAreas", "nClusters", "meanLog2FoldChange", "direction"])
    return pd.DataFrame(rows).sort_values(["nDiseaseAreas", "gene"], ascending=[False, True]).reset_index(drop=True)
