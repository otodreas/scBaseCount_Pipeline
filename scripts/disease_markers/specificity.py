"""Cluster interpretation and gene-specificity helpers for atlas discovery."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from disease_markers.candidates import study_balanced_weights
from disease_markers.concordance import as_string_series, cluster_label_purity, is_usable_label
from disease_markers.validation import same_study_case_control_profiles


def attach_gene_symbols(frame: pd.DataFrame, geneMap: pd.DataFrame, *, geneCol: str = "gene") -> pd.DataFrame:
    out = frame.copy()
    if geneCol not in out.columns or geneMap.empty:
        return out
    mapped = geneMap.rename(columns={"ensembl_id": geneCol})
    return out.merge(mapped, on=geneCol, how="left")


def cluster_interpretation_table(
    obs: pd.DataFrame,
    *,
    clusterKey: str = "leiden_atlas",
    labelKey: str = "cell_type",
    ontologyKey: str = "cell_ontology_term_id",
    studyKey: str = "study_accession",
    sampleKey: str = "SRX_accession",
    highPurity: float = 0.7,
    resolvedMinStudies: int = 3,
) -> pd.DataFrame:
    purity = cluster_label_purity(obs, clusterKey=clusterKey, labelKey=labelKey, studyKey=studyKey)
    support = (
        obs.groupby(as_string_series(obs[clusterKey]), observed=True)
        .agg(
            nCells=(sampleKey, "size"),
            nSamples=(sampleKey, "nunique"),
            nStudiesAll=(studyKey, "nunique"),
            nDiseaseAreas=("diseaseArea", "nunique"),
        )
        .rename_axis("cluster")
        .reset_index()
    )
    ontology_rows: list[dict[str, object]] = []
    for cluster, group in obs.groupby(as_string_series(obs[clusterKey]), observed=True):
        usable = (
            group.loc[is_usable_label(group[ontologyKey]), ontologyKey]
            if ontologyKey in group.columns
            else pd.Series(dtype=object)
        )
        if usable.empty:
            ontology_rows.append(
                {
                    "cluster": str(cluster),
                    "topOntology": None,
                    "topOntologyFraction": float("nan"),
                    "nOntologies": 0,
                }
            )
            continue
        counts = as_string_series(usable).value_counts()
        ontology_rows.append(
            {
                "cluster": str(cluster),
                "topOntology": str(counts.index[0]),
                "topOntologyFraction": float(counts.iloc[0] / counts.sum()),
                "nOntologies": int(counts.size),
            }
        )
    ont = pd.DataFrame(ontology_rows)
    out = support.merge(purity, on="cluster", how="left").merge(ont, on="cluster", how="left")
    out["interpretationStatus"] = np.where(
        (out["topLabelFraction"].fillna(0) >= highPurity)
        & (out["nStudies"].fillna(0) >= resolvedMinStudies)
        & out["topLabel"].notna(),
        "resolved",
        "unresolved",
    )
    out["interpretedCellType"] = np.where(
        out["interpretationStatus"].eq("resolved"),
        out["topLabel"],
        "unresolved",
    )
    return out.sort_values("cluster").reset_index(drop=True)


def tau_specificity(values: np.ndarray) -> float:
    """Yanai tau for one gene across cluster means. Near 1 means one-cluster specific."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0 or not np.isfinite(arr).any():
        return float("nan")
    arr = np.clip(arr, a_min=0.0, a_max=None)
    maximum = float(arr.max())
    if maximum <= 0:
        return float("nan")
    return float(np.sum(1.0 - (arr / maximum)) / (arr.size - 1))


def _matrix_block_to_array(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return np.asarray(matrix.toarray(), dtype=np.float64)
    return np.asarray(matrix, dtype=np.float64)


def gene_specificity_table(
    pdata: AnnData,
    *,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
    minStudies: int = 3,
    minProfilesForGene: int = 5,
    minTotalCounts: int = 20,
    geneChunkSize: int = 2000,
) -> pd.DataFrame:
    """Compute tau and detection specificity using gene-chunked dense blocks."""
    _ = minStudies
    obs = pdata.obs.copy()
    obs[clusterKey] = as_string_series(obs[clusterKey])
    obs[studyKey] = as_string_series(obs[studyKey])
    if "psbulk_props" not in pdata.layers:
        raise KeyError("pseudobulk object is missing layers['psbulk_props']")

    profile_support = np.asarray((pdata.X > 0).sum(axis=0)).ravel()
    total_counts = np.asarray(pdata.X.sum(axis=0)).ravel()
    keep_genes = (profile_support >= minProfilesForGene) & (total_counts >= minTotalCounts)
    gene_names = pdata.var_names.to_numpy()[keep_genes]
    keep_idx = np.flatnonzero(keep_genes)
    if keep_idx.size == 0:
        return pd.DataFrame()

    clusters = sorted(obs[clusterKey].unique(), key=lambda x: int(x) if str(x).isdigit() else str(x))
    cluster_index = {cluster: i for i, cluster in enumerate(clusters)}
    n_clusters = len(clusters)
    weights = study_balanced_weights(obs[studyKey]).to_numpy(dtype=np.float64)
    studies = sorted(obs[studyKey].unique())

    rows: list[dict[str, object]] = []
    for start in range(0, keep_idx.size, geneChunkSize):
        block_idx = keep_idx[start : start + geneChunkSize]
        block_genes = gene_names[start : start + geneChunkSize]
        counts = _matrix_block_to_array(pdata.X[:, block_idx])
        props = _matrix_block_to_array(pdata.layers["psbulk_props"][:, block_idx])
        n_genes = counts.shape[1]

        libsize = counts.sum(axis=1)
        libsize = np.where(libsize > 0, libsize, 1.0)
        cpm = np.log1p(counts / libsize[:, None] * 1e4)

        mean_expr = np.zeros((n_clusters, n_genes), dtype=np.float64)
        mean_det = np.zeros_like(mean_expr)
        for cluster, group in obs.groupby(clusterKey, observed=True):
            positions = obs.index.get_indexer(group.index)
            w = weights[positions]
            w = w / w.sum()
            mean_expr[cluster_index[str(cluster)]] = np.average(cpm[positions], axis=0, weights=w)
            mean_det[cluster_index[str(cluster)]] = np.average(props[positions], axis=0, weights=w)

        study_top_votes = np.zeros((n_clusters, n_genes), dtype=np.int32)
        n_studies_scored = np.zeros(n_genes, dtype=np.int32)
        for study in studies:
            study_obs = obs[obs[studyKey] == study]
            if study_obs.empty:
                continue
            study_means = np.full((n_clusters, n_genes), np.nan, dtype=np.float64)
            scored = False
            for cluster, group in study_obs.groupby(clusterKey, observed=True):
                positions = obs.index.get_indexer(group.index)
                if positions.size == 0:
                    continue
                study_means[cluster_index[str(cluster)]] = cpm[positions].mean(axis=0)
                scored = True
            if not scored:
                continue
            present = (~np.isnan(study_means)).sum(axis=0) >= 2
            if not present.any():
                continue
            filled = np.where(np.isnan(study_means), -np.inf, study_means)
            top_idx_study = np.argmax(filled, axis=0)
            gene_ids = np.where(present)[0]
            for gene_i in gene_ids:
                study_top_votes[top_idx_study[gene_i], gene_i] += 1
                n_studies_scored[gene_i] += 1

        top_idx = np.argmax(mean_expr, axis=0)
        for gene_i, gene in enumerate(block_genes):
            expr = mean_expr[:, gene_i]
            det = mean_det[:, gene_i]
            top_i = int(top_idx[gene_i])
            top_cluster = clusters[top_i]
            background_det = np.delete(det, top_i)
            background_expr = np.delete(expr, top_i)
            rows.append(
                {
                    "gene": str(gene),
                    "topCluster": top_cluster,
                    "tau": tau_specificity(expr),
                    "meanExprTop": float(expr[top_i]),
                    "meanDetectionTop": float(det[top_i]),
                    "meanDetectionBackground": float(background_det.mean()) if background_det.size else float("nan"),
                    "maxDetectionBackground": float(background_det.max()) if background_det.size else float("nan"),
                    "meanExprBackground": float(background_expr.mean()) if background_expr.size else float("nan"),
                    "nStudiesAgreeTop": int(study_top_votes[top_i, gene_i]),
                    "nStudiesScored": int(n_studies_scored[gene_i]),
                    "detectionDifference": float(det[top_i] - (background_det.mean() if background_det.size else 0.0)),
                }
            )
    return pd.DataFrame(rows)


def prevalence_by_arm(
    pdata: AnnData,
    *,
    area: str,
    cluster: str,
    gene: str,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
) -> dict[str, float | int]:
    selected = same_study_case_control_profiles(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=clusterKey,
        studyKey=studyKey,
    )
    if not bool(selected.any()) or gene not in pdata.var_names:
        return {
            "caseDetection": float("nan"),
            "controlDetection": float("nan"),
            "detectionDelta": float("nan"),
            "nCaseProfiles": 0,
            "nControlProfiles": 0,
        }
    sub = pdata[selected.to_numpy(), [gene]].copy()
    diseased = sub.obs["diseased"].astype("boolean").eq(True).fillna(False)
    props = _matrix_block_to_array(sub.layers["psbulk_props"]).ravel()
    case = props[diseased.to_numpy()]
    control = props[~diseased.to_numpy()]
    case_det = float(case.mean()) if case.size else float("nan")
    control_det = float(control.mean()) if control.size else float("nan")
    return {
        "caseDetection": case_det,
        "controlDetection": control_det,
        "detectionDelta": case_det - control_det
        if np.isfinite(case_det) and np.isfinite(control_det)
        else float("nan"),
        "nCaseProfiles": int(case.size),
        "nControlProfiles": int(control.size),
    }


def control_home_cluster(
    pdata: AnnData,
    *,
    gene: str,
    studyKey: str = "study_accession",
    clusterKey: str = "leiden_atlas",
) -> tuple[str | None, float]:
    """Cluster with highest study-balanced mean detection among nondiseased profiles."""
    if gene not in pdata.var_names:
        return None, float("nan")
    controls = pdata.obs["diseased"].astype("boolean").eq(False).fillna(False)
    if not bool(controls.any()):
        return None, float("nan")
    sub = pdata[controls.to_numpy(), [gene]]
    props = _matrix_block_to_array(sub.layers["psbulk_props"]).ravel()
    obs = sub.obs.copy()
    obs[clusterKey] = as_string_series(obs[clusterKey])
    weights = study_balanced_weights(obs[studyKey]).to_numpy(dtype=np.float64)
    best_cluster = None
    best_value = -1.0
    for cluster, group in obs.groupby(clusterKey, observed=True):
        positions = obs.index.get_indexer(group.index)
        w = weights[positions]
        value = float(np.average(props[positions], weights=w))
        if value > best_value:
            best_value = value
            best_cluster = str(cluster)
    return best_cluster, best_value


def study_direction_agreement(
    pdata: AnnData,
    *,
    gene: str,
    area: str,
    cluster: str,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
) -> dict[str, object]:
    """Per-study sign agreement from case vs control mean log1p-CPM differences."""
    selected = same_study_case_control_profiles(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=clusterKey,
        studyKey=studyKey,
    )
    if not bool(selected.any()) or gene not in pdata.var_names:
        return {
            "nStudiesAgree": 0,
            "nStudiesScored": 0,
            "studyDirectionAgreement": "0/0",
            "majorityDirection": None,
        }
    sub = pdata[selected.to_numpy(), [gene]].copy()
    obs = sub.obs.copy()
    counts = _matrix_block_to_array(sub.X).ravel()
    libsize = np.asarray(sub.X.sum(axis=1)).ravel() if sub.n_vars else np.ones(sub.n_obs)
    # Use the selected gene only for expression, but library size from full profile when available.
    if sub.n_obs and gene in pdata.var_names:
        full = pdata[selected.to_numpy()]
        libsize = np.asarray(full.X.sum(axis=1)).ravel()
    libsize = np.where(libsize > 0, libsize, 1.0)
    expr = np.log1p(counts / libsize * 1e4)
    diseased = obs["diseased"].astype("boolean").eq(True).fillna(False).to_numpy()
    studies = as_string_series(obs[studyKey])
    signs: list[int] = []
    for study in sorted(studies.unique()):
        mask = studies.eq(study).to_numpy()
        case_vals = expr[mask & diseased]
        control_vals = expr[mask & ~diseased]
        if case_vals.size == 0 or control_vals.size == 0:
            continue
        delta = float(case_vals.mean() - control_vals.mean())
        if not np.isfinite(delta) or delta == 0:
            continue
        signs.append(1 if delta > 0 else -1)
    if not signs:
        return {
            "nStudiesAgree": 0,
            "nStudiesScored": 0,
            "studyDirectionAgreement": "0/0",
            "majorityDirection": None,
        }
    majority = 1 if sum(signs) >= 0 else -1
    n_agree = int(sum(sign == majority for sign in signs))
    return {
        "nStudiesAgree": n_agree,
        "nStudiesScored": int(len(signs)),
        "studyDirectionAgreement": f"{n_agree}/{len(signs)}",
        "majorityDirection": "up" if majority > 0 else "down",
    }
