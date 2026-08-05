"""Screen integrated atlas clusters for shared and disease-associated candidates."""

import numpy as np
import pandas as pd

from disease_markers.concordance import as_string_series, is_usable_label

OTHER_AREA = "Other"


def annotate_obs_with_labels(
    obs: pd.DataFrame,
    labelTable: pd.DataFrame,
    *,
    sampleKey: str = "SRX_accession",
    eligibleOnly: bool = True,
) -> pd.DataFrame:
    """Join disease-marker labels onto cell-level obs and optionally keep eligible cells."""
    if sampleKey not in obs.columns:
        raise KeyError(f"obs is missing required column {sampleKey!r}")
    required = {
        "srxAccession",
        "studyAccession",
        "diseaseArea",
        "diseased",
        "isBiologicalControl",
        "controlType",
        "eligible",
    }
    missing = required - set(labelTable.columns)
    if missing:
        raise KeyError(f"label table missing required columns: {sorted(missing)}")

    labels = labelTable.set_index("srxAccession")
    out = obs.copy()
    for col in (
        "diseaseArea",
        "diseaseAreaSource",
        "diseaseName",
        "diseaseOntologyTermId",
        "diseased",
        "isBiologicalControl",
        "controlType",
        "eligible",
        "excludeReason",
    ):
        if col in labels.columns:
            out[col] = out[sampleKey].map(labels[col])
    if "diseaseOntologyTermId" in labels.columns and "disease_ontology_term_id" not in out.columns:
        out["disease_ontology_term_id"] = out[sampleKey].map(labels["diseaseOntologyTermId"])
    if "studyAccession" in labels.columns and "study_accession" not in out.columns:
        out["study_accession"] = out[sampleKey].map(labels["studyAccession"])
    if eligibleOnly:
        eligible = out["eligible"].fillna(False).astype(bool)
        out = out.loc[eligible.to_numpy()].copy()
    out["diseased"] = out["diseased"].astype("boolean")
    return out


def cluster_support_table(
    obs: pd.DataFrame,
    *,
    clusterKey: str = "leiden_atlas",
    sampleKey: str = "SRX_accession",
    studyKey: str = "study_accession",
    labelKey: str = "cell_type",
    ontologyKey: str = "cell_ontology_term_id",
    diseaseAreaKey: str = "diseaseArea",
    diseasedKey: str = "diseased",
) -> pd.DataFrame:
    """Summarize cross-study and disease support for each integrated cluster."""
    required = {clusterKey, sampleKey, studyKey, diseaseAreaKey, diseasedKey}
    missing = required - set(obs.columns)
    if missing:
        raise KeyError(f"obs is missing required columns: {sorted(missing)}")

    frame = obs.copy()
    frame[clusterKey] = as_string_series(frame[clusterKey])
    frame[sampleKey] = as_string_series(frame[sampleKey])
    frame[studyKey] = as_string_series(frame[studyKey])
    frame[diseaseAreaKey] = as_string_series(frame[diseaseAreaKey])
    diseased = frame[diseasedKey].astype("boolean")

    rows: list[dict[str, object]] = []
    for cluster, group in frame.groupby(clusterKey, observed=True):
        n_cells = int(len(group))
        sample_counts = group[sampleKey].value_counts()
        study_counts = group[studyKey].value_counts()
        area_counts = group[diseaseAreaKey].value_counts()
        mapped_areas = area_counts.drop(labels=[OTHER_AREA], errors="ignore")

        diseased_cells = int(diseased.loc[group.index].eq(True).fillna(False).sum())
        nondiseased_cells = int(diseased.loc[group.index].eq(False).fillna(False).sum())
        unknown_cells = n_cells - diseased_cells - nondiseased_cells

        top_study_fraction = float(study_counts.iloc[0] / n_cells) if n_cells else float("nan")
        top_label = None
        top_label_fraction = float("nan")
        n_labels = 0
        if labelKey in group.columns:
            usable = group.loc[is_usable_label(group[labelKey]), labelKey]
            if not usable.empty:
                label_counts = as_string_series(usable).value_counts()
                n_labels = int(label_counts.size)
                top_label = str(label_counts.index[0])
                top_label_fraction = float(label_counts.iloc[0] / label_counts.sum())

        top_ontology = None
        top_ontology_fraction = float("nan")
        n_ontologies = 0
        if ontologyKey in group.columns:
            usable_ont = group.loc[is_usable_label(group[ontologyKey]), ontologyKey]
            if not usable_ont.empty:
                ontology_counts = as_string_series(usable_ont).value_counts()
                n_ontologies = int(ontology_counts.size)
                top_ontology = str(ontology_counts.index[0])
                top_ontology_fraction = float(ontology_counts.iloc[0] / ontology_counts.sum())

        informative_areas = 0
        informative_studies = 0
        for area, area_group in group.groupby(diseaseAreaKey, observed=True):
            if str(area) == OTHER_AREA:
                continue
            area_diseased = diseased.loc[area_group.index]
            case_studies = set(area_group.loc[area_diseased.eq(True).fillna(False), studyKey])
            control_studies = set(group.loc[diseased.eq(False).fillna(False), studyKey])
            overlap = case_studies & control_studies
            if overlap:
                informative_areas += 1
                informative_studies += len(overlap)

        rows.append(
            {
                "cluster": str(cluster),
                "nCells": n_cells,
                "nSamples": int(sample_counts.size),
                "nStudies": int(study_counts.size),
                "nDiseaseAreas": int(mapped_areas.size),
                "diseaseAreas": ",".join(sorted(mapped_areas.index.astype(str))),
                "dominantStudy": str(study_counts.index[0]) if not study_counts.empty else None,
                "dominantStudyFraction": top_study_fraction,
                "nSourceLabels": n_labels,
                "topSourceLabel": top_label,
                "topSourceLabelFraction": top_label_fraction,
                "nOntologies": n_ontologies,
                "topOntology": top_ontology,
                "topOntologyFraction": top_ontology_fraction,
                "nDiseasedCells": diseased_cells,
                "nNondiseasedCells": nondiseased_cells,
                "nUnknownStatusCells": unknown_cells,
                "nInformativeDiseaseAreas": informative_areas,
                "nInformativeStudyOverlaps": informative_studies,
            }
        )
    return pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)


def same_study_contrast_support(
    obs: pd.DataFrame,
    *,
    clusterKey: str = "leiden_atlas",
    sampleKey: str = "SRX_accession",
    studyKey: str = "study_accession",
    diseaseAreaKey: str = "diseaseArea",
    diseasedKey: str = "diseased",
    minCellsPerProfile: int = 10,
    minCaseProfiles: int = 2,
    minControlProfiles: int = 2,
) -> pd.DataFrame:
    """Count same-study case and control pseudobulk profiles per cluster x disease area."""
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

    profile = (
        frame.groupby([clusterKey, studyKey, sampleKey, diseaseAreaKey, diseasedKey], observed=True)
        .size()
        .rename("nCells")
        .reset_index()
    )
    profile = profile[profile["nCells"] >= minCellsPerProfile].copy()

    rows: list[dict[str, object]] = []
    for (cluster, area), group in profile.groupby([clusterKey, diseaseAreaKey], observed=True):
        if str(area) == OTHER_AREA:
            continue
        case_mask = group[diseasedKey].eq(True).fillna(False)
        control_mask = group[diseasedKey].eq(False).fillna(False)
        case_studies = set(group.loc[case_mask, studyKey])
        control_studies = set(group.loc[control_mask, studyKey])
        overlap = sorted(case_studies & control_studies)
        if not overlap:
            rows.append(
                {
                    "cluster": str(cluster),
                    "diseaseArea": str(area),
                    "nOverlapStudies": 0,
                    "nCaseProfiles": 0,
                    "nControlProfiles": 0,
                    "nCaseCells": 0,
                    "nControlCells": 0,
                    "eligibleForContrast": False,
                }
            )
            continue

        overlap_group = group[group[studyKey].isin(overlap)]
        cases = overlap_group[overlap_group[diseasedKey].eq(True).fillna(False)]
        controls = overlap_group[overlap_group[diseasedKey].eq(False).fillna(False)]
        rows.append(
            {
                "cluster": str(cluster),
                "diseaseArea": str(area),
                "nOverlapStudies": len(overlap),
                "nCaseProfiles": int(len(cases)),
                "nControlProfiles": int(len(controls)),
                "nCaseCells": int(cases["nCells"].sum()),
                "nControlCells": int(controls["nCells"].sum()),
                "eligibleForContrast": bool(
                    len(cases) >= minCaseProfiles and len(controls) >= minControlProfiles and len(overlap) >= 1
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["cluster", "diseaseArea"]).reset_index(drop=True)


def classify_cluster_candidates(
    support: pd.DataFrame,
    contrastSupport: pd.DataFrame | None = None,
    *,
    minStudiesShared: int = 5,
    minDiseaseAreasShared: int = 3,
    maxDominantStudyFraction: float = 0.5,
    minLabelPurity: float = 0.5,
    minOverlapStudies: int = 2,
    minCaseProfiles: int = 2,
    minControlProfiles: int = 2,
) -> pd.DataFrame:
    """Flag shared, disease-associated, and unexpected-convergence candidates.

    Eligibility is evidence-based rather than a composite score:
    - shared: many studies and disease areas, not dominated by one study
    - diseaseAssociated: enough same-study case/control support for at least one area
    - unexpectedConvergence: multiple source labels mixed under broad study support
    """
    required = {
        "cluster",
        "nStudies",
        "nDiseaseAreas",
        "dominantStudyFraction",
        "topSourceLabelFraction",
        "nSourceLabels",
    }
    missing = required - set(support.columns)
    if missing:
        raise KeyError(f"support table missing required columns: {sorted(missing)}")

    out = support.copy()
    out["isSharedCandidate"] = (
        (out["nStudies"] >= minStudiesShared)
        & (out["nDiseaseAreas"] >= minDiseaseAreasShared)
        & (out["dominantStudyFraction"] <= maxDominantStudyFraction)
    )

    disease_ok = pd.Series(False, index=out.index)
    if contrastSupport is not None and not contrastSupport.empty:
        eligible = contrastSupport[
            (contrastSupport["nOverlapStudies"] >= minOverlapStudies)
            & (contrastSupport["nCaseProfiles"] >= minCaseProfiles)
            & (contrastSupport["nControlProfiles"] >= minControlProfiles)
            & contrastSupport["eligibleForContrast"].astype(bool)
        ]
        disease_clusters = set(eligible["cluster"].astype(str))
        disease_ok = out["cluster"].astype(str).isin(disease_clusters)
    out["isDiseaseAssociatedCandidate"] = disease_ok.to_numpy()

    purity = out["topSourceLabelFraction"].fillna(0.0)
    out["isUnexpectedConvergenceCandidate"] = (
        out["isSharedCandidate"] & (out["nSourceLabels"].fillna(0).astype(int) >= 2) & (purity < minLabelPurity)
    )
    out["isAnyCandidate"] = (
        out["isSharedCandidate"] | out["isDiseaseAssociatedCandidate"] | out["isUnexpectedConvergenceCandidate"]
    )
    return out.reset_index(drop=True)


def study_balanced_weights(studyLabels: pd.Series) -> pd.Series:
    """Return inverse-study-frequency weights that sum to the number of rows."""
    labels = as_string_series(studyLabels)
    if labels.empty:
        return pd.Series(dtype=np.float64)
    counts = labels.value_counts()
    weights = labels.map(lambda study: 1.0 / float(counts[study]))
    weights = weights.astype(np.float64)
    return weights * (len(weights) / weights.sum())
