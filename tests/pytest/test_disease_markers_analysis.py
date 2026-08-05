import anndata as ad
import numpy as np
import pandas as pd
import pytest
from disease_markers.candidates import (
    annotate_obs_with_labels,
    classify_cluster_candidates,
    cluster_support_table,
    same_study_contrast_support,
    study_balanced_weights,
)
from disease_markers.concordance import (
    cluster_label_purity,
    concordance_summary,
    source_label_contingency,
)
from disease_markers.validation import (
    differential_abundance_by_study,
    filter_pseudobulk_profiles,
    filter_two_sided_de,
    same_study_case_control_profiles,
    sample_cluster_proportions,
    shared_direction_genes,
)


def _toy_obs() -> pd.DataFrame:
    rows = []
    # Study A: IPF cases and controls, cluster 0 dominated by alveolar cells
    for i in range(20):
        rows.append(
            {
                "SRX_accession": "SRX_A_CASE" if i < 10 else "SRX_A_CTRL",
                "study_accession": "PRJ_A",
                "leiden_atlas": "0",
                "cell_type": "alveolar",
                "cell_ontology_term_id": "CL:0002063",
                "diseaseArea": "IPF / Pulmonary Fibrosis",
                "diseased": True if i < 10 else False,
                "isBiologicalControl": False if i < 10 else True,
                "controlType": None if i < 10 else "healthy",
                "eligible": True,
                "excludeReason": None,
            }
        )
    # Study B: COPD cases and controls, same cluster, mixed source labels
    for i in range(20):
        rows.append(
            {
                "SRX_accession": "SRX_B_CASE" if i < 10 else "SRX_B_CTRL",
                "study_accession": "PRJ_B",
                "leiden_atlas": "0",
                "cell_type": "macrophage" if i % 2 == 0 else "alveolar",
                "cell_ontology_term_id": "CL:0000235" if i % 2 == 0 else "CL:0002063",
                "diseaseArea": "COPD",
                "diseased": True if i < 10 else False,
                "isBiologicalControl": False if i < 10 else True,
                "controlType": None if i < 10 else "explicitControl",
                "eligible": True,
                "excludeReason": None,
            }
        )
    # Study C only: rare cluster dominated by one study
    for _i in range(12):
        rows.append(
            {
                "SRX_accession": "SRX_C_ONLY",
                "study_accession": "PRJ_C",
                "leiden_atlas": "1",
                "cell_type": "UNKNOWN",
                "cell_ontology_term_id": "",
                "diseaseArea": "Lung Cancer",
                "diseased": True,
                "isBiologicalControl": False,
                "controlType": None,
                "eligible": True,
                "excludeReason": None,
            }
        )
    return pd.DataFrame(rows)


def test_source_label_contingency_ignores_unknown() -> None:
    obs = _toy_obs()
    table = source_label_contingency(obs, clusterKey="leiden_atlas", labelKey="cell_type")
    assert "0" in table.index
    assert "1" not in table.index
    assert set(table.columns) == {"alveolar", "macrophage"}


def test_concordance_summary_has_global_and_macro() -> None:
    obs = _toy_obs()
    summary = concordance_summary(
        obs,
        clusterKey="leiden_atlas",
        labelKey="cell_type",
        studyKey="study_accession",
    )
    scopes = set(summary["scope"])
    assert {"global", "study", "studyMacro"} <= scopes
    purity = cluster_label_purity(
        obs,
        clusterKey="leiden_atlas",
        labelKey="cell_type",
        studyKey="study_accession",
    )
    assert float(purity.loc[purity["cluster"] == "0", "topLabelFraction"].iloc[0]) > 0.5


def test_study_balanced_weights_equalize_studies() -> None:
    labels = pd.Series(["PRJ_A"] * 9 + ["PRJ_B"] * 1)
    weights = study_balanced_weights(labels)
    assert pytest.approx(weights.sum(), rel=1e-6) == len(labels)
    assert pytest.approx(weights[labels == "PRJ_A"].sum(), rel=1e-6) == pytest.approx(
        weights[labels == "PRJ_B"].sum(), rel=1e-6
    )


def test_cluster_support_and_candidate_classification() -> None:
    obs = _toy_obs()
    support = cluster_support_table(obs)
    contrast = same_study_contrast_support(
        obs,
        minCellsPerProfile=5,
        minCaseProfiles=1,
        minControlProfiles=1,
    )
    candidates = classify_cluster_candidates(
        support,
        contrast,
        minStudiesShared=2,
        minDiseaseAreasShared=2,
        maxDominantStudyFraction=0.7,
        minLabelPurity=0.8,
        minOverlapStudies=1,
        minCaseProfiles=1,
        minControlProfiles=1,
    )
    row0 = candidates.set_index("cluster").loc["0"]
    row1 = candidates.set_index("cluster").loc["1"]
    assert bool(row0["isSharedCandidate"]) is True
    assert bool(row0["isDiseaseAssociatedCandidate"]) is True
    assert bool(row1["isSharedCandidate"]) is False
    assert (
        int(contrast.loc[(contrast.cluster == "0") & (contrast.diseaseArea == "COPD"), "nOverlapStudies"].iloc[0]) == 1
    )


def test_annotate_obs_with_labels_filters_ineligible() -> None:
    obs = pd.DataFrame(
        {
            "SRX_accession": ["SRX1", "SRX2"],
            "study_accession": ["PRJ1", "PRJ1"],
            "leiden_atlas": ["0", "0"],
            "cell_type": ["alveolar", "alveolar"],
        }
    )
    labels = pd.DataFrame(
        [
            {
                "srxAccession": "SRX1",
                "studyAccession": "PRJ1",
                "diseaseArea": "COPD",
                "diseased": True,
                "isBiologicalControl": False,
                "controlType": None,
                "eligible": True,
                "excludeReason": None,
            },
            {
                "srxAccession": "SRX2",
                "studyAccession": "PRJ1",
                "diseaseArea": "Other",
                "diseased": None,
                "isBiologicalControl": False,
                "controlType": None,
                "eligible": False,
                "excludeReason": "unmapped_disease",
            },
        ]
    )
    annotated = annotate_obs_with_labels(obs, labels)
    assert list(annotated["SRX_accession"]) == ["SRX1"]


def test_same_study_case_control_profiles_excludes_foreign_controls() -> None:
    pdata = ad.AnnData(np.ones((4, 2)))
    pdata.obs = pd.DataFrame(
        {
            "leiden_atlas": ["0", "0", "0", "0"],
            "study_accession": ["PRJ_A", "PRJ_A", "PRJ_B", "PRJ_C"],
            "diseaseArea": ["COPD", "COPD", "Other", "COPD"],
            "diseased": [True, False, False, True],
        },
        index=[f"pb{i}" for i in range(4)],
    )
    selected = same_study_case_control_profiles(pdata, area="COPD", cluster="0")
    # PRJ_A has both a COPD case and a control. PRJ_B has only a control.
    # PRJ_C has only a COPD case. Only PRJ_A profiles remain.
    assert selected.tolist() == [True, True, False, False]


def test_filter_pseudobulk_profiles_and_two_sided_de() -> None:
    pdata = ad.AnnData(np.ones((3, 2)))
    pdata.obs = pd.DataFrame({"psbulk_cells": [20, 5, 15]}, index=["a", "b", "c"])
    filtered = filter_pseudobulk_profiles(pdata, minCellsPerProfile=10)
    assert list(filtered.obs_names) == ["a", "c"]

    results = pd.DataFrame(
        {
            "gene": ["g1", "g2", "g3"],
            "padj": [0.01, 0.01, 0.2],
            "log2FoldChange": [1.5, -1.2, 2.0],
        }
    )
    hits = filter_two_sided_de(results, padj=0.05, lfc=1.0)
    assert set(hits["gene"]) == {"g1", "g2"}


def test_differential_abundance_and_shared_direction() -> None:
    obs = _toy_obs()
    proportions = sample_cluster_proportions(obs)
    abundance = differential_abundance_by_study(
        proportions,
        area="IPF / Pulmonary Fibrosis",
        cluster="0",
        minSamplesPerArm=1,
    )
    assert not abundance.empty
    assert set(abundance["studyAccession"]) == {"PRJ_A"}

    hits = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2", "GENE2"],
            "diseaseArea": ["COPD", "IPF / Pulmonary Fibrosis", "COPD", "IPF / Pulmonary Fibrosis"],
            "cluster": ["0", "0", "0", "0"],
            "log2FoldChange": [1.2, 1.5, 1.1, -1.3],
        }
    )
    shared = shared_direction_genes(hits, minDiseaseAreas=2)
    assert list(shared["gene"]) == ["GENE1"]
    assert shared.iloc[0]["direction"] == "up"
