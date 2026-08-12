"""Tests for full-atlas aggregation, memory preflight, and adaptive ranking."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from disease_markers.aggregation import (
    aggregate_sparse_pseudobulk,
    build_aggregate_fingerprint,
    reject_dense_cell_matrix,
    write_checkpoint,
)
from disease_markers.config import AtlasDeAnalysisConfig
from disease_markers.memory import assert_memory_available, estimate_raw_sparse_bytes
from disease_markers.ranking import (
    build_evidence_pools,
    empirical_percentile_scores,
    merge_duplicate_evidence,
    select_review_queue,
)
from disease_markers.specificity import gene_specificity_table, tau_specificity
from scipy import sparse


def _toy_cell_adata() -> AnnData:
    # 6 cells, 3 genes, 2 samples x 2 clusters with overlapping membership.
    x = sparse.csr_matrix(
        np.array(
            [
                [10, 0, 1],
                [12, 0, 0],
                [0, 8, 2],
                [0, 9, 0],
                [5, 0, 3],
                [7, 1, 0],
            ],
            dtype=np.float64,
        )
    )
    obs = pd.DataFrame(
        {
            "SRX_accession": ["SRX1", "SRX1", "SRX1", "SRX2", "SRX2", "SRX2"],
            "leiden_atlas": ["0", "0", "1", "0", "1", "1"],
            "study_accession": ["PRJ_A"] * 6,
            "diseaseArea": ["COPD"] * 6,
            "diseased": [True, True, True, False, False, False],
            "cell_type": ["macrophage", "macrophage", "epithelial", "macrophage", "epithelial", "epithelial"],
            "cell_ontology_term_id": ["CL:1", "CL:1", "CL:2", "CL:1", "CL:2", "CL:2"],
            "eligible": [True] * 6,
        },
        index=[f"c{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["G1", "G2", "G3"])
    return AnnData(X=x, obs=obs, var=var)


def test_aggregate_sparse_pseudobulk_sums_and_props() -> None:
    adata = _toy_cell_adata()
    pdata = aggregate_sparse_pseudobulk(
        adata.X,
        adata.obs,
        adata.var,
        sampleKey="SRX_accession",
        clusterKey="leiden_atlas",
    )
    assert pdata.n_obs == 4
    assert "psbulk_props" in pdata.layers
    assert "psbulk_cells" in pdata.obs.columns
    # SRX1 cluster 0 has cells [10,0,1] and [12,0,0]
    row = pdata.obs_names.get_loc("SRX1_0")
    assert pdata.obs.iloc[row]["psbulk_cells"] == 2
    assert float(pdata.X[row, 0]) == 22.0
    assert float(pdata.layers["psbulk_props"][row, 0]) == 1.0
    assert float(pdata.layers["psbulk_props"][row, 1]) == 0.0


def test_reject_dense_cell_matrix() -> None:
    with pytest.raises(TypeError, match="dense"):
        reject_dense_cell_matrix(np.ones((3, 3)), label="test")
    reject_dense_cell_matrix(sparse.csr_matrix(np.ones((3, 3))), label="test")


def test_tau_and_gene_specificity_chunked() -> None:
    assert tau_specificity(np.array([10.0, 0.0, 0.0])) == pytest.approx(1.0)
    counts = np.array(
        [
            [100, 0, 5],
            [120, 0, 4],
            [0, 80, 3],
            [0, 90, 2],
        ],
        dtype=float,
    )
    props = (counts > 0).astype(float)
    pdata = AnnData(
        X=counts,
        obs=pd.DataFrame(
            {
                "leiden_atlas": ["0", "0", "1", "1"],
                "study_accession": ["PRJ_A", "PRJ_B", "PRJ_A", "PRJ_B"],
                "diseased": [True, False, True, False],
            },
            index=[f"p{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
        layers={"psbulk_props": props},
    )
    table = gene_specificity_table(pdata, geneChunkSize=2, minProfilesForGene=1, minTotalCounts=1, minStudies=1)
    assert not table.empty
    g1 = table.set_index("gene").loc["G1"]
    assert g1["topCluster"] == "0"
    assert float(g1["tau"]) > 0.5


def test_empirical_percentile_and_review_budget() -> None:
    values = pd.Series([1.0, 2.0, 3.0, 4.0])
    scores = empirical_percentile_scores(values)
    assert scores.iloc[0] == pytest.approx(0.0)
    assert scores.iloc[-1] == pytest.approx(1.0)

    restricted = pd.DataFrame(
        {
            "gene": [f"g{i}" for i in range(10)],
            "geneSymbol": [f"G{i}" for i in range(10)],
            "topCluster": [str(i % 3) for i in range(10)],
            "tau": np.linspace(0.85, 0.99, 10),
            "meanDetectionTop": np.linspace(0.3, 0.8, 10),
            "maxDetectionBackground": [0.01] * 10,
            "nStudiesAgreeTop": [5] * 10,
            "nStudiesScored": [5] * 10,
            "detectionDifference": np.linspace(0.2, 0.7, 10),
            "interpretationStatus": ["resolved"] * 10,
            "interpretedCellType": ["macrophage"] * 10,
        }
    )
    de_hits = pd.DataFrame(
        {
            "gene": [f"d{i}" for i in range(12)],
            "geneSymbol": [f"D{i}" for i in range(12)],
            "cluster": [str(i % 4) for i in range(12)],
            "diseaseArea": ["IPF / Pulmonary Fibrosis" if i % 2 == 0 else "COPD" for i in range(12)],
            "log2FoldChange": [2.0 if i % 2 == 0 else -2.0 for i in range(12)],
            "padj": np.linspace(1e-6, 1e-2, 12),
            "detectionDelta": [0.3 if i % 2 == 0 else -0.3 for i in range(12)],
            "nStudies": [3] * 12,
            "interpretationStatus": ["resolved"] * 12,
            "interpretedCellType": ["macrophage"] * 12,
        }
    )
    pools = build_evidence_pools(
        restricted=restricted,
        deHits=de_hits,
        sharedGenes=pd.DataFrame(),
        geneClass=pd.DataFrame(),
        unexpected=pd.DataFrame(),
        padj=0.05,
        lfc=1.0,
        minDetectionDelta=0.15,
        minTau=0.8,
        minTargetDetection=0.2,
        maxBackgroundDetection=0.05,
        minStudiesForSpecificity=3,
    )
    primary, extended, thresholds = select_review_queue(
        pools,
        primaryBudget=5,
        extendedBudget=12,
        maxPerClassPrimary=3,
        maxPerClassExtended=6,
        maxPerGene=1,
        maxPerCluster=2,
        maxPerDiseaseArea=4,
    )
    assert len(primary) <= 5
    assert len(primary) + len(extended) <= 12
    assert not thresholds.empty
    assert primary["reviewTier"].eq("primary").all()


def test_merge_duplicate_evidence_keeps_classes() -> None:
    frame = pd.DataFrame(
        [
            {
                "gene": "g1",
                "geneSymbol": "G1",
                "cluster": "0",
                "diseaseArea": "COPD",
                "proposalClass": "replicatedDiseaseGain",
                "evidenceScore": 0.9,
            },
            {
                "gene": "g1",
                "geneSymbol": "G1",
                "cluster": "0",
                "diseaseArea": "COPD",
                "proposalClass": "unexpectedExpression",
                "evidenceScore": 0.8,
            },
        ]
    )
    merged = merge_duplicate_evidence(frame)
    assert len(merged) == 1
    assert "unexpectedExpression" in merged.iloc[0]["allEvidenceClasses"]
    assert float(merged.iloc[0]["evidenceScore"]) == pytest.approx(0.9)


def test_checkpoint_fingerprint_roundtrip(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas.h5ad"
    contexts = tmp_path / "contexts.jsonl"
    atlas_manifest = tmp_path / "atlas_result.json"
    gene_info = tmp_path / "geneInfo.tab"
    sample_meta = tmp_path / "samples.parquet"
    for path in (atlas, contexts, gene_info):
        path.write_text("x")
    atlas_manifest.write_text(
        json.dumps(
            {
                "files": [
                    {
                        "accession": "SRX1",
                        "studyAccession": "STUDY1",
                        "r2Key": "prefix/SRX1.h5ad",
                        "status": "success",
                        "skipReason": None,
                        "qc": None,
                    }
                ]
            }
        )
    )
    pd.DataFrame({"srxAccession": ["SRX1"]}).to_parquet(sample_meta)

    cfg = AtlasDeAnalysisConfig(
        atlasPath=atlas,
        outputDir=tmp_path / "out",
        contextsPath=contexts,
        atlasManifestPath=atlas_manifest,
        geneInfoPath=gene_info,
    )
    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    fingerprint = build_aggregate_fingerprint(cfg, sampleMetadataPath=sample_meta)
    pdata = _toy_cell_adata()
    pdata.layers["psbulk_props"] = (pdata.X.toarray() > 0).astype(float)
    pdata.obs["psbulk_cells"] = 2
    write_checkpoint(pdata, fingerprint, cfg)
    stored = json.loads(cfg.fingerprintPath.read_text())
    assert stored["sha256"] == fingerprint["sha256"]

    cfg_changed = cfg.model_copy(update={"minCellsPerProfile": 99})
    changed = build_aggregate_fingerprint(cfg_changed, sampleMetadataPath=sample_meta)
    assert changed["sha256"] != fingerprint["sha256"]


def test_memory_estimate_and_low_memory_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    x = sparse.csr_matrix(np.array([[1, 0], [0, 2], [3, 0]], dtype=np.float64))
    adata = AnnData(X=x, obs=pd.DataFrame(index=["a", "b", "c"]), var=pd.DataFrame(index=["g1", "g2"]))
    path = tmp_path / "toy.h5ad"
    adata.write_h5ad(path)
    estimate = estimate_raw_sparse_bytes(path, overheadFactor=1.0)
    assert estimate.nObs == 3
    assert estimate.nVars == 2
    assert estimate.nnz == 3
    assert estimate.estimatedBytes > 0

    monkeypatch.setattr("disease_markers.memory.available_ram_bytes", lambda: 1)
    with pytest.raises(MemoryError):
        assert_memory_available(estimate, reserveBytes=0)


def test_select_review_queue_does_not_pad_with_invalid_rows() -> None:
    pools = {
        "clusterRestricted": pd.DataFrame(),
        "replicatedDiseaseGain": pd.DataFrame(
            {
                "gene": ["g1"],
                "geneSymbol": ["G1"],
                "cluster": ["0"],
                "diseaseArea": ["COPD"],
                "proposalClass": ["replicatedDiseaseGain"],
                "evidenceScore": [0.9],
                "interpretationStatus": ["resolved"],
                "interpretedCellType": ["macrophage"],
            }
        ),
        "replicatedDiseaseDepletion": pd.DataFrame(),
        "sharedDiseaseProgram": pd.DataFrame(),
        "oppositeDiseaseEffect": pd.DataFrame(),
        "unexpectedExpression": pd.DataFrame(),
    }
    primary, extended, _thresholds = select_review_queue(
        pools,
        primaryBudget=20,
        extendedBudget=60,
    )
    assert len(primary) == 1
    assert extended.empty
