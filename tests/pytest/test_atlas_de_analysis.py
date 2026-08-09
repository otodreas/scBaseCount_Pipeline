"""Tests for full-atlas aggregation, memory preflight, and specificity helpers."""

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


def test_checkpoint_fingerprint_roundtrip(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas.h5ad"
    contexts = tmp_path / "contexts.jsonl"
    atlas_csv = tmp_path / "atlas.csv"
    gene_info = tmp_path / "geneInfo.tab"
    sample_meta = tmp_path / "samples.parquet"
    for path in (atlas, contexts, atlas_csv, gene_info):
        path.write_text("x")
    pd.DataFrame({"srxAccession": ["SRX1"]}).to_parquet(sample_meta)

    cfg = AtlasDeAnalysisConfig(
        atlasPath=atlas,
        outputDir=tmp_path / "out",
        contextsPath=contexts,
        atlasCsvPath=atlas_csv,
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
