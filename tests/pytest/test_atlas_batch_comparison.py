import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.batch_comparison import (
    EXPECTED_BASELINE,
    PCA_KEY,
    SRX_KEY,
    STUDY_KEY,
    STUDY_TECH_KEY,
    TECH_KEY,
    aggregate_scib_results,
    attach_tech_10x,
    ensure_shared_pca,
    harmony_obsm_key,
    load_run_json,
    make_study_tech_batch_key,
    preserve_study_harmony_embedding,
    run_harmony_variants,
    run_scib_evaluation_grid,
    validate_baseline_manifests,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig


def _write_datasets(path: Path, rows: list[dict[str, str]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _tiny_adata() -> sc.AnnData:
    adata = sc.AnnData(
        np.ones((6, 4), dtype=np.float32),
        obs=pd.DataFrame(
            {
                SRX_KEY: ["SRX1", "SRX1", "SRX2", "SRX2", "SRX3", "SRX3"],
                STUDY_KEY: ["PRJ1", "PRJ1", "PRJ1", "PRJ1", "PRJ2", "PRJ2"],
                "cell_type": ["t1", "t1", "t2", "t2", "t1", "t2"],
            },
            index=[f"c{i}" for i in range(6)],
        ),
    )
    adata.obsm[PCA_KEY] = np.random.default_rng(0).normal(size=(6, 5)).astype(np.float32)
    adata.obsm["X_pca_harmony"] = np.random.default_rng(1).normal(size=(6, 5)).astype(np.float32)
    return adata


def test_attach_tech_10x_join(tmp_path: Path) -> None:
    adata = _tiny_adata()
    datasets = _write_datasets(
        tmp_path / "datasets.csv",
        [
            {"srx_accession": "SRX1", "study_accession": "PRJ1", "tech_10x": "3_prime_gex"},
            {"srx_accession": "SRX2", "study_accession": "PRJ1", "tech_10x": "5_prime_gex"},
            {"srx_accession": "SRX3", "study_accession": "PRJ2", "tech_10x": "3_prime_gex"},
        ],
    )
    audit = attach_tech_10x(adata, datasets)
    assert TECH_KEY in adata.obs
    assert adata.obs.loc["c0", TECH_KEY] == "3_prime_gex"
    assert adata.obs.loc["c2", TECH_KEY] == "5_prime_gex"
    assert audit["missingCatalog"] == 0
    assert audit["nStudyTechBatches"] == 3


def test_attach_tech_10x_fails_on_missing_srx(tmp_path: Path) -> None:
    adata = _tiny_adata()
    datasets = _write_datasets(
        tmp_path / "datasets.csv",
        [
            {"srx_accession": "SRX1", "study_accession": "PRJ1", "tech_10x": "3_prime_gex"},
            {"srx_accession": "SRX2", "study_accession": "PRJ1", "tech_10x": "5_prime_gex"},
        ],
    )
    with pytest.raises(ValueError, match="lack a catalog row"):
        attach_tech_10x(adata, datasets)


def test_attach_tech_10x_fails_on_study_mismatch(tmp_path: Path) -> None:
    adata = _tiny_adata()
    datasets = _write_datasets(
        tmp_path / "datasets.csv",
        [
            {"srx_accession": "SRX1", "study_accession": "PRJ9", "tech_10x": "3_prime_gex"},
            {"srx_accession": "SRX2", "study_accession": "PRJ1", "tech_10x": "5_prime_gex"},
            {"srx_accession": "SRX3", "study_accession": "PRJ2", "tech_10x": "3_prime_gex"},
        ],
    )
    with pytest.raises(ValueError, match="study_accession mismatches"):
        attach_tech_10x(adata, datasets)


def test_make_study_tech_batch_key() -> None:
    adata = _tiny_adata()
    adata.obs[TECH_KEY] = ["3_prime_gex", "3_prime_gex", "5_prime_gex", "5_prime_gex", "3_prime_gex", "3_prime_gex"]
    info = make_study_tech_batch_key(adata)
    assert info["column"] == STUDY_TECH_KEY
    assert adata.obs.loc["c0", STUDY_TECH_KEY] == "PRJ1|3_prime_gex"
    assert adata.obs.loc["c2", STUDY_TECH_KEY] == "PRJ1|5_prime_gex"
    assert info["nBatches"] == 3


def test_make_study_tech_batch_key_rejects_separator() -> None:
    adata = _tiny_adata()
    adata.obs[TECH_KEY] = ["3|prime", "3|prime", "5_prime_gex", "5_prime_gex", "3_prime_gex", "3_prime_gex"]
    with pytest.raises(ValueError, match="contain separator"):
        make_study_tech_batch_key(adata)


def test_validate_baseline_manifests_ok() -> None:
    subset = dict(EXPECTED_BASELINE)
    production = dict(EXPECTED_BASELINE)
    validated = validate_baseline_manifests(subsetRun=subset, productionRun=production)
    assert validated["expected"]["nPcs"] == 50


def test_validate_baseline_manifests_mismatch() -> None:
    subset = dict(EXPECTED_BASELINE)
    production = dict(EXPECTED_BASELINE)
    production["nPcs"] = 20
    with pytest.raises(ValueError, match="baselines disagree"):
        validate_baseline_manifests(subsetRun=subset, productionRun=production)


def test_ensure_shared_pca_and_preserve_study_harmony() -> None:
    adata = _tiny_adata()
    ensure_shared_pca(adata, nPcs=5)
    target = preserve_study_harmony_embedding(adata)
    assert target == harmony_obsm_key(STUDY_KEY)
    assert target in adata.obsm
    np.testing.assert_array_equal(adata.obsm[target], adata.obsm["X_pca_harmony"])


def test_run_harmony_variants_distinct_obsm() -> None:
    adata = _tiny_adata()
    adata.obs[TECH_KEY] = ["3_prime_gex"] * 6
    make_study_tech_batch_key(adata)
    cfg = AtlasPostprocessingConfig(nPcs=3, nJobs=1)

    fake = MagicMock()
    fake.Z_corr = np.random.default_rng(2).normal(size=(6, 3)).astype(np.float32)

    with patch("atlas_postprocessing.core.harmonypy.run_harmony", return_value=fake) as run_harmony:
        embedding_map = run_harmony_variants(
            adata,
            cfg,
            batchKeys=[STUDY_TECH_KEY, SRX_KEY],
            nPcs=3,
            skipExisting=False,
        )

    assert run_harmony.call_count == 2
    assert set(embedding_map) == {STUDY_TECH_KEY, SRX_KEY}
    assert embedding_map[STUDY_TECH_KEY] in adata.obsm
    assert embedding_map[SRX_KEY] in adata.obsm
    # Shared PCA must remain untouched.
    assert adata.obsm[PCA_KEY].shape == (6, 5)


def test_run_scib_evaluation_grid(tmp_path: Path) -> None:
    adata = _tiny_adata()
    adata.obs[STUDY_TECH_KEY] = "PRJ1|3_prime_gex"
    embedding_keys = [PCA_KEY, harmony_obsm_key(STUDY_KEY)]

    with patch("atlas_postprocessing.batch_comparison.run_scib_benchmark") as scib:
        scib.return_value = tmp_path / "unused.csv"
        runs = run_scib_evaluation_grid(
            adata,
            outRoot=tmp_path / "scib",
            embeddingKeys=embedding_keys,
            evalBatchKeys=[STUDY_KEY, STUDY_TECH_KEY, SRX_KEY],
            labelKey="cell_type",
            nJobs=2,
            force=True,
        )

    assert scib.call_count == 3
    called_batch_keys = [call.kwargs["batchKey"] for call in scib.call_args_list]
    assert called_batch_keys == [STUDY_KEY, STUDY_TECH_KEY, SRX_KEY]
    for call in scib.call_args_list:
        assert call.kwargs["embeddingKeys"] == embedding_keys
        assert call.kwargs["preIntegratedKey"] == PCA_KEY
    assert len(runs) == 3


def test_aggregate_scib_results(tmp_path: Path) -> None:
    embedding_map = {STUDY_KEY: harmony_obsm_key(STUDY_KEY)}
    run_dir = tmp_path / "eval"
    run_dir.mkdir()
    csv_path = run_dir / "scib_results.csv"
    pd.DataFrame(
        {
            "Total": [0.4, 0.5],
            "Batch correction": [0.2, 0.3],
            "Bio conservation": [0.6, 0.7],
        },
        index=[PCA_KEY, harmony_obsm_key(STUDY_KEY)],
    ).to_csv(csv_path)
    runs = [
        {
            "evalBatchKey": STUDY_KEY,
            "csv": str(csv_path),
            "svg": str(run_dir / "scib_results.svg"),
            "embeddingKeys": [PCA_KEY, harmony_obsm_key(STUDY_KEY)],
        }
    ]
    out_csv = tmp_path / "matrix.csv"
    aggregate_scib_results(runs, embeddingMap=embedding_map, outCsv=out_csv)
    frame = pd.read_csv(out_csv)
    assert set(frame["harmonyBatchKey"]) == {"uncorrected", STUDY_KEY}
    assert (frame["evalBatchKey"] == STUDY_KEY).all()


def test_load_run_json(tmp_path: Path) -> None:
    path = tmp_path / "run.json"
    path.write_text(json.dumps(EXPECTED_BASELINE))
    payload = load_run_json(path)
    assert payload["nPcs"] == 50
