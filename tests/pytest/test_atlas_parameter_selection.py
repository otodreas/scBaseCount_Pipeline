import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from atlas_postprocessing.artifacts import (
    apply_parameters_to_config,
    load_approved_parameters,
    validate_approved_against_calibration,
    write_json,
    write_parameters_template,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.selection import FIXED_N_NEIGHBORS, FIXED_N_TOP_GENES, run_calibration
from cluster_validation import default_resolutions, select_resolution_on_graph
from cluster_validation.clustering import sweep_leiden_resolutions
from cluster_validation.resolution import score_resolutions


def test_parameters_template_seeds_recommendation_values(tmp_path: Path) -> None:
    calibration_dir = tmp_path / "parameter_selection"
    summary_path = calibration_dir / "calibration_summary.json"
    summary = {
        "baseline": {"nTopGenes": 2000, "nPcs": 18, "nNeighbors": 15, "resolution": 0.7},
        "candidates": {
            "hvg": [2000],
            "pc": [18],
            "neighbors": [15],
            "resolution": [0.5, 0.7, 0.9],
        },
        "recommendation": {
            "nTopGenes": 2000,
            "nPcs": 18,
            "nNeighbors": 15,
            "resolution": 0.7,
            "method": "matched_jaccard_argmax",
        },
    }
    write_json(summary_path, summary)

    cfg = AtlasPostprocessingConfig(
        calibrationDir=calibration_dir,
        nTopGenes=2000,
        nPcs=18,
        nNeighbors=15,
        resolution=0.7,
    )
    template_path = write_parameters_template(cfg, summary_path)
    assert template_path.is_file()
    assert not (calibration_dir / "approved_parameters.json").exists()

    approved_path = calibration_dir / "approved_parameters.json"
    payload = json.loads(template_path.read_text())
    payload["resolution"] = 0.5  # manual override within evaluated grid
    payload["calibrationSummary"] = str(summary_path)
    approved_path.write_text(json.dumps(payload, indent=2))

    loaded = load_approved_parameters(approved_path)
    validated = validate_approved_against_calibration(loaded, parametersPath=approved_path)
    assert validated["candidates"]["hvg"] == [2000]
    assert loaded.resolution == 0.5

    applied = apply_parameters_to_config(cfg, loaded, parametersPath=approved_path)
    assert applied.nTopGenes == FIXED_N_TOP_GENES
    assert applied.nNeighbors == FIXED_N_NEIGHBORS
    assert applied.nPcs == 18
    assert applied.resolution == 0.5
    assert applied.parametersJson == approved_path


def test_approved_graph_value_must_match_singleton_candidate(tmp_path: Path) -> None:
    summary_path = tmp_path / "calibration_summary.json"
    write_json(
        summary_path,
        {
            "candidates": {
                "hvg": [2000],
                "pc": [18],
                "neighbors": [15],
                "resolution": [0.5, 0.7],
            }
        },
    )
    approved_path = tmp_path / "approved_parameters.json"
    write_json(
        approved_path,
        {
            "nTopGenes": 4000,
            "nPcs": 18,
            "nNeighbors": 15,
            "resolution": 0.5,
            "calibrationSummary": str(summary_path),
        },
    )
    parameters = load_approved_parameters(approved_path)
    with pytest.raises(ValueError, match="nTopGenes"):
        validate_approved_against_calibration(parameters, parametersPath=approved_path)


def test_approved_resolution_outside_grid_is_rejected(tmp_path: Path) -> None:
    summary_path = tmp_path / "calibration_summary.json"
    write_json(
        summary_path,
        {
            "candidates": {
                "hvg": [2000],
                "pc": [18],
                "neighbors": [15],
                "resolution": [0.5, 0.7],
            }
        },
    )
    approved_path = tmp_path / "approved_parameters.json"
    write_json(
        approved_path,
        {
            "nTopGenes": 2000,
            "nPcs": 18,
            "nNeighbors": 15,
            "resolution": 1.3,
            "calibrationSummary": str(summary_path),
        },
    )
    parameters = load_approved_parameters(approved_path)
    with pytest.raises(ValueError, match="resolution"):
        validate_approved_against_calibration(parameters, parametersPath=approved_path)


def _toy_harmony_ready_adata(n_cells: int = 60) -> sc.AnnData:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(n_cells, 12)).astype(np.float32)
    x[: n_cells // 2] += 4.0
    adata = sc.AnnData(x)
    adata.obs = pd.DataFrame(
        {
            "study_accession": ["S1"] * (n_cells // 2) + ["S2"] * (n_cells - n_cells // 2),
            "cell_type": ["A"] * (n_cells // 2) + ["B"] * (n_cells - n_cells // 2),
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    sc.pp.neighbors(adata, n_neighbors=5, use_rep="X")
    return adata


def test_atlas_adapter_matches_shared_resolution_selection() -> None:
    adata = _toy_harmony_ready_adata()
    resolutions = [0.2, 0.4, 0.6]

    shared_adata = adata.copy()
    _, shared_sel = select_resolution_on_graph(
        shared_adata,
        resolutions=resolutions,
        weakPriorKey="cell_type",
    )

    manual_adata = adata.copy()
    manual_adata, n_clusters = sweep_leiden_resolutions(manual_adata, resolutions)
    manual_sel = score_resolutions(
        manual_adata,
        resolutions=resolutions,
        weakPriorKey="cell_type",
        nClusters=n_clusters,
    )

    assert shared_sel.bestIdx == manual_sel.bestIdx
    assert shared_sel.selectedResolution == manual_sel.selectedResolution
    assert np.allclose(shared_sel.jaccArr, manual_sel.jaccArr)
    assert np.array_equal(shared_sel.kArr, manual_sel.kArr)


def test_run_calibration_writes_singleton_candidates_and_no_approved_file(tmp_path: Path) -> None:
    n_cells = 40
    rng = np.random.default_rng(3)
    x = np.abs(rng.normal(size=(n_cells, 30))).astype(np.float32)
    x[:20] += 3.0
    adata = sc.AnnData(x)
    adata.obs = pd.DataFrame(
        {
            "study_accession": ["S1"] * 20 + ["S2"] * 20,
            "cell_type": ["A"] * 20 + ["B"] * 20,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    adata.var_names = [f"g{i}" for i in range(30)]

    cfg = AtlasPostprocessingConfig(
        inputH5ad=tmp_path / "atlas.h5ad",
        calibrationDir=tmp_path / "parameter_selection",
        nTopGenes=FIXED_N_TOP_GENES,
        nNeighbors=FIXED_N_NEIGHBORS,
        nPcsCompute=10,
        nPcsMin=5,
        nPcsCumvarTarget=0.3,
        resolutionCandidates=[0.2, 0.4, 0.6],
        writePlots=False,
        nJobs=1,
    )

    fake_harmony = type("H", (), {})()
    fake_harmony.Z_corr = rng.normal(size=(n_cells, 8)).astype(np.float32)

    with patch("atlas_postprocessing.selection.run_harmony_on_pcs", return_value=fake_harmony.Z_corr):
        # Force HVG count down so the tiny matrix works.
        with patch("atlas_postprocessing.selection.select_hvgs", side_effect=lambda ad, _cfg, nTopGenes=None: ad):
            with patch(
                "atlas_postprocessing.selection.scale_and_pca",
                side_effect=lambda ad, _cfg, nPcsCompute=None: _add_pca(ad, nPcsCompute or 10),
            ):
                with patch(
                    "atlas_postprocessing.selection.build_neighbors",
                    side_effect=lambda ad, **kwargs: sc.pp.neighbors(ad, n_neighbors=5, use_rep="X"),
                ):
                    with patch(
                        "atlas_postprocessing.selection.pick_n_pcs",
                        return_value=(8, 55.0),
                    ):
                        summary = run_calibration(cfg, adata=adata)

    assert summary["candidates"]["hvg"] == [2000]
    assert summary["candidates"]["neighbors"] == [15]
    assert summary["candidates"]["pc"] == [8]
    assert summary["candidates"]["resolution"] == [0.2, 0.4, 0.6]
    assert summary["recommendation"]["resolution"] == summary["baseline"]["resolution"]
    assert (cfg.calibrationDir / "parameters_template.json").is_file()
    assert not (cfg.calibrationDir / "approved_parameters.json").exists()
    assert (cfg.calibrationDir / "metrics" / "resolution.csv").is_file()


def _add_pca(adata: sc.AnnData, n_comps: int) -> sc.AnnData:
    sc.pp.pca(adata, n_comps=min(n_comps, adata.n_obs - 1, adata.n_vars - 1), svd_solver="arpack")
    return adata


def test_default_resolution_candidates_match_cluster_validation() -> None:
    assert AtlasPostprocessingConfig().resolutionCandidates == default_resolutions()
