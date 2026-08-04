import json
from pathlib import Path

import pytest
from atlas_postprocessing.artifacts import (
    apply_parameters_to_config,
    load_approved_parameters,
    validate_approved_against_calibration,
    write_json,
    write_parameters_template,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.selection import validate_candidate_lists


def test_validate_candidate_lists_requires_two_values() -> None:
    cfg = AtlasPostprocessingConfig(hvgCandidates=[2000])
    with pytest.raises(ValueError, match="hvgCandidates"):
        validate_candidate_lists(cfg)


def test_parameters_template_and_approved_reuse(tmp_path: Path) -> None:
    calibration_dir = tmp_path / "parameter_selection"
    summary_path = calibration_dir / "calibration_summary.json"
    summary = {
        "baseline": {"nTopGenes": 2000, "nPcs": 20, "nNeighbors": 15, "resolution": 1.0},
        "candidates": {
            "hvg": [1000, 2000],
            "pc": [10, 20],
            "neighbors": [5, 15],
            "resolution": [0.2, 1.0],
        },
    }
    write_json(summary_path, summary)

    cfg = AtlasPostprocessingConfig(
        calibrationDir=calibration_dir,
        nTopGenes=2000,
        nPcs=20,
        nNeighbors=15,
        resolution=1.0,
    )
    template_path = write_parameters_template(cfg, summary_path)
    assert template_path.is_file()

    approved_path = calibration_dir / "approved_parameters.json"
    payload = json.loads(template_path.read_text())
    payload["calibrationSummary"] = str(summary_path)
    approved_path.write_text(json.dumps(payload, indent=2))

    loaded = load_approved_parameters(approved_path)
    validated = validate_approved_against_calibration(loaded, parametersPath=approved_path)
    assert validated["candidates"]["hvg"] == [1000, 2000]

    applied = apply_parameters_to_config(cfg, loaded, parametersPath=approved_path)
    assert applied.nTopGenes == 2000
    assert applied.nNeighbors == 15
    assert applied.parametersJson == approved_path


def test_approved_value_must_be_evaluated_candidate(tmp_path: Path) -> None:
    summary_path = tmp_path / "calibration_summary.json"
    write_json(
        summary_path,
        {
            "candidates": {
                "hvg": [1000, 2000],
                "pc": [10, 20],
                "neighbors": [5, 15],
                "resolution": [0.2, 1.0],
            }
        },
    )
    approved_path = tmp_path / "approved_parameters.json"
    write_json(
        approved_path,
        {
            "nTopGenes": 8000,
            "nPcs": 20,
            "nNeighbors": 15,
            "resolution": 1.0,
            "calibrationSummary": str(summary_path),
        },
    )
    parameters = load_approved_parameters(approved_path)
    with pytest.raises(ValueError, match="nTopGenes"):
        validate_approved_against_calibration(parameters, parametersPath=approved_path)
