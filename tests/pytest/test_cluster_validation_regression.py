"""Golden regression test for cluster_validation against a committed baseline."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest
import scanpy as sc
from _clval_capture import ALL_FIELDS, capture_fields, config_snapshot

SRX = "SRX12708356"
DATA = Path(__file__).parent / "data" / f"{SRX}.h5ad"
BASELINE = Path(__file__).parent / "baselines" / f"clval_{SRX}.json"
CONFIG_SNAPSHOT = Path(__file__).parent / "baselines" / "clval_config_snapshot.json"

_FLOAT_REL_TOL = 1e-6
_FLOAT_ABS_TOL = 1e-9


def _matches(actual: Any, expected: Any) -> bool:
    """Recursively compare values, using pytest.approx for float leaves."""
    if isinstance(expected, bool):
        return actual == expected
    if isinstance(expected, float):
        return actual == pytest.approx(expected, rel=_FLOAT_REL_TOL, abs=_FLOAT_ABS_TOL)
    if isinstance(expected, list):
        return (
            isinstance(actual, list)
            and len(actual) == len(expected)
            and all(_matches(a, e) for a, e in zip(actual, expected, strict=True))
        )
    if isinstance(expected, dict):
        return (
            isinstance(actual, dict)
            and actual.keys() == expected.keys()
            and all(_matches(actual[k], expected[k]) for k in expected)
        )
    return actual == expected


def test_cluster_validation_config_matches_snapshot() -> None:
    """Assert the default ClusterValidationConfig matches the committed model_dump snapshot."""
    from cluster_validation import ClusterValidationConfig

    current = config_snapshot(ClusterValidationConfig())
    if os.environ.get("UPDATE_CLVAL_CONFIG_SNAPSHOT"):
        CONFIG_SNAPSHOT.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
        pytest.skip("snapshot regenerated")
    expected = json.loads(CONFIG_SNAPSHOT.read_text())
    drift = sorted(set(current) | set(expected))
    mismatches = [k for k in drift if current.get(k) != expected.get(k)]
    assert not mismatches, f"ClusterValidationConfig drift: {mismatches}; regen with UPDATE_CLVAL_CONFIG_SNAPSHOT=1"


@pytest.mark.skipif(
    not os.environ.get("RUN_CLVAL_REGRESSION"),
    reason="set RUN_CLVAL_REGRESSION=1 to run (slow clustering regression)",
)
def test_cluster_validation_matches_baseline(tmp_path: Path) -> None:
    """Run cluster validation on the committed fixture and compare to the baseline."""
    from cluster_validation import ClusterValidationConfig, run_cluster_validation_on_adata

    adata = sc.read(str(DATA))
    cfg = ClusterValidationConfig(figsDir=tmp_path / "figs", outputDir=tmp_path / "data")
    _, result = run_cluster_validation_on_adata(
        adata,
        cfg,
        SRX,
        SRX,
        write_outputs=False,
        plot=False,
    )
    captured = capture_fields(result)
    baseline = json.loads(BASELINE.read_text())
    mismatches = [field for field in ALL_FIELDS if not _matches(captured.get(field), baseline.get(field))]
    assert not mismatches, f"drift in fields: {mismatches}"
