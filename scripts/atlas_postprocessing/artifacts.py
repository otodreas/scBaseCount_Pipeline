import csv
import json
import logging
from pathlib import Path
from typing import Any

from pydantic import ValidationError
from shared.repo import REPO_ROOT, rel_to_repo

from atlas_postprocessing.config import AtlasPostprocessingConfig, AtlasPostprocessingParameters

log = logging.getLogger(__name__)


def write_metric_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    log.info("Wrote %s", rel_to_repo(path))
    return path


def write_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    log.info("Wrote %s", rel_to_repo(path))
    return path


def write_parameters_template(cfg: AtlasPostprocessingConfig, calibrationSummaryPath: Path) -> Path:
    template = AtlasPostprocessingParameters(
        nTopGenes=cfg.nTopGenes,
        nPcs=cfg.nPcs,
        nNeighbors=cfg.nNeighbors,
        resolution=cfg.resolution,
        calibrationSummary=rel_to_repo(calibrationSummaryPath),
    )
    path = cfg.calibrationDir / "parameters_template.json"
    return write_json(path, template.model_dump())


def reject_tuning_overrides_with_parameters_json(overrideFlags: list[str]) -> None:
    """Raise when scalar tuning CLI flags are combined with an approved JSON."""
    if not overrideFlags:
        return
    rendered = ", ".join(f"--{name.replace('_', '-')}" for name in overrideFlags)
    raise ValueError(
        "Do not combine --parameters-json with scalar tuning flags "
        f"({rendered}). The approved JSON is authoritative for those values."
    )


def load_approved_parameters(path: Path) -> AtlasPostprocessingParameters:
    try:
        payload = json.loads(path.read_text())
        return AtlasPostprocessingParameters.model_validate(payload)
    except (OSError, json.JSONDecodeError, ValidationError) as exc:
        raise ValueError(f"Failed to load approved parameters from {rel_to_repo(path)}: {exc}") from exc


def resolve_calibration_summary_path(parameters: AtlasPostprocessingParameters) -> Path | None:
    if parameters.calibrationSummary is None:
        return None
    summary_path = Path(parameters.calibrationSummary)
    if not summary_path.is_absolute():
        summary_path = REPO_ROOT / summary_path
    return summary_path


def validate_approved_against_calibration(
    parameters: AtlasPostprocessingParameters,
    *,
    parametersPath: Path,
) -> dict[str, Any]:
    """Confirm approved values were evaluated during calibration when a summary is referenced."""
    summary_path = resolve_calibration_summary_path(parameters)
    if summary_path is None:
        raise ValueError(
            f"{rel_to_repo(parametersPath)} must include calibrationSummary so approved values "
            "can be checked against evaluated candidates"
        )
    if not summary_path.is_file():
        raise FileNotFoundError(f"Calibration summary not found: {rel_to_repo(summary_path)}")

    summary = json.loads(summary_path.read_text())
    candidates = summary.get("candidates") or {}
    checks = {
        "nTopGenes": ("hvg", parameters.nTopGenes),
        "nPcs": ("pc", parameters.nPcs),
        "nNeighbors": ("neighbors", parameters.nNeighbors),
        "resolution": ("resolution", parameters.resolution),
    }
    for field, (candidate_key, value) in checks.items():
        evaluated = candidates.get(candidate_key)
        if evaluated is None:
            raise ValueError(f"Calibration summary missing candidates.{candidate_key}")
        if not _value_in_candidates(value, evaluated):
            raise ValueError(
                f"Approved {field}={value!r} was not among calibration candidates {evaluated!r} "
                f"in {rel_to_repo(summary_path)}"
            )
    return summary


def _value_in_candidates(value: Any, evaluated: list[Any]) -> bool:
    if isinstance(value, float):
        return any(abs(float(item) - value) < 1e-9 for item in evaluated)
    return value in evaluated


def apply_parameters_to_config(
    cfg: AtlasPostprocessingConfig,
    parameters: AtlasPostprocessingParameters,
    *,
    parametersPath: Path,
) -> AtlasPostprocessingConfig:
    n_pcs_compute = max(cfg.nPcsCompute, parameters.nPcs)
    return cfg.model_copy(
        update={
            "nTopGenes": parameters.nTopGenes,
            "nPcs": parameters.nPcs,
            "nNeighbors": parameters.nNeighbors,
            "resolution": parameters.resolution,
            "nPcsCompute": n_pcs_compute,
            "parametersJson": parametersPath,
        }
    )
