import json
import logging
import time
from typing import Any

import scanpy as sc
from shared.repo import rel_to_repo

from atlas_postprocessing.artifacts import (
    load_approved_parameters,
    validate_approved_against_calibration,
    write_json,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import run_postprocessing, timed
from atlas_postprocessing.sampling import sample_metadata
from atlas_postprocessing.scib import run_scib_benchmark

log = logging.getLogger(__name__)


def run_validation(
    cfg: AtlasPostprocessingConfig,
    adata: sc.AnnData | None = None,
    *,
    scibJobs: int = 6,
    forceScib: bool = False,
) -> dict[str, Any]:
    """Run approved parameters on a sample and write validation artifacts."""
    if cfg.parametersJson is None:
        raise ValueError("Validation requires cfg.parametersJson")

    cfg.validationDir.mkdir(parents=True, exist_ok=True)
    scib_dir = cfg.validationDir / "scib"

    started = time.perf_counter()
    validated = timed(
        "approved subset postprocessing",
        lambda: run_postprocessing(cfg, adata=adata, workflow="validation"),
        logger=log,
    )
    # RF merge of leiden_atlas could run here before scIB.
    timed(
        "scIB benchmark",
        lambda: run_scib_benchmark(
            validated,
            outDir=scib_dir,
            batchKey=cfg.batchKey,
            labelKey=cfg.cellTypeKey,
            nJobs=scibJobs,
            force=forceScib,
        ),
        logger=log,
    )

    parameters = load_approved_parameters(cfg.parametersJson)
    summary = validate_approved_against_calibration(parameters, parametersPath=cfg.parametersJson)
    recommendation = summary.get("recommendation") or {}
    validation_summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "outputDir": rel_to_repo(cfg.validationDir),
        "parametersJson": rel_to_repo(cfg.parametersJson),
        "calibrationSummary": parameters.calibrationSummary,
        "resolved": {
            "nTopGenes": cfg.nTopGenes,
            "nPcs": cfg.nPcs,
            "nNeighbors": cfg.nNeighbors,
            "resolution": cfg.resolution,
        },
        "recommendation": recommendation,
        "approvedVersusRecommendedResolution": {
            "approved": cfg.resolution,
            "recommended": recommendation.get("resolution"),
            "matchesRecommendation": (
                recommendation.get("resolution") is not None
                and abs(float(recommendation["resolution"]) - float(cfg.resolution)) < 1e-9
            ),
        },
        "rfMerge": None,  # not implemented at time of submission
        "sampling": sample_metadata(validated),
        "subsetH5ad": rel_to_repo(cfg.outputH5ad),
        "runJson": rel_to_repo(cfg.outputH5ad.with_name(f"{cfg.outputH5ad.stem}_run.json")),
        "figuresDir": rel_to_repo(cfg.figsDir),
        "scib": {
            "csv": rel_to_repo(scib_dir / "scib_results.csv"),
            "svg": rel_to_repo(scib_dir / "scib_results.svg"),
        },
        "timingsSeconds": round(time.perf_counter() - started, 3),
        "note": (
            "Review the full scIB metric table before launching full-atlas production. "
            "There is no automatic pass/fail threshold."
        ),
    }
    write_json(cfg.validationDir / "subset_validation_summary.json", validation_summary)
    log.info("Validation summary: %s", json.dumps(validation_summary["resolved"]))
    return validation_summary
