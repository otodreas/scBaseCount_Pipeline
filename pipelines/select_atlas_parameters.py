import argparse
import datetime
from pathlib import Path

import scanpy as sc
from atlas_postprocessing.artifacts import (
    apply_parameters_to_config,
    load_approved_parameters,
    validate_approved_against_calibration,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import apply_thread_settings, timed
from atlas_postprocessing.sampling import sample_study_proportional
from atlas_postprocessing.selection import FIXED_N_NEIGHBORS, FIXED_N_TOP_GENES, run_calibration
from atlas_postprocessing.validation import run_validation
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

_LOG_FILENAME = "select_atlas_parameters.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
configure_file_logger(_LOG_FILENAME, "atlas_postprocessing")
add_stdout_handler()

_DEFAULT_CFG = AtlasPostprocessingConfig()


def _add_shared_args(parser: argparse.ArgumentParser) -> None:
    d = _DEFAULT_CFG
    parser.add_argument("--input", type=Path, required=True, metavar="PATH", help="Full atlas h5ad")
    parser.add_argument(
        "--sample-cells",
        type=int,
        required=True,
        metavar="N",
        help="Exact number of cells in the in-memory representative sample",
    )
    parser.add_argument("--batch-key", type=str, default=d.batchKey, metavar="COL", help="obs batch column")
    parser.add_argument("--cell-type-key", type=str, default=d.cellTypeKey, metavar="COL", help="obs cell type column")
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        metavar="N",
        help="Thread budget for Scanpy and Harmony (0 leaves library defaults)",
    )


def _parse_args() -> argparse.Namespace:
    d = _DEFAULT_CFG
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate atlas postprocessing parameters with cluster-validation selection on a "
            "Harmony graph, or validate an approved set on the same sampling policy."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)

    calibrate = sub.add_parser(
        "calibrate",
        help="Fix HVGs/neighbors, choose PCs by cumvar, sweep resolutions, write advisory recommendation",
    )
    _add_shared_args(calibrate)
    calibrate.add_argument(
        "--output-dir",
        type=Path,
        default=d.calibrationDir,
        metavar="PATH",
        help="Calibration output root (metrics/, figures/, JSON manifests)",
    )
    calibrate.add_argument(
        "--n-pcs-compute",
        type=int,
        default=d.nPcsCompute,
        metavar="N",
        help="PCs computed by PCA before the adaptive chooser",
    )
    calibrate.add_argument(
        "--resolution-candidates",
        type=float,
        nargs="+",
        default=None,
        metavar="R",
        help="Leiden resolution candidate list (default: cluster-validation 0.1..1.9 step 0.1)",
    )

    validate = sub.add_parser("validate", help="Run approved parameters on the sample and scIB-benchmark")
    _add_shared_args(validate)
    validate.add_argument(
        "--parameters-json",
        type=Path,
        required=True,
        metavar="PATH",
        help="Approved parameter JSON from calibration review",
    )
    validate.add_argument(
        "--output-dir",
        type=Path,
        default=d.validationDir,
        metavar="PATH",
        help="Subset validation output root",
    )
    validate.add_argument("--no-plots", action="store_true", help="Skip writing UMAP PNGs")
    validate.add_argument("--scib-jobs", type=int, default=6, metavar="N", help="scIB n_jobs")
    validate.add_argument(
        "--force-scib",
        action="store_true",
        help="Re-run scIB even when artifacts already exist",
    )
    return parser.parse_args()


def _cfg_from_calibrate_args(args: argparse.Namespace) -> AtlasPostprocessingConfig:
    update: dict = {
        "inputH5ad": args.input,
        "calibrationDir": args.output_dir,
        "batchKey": args.batch_key,
        "cellTypeKey": args.cell_type_key,
        "nTopGenes": FIXED_N_TOP_GENES,
        "nNeighbors": FIXED_N_NEIGHBORS,
        "nPcsCompute": args.n_pcs_compute,
        "writePlots": False,
        "nJobs": args.threads,
    }
    if args.resolution_candidates is not None:
        update["resolutionCandidates"] = list(args.resolution_candidates)
    return _DEFAULT_CFG.model_copy(update=update)


def _cfg_from_validate_args(args: argparse.Namespace) -> AtlasPostprocessingConfig:
    parameters = load_approved_parameters(args.parameters_json)
    summary = validate_approved_against_calibration(parameters, parametersPath=args.parameters_json)
    output_h5ad = args.output_dir / "atlas_pp_subset.h5ad"
    figs_dir = args.output_dir / "figures"
    cfg = _DEFAULT_CFG.model_copy(
        update={
            "inputH5ad": args.input,
            "outputH5ad": output_h5ad,
            "figsDir": figs_dir,
            "validationDir": args.output_dir,
            "batchKey": args.batch_key,
            "cellTypeKey": args.cell_type_key,
            "writePlots": not args.no_plots,
            "r2Key": None,
            "nJobs": args.threads,
        }
    )
    cfg = apply_parameters_to_config(cfg, parameters, parametersPath=args.parameters_json)
    log.info("Imported parameters from %s", rel_to_repo(args.parameters_json))
    log.info(
        "Resolved nTopGenes=%s nPcs=%s nNeighbors=%s resolution=%s",
        cfg.nTopGenes,
        cfg.nPcs,
        cfg.nNeighbors,
        cfg.resolution,
    )
    recommendation = summary.get("recommendation") or {}
    if recommendation:
        log.info("Calibration recommendation: %s", recommendation)
    log.info("Calibration summary baseline: %s", summary.get("baseline"))
    return cfg


def _load_and_sample(cfg: AtlasPostprocessingConfig, sampleCells: int, command: str) -> sc.AnnData:
    """Load the full atlas once and overwrite the in-memory object with a representative sample."""
    adata = timed("load full atlas", lambda: sc.read_h5ad(cfg.inputH5ad), logger=log)

    # Use a different random seed for calibrate and validate.
    seed = 1234 if command == "calibrate" else 6789

    log.info(
        "Full atlas loaded: %s cells x %s genes (%s studies) with random seed %s",
        f"{adata.n_obs:,}",
        f"{adata.n_vars:,}",
        adata.obs[cfg.batchKey].nunique() if cfg.batchKey in adata.obs else "?",
        seed,
    )

    adata = timed(
        "study-proportional sample",
        lambda: sample_study_proportional(
            adata,
            n=sampleCells,
            stratifyKey=cfg.batchKey,
            seed=seed,
        ),
        logger=log,
    )
    return adata


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info("select_atlas_parameters %s started", args.command)

    started = datetime.datetime.now()
    if args.command == "calibrate":
        cfg = _cfg_from_calibrate_args(args)
    elif args.command == "validate":
        cfg = _cfg_from_validate_args(args)
    else:
        raise ValueError(f"Unknown command {args.command!r}")

    apply_thread_settings(cfg)
    log.info("config: %s", cfg.model_dump_json())
    sampled = _load_and_sample(cfg, args.sample_cells, args.command)

    if args.command == "calibrate":
        run_calibration(cfg, adata=sampled)
    elif args.command == "validate":
        run_validation(
            cfg,
            adata=sampled,
            scibJobs=args.scib_jobs,
            forceScib=args.force_scib,
        )
    else:
        raise ValueError(f"Unknown command {args.command!r}")

    log.info("select_atlas_parameters %s complete in %s", args.command, datetime.datetime.now() - started)


if __name__ == "__main__":
    main()
