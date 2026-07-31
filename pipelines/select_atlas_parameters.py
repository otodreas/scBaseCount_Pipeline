import argparse
import datetime
import json
import time
from pathlib import Path

import scanpy as sc
from atlas_postprocessing.artifacts import (
    apply_parameters_to_config,
    load_approved_parameters,
    validate_approved_against_calibration,
    write_json,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import run_postprocessing, timed
from atlas_postprocessing.sampling import SAMPLE_SEED, sample_metadata, sample_study_proportional
from atlas_postprocessing.scib import run_scib_benchmark
from atlas_postprocessing.selection import run_calibration
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

_LOG_FILENAME = "select_atlas_parameters.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
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
    parser.add_argument("--threads", type=int, default=0, metavar="N", help="scanpy n_jobs (0 leaves the default)")


def _parse_args() -> argparse.Namespace:
    d = _DEFAULT_CFG
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate atlas postprocessing parameters on a study-proportional sample of the full atlas, "
            "or validate an approved set on the same sampling policy."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)

    calibrate = sub.add_parser("calibrate", help="Sweep HVG/PC/neighbor/resolution and write diagnostics")
    _add_shared_args(calibrate)
    calibrate.add_argument(
        "--output-dir",
        type=Path,
        default=d.calibrationDir,
        metavar="PATH",
        help="Calibration output root (metrics/, figures/, JSON manifests)",
    )
    calibrate.add_argument("--n-top-genes", type=int, default=d.nTopGenes, metavar="N", help="Baseline HVGs")
    calibrate.add_argument("--n-pcs", type=int, default=d.nPcs, metavar="N", help="Baseline PCs")
    calibrate.add_argument("--n-pcs-compute", type=int, default=d.nPcsCompute, metavar="N", help="PCs computed by PCA")
    calibrate.add_argument("--n-neighbors", type=int, default=d.nNeighbors, metavar="N", help="Baseline neighbors")
    calibrate.add_argument(
        "--resolution", type=float, default=d.resolution, metavar="R", help="Baseline Leiden resolution"
    )
    calibrate.add_argument(
        "--hvg-candidates",
        type=int,
        nargs="+",
        default=d.hvgCandidates,
        metavar="N",
        help="HVG candidate list",
    )
    calibrate.add_argument(
        "--pc-candidates",
        type=int,
        nargs="+",
        default=d.pcCandidates,
        metavar="N",
        help="PC candidate list",
    )
    calibrate.add_argument(
        "--neighbor-candidates",
        type=int,
        nargs="+",
        default=d.neighborCandidates,
        metavar="N",
        help="Neighbor candidate list",
    )
    calibrate.add_argument(
        "--resolution-candidates",
        type=float,
        nargs="+",
        default=d.resolutionCandidates,
        metavar="R",
        help="Leiden resolution candidate list",
    )

    validate = sub.add_parser("validate", help="Run approved parameters on the sample and scIB-benchmark Harmony")
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
    return _DEFAULT_CFG.model_copy(
        update={
            "inputH5ad": args.input,
            "calibrationDir": args.output_dir,
            "batchKey": args.batch_key,
            "cellTypeKey": args.cell_type_key,
            "nTopGenes": args.n_top_genes,
            "nPcs": args.n_pcs,
            "nPcsCompute": args.n_pcs_compute,
            "nNeighbors": args.n_neighbors,
            "resolution": args.resolution,
            "hvgCandidates": list(args.hvg_candidates),
            "pcCandidates": list(args.pc_candidates),
            "neighborCandidates": list(args.neighbor_candidates),
            "resolutionCandidates": list(args.resolution_candidates),
            "writePlots": False,
        }
    )


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
    log.info("Calibration summary baseline: %s", summary.get("baseline"))
    return cfg


def _load_and_sample(cfg: AtlasPostprocessingConfig, sampleCells: int) -> sc.AnnData:
    """Load the full atlas once and overwrite the in-memory object with a representative sample."""
    adata = timed("load full atlas", lambda: sc.read_h5ad(cfg.inputH5ad), logger=log)
    log.info(
        "Full atlas loaded: %s cells x %s genes (%s studies)",
        f"{adata.n_obs:,}",
        f"{adata.n_vars:,}",
        adata.obs[cfg.batchKey].nunique() if cfg.batchKey in adata.obs else "?",
    )
    adata = timed(
        "study-proportional sample",
        lambda: sample_study_proportional(
            adata,
            n=sampleCells,
            stratifyKey=cfg.batchKey,
            seed=SAMPLE_SEED,
        ),
        logger=log,
    )
    return adata


def _run_validate(args: argparse.Namespace) -> None:
    cfg = _cfg_from_validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scib_dir = args.output_dir / "scib"

    started = time.perf_counter()
    sampled = _load_and_sample(cfg, args.sample_cells)
    adata = timed(
        "approved subset postprocessing",
        lambda: run_postprocessing(cfg, adata=sampled),
        logger=log,
    )
    timed(
        "scIB benchmark",
        lambda: run_scib_benchmark(
            adata,
            outDir=scib_dir,
            batchKey=cfg.batchKey,
            labelKey=cfg.cellTypeKey,
            nJobs=args.scib_jobs,
            force=args.force_scib,
        ),
        logger=log,
    )

    parameters = load_approved_parameters(args.parameters_json)
    validation_summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "outputDir": rel_to_repo(args.output_dir),
        "parametersJson": rel_to_repo(args.parameters_json),
        "calibrationSummary": parameters.calibrationSummary,
        "resolved": {
            "nTopGenes": cfg.nTopGenes,
            "nPcs": cfg.nPcs,
            "nNeighbors": cfg.nNeighbors,
            "resolution": cfg.resolution,
        },
        "sampling": sample_metadata(adata),
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
    write_json(args.output_dir / "subset_validation_summary.json", validation_summary)
    log.info("Validation summary: %s", json.dumps(validation_summary["resolved"]))


def main() -> None:
    args = _parse_args()

    if args.threads > 0:
        sc.settings.n_jobs = args.threads

    log_run_separator(log)
    log.info("select_atlas_parameters %s started", args.command)

    started = datetime.datetime.now()
    if args.command == "calibrate":
        cfg = _cfg_from_calibrate_args(args)
        log.info("config: %s", cfg.model_dump_json())
        sampled = _load_and_sample(cfg, args.sample_cells)
        run_calibration(cfg, adata=sampled)
    elif args.command == "validate":
        _run_validate(args)
    else:
        raise ValueError(f"Unknown command {args.command!r}")

    log.info("select_atlas_parameters %s complete in %s", args.command, datetime.datetime.now() - started)


if __name__ == "__main__":
    main()
