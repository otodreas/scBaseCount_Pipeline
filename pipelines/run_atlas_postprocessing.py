import argparse
import datetime
from pathlib import Path

from atlas_postprocessing.artifacts import (
    apply_parameters_to_config,
    load_approved_parameters,
    reject_tuning_overrides_with_parameters_json,
    validate_approved_against_calibration,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import run_postprocessing
from dotenv import load_dotenv
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

load_dotenv()

_LOG_FILENAME = "atlas_postprocessing.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
configure_file_logger(_LOG_FILENAME, "atlas_postprocessing")
add_stdout_handler()

_DEFAULT_CFG = AtlasPostprocessingConfig()
_TUNING_FLAGS = ("n_top_genes", "n_pcs", "n_neighbors", "resolution")


def _parse_args() -> argparse.Namespace:
    d = _DEFAULT_CFG
    parser = argparse.ArgumentParser(
        description="Run atlas postprocessing (HVG, PCA, Harmony, neighbors, UMAP, Leiden) with one parameter set."
    )
    parser.add_argument("--input", type=Path, default=d.inputH5ad, metavar="PATH", help="Input atlas h5ad")
    parser.add_argument("--output", type=Path, default=d.outputH5ad, metavar="PATH", help="Output atlas h5ad")
    parser.add_argument("--figs-dir", type=Path, default=d.figsDir, metavar="PATH", help="Directory for UMAP PNGs")
    parser.add_argument("--batch-key", type=str, default=d.batchKey, metavar="COL", help="obs batch column")
    parser.add_argument("--cell-type-key", type=str, default=d.cellTypeKey, metavar="COL", help="obs cell type column")
    parser.add_argument("--n-top-genes", type=int, default=None, metavar="N", help="Number of HVGs")
    parser.add_argument("--n-pcs", type=int, default=None, metavar="N", help="PCs used for the neighbor graph")
    parser.add_argument("--n-pcs-compute", type=int, default=d.nPcsCompute, metavar="N", help="PCs computed by PCA")
    parser.add_argument("--n-neighbors", type=int, default=None, metavar="N", help="Neighbors for the graph")
    parser.add_argument("--resolution", type=float, default=None, metavar="R", help="Leiden resolution")
    parser.add_argument(
        "--parameters-json",
        type=Path,
        default=None,
        metavar="PATH",
        help="Approved parameter JSON from calibration (authoritative for the four tuning knobs)",
    )
    parser.add_argument("--no-plots", action="store_true", help="Skip writing UMAP PNGs")
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        metavar="N",
        help="Thread budget for Scanpy, Harmony, and parallel UMAP (0 leaves library defaults)",
    )
    parser.add_argument("--r2-key", type=str, default=d.r2Key, metavar="KEY", help="R2 key")
    return parser.parse_args()


def _tuning_override_requested(args: argparse.Namespace) -> list[str]:
    return [name for name in _TUNING_FLAGS if getattr(args, name) is not None]


def build_config(args: argparse.Namespace) -> AtlasPostprocessingConfig:
    overrides = _tuning_override_requested(args)
    if args.parameters_json is not None:
        reject_tuning_overrides_with_parameters_json(overrides)

    cfg = _DEFAULT_CFG.model_copy(
        update={
            "inputH5ad": args.input,
            "outputH5ad": args.output,
            "figsDir": args.figs_dir,
            "batchKey": args.batch_key,
            "cellTypeKey": args.cell_type_key,
            "nPcsCompute": args.n_pcs_compute,
            "writePlots": not args.no_plots,
            "r2Key": args.r2_key,
            "nJobs": args.threads,
        }
    )

    if args.parameters_json is not None:
        parameters = load_approved_parameters(args.parameters_json)
        summary = validate_approved_against_calibration(parameters, parametersPath=args.parameters_json)
        cfg = apply_parameters_to_config(cfg, parameters, parametersPath=args.parameters_json)
        log.info("Imported parameters from %s", rel_to_repo(args.parameters_json))
        log.info(
            "Resolved nTopGenes=%s nPcs=%s nNeighbors=%s resolution=%s",
            cfg.nTopGenes,
            cfg.nPcs,
            cfg.nNeighbors,
            cfg.resolution,
        )
        if parameters.calibrationSummary:
            log.info("Calibration summary: %s", parameters.calibrationSummary)
        log.info("Calibration baseline recorded in summary: %s", summary.get("baseline"))
        return cfg

    return cfg.model_copy(
        update={
            "nTopGenes": args.n_top_genes if args.n_top_genes is not None else _DEFAULT_CFG.nTopGenes,
            "nPcs": args.n_pcs if args.n_pcs is not None else _DEFAULT_CFG.nPcs,
            "nNeighbors": args.n_neighbors if args.n_neighbors is not None else _DEFAULT_CFG.nNeighbors,
            "resolution": args.resolution if args.resolution is not None else _DEFAULT_CFG.resolution,
        }
    )


def main() -> None:
    args = _parse_args()
    cfg = build_config(args)

    log_run_separator(log)
    log.info("new atlas postprocessing run started")
    log.info("config: %s", cfg.model_dump_json())

    started = datetime.datetime.now()
    # RF merge of leiden_atlas could slot inside this production run, after Harmony Leiden.
    run_postprocessing(cfg, workflow="production")
    log.info("atlas postprocessing run complete in %s", datetime.datetime.now() - started)


if __name__ == "__main__":
    main()
