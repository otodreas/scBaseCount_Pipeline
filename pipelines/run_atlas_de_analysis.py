"""Checkpointed full-atlas disease DE and noteworthy-gene discovery runner."""

import argparse
import datetime
import json
from pathlib import Path

import pandas as pd
from disease_markers.aggregation import aggregate_atlas, build_aggregate_fingerprint, load_checkpoint_if_valid
from disease_markers.analysis import analyze_from_pseudobulk
from disease_markers.config import AtlasDeAnalysisConfig
from metadata.config import MetadataConfig
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

_LOG_FILENAME = "atlas_de_analysis.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
configure_file_logger(_LOG_FILENAME, "disease_markers")
add_stdout_handler()

_DEFAULT_CFG = AtlasDeAnalysisConfig()


def _parse_args() -> argparse.Namespace:
    d = _DEFAULT_CFG
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate atlas raw counts into sample x Leiden pseudobulks, run study-aware "
            "disease DE, and write an adaptive noteworthy-gene review shortlist."
        )
    )
    parser.add_argument(
        "--stage",
        choices=("aggregate", "analyze", "all"),
        default="all",
        help="Pipeline stage to run",
    )
    parser.add_argument("--atlas", type=Path, default=d.atlasPath, metavar="PATH", help="Postprocessed atlas h5ad")
    parser.add_argument(
        "--output-dir", type=Path, default=d.outputDir, metavar="PATH", help="Analysis output directory"
    )
    parser.add_argument("--contexts", type=Path, default=d.contextsPath, metavar="PATH", help="contexts JSONL")
    parser.add_argument("--atlas-csv", type=Path, default=d.atlasCsvPath, metavar="PATH", help="Atlas accession CSV")
    parser.add_argument("--gene-info", type=Path, default=d.geneInfoPath, metavar="PATH", help="STAR geneInfo.tab")
    parser.add_argument(
        "--memory-reserve-gib",
        type=float,
        default=d.memoryReserveBytes / (1024**3),
        metavar="GIB",
        help="RAM reserved beyond the estimated atlas sparse matrix",
    )
    parser.add_argument("--primary-budget", type=int, default=d.primaryBudget, metavar="N")
    parser.add_argument("--extended-budget", type=int, default=d.extendedBudget, metavar="N")
    parser.add_argument("--min-cells-per-profile", type=int, default=d.minCellsPerProfile, metavar="N")
    parser.add_argument("--force-aggregate", action="store_true", help="Ignore an existing pseudobulk checkpoint")
    return parser.parse_args()


def build_config(args: argparse.Namespace) -> AtlasDeAnalysisConfig:
    return _DEFAULT_CFG.model_copy(
        update={
            "atlasPath": args.atlas,
            "outputDir": args.output_dir,
            "contextsPath": args.contexts,
            "atlasCsvPath": args.atlas_csv,
            "geneInfoPath": args.gene_info,
            "memoryReserveBytes": int(args.memory_reserve_gib * 1024**3),
            "primaryBudget": args.primary_budget,
            "extendedBudget": args.extended_budget,
            "minCellsPerProfile": args.min_cells_per_profile,
        }
    )


def main() -> None:
    args = _parse_args()
    cfg = build_config(args)
    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    cfg.figuresDir.mkdir(parents=True, exist_ok=True)
    cfg.checkpointsDir.mkdir(parents=True, exist_ok=True)

    log_run_separator(log)
    log.info("new atlas DE analysis run started")
    log.info("stage=%s config=%s", args.stage, cfg.model_dump_json())
    log.info("sample metadata=%s", rel_to_repo(MetadataConfig().sampleParquetPath))
    started = datetime.datetime.now()

    fingerprint = build_aggregate_fingerprint(cfg)
    summary: dict = {
        "stage": args.stage,
        "atlasPath": rel_to_repo(cfg.atlasPath),
        "outputDir": rel_to_repo(cfg.outputDir),
        "fingerprintSha": fingerprint["sha256"],
        "memoryReserveBytes": cfg.memoryReserveBytes,
        "startedAt": started.isoformat(timespec="seconds"),
    }

    if args.stage in {"aggregate", "all"}:
        pdata = aggregate_atlas(cfg, reuseCheckpoint=not args.force_aggregate)
        summary["pseudobulkProfiles"] = int(pdata.n_obs)
        summary["pseudobulkGenes"] = int(pdata.n_vars)
    else:
        pdata = load_checkpoint_if_valid(cfg)
        if pdata is None:
            raise SystemExit(
                f"No valid pseudobulk checkpoint at {rel_to_repo(cfg.pseudobulkPath)}. Run --stage aggregate first."
            )

    if args.stage in {"analyze", "all"}:
        analysis_summary = analyze_from_pseudobulk(
            cfg,
            pdata,
            clusterObs=pd.DataFrame(pdata.obs),
            fingerprintSha=str(fingerprint["sha256"]),
        )
        summary.update(analysis_summary)

    summary["finishedAt"] = datetime.datetime.now().isoformat(timespec="seconds")
    summary["elapsedSeconds"] = (datetime.datetime.now() - started).total_seconds()
    run_path = cfg.outputDir / "run_summary.json"
    run_path.write_text(json.dumps(summary, indent=2) + "\n")
    log.info("atlas DE analysis complete in %s", datetime.datetime.now() - started)
    log.info("Wrote %s", rel_to_repo(run_path))


if __name__ == "__main__":
    main()
