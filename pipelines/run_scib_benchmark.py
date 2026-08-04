import argparse
from pathlib import Path

from atlas_postprocessing.scib import run_scib_benchmark
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

log = configure_file_logger("scib_benchmark.log", __name__)
add_stdout_handler()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scib_metrics benchmark on atlas embeddings.")
    parser.add_argument("--input", type=Path, required=True, metavar="PATH", help="Input atlas h5ad")
    parser.add_argument("--out-dir", type=Path, required=True, metavar="PATH", help="Output directory")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run benchmark even when scib_results.csv and scib_results.svg exist",
    )
    parser.add_argument("--batch-key", type=str, default="study_accession", metavar="COL", help="obs batch column")
    parser.add_argument("--label-key", type=str, default="cell_type", metavar="COL", help="obs label column")
    parser.add_argument("--jobs", type=int, default=6, metavar="N", help="scIB n_jobs")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    log.info("scib benchmark run started")
    log.info(
        "input=%s out_dir=%s force=%s",
        rel_to_repo(args.input),
        rel_to_repo(args.out_dir),
        args.force,
    )
    run_scib_benchmark(
        args.input,
        outDir=args.out_dir,
        batchKey=args.batch_key,
        labelKey=args.label_key,
        nJobs=args.jobs,
        force=args.force,
    )
    log.info("scib benchmark complete")


if __name__ == "__main__":
    main()
