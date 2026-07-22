import argparse
from pathlib import Path

from scib_benchmark import DEFAULT_EMBEDDING_KEYS, run_scib_benchmark
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

log = configure_file_logger("scib_benchmark.log", __name__)
add_stdout_handler()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scib_metrics benchmark on atlas embeddings.")
    parser.add_argument("--input", type=Path, required=True, metavar="PATH", help="Input atlas h5ad")
    parser.add_argument("--out-dir", type=Path, required=True, metavar="PATH", help="Output directory")
    parser.add_argument("--batch-key", type=str, default="study_accession", metavar="COL", help="obs batch column")
    parser.add_argument("--label-key", type=str, default="cell_type", metavar="COL", help="obs label column")
    parser.add_argument(
        "--embeddings",
        nargs="+",
        default=DEFAULT_EMBEDDING_KEYS,
        metavar="KEY",
        help="obsm keys to benchmark",
    )
    parser.add_argument("--n-jobs", type=int, default=6, metavar="N", help="Parallel jobs for neighbor search")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run benchmark even when scib_results.csv and scib_results.svg exist",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    log.info("scib benchmark run started")
    log.info(
        "input=%s out_dir=%s batch_key=%s label_key=%s embeddings=%s n_jobs=%d force=%s",
        rel_to_repo(args.input),
        rel_to_repo(args.out_dir),
        args.batch_key,
        args.label_key,
        args.embeddings,
        args.n_jobs,
        args.force,
    )
    run_scib_benchmark(
        input_h5ad=args.input,
        out_dir=args.out_dir,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embedding_keys=args.embeddings,
        n_jobs=args.n_jobs,
        force=args.force,
        log=log,
    )
    log.info("scib benchmark run complete")


if __name__ == "__main__":
    main()
