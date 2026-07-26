import argparse
from pathlib import Path

import scanpy as sc
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation
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
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    log.info("scib benchmark run started")
    log.info(
        "input=%s out_dir=%s, force=%s",
        rel_to_repo(args.input),
        rel_to_repo(args.out_dir),
        args.force,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    bm = Benchmarker(
        adata=sc.read_h5ad(args.input, backed="r"),
        batch_key="study_accession",
        label_key="cell_type",
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=["X_pca", "X_pca_harmony"],
        pre_integrated_embedding_obsm_key="X_pca",
        n_jobs=6,
    )
    bm.benchmark()
    log.info("scib benchmark metrics computed")
    df = bm.get_results()
    df.to_csv(args.out_dir / "scib_results.csv")
    bm.plot_results_table(show=False, save_dir=str(args.out_dir))
    log.info("scib benchmark plot & csv saved")


if __name__ == "__main__":
    main()
