"""Compare a baseline atlas against a ribosomal-filtered candidate atlas."""

import argparse
from pathlib import Path

from dotenv import load_dotenv
from h5ad_concat.compare import AtlasCompareConfig, compare_atlases
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from storage import download_from_r2

load_dotenv()

_LOG_FILENAME = "atlas_compare.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


def _ensure_local_h5ad(local_path: Path, r2_key: str) -> None:
    if local_path.exists():
        log.info("Using local atlas %s", rel_to_repo(local_path))
        return
    log.info("Local atlas missing at %s; downloading r2 key %s", rel_to_repo(local_path), r2_key)
    download_from_r2(r2_key, local_path, verify_md5=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare baseline and candidate atlas h5ad files")
    parser.add_argument(
        "--baseline",
        type=Path,
        default=REPO_ROOT / "output/atlas/v2/atlas_v2.h5ad",
        help="Baseline atlas h5ad",
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        default=REPO_ROOT / "output/atlas/2026-08-12/atlas.h5ad",
        help="Candidate atlas h5ad",
    )
    parser.add_argument(
        "--baseline-r2-fallback",
        required=True,
        help="R2 object key used when --baseline is missing locally",
    )
    parser.add_argument(
        "--candidate-r2-fallback",
        required=True,
        help="R2 object key used when --candidate is missing locally",
    )
    parser.add_argument(
        "--baseline-manifest",
        type=Path,
        default=REPO_ROOT / "output/atlas/v2/atlas_v2_result.json",
        help="Baseline result manifest JSON",
    )
    parser.add_argument(
        "--candidate-manifest",
        type=Path,
        default=REPO_ROOT / "output/atlas/2026-08-12/atlas_result.json",
        help="Candidate result manifest JSON",
    )
    parser.add_argument(
        "--baseline-status-csv",
        type=Path,
        default=REPO_ROOT / "output/atlas/v2/atlas_v2.csv",
        help="Baseline status CSV used when the baseline manifest lacks files[]",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=REPO_ROOT / "output/atlas/2026-08-12/atlas_compare_report.json",
        help="Repository-relative JSON comparison report path",
    )
    parser.add_argument("--chunk-size", type=int, default=2048, help="Obs chunk size for matrix compares")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    _ensure_local_h5ad(args.baseline, args.baseline_r2_fallback)
    _ensure_local_h5ad(args.candidate, args.candidate_r2_fallback)
    cfg = AtlasCompareConfig(
        baselinePath=args.baseline,
        candidatePath=args.candidate,
        baselineManifestPath=args.baseline_manifest,
        candidateManifestPath=args.candidate_manifest,
        baselineStatusCsvPath=args.baseline_status_csv,
        reportPath=args.report,
        chunkSize=args.chunk_size,
    )
    log.info(
        "Comparing baseline=%s candidate=%s",
        rel_to_repo(cfg.baselinePath),
        rel_to_repo(cfg.candidatePath),
    )
    report = compare_atlases(cfg)
    log.info("Wrote comparison report to %s", rel_to_repo(cfg.reportPath))
    log.info(
        "byteIdentical=%s fullLogicalIdentical=%s candidateOnly=%d baselineOnly=%d",
        report.byteIdentical,
        report.fullLogicalIdentical,
        report.nCandidateOnlyCells,
        report.nBaselineOnlyCells,
    )


if __name__ == "__main__":
    main()
