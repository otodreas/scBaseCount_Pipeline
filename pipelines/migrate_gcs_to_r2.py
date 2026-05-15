from __future__ import annotations

import argparse
import datetime
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from gcs import download_from_gcs, gcs_blob_md5, gcs_local_path
from r2 import gcs_uri_to_r2_raw_key, r2_raw_matches_gcs, upload_to_r2
from r2.client import _MD5_METADATA_KEY, _local_md5_b64
from shared.csv_writer import append_csv_row
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "migration"

log = configure_file_logger("migrate_gcs_to_r2.log", __name__)
add_stdout_handler()


_CSV_COLUMNS = ["position", "srx", "status", "gs_uri", "r2_key", "md5", "timestamp", "error"]


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    gs_uri: str = "",
    r2_key: str = "",
    md5: str = "",
    error: str = "",
) -> None:
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [position, srx, status, gs_uri, r2_key, md5, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), error],
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy raw h5ad files from GCS to R2. Safe to re-run: files already in R2 with a matching MD5 are skipped."
    )
    parser.add_argument(
        "--datasets",
        type=Path,
        default=_DEFAULT_DATASETS_CSV,
        metavar="PATH",
        help=f"Path to datasets CSV (default: {_DEFAULT_DATASETS_CSV})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be uploaded without doing anything.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info("new migration run started (datasets: %s)", args.datasets)

    run_dir = RUN_OUTPUT_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / "run.csv"
    log.info("run summary output directory: %s", run_dir)

    datasets = pd.read_csv(args.datasets)
    log.info("Loaded %d accession(s) from %s", len(datasets), args.datasets)

    total = len(datasets)
    skipped = 0
    failed = 0

    for n, (_, row) in enumerate(datasets.iterrows(), start=1):
        srx = row["srx_accession"]
        gs_uri = row["file_path"]
        r2_key = gcs_uri_to_r2_raw_key(gs_uri)
        position = f"{n}/{total}"

        gcs_md5 = gcs_blob_md5(gs_uri)

        if r2_raw_matches_gcs(r2_key, gcs_md5):
            log.info("%s (%s): already in R2 with matching MD5, skipping", srx, position)
            _append_summary_row(csv_summary_path, srx, "skipped", position, gs_uri, r2_key, gcs_md5)
            skipped += 1
            continue

        if args.dry_run:
            log.info("%s (%s): would upload %s -> %s", srx, position, gs_uri, r2_key)
            _append_summary_row(csv_summary_path, srx, "dry_run", position, gs_uri, r2_key, gcs_md5)
            continue

        local_path = gcs_local_path(gs_uri, REPO_ROOT / "data")
        downloaded = False
        try:
            if local_path.exists():
                log.info("%s (%s): found local file at %s, skipping GCS download", srx, position, local_path)
            else:
                log.info("%s (%s): downloading from GCS", srx, position)
                download_from_gcs(gs_uri, REPO_ROOT / "data")
                downloaded = True

            local_md5 = _local_md5_b64(local_path)
            if local_md5 != gcs_md5:
                raise RuntimeError(f"Download integrity check failed: local MD5 {local_md5} != GCS MD5 {gcs_md5}")
            log.info("%s (%s): MD5 verified (%s)", srx, position, local_md5)

            log.info("%s (%s): uploading to R2 at %s", srx, position, r2_key)
            upload_to_r2(local_path, r2_key, extra_metadata={_MD5_METADATA_KEY: gcs_md5})
            if not r2_raw_matches_gcs(r2_key, gcs_md5):
                raise RuntimeError(
                    f"Post-upload metadata check failed: R2 object at {r2_key} does not reflect expected MD5"
                )
            log.info("%s (%s): post-upload metadata verified", srx, position)
            _append_summary_row(csv_summary_path, srx, "uploaded", position, gs_uri, r2_key, gcs_md5)
            log.info("%s (%s): done", srx, position)
        except Exception as exc:
            log.exception("%s (%s): failed", srx, position)
            _append_summary_row(
                csv_summary_path, srx, "failed", position, gs_uri, r2_key, gcs_md5, error=f"{type(exc).__name__}: {exc}"
            )
            failed += 1
        finally:
            if downloaded and local_path.exists():
                local_path.unlink()
                log.debug("%s: deleted local copy %s", srx, local_path)

    log.info(
        "Migration complete. %d/%d skipped (already in R2), %d failed.",
        skipped,
        total,
        failed,
    )


if __name__ == "__main__":
    main()
