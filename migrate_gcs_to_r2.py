from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from gcs import download_from_gcs, gcs_blob_md5, gcs_local_path
from r2 import gcs_uri_to_r2_raw_key, r2_raw_matches_gcs, upload_to_r2
from r2.client import _MD5_METADATA_KEY, _local_md5_b64
from shared.logger import configure_file_logger
from shared.repo import REPO_ROOT

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"

log = configure_file_logger("migrate_gcs_to_r2.log", __name__)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


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

    for handler in log.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write("\n")
    log.info("new migration run started (datasets: %s)", args.datasets)

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
            skipped += 1
            continue

        if args.dry_run:
            log.info("%s (%s): would upload %s -> %s", srx, position, gs_uri, r2_key)
            continue

        local_path = gcs_local_path(gs_uri, REPO_ROOT / "data")
        downloaded = False
        try:
            log.info("%s (%s): downloading from GCS", srx, position)
            download_from_gcs(gs_uri, REPO_ROOT / "data")
            downloaded = True

            local_md5 = _local_md5_b64(local_path)
            if local_md5 != gcs_md5:
                raise RuntimeError(
                    f"Download integrity check failed: local MD5 {local_md5} != GCS MD5 {gcs_md5}"
                )
            log.info("%s (%s): MD5 verified (%s)", srx, position, local_md5)

            log.info("%s (%s): uploading to R2 at %s", srx, position, r2_key)
            upload_to_r2(local_path, r2_key, extra_metadata={_MD5_METADATA_KEY: gcs_md5})
            log.info("%s (%s): done", srx, position)
        except Exception:
            log.exception("%s (%s): failed", srx, position)
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
