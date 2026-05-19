from __future__ import annotations

import argparse
import datetime
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from cluster_validation import ClusterValidationConfig, ClusterValidationResult, run_cluster_validation
from dotenv import load_dotenv
from gcs import download_from_gcs, gcs_local_path, verify_download
from r2 import (
    download_from_r2,
    fetch_uploaded_r2_keys,
    gcs_uri_to_r2_raw_key,
    r2_key_exists,
    r2_object_md5,
    upload_to_r2,
    verify_upload,
)
from r2.client import _local_md5_b64
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
GCS_LOCAL_ROOT = REPO_ROOT / "data"
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "clustering_pipeline"
_LOG_FILENAME = "clustering_pipeline.log"

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


def read_datasets(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
    run_dir: Path,
    figs_dir: Path,
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "run_dir": rel_to_repo(run_dir),
        "run_csv": rel_to_repo(run_dir / "run.csv"),
        "figs_dir": rel_to_repo(figs_dir),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "r2_prefix": args.r2_prefix,
        "datasets_path": str(args.datasets),
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


_CSV_COLUMNS = [
    "position",
    "srx",
    "status",
    "r2_file",
    "timestamp",
    "error",
    "selectedResolution",
    "nPcs",
    "cumvar",
    "nClustersPreMerge",
    "nClustersPostMerge",
    "nMerges",
    "kPrior",
    "kFiltered",
    "nCellsDropped",
    "nCellsFinal",
    "jaccAtSelected",
]


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    r2_key: str = "",
    error: str = "",
    result: ClusterValidationResult | None = None,
) -> None:
    if result is None:
        stats_cells: list[str] = [""] * 11
    else:
        n_merges = result.nClustersPreMerge - result.nClustersPostMerge
        jacc_at_selected = result.jaccArr[result.resolutions.index(result.selectedResolution)]
        stats_cells = [
            str(result.selectedResolution),
            str(result.nPcs),
            str(result.cumvar),
            str(result.nClustersPreMerge),
            str(result.nClustersPostMerge),
            str(n_merges),
            str(result.kPrior),
            str(result.kFiltered),
            str(result.nCellsDropped),
            str(result.nCellsFinal),
            str(jacc_at_selected),
        ]
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [
            position,
            srx,
            status,
            r2_key,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            error,
            *stats_cells,
        ],
    )


def _fetch_raw_h5ad(srx: str, gs_uri: str, raw_h5ad: Path) -> None:
    r2_raw_key = gcs_uri_to_r2_raw_key(gs_uri)
    if r2_key_exists(r2_raw_key):
        log.info("%s: downloading raw h5ad from R2 (%s)", srx, r2_raw_key)
        download_from_r2(r2_raw_key, raw_h5ad)
        stored_md5 = r2_object_md5(r2_raw_key)
        if stored_md5:
            local_md5 = _local_md5_b64(raw_h5ad)
            if local_md5 != stored_md5:
                raise RuntimeError(
                    f"{srx}: R2 download integrity check failed: local MD5 {local_md5} != stored MD5 {stored_md5}"
                )
            log.info("%s: R2 download MD5 verified (%s)", srx, local_md5)
        return
    log.info("%s: raw h5ad not in R2, falling back to GCS", srx)
    download_from_gcs(gs_uri, GCS_LOCAL_ROOT)
    if not verify_download(gs_uri, GCS_LOCAL_ROOT):
        raise RuntimeError(f"Download verification failed for {srx}: file not found at {raw_h5ad}")


def process_accession(
    srx: str,
    gs_uri: str,
    r2_key: str,
    datasets_path: Path,
    figs_root: Path,
) -> tuple[ClusterValidationResult | None, Exception | None]:
    cfg = ClusterValidationConfig(srxAccession=srx, summaryPath=datasets_path, figsDir=figs_root)
    raw_h5ad = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)

    downloaded = False
    try:
        if raw_h5ad.exists():
            log.info("%s: found local h5ad at %s, skipping download", srx, raw_h5ad)
        else:
            _fetch_raw_h5ad(srx, gs_uri, raw_h5ad)
            downloaded = True

        _, result = run_cluster_validation(cfg)

        if downloaded:
            safe_delete(raw_h5ad, log)

        upload_to_r2(result.adataPath, r2_key)
        verify_upload(r2_key)
        safe_delete(result.adataPath, log)

        log.info("%s: done", srx)
        return result, None

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        safe_delete(raw_h5ad, log)
        safe_delete(cfg.outputDir / f"{srx}_clustered.h5ad", log)
        return None, exc


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the clustering pipeline (clustering only, no annotation).")
    parser.add_argument(
        "--datasets",
        type=Path,
        default=_DEFAULT_DATASETS_CSV,
        metavar="PATH",
        help=f"Path to datasets CSV (default: {_DEFAULT_DATASETS_CSV})",
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        default=f"clustering_pipeline_{RUN_TIMESTAMP}",
        metavar="PREFIX",
        help="R2 prefix (default: clustering_pipeline_{RUN_TIMESTAMP})",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Write a metadata JSON file next to the run CSV with this note.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Number of accessions to process in parallel (default: 1)",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info("new clustering pipeline run started (r2 prefix: %s)", args.r2_prefix)

    datasets = read_datasets(args.datasets)
    log.info("Loaded %d accession(s) from %s", len(datasets), args.datasets)
    uploaded = fetch_uploaded_r2_keys()

    run_dir = RUN_OUTPUT_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / "run.csv"
    metadata_path = run_dir / "metadata.json"
    figs_root = run_dir / "figs"
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir, figs_root)
    log.info("run summary output directory: %s", run_dir)

    total = len(datasets)
    skipped = 0
    r2_suffix = "_clustered.h5ad"

    work_items: list[tuple[str, str, str, str]] = []
    for n, (_, row) in enumerate(datasets.iterrows(), start=1):
        srx = row["srx_accession"]
        gs_uri = row["file_path"]
        r2_key = f"{args.r2_prefix}/{srx}{r2_suffix}"
        position = f"{n}/{total}"
        if r2_key in uploaded:
            log.info("%s: already uploaded, skipping", srx)
            _append_summary_row(csv_summary_path, srx, "skipped", position, r2_key)
            skipped += 1
            continue
        work_items.append((srx, gs_uri, r2_key, position))

    log.info("Submitting %d accession(s) to %d worker(s)", len(work_items), args.workers)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_accession, srx, gs_uri, r2_key, args.datasets, figs_root): (srx, r2_key, position)
            for srx, gs_uri, r2_key, position in work_items
        }
        for future in as_completed(futures):
            srx, r2_key, position = futures[future]
            result, exc = future.result()
            if exc is not None:
                log.warning("%s: failed", srx)
                _append_summary_row(csv_summary_path, srx, "failed", position, error=f"{type(exc).__name__}: {exc}")
                skipped += 1
            else:
                _append_summary_row(csv_summary_path, srx, "success", position, r2_key, result=result)

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)


if __name__ == "__main__":
    main()
