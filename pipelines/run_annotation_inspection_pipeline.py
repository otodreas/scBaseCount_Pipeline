from __future__ import annotations

import argparse
import datetime
import json
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from annotation_inspector.extremes import write_extremes_csv
from annotation_inspector.inspect import PAIR_COLUMNS, inspect_accession
from dotenv import load_dotenv
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from storage import download_from_r2, fetch_uploaded_r2_keys, r2_key_exists

load_dotenv()

RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
DOWNLOAD_ROOT = REPO_ROOT / "data" / "annotation_inspection"
RUN_OUTPUT_ROOT = REPO_ROOT / "output" / "annotation_inspection_pipeline"
_LOG_FILENAME = "annotation_inspection_pipeline.log"
_RUN_CSV_FILENAME = "run.csv"
_SUMMARY_FILENAME = "summary.csv"
_EXTREMES_FILENAME = "extremes.csv"
_ANNOTATED_H5AD_PATTERN = re.compile(r"^(?P<srx>SRX\d+)_annotated\.h5ad$")

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


_RUN_CSV_COLUMNS = [
    "position",
    "srx",
    "status",
    "input_r2_file",
    "cyteonto_r2_file",
    "cyteonto_found",
    "n_pairs",
    "n_cells",
    "timestamp",
    "duration_seconds",
    "error",
]


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
    run_dir: Path,
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "run_dir": rel_to_repo(run_dir),
        "run_csv": rel_to_repo(run_dir / _RUN_CSV_FILENAME),
        "summary_csv": rel_to_repo(run_dir / _SUMMARY_FILENAME),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "dry_run": bool(args.dry_run),
        "workers": args.workers,
        "top_n": args.top_n,
        "emit_extremes": args.emit_extremes,
    }
    if args.from_summary is not None:
        payload["from_summary"] = rel_to_repo(args.from_summary)
    else:
        payload["input_prefix"] = args.input_prefix
        payload["cyteonto_prefix"] = args.cyteonto_prefix
    if args.metadata is not None:
        payload["notes"] = args.metadata
    if args.emit_extremes:
        payload["extremes_csv"] = rel_to_repo(run_dir / _EXTREMES_FILENAME)
    metadata_path.write_text(json.dumps(payload, indent=2))


def _record_accession(
    csv_path: Path,
    srx: str,
    status: str,
    position: str = "",
    input_r2_key: str = "",
    cyteonto_r2_key: str = "",
    cyteonto_found: bool | str = "",
    n_pairs: int | str = "",
    n_cells: int | str = "",
    duration_seconds: float | None = None,
    error: str = "",
) -> None:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    duration_str = "" if duration_seconds is None else f"{duration_seconds:.2f}"
    append_csv_row(
        csv_path,
        _RUN_CSV_COLUMNS,
        [
            position,
            srx,
            status,
            input_r2_key,
            cyteonto_r2_key,
            cyteonto_found,
            n_pairs,
            n_cells,
            timestamp,
            duration_str,
            error,
        ],
    )


def _srx_from_input_key(input_prefix: str, r2_key: str) -> str | None:
    prefix = f"{input_prefix}/"
    if not r2_key.startswith(prefix):
        return None
    filename = r2_key[len(prefix) :]
    match = _ANNOTATED_H5AD_PATTERN.match(filename)
    if match is None:
        return None
    return match.group("srx")


def _discover_work_items(input_prefix: str) -> list[tuple[str, str, str]]:
    uploaded = fetch_uploaded_r2_keys(prefix=input_prefix)
    input_keys = sorted(key for key in uploaded if key.startswith(f"{input_prefix}/"))
    total = len(input_keys)
    work_items: list[tuple[str, str, str]] = []

    for n, input_r2_key in enumerate(input_keys, start=1):
        srx = _srx_from_input_key(input_prefix, input_r2_key)
        if srx is None:
            log.warning("Skipping unrecognized input key: %s", input_r2_key)
            continue
        work_items.append((srx, input_r2_key, f"{n}/{total}"))

    return work_items


def _append_summary_rows(summary_path: Path, pair_df: pd.DataFrame, lock: threading.Lock) -> None:
    with lock:
        write_header = not summary_path.exists()
        pair_df.to_csv(summary_path, mode="a", header=write_header, index=False)


def process_accession(
    srx: str,
    input_r2_key: str,
    cyteonto_prefix: str,
) -> tuple[pd.DataFrame | None, bool, float, Exception | None]:
    """Download, inspect, and clean up one accession; returns pair_df, cyteonto_found, elapsed, error."""
    run_start = time.monotonic()
    local_h5ad = DOWNLOAD_ROOT / f"{srx}_annotated.h5ad"
    cyteonto_r2_key = f"{cyteonto_prefix}/{srx}_cyteonto.csv"
    local_cyteonto = DOWNLOAD_ROOT / f"{srx}_cyteonto.csv"
    cyteonto_found = False

    try:
        if not r2_key_exists(input_r2_key):
            raise FileNotFoundError(f"annotated h5ad not found in R2 at {input_r2_key}")

        log.info("%s: downloading annotated h5ad from R2 (%s)", srx, input_r2_key)
        download_from_r2(input_r2_key, local_h5ad)

        cyteonto_path: Path | None = None
        if r2_key_exists(cyteonto_r2_key):
            log.info("%s: downloading CyteOnto CSV from R2 (%s)", srx, cyteonto_r2_key)
            download_from_r2(cyteonto_r2_key, local_cyteonto)
            cyteonto_path = local_cyteonto
            cyteonto_found = True
        else:
            log.warning("%s: CyteOnto CSV not found at %s; cytescore will be NaN", srx, cyteonto_r2_key)

        pair_df = inspect_accession(srx, local_h5ad, cyteonto_path)
        return pair_df, cyteonto_found, time.monotonic() - run_start, None

    except Exception as exc:
        log.exception("%s: inspection failed", srx)
        return None, cyteonto_found, time.monotonic() - run_start, exc

    finally:
        safe_delete(local_h5ad, log)
        safe_delete(local_cyteonto, log)


def run_from_summary(args: argparse.Namespace) -> None:
    summary_path = args.from_summary.resolve()
    if not summary_path.is_file():
        raise FileNotFoundError(f"summary file not found: {summary_path}")

    run_dir = args.output_dir.resolve() if args.output_dir is not None else summary_path.parent
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = run_dir / "metadata.json"
    extremes_path = run_dir / _EXTREMES_FILENAME

    log_run_separator(log)
    log.info("standalone extremes run from %s", summary_path)

    pair_df = pd.read_csv(summary_path)
    missing = [col for col in PAIR_COLUMNS if col not in pair_df.columns]
    if missing:
        raise ValueError(f"summary.csv missing required column(s): {missing}")

    write_extremes_csv(pair_df, args.top_n, extremes_path)
    log.info("wrote extremes to %s", extremes_path)

    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect CyteType-annotated h5ads from R2, join CyteOnto cytescores, "
            "and emit summary and optional extremes CSV."
        )
    )
    parser.add_argument(
        "--from-summary",
        type=Path,
        default=None,
        metavar="PATH",
        help="Build extremes.csv from an existing summary.csv without R2 fetch.",
    )
    parser.add_argument(
        "--input-prefix",
        type=str,
        default=None,
        metavar="PREFIX",
        help="R2 prefix containing annotated h5ads (required unless --from-summary).",
    )
    parser.add_argument(
        "--cyteonto-prefix",
        type=str,
        default=None,
        metavar="PREFIX",
        help="R2 prefix containing {srx}_cyteonto.csv files (required unless --from-summary).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Output directory for --from-summary (default: directory containing summary.csv).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Concurrent R2 fetch + inspect workers (default: 1).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        metavar="N",
        help="Top/bottom STATE cell types per CyteType label for extremes.csv (default: 10).",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Free-form note recorded as 'notes' in metadata.json.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan the run without R2 downloads or inspection.",
    )
    parser.add_argument(
        "--no-extremes",
        dest="emit_extremes",
        action="store_false",
        help="Skip extremes.csv generation (default: write extremes.csv from summary).",
    )
    args = parser.parse_args()

    if args.from_summary is not None:
        if args.input_prefix is not None or args.cyteonto_prefix is not None:
            parser.error("--from-summary cannot be combined with --input-prefix or --cyteonto-prefix")
        return args

    if args.input_prefix is None or args.cyteonto_prefix is None:
        parser.error("--input-prefix and --cyteonto-prefix are required unless --from-summary is set")
    return args


def main() -> None:
    args = _parse_args()

    if args.from_summary is not None:
        run_from_summary(args)
        return

    log_run_separator(log)
    log.info(
        "new annotation inspection run (dry_run=%s, input=%s, cyteonto=%s, workers=%d, emit_extremes=%s)",
        args.dry_run,
        args.input_prefix,
        args.cyteonto_prefix,
        args.workers,
        args.emit_extremes,
    )
    if args.metadata is not None:
        log.info("run-level --metadata note: %s", args.metadata)

    run_dir = RUN_OUTPUT_ROOT / (f"dry_run_{RUN_TIMESTAMP}" if args.dry_run else RUN_TIMESTAMP)
    run_dir.mkdir(parents=True, exist_ok=True)
    run_csv_path = run_dir / _RUN_CSV_FILENAME
    summary_path = run_dir / _SUMMARY_FILENAME
    metadata_path = run_dir / "metadata.json"
    extremes_path = run_dir / _EXTREMES_FILENAME
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir)
    log.info("run output directory: %s", run_dir)

    work_items = _discover_work_items(args.input_prefix)
    if args.dry_run:
        log.info("Dry-run: would inspect %d accession(s)", len(work_items))
        for srx, input_r2_key, position in work_items:
            cyteonto_r2_key = f"{args.cyteonto_prefix}/{srx}_cyteonto.csv"
            log.info(
                "%s (%s): would inspect (input=%s, cyteonto=%s)",
                srx,
                position,
                input_r2_key,
                cyteonto_r2_key,
            )
            _record_accession(
                run_csv_path,
                srx,
                "dry_run",
                position,
                input_r2_key,
                cyteonto_r2_key,
            )
        log.info("Dry-run complete. Plan written to %s", run_dir)
        return

    summary_lock = threading.Lock()
    summary_frames: list[pd.DataFrame] = []
    processed = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                process_accession,
                srx,
                input_r2_key,
                args.cyteonto_prefix,
            ): (srx, input_r2_key, position)
            for srx, input_r2_key, position in work_items
        }

        for future in as_completed(futures):
            srx, input_r2_key, position = futures[future]
            cyteonto_r2_key = f"{args.cyteonto_prefix}/{srx}_cyteonto.csv"
            pair_df, cyteonto_found, elapsed, exc = future.result()

            if exc is not None:
                failed += 1
                _record_accession(
                    run_csv_path,
                    srx,
                    "failed",
                    position,
                    input_r2_key,
                    cyteonto_r2_key,
                    cyteonto_found=cyteonto_found,
                    duration_seconds=elapsed,
                    error=f"{type(exc).__name__}: {exc}",
                )
                continue

            assert pair_df is not None
            _append_summary_rows(summary_path, pair_df, summary_lock)
            summary_frames.append(pair_df)

            n_cells = int(pair_df["n_cells"].sum())
            _record_accession(
                run_csv_path,
                srx,
                "success",
                position,
                input_r2_key,
                cyteonto_r2_key,
                cyteonto_found=cyteonto_found,
                n_pairs=len(pair_df),
                n_cells=n_cells,
                duration_seconds=elapsed,
            )
            processed += 1
            log.info("%s (%s): done in %.2fs (%d pairs)", srx, position, elapsed, len(pair_df))

    if args.emit_extremes and summary_frames:
        combined = pd.concat(summary_frames, ignore_index=True)
        write_extremes_csv(combined, args.top_n, extremes_path)
        log.info("wrote extremes to %s", extremes_path)

    log.info("Pipeline complete. %d processed, %d failed.", processed, failed)


if __name__ == "__main__":
    main()
