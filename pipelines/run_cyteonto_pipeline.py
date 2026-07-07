from __future__ import annotations

import argparse
import datetime
import json
import re
import shutil
import time
from pathlib import Path

import pandas as pd
from cyteonto import CyteOntoConfig, run_cyteonto
from dotenv import load_dotenv
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from storage import download_from_r2, fetch_uploaded_r2_keys, r2_key_exists, upload_to_r2, verify_upload

load_dotenv()

RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
ANNOTATED_DOWNLOAD_ROOT = REPO_ROOT / "data" / "annotated_h5ad"
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "cyteonto_pipeline"
_CYTEONTO_RUNS_DIR = REPO_ROOT / "output" / "cyteonto" / "runs"
_LOG_FILENAME = "cyteonto_pipeline.log"
_RUN_CSV_FILENAME = "run.csv"
_RESULTS_DIRNAME = "results"
_ANNOTATED_H5AD_PATTERN = re.compile(r"^(?P<srx>SRX\d+)_annotated\.h5ad$")

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
    run_dir: Path,
    results_dir: Path,
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "run_dir": rel_to_repo(run_dir),
        "run_csv": rel_to_repo(run_dir / _RUN_CSV_FILENAME),
        "results_dir": rel_to_repo(results_dir),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "dry_run": bool(args.dry_run),
        "input_prefix": args.input_prefix,
        "r2_prefix": args.r2_prefix,
        "poll_interval_s": args.poll_interval_s,
        "poll_timeout_s": args.poll_timeout_s,
        "min_interval_seconds": args.min_interval,
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


_CSV_COLUMNS = [
    "position",
    "srx",
    "status",
    "input_r2_file",
    "output_r2_file",
    "run_id",
    "local_csv",
    "timestamp",
    "duration_seconds",
    "error",
]


def _record_accession(
    csv_path: Path,
    srx: str,
    status: str,
    position: str = "",
    input_r2_key: str = "",
    output_r2_key: str = "",
    run_id: str = "",
    local_csv: str = "",
    duration_seconds: float | None = None,
    error: str = "",
) -> None:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    duration_str = "" if duration_seconds is None else f"{duration_seconds:.2f}"
    append_csv_row(
        csv_path,
        _CSV_COLUMNS,
        [
            position,
            srx,
            status,
            input_r2_key,
            output_r2_key,
            run_id,
            local_csv,
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


def _discover_work_items(
    input_prefix: str,
    output_prefix: str,
    uploaded_keys: set[str],
    csv_summary_path: Path | None = None,
) -> tuple[list[tuple[str, str, str, str]], int]:
    work_items: list[tuple[str, str, str, str]] = []
    skipped = 0
    input_keys = sorted(key for key in uploaded_keys if key.startswith(f"{input_prefix}/"))
    total = len(input_keys)

    for n, input_r2_key in enumerate(input_keys, start=1):
        srx = _srx_from_input_key(input_prefix, input_r2_key)
        position = f"{n}/{total}"
        if srx is None:
            log.warning("Skipping unrecognized input key: %s", input_r2_key)
            skipped += 1
            continue

        output_r2_key = f"{output_prefix}/{srx}_cyteonto.csv"
        if output_r2_key in uploaded_keys:
            log.info("%s: already uploaded, skipping", srx)
            if csv_summary_path is not None:
                _record_accession(
                    csv_summary_path,
                    srx,
                    "skipped",
                    position,
                    input_r2_key,
                    output_r2_key,
                )
            skipped += 1
            continue

        work_items.append((srx, input_r2_key, output_r2_key, position))

    return work_items, skipped


def process_accession(
    srx: str,
    input_r2_key: str,
    output_r2_key: str,
    results_dir: Path,
    poll_interval_s: int,
    poll_timeout_s: int,
) -> tuple[Exception | None, str, str, pd.DataFrame | None]:
    """Run CyteOnto for one accession; returns (error_or_None, run_id, local_csv_rel, dataframe_or_None)."""
    local_h5ad = ANNOTATED_DOWNLOAD_ROOT / f"{srx}_annotated.h5ad"
    run_id = ""
    local_csv_rel = ""
    try:
        if not r2_key_exists(input_r2_key):
            raise FileNotFoundError(f"input annotated h5ad not found in R2 at {input_r2_key}")

        log.info("%s: downloading annotated h5ad from R2 (%s)", srx, input_r2_key)
        download_from_r2(input_r2_key, local_h5ad)

        cfg = CyteOntoConfig(
            h5adPath=local_h5ad,
            pollIntervalS=poll_interval_s,
            pollTimeoutS=poll_timeout_s,
        )
        df = run_cyteonto(cfg)
        if df is None:
            return KeyboardInterrupt("CyteOnto polling interrupted"), run_id, local_csv_rel, None

        run_id = str(df["run_id"].iloc[0])
        run_id_csv = _CYTEONTO_RUNS_DIR / f"{run_id}.csv"
        if not run_id_csv.exists():
            raise FileNotFoundError(f"expected CyteOnto result CSV not found at {run_id_csv}")

        results_dir.mkdir(parents=True, exist_ok=True)
        local_srx_csv = results_dir / f"{srx}_cyteonto.csv"
        shutil.move(str(run_id_csv), str(local_srx_csv))
        local_csv_rel = rel_to_repo(local_srx_csv)

        upload_to_r2(local_srx_csv, output_r2_key)
        verify_upload(output_r2_key)
        safe_delete(local_h5ad, log)

        log.info("%s: done (run_id=%s)", srx, run_id)
        return None, run_id, local_csv_rel, df

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        safe_delete(local_h5ad, log)
        return exc, run_id, local_csv_rel, None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run CyteOnto on annotated h5ads stored under an R2 prefix, one accession at a time."
    )
    parser.add_argument(
        "--input-prefix",
        type=str,
        required=True,
        metavar="PREFIX",
        help="R2 prefix containing annotated h5ads (e.g. cytetype_pipeline_20260522_175813).",
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        default=f"cyteonto_pipeline_{RUN_TIMESTAMP}",
        metavar="PREFIX",
        help="R2 prefix for CyteOnto result CSVs (default: cyteonto_pipeline_{RUN_TIMESTAMP})",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Free-form note recorded as 'notes' in the run's metadata.json.",
    )
    parser.add_argument(
        "--poll-interval-s",
        type=int,
        default=10,
        metavar="SECONDS",
        help="Seconds between CyteOnto result polls (default: 10).",
    )
    parser.add_argument(
        "--poll-timeout-s",
        type=int,
        default=3600,
        metavar="SECONDS",
        help="Total seconds to wait for a CyteOnto run before raising TimeoutError (default: 3600).",
    )
    parser.add_argument(
        "--min-interval",
        type=int,
        default=0,
        metavar="SECONDS",
        help=(
            "Minimum seconds between the START of consecutive accession runs "
            "(default: 0, no spacing). If a run takes longer than this, the next run starts "
            "immediately when it finishes. Accessions always run serially."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Plan the run without performing R2 downloads, CyteOnto API calls, or R2 uploads. "
            "Writes run.csv and metadata.json under output/cyteonto_pipeline/dry_run_{timestamp}/."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info(
        "new cyteonto pipeline run started (dry_run=%s, input prefix: %s, r2 prefix: %s, min_interval: %ds)",
        args.dry_run,
        args.input_prefix,
        args.r2_prefix,
        args.min_interval,
    )
    if args.metadata is not None:
        log.info("run-level --metadata note (recorded as 'notes' in metadata.json): %s", args.metadata)

    uploaded = fetch_uploaded_r2_keys(prefix=args.input_prefix) | fetch_uploaded_r2_keys(prefix=args.r2_prefix)

    run_dir = RUN_OUTPUT_DIR / (f"dry_run_{RUN_TIMESTAMP}" if args.dry_run else RUN_TIMESTAMP)
    results_dir = run_dir / _RESULTS_DIRNAME
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / _RUN_CSV_FILENAME
    metadata_path = run_dir / "metadata.json"
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir, results_dir)
    log.info("run summary output directory: %s", run_dir)

    work_items, skipped = _discover_work_items(
        args.input_prefix,
        args.r2_prefix,
        uploaded,
        csv_summary_path=csv_summary_path,
    )

    if args.dry_run:
        log.info("Dry-run: no R2 downloads, CyteOnto API calls, or R2 uploads will be performed")
        if args.min_interval > 0:
            log.info(
                "Would run %d accession(s) serially with at least %d seconds between run starts",
                len(work_items),
                args.min_interval,
            )
        else:
            log.info("Would run %d accession(s) serially", len(work_items))
        for srx, input_r2_key, output_r2_key, position in work_items:
            log.info(
                "%s (%s): would submit to CyteOnto (input=%s, output=%s)",
                srx,
                position,
                input_r2_key,
                output_r2_key,
            )
            _record_accession(
                csv_summary_path,
                srx,
                "dry_run",
                position,
                input_r2_key,
                output_r2_key,
            )
        log.info(
            "Dry-run complete. %d would-run, %d skipped or unrecognized. Plan written to %s",
            len(work_items),
            skipped,
            run_dir,
        )
        return

    if args.min_interval > 0:
        log.info(
            "Running %d accession(s) serially with at least %d seconds between run starts",
            len(work_items),
            args.min_interval,
        )
    else:
        log.info("Running %d accession(s) serially", len(work_items))

    processed = 0
    interrupted = False

    for i, (srx, input_r2_key, output_r2_key, position) in enumerate(work_items):
        run_start = time.monotonic()
        log.info("%s (%s): run starting", srx, position)
        exc, run_id, local_csv_rel, _df = process_accession(
            srx,
            input_r2_key,
            output_r2_key,
            results_dir,
            args.poll_interval_s,
            args.poll_timeout_s,
        )
        elapsed = time.monotonic() - run_start
        log.info("%s (%s): run finished in %.2fs", srx, position, elapsed)

        if exc is not None and isinstance(exc, KeyboardInterrupt):
            log.warning("%s: interrupted", srx)
            _record_accession(
                csv_summary_path,
                srx,
                "interrupted",
                position,
                input_r2_key,
                output_r2_key,
                run_id=run_id,
                local_csv=local_csv_rel,
                duration_seconds=elapsed,
                error=str(exc),
            )
            interrupted = True
            break

        if exc is not None:
            log.warning("%s: failed", srx)
            _record_accession(
                csv_summary_path,
                srx,
                "failed",
                position,
                input_r2_key,
                output_r2_key,
                run_id=run_id,
                local_csv=local_csv_rel,
                duration_seconds=elapsed,
                error=f"{type(exc).__name__}: {exc}",
            )
            skipped += 1
        else:
            _record_accession(
                csv_summary_path,
                srx,
                "success",
                position,
                input_r2_key,
                output_r2_key,
                run_id=run_id,
                local_csv=local_csv_rel,
                duration_seconds=elapsed,
            )
            processed += 1

        if args.min_interval > 0 and i < len(work_items) - 1 and not interrupted:
            next_srx = work_items[i + 1][0]
            remaining = len(work_items) - (i + 1)
            wait = max(0.0, args.min_interval - elapsed)
            if wait > 0:
                log.info(
                    "Sleeping %.0fs before next run (%s); %s took %.0fs of %ds spacing window; %d accession(s) remaining",
                    wait,
                    next_srx,
                    srx,
                    elapsed,
                    args.min_interval,
                    remaining,
                )
                time.sleep(wait)
                log.info("Resuming after %.0f second sleep", wait)
            else:
                log.info(
                    "%s took %.0fs, exceeding %ds spacing window by %.0fs; starting next run (%s) immediately; %d accession(s) remaining",
                    srx,
                    elapsed,
                    args.min_interval,
                    elapsed - args.min_interval,
                    next_srx,
                    remaining,
                )

    if interrupted:
        remaining = len(work_items) - processed - 1
        log.info(
            "Pipeline interrupted. %d skipped, %d processed, %d remaining.",
            skipped,
            processed,
            remaining,
        )
    else:
        log.info("Pipeline complete. %d skipped, %d processed.", skipped, processed)


if __name__ == "__main__":
    main()
