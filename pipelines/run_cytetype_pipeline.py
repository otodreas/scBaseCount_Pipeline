from __future__ import annotations

import argparse
import datetime
import json
import time
from pathlib import Path

import pandas as pd
from cluster_validation.merge import MERGED_CLUSTER_KEY
from cytetype_runner import CyteTypeRunnerConfig, require_api_key, run_cytetype
from dotenv import load_dotenv
from r2 import download_from_r2, fetch_uploaded_r2_keys, r2_key_exists, upload_to_r2, verify_upload
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from study_context import ExperimentContext, experiment_context_summary

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
_DEFAULT_CONTEXTS_JSONL = REPO_ROOT / "output" / "context" / "contexts.jsonl"
RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
CLUSTERED_DOWNLOAD_ROOT = REPO_ROOT / "data" / "clustered_h5ad"
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "cytetype_pipeline"
_LOG_FILENAME = "cytetype_pipeline.log"
_RUN_CSV_FILENAME = "run.csv"
_JOB_DETAILS_FILENAME = "job_details.csv"

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


def read_datasets(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def load_contexts(path: Path) -> dict[str, ExperimentContext]:
    if not path.exists():
        log.warning("contexts.jsonl not found at %s; study context will be empty for all samples", path)
        return {}
    contexts: dict[str, ExperimentContext] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ctx = ExperimentContext.model_validate_json(line)
            contexts[ctx.accession] = ctx
    return contexts


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
    run_dir: Path,
    metadata_columns: list[str],
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "run_dir": rel_to_repo(run_dir),
        "run_csv": rel_to_repo(run_dir / _RUN_CSV_FILENAME),
        "job_details_csv": rel_to_repo(run_dir / _JOB_DETAILS_FILENAME),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "dry_run": bool(args.dry_run),
        "clustering_prefix": args.clustering_prefix,
        "r2_prefix": args.r2_prefix,
        "datasets_path": str(args.datasets),
        "contexts_path": str(args.contexts) if args.contexts is not None else str(_DEFAULT_CONTEXTS_JSONL),
        "min_interval_seconds": args.min_interval,
        "cytetype_metadata_columns": metadata_columns,
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
    "timestamp",
    "duration_seconds",
    "error",
]
_JOB_DETAILS_COLUMNS = [
    "srx",
    "status",
    "position",
    "input_r2_file",
    "output_r2_file",
    "job_id",
    "report_url",
    "api_url",
    "timestamp",
    "duration_seconds",
    "error",
]


def _record_accession(
    csv_path: Path,
    job_details_path: Path,
    srx: str,
    status: str,
    position: str = "",
    input_r2_key: str = "",
    output_r2_key: str = "",
    job_id: str = "",
    report_url: str = "",
    api_url: str = "",
    duration_seconds: float | None = None,
    error: str = "",
) -> None:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    duration_str = "" if duration_seconds is None else f"{duration_seconds:.2f}"
    append_csv_row(
        csv_path,
        _CSV_COLUMNS,
        [position, srx, status, input_r2_key, output_r2_key, timestamp, duration_str, error],
    )
    append_csv_row(
        job_details_path,
        _JOB_DETAILS_COLUMNS,
        [
            srx,
            status,
            position,
            input_r2_key,
            output_r2_key,
            job_id,
            report_url,
            api_url,
            timestamp,
            duration_str,
            error,
        ],
    )


def _row_to_metadata(row: pd.Series) -> dict[str, str]:
    """Convert a datasets-CSV row into a CyteType metadata dict, stringifying values and skipping NaN."""
    return {str(col): str(val) for col, val in row.items() if not pd.isna(val)}


def process_accession(
    srx: str,
    input_r2_key: str,
    output_r2_key: str,
    contexts: dict[str, ExperimentContext],
    metadata: dict[str, str],
) -> tuple[Exception | None, str, str, str]:
    """Run CyteType for one accession; returns (error_or_None, job_id, report_url, api_url)."""
    local_clustered = CLUSTERED_DOWNLOAD_ROOT / f"{srx}_clustered.h5ad"
    cytetype_h5ad: Path | None = None
    job_id = ""
    report_url = ""
    api_url = ""
    try:
        if not r2_key_exists(input_r2_key):
            raise FileNotFoundError(f"input clustered h5ad not found in R2 at {input_r2_key}")

        log.info("%s: downloading clustered h5ad from R2 (%s)", srx, input_r2_key)
        download_from_r2(input_r2_key, local_clustered)

        cytetype_cfg = CyteTypeRunnerConfig(srxAccession=srx)
        ctx = contexts.get(srx)
        study_context = experiment_context_summary(ctx) if ctx else ""
        if not ctx:
            log.warning("%s: no study context found in contexts.jsonl; proceeding with empty context", srx)

        result = run_cytetype(
            cytetype_cfg, local_clustered, MERGED_CLUSTER_KEY, study_context, metadata=metadata
        )
        cytetype_h5ad = result.outputPath
        job_id = result.jobId
        report_url = result.reportUrl
        api_url = result.apiUrl
        if report_url:
            log.info("%s: report URL: %s", srx, report_url)
        else:
            log.warning("%s: report URL not found in adata.uns['cytetype_jobDetails']", srx)
        safe_delete(local_clustered, log)

        upload_to_r2(cytetype_h5ad, output_r2_key)
        verify_upload(output_r2_key)
        safe_delete(cytetype_h5ad, log)

        log.info("%s: done", srx)
        return None, job_id, report_url, api_url

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        safe_delete(local_clustered, log)
        if cytetype_h5ad is not None:
            safe_delete(cytetype_h5ad, log)
        return exc, job_id, report_url, api_url


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run CyteType annotation on clustered h5ads stored under an R2 prefix."
    )
    parser.add_argument(
        "--datasets",
        type=Path,
        default=_DEFAULT_DATASETS_CSV,
        metavar="PATH",
        help=f"Path to datasets CSV (default: {_DEFAULT_DATASETS_CSV})",
    )
    parser.add_argument(
        "--contexts",
        type=Path,
        default=None,
        metavar="PATH",
        help=f"Path to contexts JSONL (default: {_DEFAULT_CONTEXTS_JSONL})",
    )
    parser.add_argument(
        "--clustering-prefix",
        type=str,
        required=True,
        metavar="PREFIX",
        help="R2 prefix containing the clustered h5ads to annotate (e.g. clustering_pipeline_20260511_140000).",
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        default=f"cytetype_pipeline_{RUN_TIMESTAMP}",
        metavar="PREFIX",
        help="R2 prefix for annotated outputs (default: cytetype_pipeline_{RUN_TIMESTAMP})",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help=(
            "Free-form note recorded as 'notes' in the run's metadata.json. "
            "Does NOT control the per-accession metadata sent to CyteType; "
            "that is always derived from the --datasets CSV row."
        ),
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
            "Plan the run without performing R2 downloads, CyteType API calls, or R2 uploads. "
            "Writes run.csv, job_details.csv, and metadata.json under output/cytetype_pipeline/dry_run_{timestamp}/."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info(
        "new cytetype pipeline run started (dry_run=%s, clustering prefix: %s, r2 prefix: %s, min_interval: %ds)",
        args.dry_run,
        args.clustering_prefix,
        args.r2_prefix,
        args.min_interval,
    )

    if not args.dry_run:
        require_api_key()

    datasets = read_datasets(args.datasets)
    log.info("Loaded %d accession(s) from %s", len(datasets), args.datasets)
    metadata_columns = [str(c) for c in datasets.columns]
    log.info(
        "cytetype per-accession metadata will be derived from %d CSV column(s): %s",
        len(metadata_columns),
        ", ".join(metadata_columns),
    )
    if args.metadata is not None:
        log.info("run-level --metadata note (recorded as 'notes' in metadata.json): %s", args.metadata)
    uploaded = fetch_uploaded_r2_keys()

    contexts_path = args.contexts if args.contexts is not None else _DEFAULT_CONTEXTS_JSONL
    contexts = load_contexts(contexts_path)

    run_dir = RUN_OUTPUT_DIR / (f"dry_run_{RUN_TIMESTAMP}" if args.dry_run else RUN_TIMESTAMP)
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / _RUN_CSV_FILENAME
    job_details_path = run_dir / _JOB_DETAILS_FILENAME
    metadata_path = run_dir / "metadata.json"
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir, metadata_columns)
    log.info("run summary output directory: %s", run_dir)

    total = len(datasets)
    skipped = 0

    work_items: list[tuple[str, str, str, str, dict[str, str]]] = []
    for n, (_, row) in enumerate(datasets.iterrows(), start=1):
        srx = row["srx_accession"]
        input_r2_key = f"{args.clustering_prefix}/{srx}_clustered.h5ad"
        output_r2_key = f"{args.r2_prefix}/{srx}_annotated.h5ad"
        position = f"{n}/{total}"
        if output_r2_key in uploaded:
            log.info("%s: already uploaded, skipping", srx)
            _record_accession(csv_summary_path, job_details_path, srx, "skipped", position, input_r2_key, output_r2_key)
            skipped += 1
            continue
        if input_r2_key not in uploaded:
            log.warning("%s: clustered input not found at %s, skipping", srx, input_r2_key)
            _record_accession(
                csv_summary_path,
                job_details_path,
                srx,
                "missing_input",
                position,
                input_r2_key,
                output_r2_key,
                error=f"input clustered h5ad not found in R2 at {input_r2_key}",
            )
            skipped += 1
            continue
        work_items.append((srx, input_r2_key, output_r2_key, position, _row_to_metadata(row)))

    if args.dry_run:
        log.info("Dry-run: no R2 downloads, CyteType API calls, or R2 uploads will be performed")
        if args.min_interval > 0:
            log.info(
                "Would run %d accession(s) serially with at least %d seconds between run starts",
                len(work_items),
                args.min_interval,
            )
        else:
            log.info("Would run %d accession(s) serially", len(work_items))
        for srx, input_r2_key, output_r2_key, position, _metadata in work_items:
            log.info(
                "%s (%s): would submit to CyteType (input=%s, output=%s)",
                srx,
                position,
                input_r2_key,
                output_r2_key,
            )
            _record_accession(csv_summary_path, job_details_path, srx, "dry_run", position, input_r2_key, output_r2_key)
        log.info(
            "Dry-run complete. %d would-run, %d skipped or missing input. Plan written to %s",
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

    for i, (srx, input_r2_key, output_r2_key, position, metadata) in enumerate(work_items):
        run_start = time.monotonic()
        log.info("%s (%s): run starting", srx, position)
        exc, job_id, report_url, api_url = process_accession(
            srx, input_r2_key, output_r2_key, contexts, metadata
        )
        elapsed = time.monotonic() - run_start
        log.info("%s (%s): run finished in %.2fs", srx, position, elapsed)
        if exc is not None:
            log.warning("%s: failed", srx)
            _record_accession(
                csv_summary_path,
                job_details_path,
                srx,
                "failed",
                position,
                input_r2_key,
                output_r2_key,
                job_id=job_id,
                report_url=report_url,
                api_url=api_url,
                duration_seconds=elapsed,
                error=f"{type(exc).__name__}: {exc}",
            )
            skipped += 1
        else:
            _record_accession(
                csv_summary_path,
                job_details_path,
                srx,
                "success",
                position,
                input_r2_key,
                output_r2_key,
                job_id=job_id,
                report_url=report_url,
                api_url=api_url,
                duration_seconds=elapsed,
            )
        if args.min_interval > 0 and i < len(work_items) - 1:
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

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)


if __name__ == "__main__":
    main()
