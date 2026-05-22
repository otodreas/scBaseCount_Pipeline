from __future__ import annotations

import argparse
import datetime
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from cluster_validation.merge import MERGED_CLUSTER_KEY
from cytetype_runner import CyteTypeRunnerConfig, run_cytetype
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
        "run_csv": rel_to_repo(run_dir / "run.csv"),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "clustering_prefix": args.clustering_prefix,
        "r2_prefix": args.r2_prefix,
        "datasets_path": str(args.datasets),
        "contexts_path": str(args.contexts) if args.contexts is not None else str(_DEFAULT_CONTEXTS_JSONL),
        "timeout_seconds": args.timeout,
        "cytetype_metadata_columns": metadata_columns,
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


_CSV_COLUMNS = ["position", "srx", "status", "input_r2_file", "output_r2_file", "timestamp", "error"]


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    input_r2_key: str = "",
    output_r2_key: str = "",
    error: str = "",
) -> None:
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [
            position,
            srx,
            status,
            input_r2_key,
            output_r2_key,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
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
) -> Exception | None:
    local_clustered = CLUSTERED_DOWNLOAD_ROOT / f"{srx}_clustered.h5ad"
    cytetype_h5ad: Path | None = None
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

        cytetype_h5ad = run_cytetype(
            cytetype_cfg, local_clustered, MERGED_CLUSTER_KEY, study_context, metadata=metadata
        )
        safe_delete(local_clustered, log)

        upload_to_r2(cytetype_h5ad, output_r2_key)
        verify_upload(output_r2_key)
        safe_delete(cytetype_h5ad, log)

        log.info("%s: done", srx)
        return None

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        safe_delete(local_clustered, log)
        if cytetype_h5ad is not None:
            safe_delete(cytetype_h5ad, log)
        return exc


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
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Number of accessions to process in parallel (default: 1)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=0,
        metavar="SECONDS",
        help="Seconds to sleep between accession runs. When > 0, runs are forced serial (workers is ignored).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info(
        "new cytetype pipeline run started (clustering prefix: %s, r2 prefix: %s, timeout: %ds)",
        args.clustering_prefix,
        args.r2_prefix,
        args.timeout,
    )

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

    run_dir = RUN_OUTPUT_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / "run.csv"
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
            _append_summary_row(csv_summary_path, srx, "skipped", position, input_r2_key, output_r2_key)
            skipped += 1
            continue
        if input_r2_key not in uploaded:
            log.warning("%s: clustered input not found at %s, skipping", srx, input_r2_key)
            _append_summary_row(
                csv_summary_path,
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

    if args.timeout > 0:
        if args.workers != 1:
            log.info("--timeout set; forcing serial execution (ignoring --workers=%d)", args.workers)
        log.info(
            "Running %d accession(s) serially with %d seconds sleep between runs",
            len(work_items),
            args.timeout,
        )
        for i, (srx, input_r2_key, output_r2_key, position, metadata) in enumerate(work_items):
            exc = process_accession(srx, input_r2_key, output_r2_key, contexts, metadata)
            if exc is not None:
                log.warning("%s: failed", srx)
                _append_summary_row(
                    csv_summary_path,
                    srx,
                    "failed",
                    position,
                    input_r2_key,
                    output_r2_key,
                    error=f"{type(exc).__name__}: {exc}",
                )
                skipped += 1
            else:
                _append_summary_row(csv_summary_path, srx, "success", position, input_r2_key, output_r2_key)
            if i < len(work_items) - 1:
                next_srx = work_items[i + 1][0]
                remaining = len(work_items) - (i + 1)
                log.info(
                    "Sleeping %d seconds after %s before next run (%s); %d accession(s) remaining",
                    args.timeout,
                    srx,
                    next_srx,
                    remaining,
                )
                time.sleep(args.timeout)
                log.info("Resuming after %d second sleep", args.timeout)
    else:
        log.info("Submitting %d accession(s) to %d worker(s)", len(work_items), args.workers)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(process_accession, srx, input_r2_key, output_r2_key, contexts, metadata): (
                    srx,
                    input_r2_key,
                    output_r2_key,
                    position,
                )
                for srx, input_r2_key, output_r2_key, position, metadata in work_items
            }
            for future in as_completed(futures):
                srx, input_r2_key, output_r2_key, position = futures[future]
                exc = future.result()
                if exc is not None:
                    log.warning("%s: failed", srx)
                    _append_summary_row(
                        csv_summary_path,
                        srx,
                        "failed",
                        position,
                        input_r2_key,
                        output_r2_key,
                        error=f"{type(exc).__name__}: {exc}",
                    )
                    skipped += 1
                else:
                    _append_summary_row(csv_summary_path, srx, "success", position, input_r2_key, output_r2_key)

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)


if __name__ == "__main__":
    main()
