from __future__ import annotations

import argparse
import datetime
import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import scanpy as sc
from celltypist_runner import CellTypistRunnerConfig, annotate_celltypist
from cluster_validation import ClusterValidationConfig, run_cluster_validation_on_adata
from cyteonto import CyteOntoConfig, run_cyteonto
from cytetype_runner import CyteTypeRunnerConfig, require_api_key, run_cytetype
from dotenv import load_dotenv
from shared.csv_writer import append_csv_row
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from study_context import ExperimentContext, experiment_context_summary

load_dotenv()

CELLTYPIST_COL = "predicted_labels"
CXG_COL = "cell_type"
CYTETYPE_COL = "cytetype_annotation_leiden_merged"
AUTHOR_COL = CELLTYPIST_COL
ALGORITHM_COLS = {"cxg": CXG_COL, "cytetype": CYTETYPE_COL}

OUTPUT_ROOT = REPO_ROOT / "output" / "celltypist_vs_cxg"
DATA_DIR = OUTPUT_ROOT / "data"
CYTEONTO_RESULTS_DIR = OUTPUT_ROOT / "cyteonto_results"
CYTEONTO_RUNS_DIR = OUTPUT_ROOT / "cyteonto_runs"
CYTEONTO_PAYLOADS_DIR = OUTPUT_ROOT / "cyteonto_payloads"
FIGS_DIR = OUTPUT_ROOT / "figs"
CONTEXTS_JSONL = REPO_ROOT / "output" / "context" / "contexts.jsonl"
LOCAL_H5AD_ROOT = REPO_ROOT / "data" / "scbasecount" / "2026-01-12" / "h5ad" / "GeneFull" / "Homo_sapiens"
RUNS_DIR = OUTPUT_ROOT / "runs"

RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
_LOG_FILENAME = "celltypist_vs_cxg_pipeline.log"
_RUN_CSV_FILENAME = "run.csv"

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()

_CSV_COLUMNS = [
    "position",
    "srx",
    "status",
    "clustered_h5ad",
    "annotated_h5ad",
    "cyteonto_csv",
    "timestamp",
    "duration_seconds",
    "error",
]


@dataclass(frozen=True)
class CachePaths:
    srx: str
    clustered: Path
    annotated: Path
    cyteonto_csv: Path

    @classmethod
    def for_srx(cls, srx: str) -> CachePaths:
        return cls(
            srx=srx,
            clustered=DATA_DIR / f"{srx}_celltypist_prior_clustered.h5ad",
            annotated=DATA_DIR / f"{srx}_cytetype_annotated.h5ad",
            cyteonto_csv=CYTEONTO_RESULTS_DIR / f"{srx}_cyteonto.csv",
        )


def load_context(accession: str, contexts_path: Path) -> str:
    """Return the study context summary for an accession, or an empty string if not found."""
    if not contexts_path.is_file():
        return ""
    for line in contexts_path.read_text().splitlines():
        if not line.strip():
            continue
        ctx = ExperimentContext.model_validate_json(line)
        if ctx.accession == accession:
            return experiment_context_summary(ctx)
    return ""


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
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "model_name": args.model_name,
        "force": bool(args.force),
        "poll_interval_s": args.poll_interval_s,
        "poll_timeout_s": args.poll_timeout_s,
        "srx_accessions": args.srx,
        "output_root": rel_to_repo(OUTPUT_ROOT),
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


def _record_accession(
    csv_path: Path,
    srx: str,
    status: str,
    paths: CachePaths,
    position: str = "",
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
            rel_to_repo(paths.clustered) if paths.clustered.is_file() else "",
            rel_to_repo(paths.annotated) if paths.annotated.is_file() else "",
            rel_to_repo(paths.cyteonto_csv) if paths.cyteonto_csv.is_file() else "",
            timestamp,
            duration_str,
            error,
        ],
    )


def _run_celltypist_and_cluster(srx: str, model_name: str, paths: CachePaths) -> None:
    raw_path = LOCAL_H5AD_ROOT / f"{srx}.h5ad"
    if not raw_path.is_file():
        raise FileNotFoundError(f"raw h5ad not found at {raw_path}")

    log.info("%s: loading raw h5ad from %s", srx, raw_path)
    adata_raw = sc.read(raw_path)

    celltypist_cfg = CellTypistRunnerConfig(modelName=model_name)
    adata_raw = annotate_celltypist(adata_raw, celltypist_cfg)

    cluster_cfg = ClusterValidationConfig(
        weakPriorKey=CELLTYPIST_COL,
        runLabel=f"{srx}_celltypist_prior",
        outputDir=DATA_DIR,
        figsDir=FIGS_DIR,
    )
    adata_clustered, cluster_result = run_cluster_validation_on_adata(
        adata_raw.copy(),
        cluster_cfg,
        srx,
        plot=False,
    )
    paths.clustered.parent.mkdir(parents=True, exist_ok=True)
    adata_clustered.write(paths.clustered)
    log.info(
        "%s: clustered h5ad written  resolution=%s  n_clusters=%d",
        srx,
        cluster_result.selectedResolution,
        cluster_result.nClustersPostMerge,
    )


def _run_cytetype(srx: str, paths: CachePaths) -> None:
    if not paths.clustered.is_file():
        raise FileNotFoundError(f"clustered h5ad not found at {paths.clustered}")

    study_context = load_context(srx, CONTEXTS_JSONL)
    if not study_context:
        log.warning("%s: no study context found in contexts.jsonl; proceeding with empty context", srx)

    cytetype_cfg = CyteTypeRunnerConfig(srxAccession=srx, outputDir=DATA_DIR)
    cytetype_result = run_cytetype(
        cytetype_cfg,
        paths.clustered,
        group_key="leiden_merged",
        study_context=study_context,
    )
    if cytetype_result.outputPath != paths.annotated:
        raise RuntimeError(f"{srx}: expected annotated h5ad at {paths.annotated}, got {cytetype_result.outputPath}")
    log.info("%s: annotated h5ad written to %s", srx, paths.annotated)


def _run_cyteonto(srx: str, paths: CachePaths, poll_interval_s: int, poll_timeout_s: int) -> bool:
    """Submit CyteOnto and cache the result CSV; return False if polling was interrupted."""
    if not paths.annotated.is_file():
        raise FileNotFoundError(f"annotated h5ad not found at {paths.annotated}")

    cyteonto_cfg = CyteOntoConfig(
        h5adPath=paths.annotated,
        authorCol=AUTHOR_COL,
        algorithmCols=ALGORITHM_COLS,
        runsDir=CYTEONTO_RUNS_DIR,
        payloadDir=CYTEONTO_PAYLOADS_DIR,
        pollIntervalS=poll_interval_s,
        pollTimeoutS=poll_timeout_s,
    )
    cyteonto_df = run_cyteonto(cyteonto_cfg)
    if cyteonto_df is None:
        log.warning(
            "%s: CyteOnto polling interrupted; call cyteonto.check_pending_runs() and rerun",
            srx,
        )
        return False

    run_id = str(cyteonto_df["run_id"].iloc[0])
    run_id_csv = CYTEONTO_RUNS_DIR / f"{run_id}.csv"
    if not run_id_csv.is_file():
        raise FileNotFoundError(f"expected CyteOnto result CSV not found at {run_id_csv}")

    paths.cyteonto_csv.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(run_id_csv), str(paths.cyteonto_csv))
    log.info("%s: cyteonto CSV written to %s (run_id=%s)", srx, paths.cyteonto_csv, run_id)
    return True


def process_accession(
    srx: str,
    model_name: str,
    *,
    force: bool,
    poll_interval_s: int,
    poll_timeout_s: int,
) -> str:
    """Run the pipeline for one accession and return its completion status."""
    paths = CachePaths.for_srx(srx)

    if force:
        for path in (paths.cyteonto_csv, paths.annotated, paths.clustered):
            if path.is_file():
                path.unlink()
                log.info("%s: removed cached file %s (--force)", srx, path)

    if paths.cyteonto_csv.is_file() and paths.annotated.is_file():
        log.info("%s: cyteonto CSV and annotated h5ad already cached, skipping", srx)
        return "skipped"

    if not paths.clustered.is_file():
        _run_celltypist_and_cluster(srx, model_name, paths)

    if not paths.annotated.is_file():
        _run_cytetype(srx, paths)

    if not paths.cyteonto_csv.is_file():
        if not _run_cyteonto(srx, paths, poll_interval_s, poll_timeout_s):
            return "interrupted"

    return "success"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CellTypist -> cluster validation -> CyteType -> CyteOnto for one or more "
            "accessions, caching outputs under output/celltypist_vs_cxg/."
        )
    )
    parser.add_argument(
        "--srx",
        nargs="+",
        default=["SRX12708356"],
        metavar="ACCESSION",
        help="One or more SRX accessions (default: SRX12708356)",
    )
    parser.add_argument(
        "--model-name",
        type=str,
        default="Adult_COVID19_PBMC.pkl",
        metavar="MODEL",
        help="CellTypist model name (default: Adult_COVID19_PBMC.pkl)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore cached outputs and recompute from scratch",
    )
    parser.add_argument(
        "--poll-interval-s",
        type=int,
        default=10,
        metavar="SECONDS",
        help="Seconds between CyteOnto result polls (default: 10)",
    )
    parser.add_argument(
        "--poll-timeout-s",
        type=int,
        default=3600,
        metavar="SECONDS",
        help="Total seconds to wait for a CyteOnto run (default: 3600)",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Free-form note recorded as 'notes' in the run metadata.json",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info(
        "celltypist vs cxg pipeline started (srx=%s, model=%s, force=%s)",
        args.srx,
        args.model_name,
        args.force,
    )

    require_api_key()

    for directory in (DATA_DIR, CYTEONTO_RESULTS_DIR, CYTEONTO_RUNS_DIR, CYTEONTO_PAYLOADS_DIR, FIGS_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    run_dir = RUNS_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / _RUN_CSV_FILENAME
    metadata_path = run_dir / "metadata.json"
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, run_dir)
    log.info("run summary output directory: %s", run_dir)

    total = len(args.srx)
    skipped = 0
    failed = 0

    for i, srx in enumerate(args.srx, start=1):
        position = f"{i}/{total}"
        run_start = time.monotonic()
        log.info("%s (%s): run starting", srx, position)
        try:
            status = process_accession(
                srx,
                args.model_name,
                force=args.force,
                poll_interval_s=args.poll_interval_s,
                poll_timeout_s=args.poll_timeout_s,
            )
            elapsed = time.monotonic() - run_start
            paths = CachePaths.for_srx(srx)
            log.info("%s (%s): %s in %.2fs", srx, position, status, elapsed)
            _record_accession(
                csv_summary_path,
                srx,
                status,
                paths,
                position,
                duration_seconds=elapsed,
            )
            if status == "skipped":
                skipped += 1
            elif status == "interrupted":
                log.info("Pipeline interrupted after %s", srx)
                break
        except KeyboardInterrupt as exc:
            elapsed = time.monotonic() - run_start
            paths = CachePaths.for_srx(srx)
            log.warning("%s (%s): interrupted after %.2fs", srx, position, elapsed)
            _record_accession(
                csv_summary_path,
                srx,
                "interrupted",
                paths,
                position,
                duration_seconds=elapsed,
                error=f"{type(exc).__name__}: {exc}",
            )
            log.info("Pipeline interrupted after %s", srx)
            break
        except Exception as exc:
            elapsed = time.monotonic() - run_start
            paths = CachePaths.for_srx(srx)
            log.exception("%s (%s): pipeline failed after %.2fs", srx, position, elapsed)
            _record_accession(
                csv_summary_path,
                srx,
                "failed",
                paths,
                position,
                duration_seconds=elapsed,
                error=f"{type(exc).__name__}: {exc}",
            )
            failed += 1

    processed = total - skipped - failed
    log.info("Pipeline complete. %d skipped, %d processed, %d failed.", skipped, processed, failed)


if __name__ == "__main__":
    main()
