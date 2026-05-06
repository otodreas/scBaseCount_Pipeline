from __future__ import annotations

import argparse
import csv
import datetime
import logging
import sys
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from cluster_validation import ClusterValidationConfig, run_cluster_validation
from cytetype_runner import CyteTypeRunnerConfig, run_cytetype
from gcs import download_from_gcs, gcs_local_path, verify_download
from r2 import fetch_uploaded_r2_keys, upload_to_r2, verify_upload
from shared.logger import configure_file_logger
from shared.repo import REPO_ROOT
from study_context import ExperimentContext, experiment_context_summary

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
_DEFAULT_CONTEXTS_JSONL = REPO_ROOT / "output" / "context" / "contexts.jsonl"
RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
GCS_LOCAL_ROOT = REPO_ROOT / "data"
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "annotation_pipeline"

log = configure_file_logger("annotation_pipeline.log", __name__)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def read_datasets(path: Path) -> list[dict]:
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


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    r2_key: str = "",
    error: str = "",
) -> None:
    write_header = not summary_path.exists()
    with open(summary_path, "a", newline="") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(["position", "srx", "status", "r2_file", "timestamp", "error"])
        writer.writerow([position, srx, status, r2_key, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), error])


def _safe_delete(path: Path) -> None:
    try:
        if path.exists():
            path.unlink()
            log.debug("Deleted %s", path)
    except OSError as exc:
        log.warning("Could not delete %s: %s", path, exc)


def process_accession(
    srx: str,
    gs_uri: str,
    r2_key: str,
    contexts: dict[str, ExperimentContext],
    datasets_path: Path,
) -> Exception | None:
    cfg = ClusterValidationConfig(srxAccession=srx, summaryPath=datasets_path)
    cytetype_cfg = CyteTypeRunnerConfig(srxAccession=srx)
    raw_h5ad = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)
    cytetype_h5ad = cytetype_cfg.outputDir / f"{srx}_cytetype_annotated.h5ad"

    downloaded = False
    try:
        if raw_h5ad.exists():
            log.info("%s: found local h5ad at %s, skipping GCS download", srx, raw_h5ad)
        else:
            download_from_gcs(gs_uri, GCS_LOCAL_ROOT)
            downloaded = True
            if not verify_download(gs_uri, GCS_LOCAL_ROOT):
                raise RuntimeError(f"Download verification failed for {srx}: file not found at {raw_h5ad}")

        _, result = run_cluster_validation(cfg)
        if downloaded:
            _safe_delete(raw_h5ad)

        ctx = contexts.get(srx)
        study_context = experiment_context_summary(ctx) if ctx else ""
        if not ctx:
            log.warning("%s: no study context found in contexts.jsonl; proceeding with empty context", srx)

        cytetype_h5ad = run_cytetype(cytetype_cfg, result.adataPath, result.mergedKey, study_context)
        _safe_delete(result.adataPath)

        upload_to_r2(cytetype_h5ad, r2_key)
        verify_upload(r2_key)
        _safe_delete(cytetype_h5ad)
        log.info("%s: done", srx)
        return None

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        if downloaded:
            _safe_delete(raw_h5ad)
        _safe_delete(cfg.outputDir / f"{srx}_clustered.h5ad")
        _safe_delete(cytetype_h5ad)
        return exc


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the CyteType annotation pipeline.")
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
        default=_DEFAULT_CONTEXTS_JSONL,
        metavar="PATH",
        help=f"Path to contexts JSONL (default: {_DEFAULT_CONTEXTS_JSONL})",
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        default=f"annotation_pipeline_{RUN_TIMESTAMP}",
        metavar="PREFIX",
        help=f"R2 prefix (default: {f"annotation_pipeline_{RUN_TIMESTAMP}"})",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    datasets = read_datasets(args.datasets)
    log.info("Loaded %d accession(s) from %s", len(datasets), args.datasets)
    datasets = datasets.sort_values("obs_count").iloc[0:1]  # TODO: remove this
    uploaded = fetch_uploaded_r2_keys()
    contexts = load_contexts(args.contexts)

    RUN_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = RUN_OUTPUT_DIR / f"run_{args.r2_prefix}.csv"
    log.info("Writing run summary to %s", summary_path)

    total = len(datasets)
    skipped = 0
    for n, (_, row) in enumerate(datasets.iterrows(), start=1):
        srx = row["srx_accession"]
        gs_uri = row["file_path"]
        r2_key = f"{args.r2_prefix}/{srx}_annotated.h5ad"
        position = f"{n}/{total}"

        if r2_key in uploaded:
            log.info("%s: already uploaded, skipping", srx)
            _append_summary_row(summary_path, srx, "skipped", position, r2_key)
            skipped += 1
            continue

        log.info("%s (%s): starting", srx, position)
        exc = process_accession(srx, gs_uri, r2_key, contexts, args.datasets)
        if exc is not None:
            log.warning("%s: failed, skipping", srx)
            _append_summary_row(summary_path, srx, "failed", position, error=f"{type(exc).__name__}: {exc}")
            skipped += 1
        else:
            _append_summary_row(summary_path, srx, "success", position, r2_key)

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)


if __name__ == "__main__":
    main()
