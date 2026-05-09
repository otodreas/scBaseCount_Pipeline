from __future__ import annotations

import argparse
import datetime
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from cluster_validation import ClusterValidationConfig, compute_cell_type_entropy_row, run_cluster_validation
from cytetype_runner import CyteTypeRunnerConfig, run_cytetype
from gcs import download_from_gcs, gcs_local_path, verify_download
from r2 import download_from_r2, fetch_uploaded_r2_keys, gcs_uri_to_r2_raw_key, r2_key_exists, r2_object_md5, upload_to_r2, verify_upload
from r2.client import _local_md5_b64
from shared.csv_writer import append_csv_row, append_jsonl_row
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


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "r2_prefix": args.r2_prefix,
        "datasets_path": str(args.datasets),
        "contexts_path": str(args.contexts),
        "skip_cytetype": args.skip_cytetype,
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


_CSV_COLUMNS = ["position", "srx", "status", "r2_file", "timestamp", "error"]


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    r2_key: str = "",
    error: str = "",
) -> None:
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [position, srx, status, r2_key, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), error],
    )


def _safe_delete(path: Path) -> None:
    try:
        if path.exists():
            path.unlink()
            log.debug("Deleted %s", path)
    except OSError as exc:
        log.warning("Could not delete %s: %s", path, exc)


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
    contexts: dict[str, ExperimentContext],
    datasets_path: Path,
    skip_cytetype: bool = False,
) -> tuple[dict[str, float] | None, Exception | None]:
    cfg = ClusterValidationConfig(srxAccession=srx, summaryPath=datasets_path)
    raw_h5ad = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)

    downloaded = False
    cytetype_h5ad: Path | None = None
    try:
        if raw_h5ad.exists():
            log.info("%s: found local h5ad at %s, skipping download", srx, raw_h5ad)
        else:
            _fetch_raw_h5ad(srx, gs_uri, raw_h5ad)
            downloaded = True

        adata, result = run_cluster_validation(cfg)
        entropy_row = compute_cell_type_entropy_row(adata, result.mergedKey)
        if downloaded:
            _safe_delete(raw_h5ad)

        if skip_cytetype:
            upload_to_r2(result.adataPath, r2_key)
            verify_upload(r2_key)
            _safe_delete(result.adataPath)
        else:
            cytetype_cfg = CyteTypeRunnerConfig(srxAccession=srx)
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
        return entropy_row, None

    except Exception as exc:
        log.exception("%s: pipeline failed", srx)
        _safe_delete(raw_h5ad)
        _safe_delete(cfg.outputDir / f"{srx}_clustered.h5ad")
        if cytetype_h5ad is not None:
            _safe_delete(cytetype_h5ad)
        return None, exc


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
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Write a metadata JSON file next to the run CSV with this note.",
    )
    parser.add_argument(
        "--skip-cytetype",
        action="store_true",
        default=False,
        help="Skip the CyteType step and upload the clustering output directly.",
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

    for handler in log.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write("\n")
    log.info("new annotation pipeline run started (r2 prefix: %s)", args.r2_prefix)

    datasets = read_datasets(args.datasets)
    log.info("Loaded %d accession(s) from %s", len(datasets), args.datasets)
    # datasets = datasets.sort_values("obs_count").iloc[0:1]  # TODO: remove this
    uploaded = fetch_uploaded_r2_keys()
    contexts = load_contexts(args.contexts)

    run_dir = RUN_OUTPUT_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_summary_path = run_dir / "run.csv"
    metadata_path = run_dir / "metadata.json"
    entropy_jsonl_path = run_dir / "entropy_matrix.jsonl"
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP)
    log.info("run summary output directory: %s", run_dir)

    total = len(datasets)
    skipped = 0
    r2_suffix = "_clustered.h5ad" if args.skip_cytetype else "_annotated.h5ad"

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
            pool.submit(process_accession, srx, gs_uri, r2_key, contexts, args.datasets, args.skip_cytetype): (srx, r2_key, position)
            for srx, gs_uri, r2_key, position in work_items
        }
        for future in as_completed(futures):
            srx, r2_key, position = futures[future]
            entropy_row, exc = future.result()
            if exc is not None:
                log.warning("%s: failed", srx)
                _append_summary_row(csv_summary_path, srx, "failed", position, error=f"{type(exc).__name__}: {exc}")
                skipped += 1
            else:
                append_jsonl_row(entropy_jsonl_path, {"srx": srx, "cell_types": entropy_row})
                _append_summary_row(csv_summary_path, srx, "success", position, r2_key)

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)

    if entropy_jsonl_path.exists():
        rows = [json.loads(line) for line in entropy_jsonl_path.read_text().splitlines() if line.strip()]
        entropy_df = pd.DataFrame({r["srx"]: r["cell_types"] for r in rows}).T
        entropy_df.index.name = "srx"
        entropy_csv_path = run_dir / "entropy_matrix.csv"
        entropy_df.to_csv(entropy_csv_path)
        log.info("Entropy matrix written to %s", entropy_csv_path)
        entropy_jsonl_path.unlink()
        log.debug("Deleted %s", entropy_jsonl_path)


if __name__ == "__main__":
    main()
