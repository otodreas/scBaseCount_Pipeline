from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path

from atlas_integration import (
    AtlasIntegrationConfig,
    build_accession_study_map,
    concat_accession_adatas,
    load_accession_h5ad,
    prepare_accession_adata,
    read_datasets_csv,
    run_atlas_integration,
)
from atlas_integration.models import MergeStats
from dotenv import load_dotenv
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from storage import (
    download_from_gcs,
    download_from_r2,
    gcs_local_path,
    gcs_uri_to_r2_raw_key,
    r2_key_exists,
    r2_object_md5,
    upload_to_r2,
    verify_download,
    verify_upload,
)
from storage.transfer import _local_md5_b64

load_dotenv()

_DEFAULT_DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
GCS_LOCAL_ROOT = REPO_ROOT / "data"
RUN_OUTPUT_DIR = REPO_ROOT / "output" / "atlas_pipeline"
_LOG_FILENAME = "atlas_integration_pipeline.log"

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()

_CSV_COLUMNS = [
    "position",
    "srx",
    "status",
    "nCellsKept",
    "study",
    "timestamp",
    "error",
]


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


def _append_merge_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str,
    n_cells_kept: int = 0,
    study: str = "",
    error: str = "",
) -> None:
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [
            position,
            srx,
            status,
            str(n_cells_kept),
            study,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            error,
        ],
    )


def merge_datasets_incrementally(cfg: AtlasIntegrationConfig, run_dir: Path):
    """Download, QC, and concatenate accessions one at a time to limit peak disk use."""
    datasets = read_datasets_csv(cfg.datasetsCsvPath)
    accessions = datasets["srx_accession"].astype(str).tolist()
    study_map = build_accession_study_map(accessions, cfg.contextsPath)
    merge_csv = run_dir / "merge.csv"

    merged_parts = []
    skipped: list[str] = []
    total = len(datasets)

    for n, (_, row) in enumerate(datasets.iterrows(), start=1):
        srx = str(row["srx_accession"])
        gs_uri = str(row["file_path"])
        position = f"{n}/{total}"
        raw_h5ad = gcs_local_path(gs_uri, GCS_LOCAL_ROOT)
        downloaded = False

        try:
            if raw_h5ad.exists():
                log.info("%s: found local h5ad at %s, skipping download", srx, raw_h5ad)
            else:
                _fetch_raw_h5ad(srx, gs_uri, raw_h5ad)
                downloaded = True

            adata = load_accession_h5ad(raw_h5ad, srx, gs_uri)
            prepared = prepare_accession_adata(adata, srx, study_map[srx], cfg)
            if prepared is None:
                skipped.append(srx)
                _append_merge_row(merge_csv, srx, "skipped", position, study=study_map[srx])
                continue

            merged_parts.append(prepared)
            _append_merge_row(
                merge_csv,
                srx,
                "merged",
                position,
                n_cells_kept=prepared.n_obs,
                study=study_map[srx],
            )
        except Exception as exc:
            log.exception("%s: merge step failed", srx)
            skipped.append(srx)
            _append_merge_row(
                merge_csv,
                srx,
                "failed",
                position,
                study=study_map.get(srx, ""),
                error=f"{type(exc).__name__}: {exc}",
            )
        finally:
            if downloaded:
                safe_delete(raw_h5ad, log)

    if not merged_parts:
        raise RuntimeError("No accessions passed QC during incremental merge")

    merged = concat_accession_adatas(merged_parts)
    stats = MergeStats(
        nAccessionsRequested=len(accessions),
        nAccessionsMerged=len(merged_parts),
        nAccessionsSkipped=len(skipped),
        nCellsFinal=merged.n_obs,
        nGenesFinal=merged.n_vars,
        nStudies=int(merged.obs[cfg.batchKey].nunique()),
        skippedAccessions=skipped,
    )
    return merged, stats


def _write_run_metadata(metadata_path: Path, args: argparse.Namespace, run_dir: Path) -> None:
    payload: dict = {
        "run_timestamp": RUN_TIMESTAMP,
        "run_dir": rel_to_repo(run_dir),
        "run_csv": rel_to_repo(run_dir / "run.csv"),
        "merge_csv": rel_to_repo(run_dir / "merge.csv"),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "r2_prefix": args.r2_prefix,
        "datasets_path": str(args.datasets),
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and integrate the lung atlas from raw scBaseCount h5ad files.")
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
        default=REPO_ROOT / "output" / "context" / "contexts.jsonl",
        metavar="PATH",
        help="Path to contexts.jsonl",
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        default=f"atlas_pipeline_{RUN_TIMESTAMP}",
        metavar="PREFIX",
        help="R2 prefix (default: atlas_pipeline_{RUN_TIMESTAMP})",
    )
    parser.add_argument(
        "--subsample-n",
        type=int,
        default=None,
        metavar="N",
        help="Optional per-accession cell cap before concat",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Write a metadata JSON file next to the run CSV with this note.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    log.info("new atlas integration pipeline run started (r2 prefix: %s)", args.r2_prefix)

    run_dir = RUN_OUTPUT_DIR / RUN_TIMESTAMP
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = run_dir / "metadata.json"
    _write_run_metadata(metadata_path, args, run_dir)

    cfg = AtlasIntegrationConfig(
        datasetsCsvPath=args.datasets,
        contextsPath=args.contexts,
        outputDir=run_dir / "atlas",
        figsDir=run_dir / "figs",
        subsampleN=args.subsample_n,
    )

    merged, merge_stats = merge_datasets_incrementally(cfg, run_dir)
    adata, result = run_atlas_integration(cfg, adata=merged, merge_stats=merge_stats)

    atlas_r2_key = f"{args.r2_prefix}/{cfg.atlasH5adName}"
    upload_to_r2(result.atlasPath, atlas_r2_key)
    verify_upload(atlas_r2_key)

    run_csv = run_dir / "run.csv"
    append_csv_row(
        run_csv,
        ["status", "atlasPath", "r2Key", "nCellsFinal", "nStudies", "timestamp"],
        [
            "success",
            rel_to_repo(result.atlasPath),
            atlas_r2_key,
            str(adata.n_obs),
            str(adata.obs[cfg.batchKey].nunique()),
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        ],
    )
    log.info("Atlas pipeline complete: %s uploaded to %s", result.atlasPath, atlas_r2_key)


if __name__ == "__main__":
    main()
