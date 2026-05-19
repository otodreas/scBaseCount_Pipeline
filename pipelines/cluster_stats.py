from __future__ import annotations

import argparse
import datetime
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import scanpy as sc
from cluster_validation import build_metric_dataframes, compute_nse_kld_row, save_metric_plot
from cluster_validation.merge import MERGED_CLUSTER_KEY
from dotenv import load_dotenv
from r2 import download_from_r2, fetch_uploaded_r2_keys
from shared.csv_writer import append_csv_row
from shared.files import safe_delete
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo

load_dotenv()

RUN_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "output" / "cluster_stats"
CLUSTERED_SUFFIX = "_clustered.h5ad"
_LOG_FILENAME = "cluster_stats.log"

log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


_CSV_COLUMNS = ["position", "srx", "status", "r2_file", "n_cell_types", "n_clusters", "timestamp", "error"]


def _append_summary_row(
    summary_path: Path,
    srx: str,
    status: str,
    position: str = "",
    r2_key: str = "",
    n_cell_types: int | str = "",
    n_clusters: int | str = "",
    error: str = "",
) -> None:
    append_csv_row(
        summary_path,
        _CSV_COLUMNS,
        [
            position,
            srx,
            status,
            r2_key,
            n_cell_types,
            n_clusters,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            error,
        ],
    )


def _write_run_metadata(
    metadata_path: Path,
    args: argparse.Namespace,
    run_ts: str,
    n_files: int,
    output_dir: Path,
) -> None:
    payload: dict = {
        "run_timestamp": run_ts,
        "run_dir": rel_to_repo(output_dir),
        "run_csv": rel_to_repo(output_dir / "run.csv"),
        "log_path": rel_to_repo(REPO_ROOT / "logs" / _LOG_FILENAME),
        "r2_prefix": args.r2_prefix,
        "n_files_matched": n_files,
    }
    if args.metadata is not None:
        payload["notes"] = args.metadata
    metadata_path.write_text(json.dumps(payload, indent=2))


def _srx_from_key(r2_key: str) -> str:
    filename = Path(r2_key).name
    return filename[: -len(CLUSTERED_SUFFIX)] if filename.endswith(CLUSTERED_SUFFIX) else filename


def _compute_count_matrix(adata: sc.AnnData) -> dict[str, dict[str, int]]:
    counts = pd.crosstab(adata.obs["cell_type"], adata.obs[MERGED_CLUSTER_KEY])
    matrix: dict[str, dict[str, int]] = {}
    for cell_type, row in counts.iterrows():
        matrix[str(cell_type)] = {str(cluster): int(value) for cluster, value in row.items()}
    return matrix


def process_accession(
    srx: str,
    r2_key: str,
    download_root: Path,
) -> tuple[dict[str, dict[str, int]] | None, dict[str, float] | None, dict[str, float] | None, Exception | None]:
    local_path = download_root / f"{srx}{CLUSTERED_SUFFIX}"
    try:
        log.info("%s: downloading %s", srx, r2_key)
        download_from_r2(r2_key, local_path)

        adata = sc.read_h5ad(local_path, backed="r")
        try:
            count_matrix = _compute_count_matrix(adata)
            nse_row, kld_row = compute_nse_kld_row(adata, MERGED_CLUSTER_KEY)
        finally:
            try:
                adata.file.close()
            except Exception:
                pass

        safe_delete(local_path, log)
        log.info(
            "%s: done (%d cell types, %d clusters)", srx, len(count_matrix), len(next(iter(count_matrix.values()), {}))
        )
        return count_matrix, nse_row, kld_row, None

    except Exception as exc:
        log.exception("%s: failed", srx)
        safe_delete(local_path, log)
        return None, None, None, exc


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute per-SRX cell_type x leiden_merged count matrices and NSE/KLD analytics from clustered h5ads stored under an R2 prefix."
    )
    parser.add_argument(
        "--r2-prefix",
        type=str,
        required=True,
        metavar="PREFIX",
        help="R2 prefix produced by a clustering run (e.g. clustering_pipeline_20260511_140000).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        metavar="PATH",
        help="Output directory (default: output/cluster_stats/{r2_prefix}).",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="TEXT",
        help="Write a metadata JSON file with this note next to the outputs.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Number of accessions to process in parallel (default: 1).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    log_run_separator(log)
    log.info("new cluster_stats run started (r2 prefix: %s)", args.r2_prefix)

    output_dir = args.output_dir if args.output_dir is not None else DEFAULT_OUTPUT_ROOT / args.r2_prefix
    output_dir.mkdir(parents=True, exist_ok=True)
    download_dir = output_dir / "data"
    download_dir.mkdir(parents=True, exist_ok=True)

    csv_summary_path = output_dir / "run.csv"
    metadata_path = output_dir / "metadata.json"
    cluster_stats_json_path = output_dir / "cluster_stats.json"
    nse_csv_path = output_dir / "nse_matrix.csv"
    kld_csv_path = output_dir / "kld_matrix.csv"
    summary_csv_path = output_dir / "cell_type_summary.csv"
    plot_path = output_dir / "cell_type_metrics.png"

    log.info("output directory: %s", output_dir)

    prefix_with_slash = args.r2_prefix.rstrip("/") + "/"
    all_keys = fetch_uploaded_r2_keys()
    matched_keys = sorted(
        key for key in all_keys if key.startswith(prefix_with_slash) and key.endswith(CLUSTERED_SUFFIX)
    )
    log.info("Matched %d clustered h5ad(s) under prefix %s", len(matched_keys), args.r2_prefix)
    _write_run_metadata(metadata_path, args, RUN_TIMESTAMP, len(matched_keys), output_dir)

    if not matched_keys:
        log.warning("No clustered h5ad files matched the prefix; nothing to do.")
        return

    total = len(matched_keys)
    work_items: list[tuple[str, str, str]] = []
    for n, key in enumerate(matched_keys, start=1):
        srx = _srx_from_key(key)
        position = f"{n}/{total}"
        work_items.append((srx, key, position))

    nested: dict[str, dict[str, dict[str, int]]] = {}
    metric_rows: list[dict] = []

    log.info("Submitting %d accession(s) to %d worker(s)", len(work_items), args.workers)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_accession, srx, key, download_dir): (srx, key, position)
            for srx, key, position in work_items
        }
        for future in as_completed(futures):
            srx, key, position = futures[future]
            count_matrix, nse_row, kld_row, exc = future.result()
            if exc is not None or count_matrix is None or nse_row is None or kld_row is None:
                _append_summary_row(
                    csv_summary_path,
                    srx,
                    "failed",
                    position,
                    key,
                    error=f"{type(exc).__name__}: {exc}" if exc else "empty result",
                )
                continue
            nested[srx] = count_matrix
            metric_rows.append({"srx": srx, "nse": nse_row, "kld": kld_row})
            n_clusters = len(next(iter(count_matrix.values()), {}))
            _append_summary_row(
                csv_summary_path,
                srx,
                "success",
                position,
                key,
                n_cell_types=len(count_matrix),
                n_clusters=n_clusters,
            )

    try:
        download_dir.rmdir()
        log.debug("Removed empty download directory %s", download_dir)
    except OSError:
        log.warning("Could not remove download directory %s (may contain leftover files)", download_dir)

    cluster_stats_json_path.write_text(json.dumps(nested, indent=2, sort_keys=True))
    log.info("cluster stats JSON written to %s (%d srx)", cluster_stats_json_path, len(nested))

    if metric_rows:
        nse_df, kld_df, summary_df = build_metric_dataframes(metric_rows)
        nse_df.to_csv(nse_csv_path)
        log.info("NSE matrix written to %s", nse_csv_path)
        kld_df.to_csv(kld_csv_path)
        log.info("KLD matrix written to %s", kld_csv_path)
        summary_df.to_csv(summary_csv_path)
        log.info("Cell type summary written to %s", summary_csv_path)
        save_metric_plot(summary_df, plot_path)
        log.info("Cell type metrics plot saved to %s", plot_path)
    else:
        log.warning("No successful accessions; skipping NSE/KLD aggregation outputs.")


if __name__ == "__main__":
    main()
