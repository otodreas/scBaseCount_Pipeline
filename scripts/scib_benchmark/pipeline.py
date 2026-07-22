import time
from collections.abc import Callable
from logging import Logger
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import scanpy as sc
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation
from shared.logger import configure_file_logger
from shared.repo import rel_to_repo

_log = configure_file_logger("scib_benchmark.log", __name__)

DEFAULT_EMBEDDING_KEYS = ["X_pca", "X_pca_harmony", "X_umap", "X_umap_uncorrected"]
RESULTS_CSV = "scib_results.csv"
RESULTS_SVG = "scib_results.svg"


def _timed[T](log: Logger, name: str, action: Callable[[], T]) -> T:
    """Run action, logging START/DONE with elapsed seconds, and return its result."""
    started = time.perf_counter()
    log.info("START %s", name)
    result = action()
    log.info("DONE %s in %.1fs", name, time.perf_counter() - started)
    return result


def _validate_embedding_keys(adata: sc.AnnData, embedding_keys: list[str]) -> None:
    """Raise if any embedding key is missing from adata.obsm."""
    missing = [key for key in embedding_keys if key not in adata.obsm]
    if missing:
        msg = f"Missing embedding keys in adata.obsm: {missing}"
        raise ValueError(msg)


def _results_paths(out_dir: Path) -> tuple[Path, Path]:
    """Return the cached CSV and SVG paths under out_dir."""
    return out_dir / RESULTS_CSV, out_dir / RESULTS_SVG


def run_scib_benchmark(
    input_h5ad: Path,
    out_dir: Path,
    batch_key: str,
    label_key: str,
    embedding_keys: list[str],
    n_jobs: int,
    force: bool,
    log: Logger | None = None,
) -> None:
    """Run scib_metrics benchmark unless cached CSV and SVG already exist."""
    log = log or _log
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path, svg_path = _results_paths(out_dir)

    if not force:
        if csv_path.is_file() and svg_path.is_file():
            log.info("Using cached results at %s and %s", rel_to_repo(csv_path), rel_to_repo(svg_path))
            return
        if csv_path.is_file() or svg_path.is_file():
            log.warning(
                "Incomplete cache in %s (expected %s and %s); re-running benchmark",
                rel_to_repo(out_dir),
                RESULTS_CSV,
                RESULTS_SVG,
            )

    adata = sc.read_h5ad(input_h5ad, backed="r")
    log.info("Loaded %s cells x %s genes from %s", f"{adata.n_obs:,}", f"{adata.n_vars:,}", rel_to_repo(input_h5ad))
    _validate_embedding_keys(adata, embedding_keys)

    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=embedding_keys,
        # pick uncorrected pca as benchmarked embedding
        pre_integrated_embedding_obsm_key="X_pca",
        n_jobs=n_jobs,
    )
    _timed(log, "benchmark", bm.benchmark)

    df = bm.get_results()
    df.to_csv(csv_path)
    log.info("Wrote %s", rel_to_repo(csv_path))

    bm.plot_results_table(show=False, save_dir=str(out_dir))
    log.info("Wrote %s", rel_to_repo(svg_path))
