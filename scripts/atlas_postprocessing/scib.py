import logging
from pathlib import Path

import scanpy as sc
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation
from shared.repo import rel_to_repo

log = logging.getLogger(__name__)


def run_scib_benchmark(
    adata: sc.AnnData | Path,
    *,
    outDir: Path,
    batchKey: str = "study_accession",
    labelKey: str = "cell_type",
    embeddingKeys: list[str] | None = None,
    preIntegratedKey: str = "X_pca",
    nJobs: int = 6,
    force: bool = False,
) -> Path:
    """Run scIB bio-conservation and batch-correction metrics and write CSV/SVG artifacts."""
    outDir.mkdir(parents=True, exist_ok=True)
    csv_path = outDir / "scib_results.csv"
    svg_path = outDir / "scib_results.svg"
    if not force and csv_path.is_file() and svg_path.is_file():
        log.info("scIB artifacts already exist under %s; skipping (pass force=True to re-run)", rel_to_repo(outDir))
        return csv_path

    if isinstance(adata, Path):
        log.info("START scIB load %s", rel_to_repo(adata))
        loaded = sc.read_h5ad(adata, backed="r")
        log.info("DONE scIB load")
    else:
        loaded = adata

    keys = embeddingKeys or ["X_pca", "X_pca_harmony"]
    log.info("START scIB Benchmarker setup embeddings=%s", keys)
    bm = Benchmarker(
        adata=loaded,
        batch_key=batchKey,
        label_key=labelKey,
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=keys,
        pre_integrated_embedding_obsm_key=preIntegratedKey,
        n_jobs=nJobs,
    )
    log.info("DONE scIB Benchmarker setup")

    log.info("START scIB bio conservation + batch correction metrics")
    bm.benchmark()
    log.info("DONE scIB metrics")

    log.info("START scIB write artifacts")
    df = bm.get_results()
    df.to_csv(csv_path)
    bm.plot_results_table(show=False, save_dir=str(outDir))
    log.info("DONE scIB write artifacts csv=%s svg=%s", rel_to_repo(csv_path), rel_to_repo(svg_path))
    return csv_path
