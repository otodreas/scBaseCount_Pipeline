import argparse
import datetime
import time
from pathlib import Path

import scanpy as sc
from atlas_postprocessing.batch_comparison import (
    EXPECTED_BASELINE,
    PCA_KEY,
    SRX_KEY,
    STUDY_KEY,
    STUDY_TECH_KEY,
    TECH_KEY,
    aggregate_scib_results,
    attach_tech_10x,
    batch_cardinalities,
    ensure_shared_pca,
    harmony_obsm_key,
    load_run_json,
    make_study_tech_batch_key,
    preserve_study_harmony_embedding,
    run_harmony_variants,
    run_scib_evaluation_grid,
    validate_baseline_manifests,
    write_comparison_summary,
)
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import apply_thread_settings, timed
from atlas_postprocessing.sampling import sample_metadata
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo

_LOG_FILENAME = "compare_atlas_batch_keys.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
configure_file_logger(_LOG_FILENAME, "atlas_postprocessing")
add_stdout_handler()

_DEFAULT_SUBSET = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "subset_validation" / "atlas_pp_subset.h5ad"
_DEFAULT_SUBSET_RUN = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "subset_validation" / "atlas_pp_subset_run.json"
_DEFAULT_PRODUCTION_RUN = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "production" / "atlas_v2_post_run.json"
_DEFAULT_DATASETS = REPO_ROOT / "output" / "metadata" / "datasets_v2.csv"
_DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / "atlas" / "v2" / "post" / "batch_key_comparison"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Harmony batch keys on a frozen atlas subset with shared PCA. "
            "Joins tech_10x metadata, builds study×technology and SRX Harmony embeddings, "
            "and runs scIB once per evaluation batch key across all embeddings."
        )
    )
    parser.add_argument(
        "--subset-h5ad",
        type=Path,
        default=_DEFAULT_SUBSET,
        metavar="PATH",
        help="Frozen subset validation h5ad with shared X_pca and study X_pca_harmony",
    )
    parser.add_argument(
        "--subset-run-json",
        type=Path,
        default=_DEFAULT_SUBSET_RUN,
        metavar="PATH",
        help="Subset validation run JSON used to pin the baseline",
    )
    parser.add_argument(
        "--production-run-json",
        type=Path,
        default=_DEFAULT_PRODUCTION_RUN,
        metavar="PATH",
        help="Production run JSON used to confirm the same baseline",
    )
    parser.add_argument(
        "--datasets",
        type=Path,
        default=_DEFAULT_DATASETS,
        metavar="PATH",
        help="datasets CSV containing srx_accession and tech_10x",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTPUT_DIR,
        metavar="PATH",
        help="Comparison output root",
    )
    parser.add_argument(
        "--cell-type-key",
        type=str,
        default="cell_type",
        metavar="COL",
        help="obs label column for scIB bio-conservation metrics",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        metavar="N",
        help="Thread budget for Harmony (0 leaves library defaults)",
    )
    parser.add_argument("--scib-jobs", type=int, default=6, metavar="N", help="scIB n_jobs")
    parser.add_argument(
        "--force-scib",
        action="store_true",
        help="Re-run scIB even when evaluation artifacts already exist",
    )
    parser.add_argument(
        "--skip-scib",
        action="store_true",
        help="Stop after writing Harmony embeddings and the comparison h5ad",
    )
    parser.add_argument(
        "--reuse-comparison-h5ad",
        action="store_true",
        help="Reuse an existing comparison h5ad in the output directory when present",
    )
    return parser.parse_args()


def _comparison_h5ad_path(outputDir: Path) -> Path:
    return outputDir / "atlas_pp_batch_comparison.h5ad"


def _build_cfg(args: argparse.Namespace) -> AtlasPostprocessingConfig:
    return AtlasPostprocessingConfig(
        batchKey=STUDY_KEY,
        cellTypeKey=args.cell_type_key,
        nTopGenes=int(EXPECTED_BASELINE["nTopGenes"]),
        nPcs=int(EXPECTED_BASELINE["nPcs"]),
        nPcsCompute=int(EXPECTED_BASELINE["nPcs"]),
        nNeighbors=int(EXPECTED_BASELINE["nNeighbors"]),
        resolution=float(EXPECTED_BASELINE["resolution"]),
        writePlots=False,
        nJobs=args.threads,
    )


def _prepare_adata(args: argparse.Namespace, cfg: AtlasPostprocessingConfig) -> tuple[sc.AnnData, dict]:
    comparison_h5ad = _comparison_h5ad_path(args.output_dir)
    if args.reuse_comparison_h5ad and comparison_h5ad.is_file():
        adata = timed("load comparison h5ad", lambda: sc.read_h5ad(comparison_h5ad), logger=log)
        join_audit = {
            "datasetsPath": rel_to_repo(args.datasets),
            "reusedComparisonH5ad": True,
            "nCells": int(adata.n_obs),
            "nExperiments": int(adata.obs[SRX_KEY].nunique()) if SRX_KEY in adata.obs else None,
            "nStudyTechBatches": int(adata.obs[STUDY_TECH_KEY].nunique()) if STUDY_TECH_KEY in adata.obs else None,
            "techValues": sorted(adata.obs[TECH_KEY].astype(str).unique().tolist()) if TECH_KEY in adata.obs else [],
        }
        return adata, join_audit

    adata = timed("load frozen subset", lambda: sc.read_h5ad(args.subset_h5ad), logger=log)
    ensure_shared_pca(adata, nPcs=cfg.nPcs)
    join_audit = attach_tech_10x(adata, args.datasets)
    make_study_tech_batch_key(adata)
    preserve_study_harmony_embedding(adata, studyKey=STUDY_KEY)
    return adata, join_audit


def main() -> None:
    args = _parse_args()
    log_run_separator(log)
    log.info("compare_atlas_batch_keys started")
    started = time.perf_counter()
    started_at = datetime.datetime.now()

    subset_run = load_run_json(args.subset_run_json)
    production_run = load_run_json(args.production_run_json)
    baseline = validate_baseline_manifests(subsetRun=subset_run, productionRun=production_run)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cfg = _build_cfg(args)
    apply_thread_settings(cfg)

    adata, join_audit = _prepare_adata(args, cfg)
    harmony_keys = [STUDY_KEY, STUDY_TECH_KEY, SRX_KEY]
    for key in harmony_keys:
        if key not in adata.obs:
            raise ValueError(f"Prepared AnnData is missing required obs column {key!r}")

    embedding_map = run_harmony_variants(
        adata,
        cfg,
        batchKeys=[STUDY_TECH_KEY, SRX_KEY],
        nPcs=cfg.nPcs,
        skipExisting=True,
    )
    # Study Harmony is preserved from the frozen subset; keep it in the map for scIB.
    embedding_map[STUDY_KEY] = harmony_obsm_key(STUDY_KEY)
    if embedding_map[STUDY_KEY] not in adata.obsm:
        preserve_study_harmony_embedding(adata, studyKey=STUDY_KEY)

    comparison_h5ad = _comparison_h5ad_path(args.output_dir)
    timed(
        "write comparison h5ad",
        lambda: adata.write_h5ad(comparison_h5ad, compression=cfg.compression),
        logger=log,
    )

    embedding_keys = [PCA_KEY, *[embedding_map[key] for key in harmony_keys]]
    scib_root = args.output_dir / "scib"
    scib_runs: list[dict] = []
    matrix_csv = scib_root / "scib_matrix_long.csv"
    if args.skip_scib:
        log.info("Skipping scIB (--skip-scib)")
    else:
        scib_runs = run_scib_evaluation_grid(
            adata,
            outRoot=scib_root,
            embeddingKeys=embedding_keys,
            evalBatchKeys=harmony_keys,
            labelKey=cfg.cellTypeKey,
            nJobs=args.scib_jobs,
            force=args.force_scib,
        )
        aggregate_scib_results(scib_runs, embeddingMap=embedding_map, outCsv=matrix_csv)

    summary = {
        "inputSubsetH5ad": rel_to_repo(args.subset_h5ad),
        "subsetRunJson": rel_to_repo(args.subset_run_json),
        "productionRunJson": rel_to_repo(args.production_run_json),
        "datasetsPath": rel_to_repo(args.datasets),
        "outputDir": rel_to_repo(args.output_dir),
        "comparisonH5ad": rel_to_repo(comparison_h5ad),
        "resolved": dict(EXPECTED_BASELINE),
        "baselineValidation": baseline,
        "sampling": sample_metadata(adata),
        "joinAudit": join_audit,
        "batchCardinalities": batch_cardinalities(adata, harmony_keys),
        "embeddings": {
            "uncorrected": PCA_KEY,
            "harmonyByBatchKey": embedding_map,
        },
        "scib": {
            "skipped": bool(args.skip_scib),
            "matrixLong": None if args.skip_scib else rel_to_repo(matrix_csv),
            "runs": scib_runs,
        },
        "startedAt": started_at.isoformat(timespec="seconds"),
        "timingsSeconds": round(time.perf_counter() - started, 3),
        "note": (
            "Interpret SRX evaluation scores cautiously: SRX often encodes donor, sample, "
            "or condition biology in addition to technical effects."
        ),
    }
    write_comparison_summary(args.output_dir / "batch_key_comparison_summary.json", summary)
    log.info("compare_atlas_batch_keys complete in %.1fs", summary["timingsSeconds"])


if __name__ == "__main__":
    main()
