import argparse
import datetime
import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import harmonypy
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
from pydantic import BaseModel
from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import REPO_ROOT, rel_to_repo
from umap_plots import plot_umap
from umap_plots.config import UmapPlotConfig

_LOG_FILENAME = "atlas_harmony.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()

_DEFAULT_INPUT = REPO_ROOT / "output" / "atlas" / "data" / "atlas_sample20.h5ad"
_DEFAULT_OUTPUT = REPO_ROOT / "output" / "atlas" / "data" / "atlas_sample20_harmony.h5ad"
_DEFAULT_FIGS = REPO_ROOT / "output" / "atlas" / "figs"


class AtlasHarmonyConfig(BaseModel):
    inputH5ad: Path
    outputH5ad: Path
    figsDir: Path
    batchKey: str = "study_accession"
    cellTypeKey: str = "cell_type"
    nTopGenes: int = 2000
    nPcs: int = 20
    nPcsCompute: int = 25
    resolution: float = 1.0
    scaleMaxValue: float = 10.0  # clip z-scores to +scaleMaxValue SD
    writePlots: bool = True
    compression: Literal["gzip", "lzf"] | None = "gzip"


def load_and_normalize(cfg: AtlasHarmonyConfig) -> sc.AnnData:
    """Read the atlas h5ad, stash raw counts, and apply normalize_total + log1p."""
    adata = sc.read_h5ad(cfg.inputH5ad)
    log.info("Loaded %s cells x %s genes", f"{adata.n_obs:,}", f"{adata.n_vars:,}")
    log.info("Studies (%s): %d", cfg.batchKey, adata.obs[cfg.batchKey].nunique())

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def embed_uncorrected(adata: sc.AnnData, cfg: AtlasHarmonyConfig) -> sc.AnnData:
    """Select batch-aware HVGs, scale, run PCA, and build the pre-correction UMAP and leiden partition."""
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.nTopGenes, batch_key=cfg.batchKey)
    # Subset counts to only the HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.info("Retained %s HVGs", f"{adata.n_vars:,}")

    sc.pp.scale(adata, max_value=cfg.scaleMaxValue)
    sc.tl.pca(adata, n_comps=cfg.nPcsCompute, svd_solver="arpack")

    sc.pp.neighbors(adata, n_pcs=cfg.nPcs)
    sc.tl.umap(adata)
    adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()

    sc.tl.leiden(
        adata,
        resolution=cfg.resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden_uncorrected",
    )
    log.info("Uncorrected leiden clusters: %d", adata.obs["leiden_uncorrected"].nunique())
    return adata


def integrate_harmony(adata: sc.AnnData, cfg: AtlasHarmonyConfig) -> sc.AnnData:
    """Run Harmony on the PCA embedding and build the corrected UMAP and leiden partition."""
    harmony_out = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, cfg.batchKey)
    adata.obsm["X_pca_harmony"] = harmony_out.Z_corr

    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=cfg.resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden_atlas",
    )
    log.info("Harmony-corrected leiden clusters: %d", adata.obs["leiden_atlas"].nunique())
    return adata


def make_plots(adata: sc.AnnData, cfg: AtlasHarmonyConfig) -> None:
    """Write pre- and post-correction UMAP PNGs colored by batch and cell type."""
    cfg.figsDir.mkdir(parents=True, exist_ok=True)
    for color_by in (cfg.batchKey, cfg.cellTypeKey):
        fig = plot_umap(
            adata,
            colorBy=color_by,
            cfg=UmapPlotConfig(figsDir=cfg.figsDir, umapKey="X_umap_uncorrected"),
            nameSuffix="uncorrected",
        )
        plt.close(fig)
        fig = plot_umap(
            adata,
            colorBy=color_by,
            cfg=UmapPlotConfig(figsDir=cfg.figsDir, umapKey="X_umap"),
            nameSuffix="harmony",
        )
        plt.close(fig)


def save_atlas(adata: sc.AnnData, cfg: AtlasHarmonyConfig) -> None:
    """Write the processed atlas h5ad and a run summary JSON alongside it."""
    cfg.outputH5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(cfg.outputH5ad, compression=cfg.compression)

    summary = {
        "input": rel_to_repo(cfg.inputH5ad),
        "output": rel_to_repo(cfg.outputH5ad),
        "figsDir": rel_to_repo(cfg.figsDir),
        "cells": int(adata.n_obs),
        "hvgs": int(adata.n_vars),
        "studies": int(adata.obs[cfg.batchKey].nunique()),
        "clustersUncorrected": int(adata.obs["leiden_uncorrected"].nunique()),
        "clustersHarmony": int(adata.obs["leiden_atlas"].nunique()),
        "resolution": cfg.resolution,
        "nPcs": cfg.nPcs,
        "nTopGenes": cfg.nTopGenes,
    }
    summary_path = cfg.outputH5ad.with_name(f"{cfg.outputH5ad.stem}_run.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    log.info("Wrote %s", rel_to_repo(cfg.outputH5ad))
    log.info("Wrote %s", rel_to_repo(summary_path))


def _timed[T](name: str, action: Callable[[], T]) -> T:
    """Run action, logging START/DONE with elapsed seconds, and return its result."""
    started = time.perf_counter()
    log.info("START %s", name)
    result = action()
    log.info("DONE %s in %.1fs", name, time.perf_counter() - started)
    return result


def run(cfg: AtlasHarmonyConfig) -> None:
    """Run the full atlas Harmony flow: load, embed, integrate, plot, and save."""
    adata = _timed("load + normalize", lambda: load_and_normalize(cfg))
    adata = _timed("HVG + PCA + pre-correction embedding", lambda: embed_uncorrected(adata, cfg))
    adata = _timed("harmony integration", lambda: integrate_harmony(adata, cfg))

    if cfg.writePlots:
        _timed("plots", lambda: make_plots(adata, cfg))
    _timed("save", lambda: save_atlas(adata, cfg))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Harmony batch correction and clustering on a merged atlas h5ad.")
    parser.add_argument("--input", type=Path, default=_DEFAULT_INPUT, metavar="PATH", help="Input atlas h5ad")
    parser.add_argument("--output", type=Path, default=_DEFAULT_OUTPUT, metavar="PATH", help="Output atlas h5ad")
    parser.add_argument("--figs-dir", type=Path, default=_DEFAULT_FIGS, metavar="PATH", help="Directory for UMAP PNGs")
    parser.add_argument("--batch-key", type=str, default="study_accession", metavar="COL", help="obs batch column")
    parser.add_argument("--cell-type-key", type=str, default="cell_type", metavar="COL", help="obs cell type column")
    parser.add_argument("--n-top-genes", type=int, default=2000, metavar="N", help="Number of HVGs")
    parser.add_argument("--n-pcs", type=int, default=30, metavar="N", help="PCs used for the neighbor graph")
    parser.add_argument("--n-pcs-compute", type=int, default=50, metavar="N", help="PCs computed by PCA")
    parser.add_argument("--resolution", type=float, default=1.0, metavar="R", help="Leiden resolution")
    parser.add_argument("--no-plots", action="store_true", help="Skip writing UMAP PNGs")
    parser.add_argument("--threads", type=int, default=0, metavar="N", help="scanpy n_jobs (0 leaves the default)")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    if args.threads > 0:
        sc.settings.n_jobs = args.threads

    cfg = AtlasHarmonyConfig(
        inputH5ad=args.input,
        outputH5ad=args.output,
        figsDir=args.figs_dir,
        batchKey=args.batch_key,
        cellTypeKey=args.cell_type_key,
        nTopGenes=args.n_top_genes,
        nPcs=args.n_pcs,
        nPcsCompute=args.n_pcs_compute,
        resolution=args.resolution,
        writePlots=not args.no_plots,
    )

    log_run_separator(log)
    log.info("new atlas harmony run started")
    log.info("config: %s", cfg.model_dump_json())

    started = datetime.datetime.now()
    run(cfg)
    log.info("atlas harmony run complete in %s", datetime.datetime.now() - started)


if __name__ == "__main__":
    main()
