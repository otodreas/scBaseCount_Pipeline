from dataclasses import dataclass

import anndata as ad
import pandas as pd
import scanpy as sc

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason


def _gene_names(adata: ad.AnnData) -> pd.Series:
    """Return uppercased gene names from gene_symbols when present, else var_names."""
    if "gene_symbols" in adata.var.columns:
        raw = adata.var["gene_symbols"].to_numpy()
    else:
        raw = adata.var_names.to_numpy()
    return pd.Series(raw, index=adata.var_names, dtype="string").str.upper()


@dataclass
class QcStats:
    nCellsBefore: int
    nCellsAfter: int
    nGenesBefore: int
    nGenesAfter: int
    medianGenesPerCell: float
    medianPctMito: float
    medianPctRibo: float
    medianPctHb: float
    pctCellsDropped: float


def flag_qc_genes(adata: ad.AnnData) -> None:
    """Set adata.var flags mt, ribo, and hb from gene names for QC metric computation."""
    names = _gene_names(adata)
    adata.var["mt"] = names.str.startswith("MT-").to_numpy()
    adata.var["ribo"] = names.str.match(r"^RP[SL]\d").to_numpy()
    # Hemoglobin subunit genes; anchored to exclude non-globin HB* genes such as HBP1, HBS1L, HBEGF.
    adata.var["hb"] = names.str.match(r"^HB[ABDEGMQZ]\d?$").to_numpy()


def apply_qc_gate(adata: ad.AnnData, cfg: H5adConcatConfig) -> tuple[ad.AnnData, QcStats]:
    """Filter low-quality cells and genes; raise FileRejected when no cells remain."""
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    flag_qc_genes(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False)
    sc.pp.filter_cells(adata, min_genes=cfg.minGenesPerCell)
    if cfg.minCellsPerGene > 0:
        sc.pp.filter_genes(adata, min_cells=cfg.minCellsPerGene)
    adata = adata[adata.obs["pct_counts_mt"] < cfg.maxPctMito].copy()
    if cfg.maxPctHb is not None:
        adata = adata[adata.obs["pct_counts_hb"] < cfg.maxPctHb].copy()

    n_cells_after = adata.n_obs
    if n_cells_after < cfg.minCellsAfterQc:
        raise FileRejected(SkipReason.preprocess_failed)

    median_genes = float(adata.obs["n_genes_by_counts"].median()) if n_cells_after > 0 else 0.0
    median_mito = float(adata.obs["pct_counts_mt"].median()) if n_cells_after > 0 else 0.0
    median_ribo = float(adata.obs["pct_counts_ribo"].median()) if n_cells_after > 0 else 0.0
    median_hb = float(adata.obs["pct_counts_hb"].median()) if n_cells_after > 0 else 0.0
    pct_dropped = 100.0 * (n_cells_before - n_cells_after) / n_cells_before if n_cells_before > 0 else 0.0

    return adata, QcStats(
        nCellsBefore=n_cells_before,
        nCellsAfter=n_cells_after,
        nGenesBefore=n_genes_before,
        nGenesAfter=adata.n_vars,
        medianGenesPerCell=median_genes,
        medianPctMito=median_mito,
        medianPctRibo=median_ribo,
        medianPctHb=median_hb,
        pctCellsDropped=pct_dropped,
    )
