import anndata as ad
import pandas as pd
import scanpy as sc

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import CELL_FILTER_ORDER, PostFilterMedians, QcStats, SkipReason

_QC_GENE_FLAGS = ("mt", "ribo", "hb")

# var columns apply_qc_gate writes: the gene flags from flag_qc_genes plus the per-gene metrics
# scanpy adds (calculate_qc_metrics with log1p=False). Exposed so downstream steps can tell these
# per-file QC stats apart from genuine gene annotations.
QC_VAR_KEYS: frozenset[str] = frozenset(_QC_GENE_FLAGS) | {
    "n_cells_by_counts",
    "mean_counts",
    "pct_dropout_by_counts",
    "total_counts",
}


def _gene_names(adata: ad.AnnData) -> pd.Series:
    """Return uppercased gene names from gene_symbols when present, else var_names."""
    if "gene_symbols" in adata.var.columns:
        raw = adata.var["gene_symbols"].to_numpy()
    else:
        raw = adata.var_names.to_numpy()
    return pd.Series(raw, index=adata.var_names, dtype="string").str.upper()


def flag_qc_genes(adata: ad.AnnData) -> None:
    """Set adata.var flags mt, ribo, and hb from gene names for QC metric computation."""
    names = _gene_names(adata)
    adata.var["mt"] = names.str.startswith("MT-").to_numpy()
    adata.var["ribo"] = names.str.match(r"^RP[SL]\d").to_numpy()
    # Hemoglobin subunit genes; anchored to exclude non-globin HB* genes such as HBP1, HBS1L, HBEGF.
    adata.var["hb"] = names.str.match(r"^HB[ABDEGMQZ]\d?$").to_numpy()


def _drop_cells_by_fraction(adata: ad.AnnData, obs_key: str, max_fraction: float) -> tuple[ad.AnnData, int]:
    """Keep cells with obs_key strictly below max_fraction * 100; return filtered adata and drop count."""
    keep = adata.obs[obs_key] < max_fraction * 100
    n_dropped = int((~keep).sum())
    if n_dropped == 0:
        return adata, 0
    return adata[keep].copy(), n_dropped


def apply_qc_gate(adata: ad.AnnData, cfg: H5adConcatConfig) -> tuple[ad.AnnData, QcStats]:
    """Filter low-quality cells and genes; raise FileRejected when remaining cells fail file-level gates."""
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    n_cells_dropped_by_filter = {name: 0 for name in CELL_FILTER_ORDER}

    flag_qc_genes(adata)  # mentioned in methods--flag ribo/mito/hb with scanpy advised regex
    sc.pp.calculate_qc_metrics(adata, qc_vars=list(_QC_GENE_FLAGS), inplace=True, log1p=False)

    before_min_genes = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=cfg.minGenesPerCell)  # mentioned in methods
    n_cells_dropped_by_filter["minGenesPerCell"] = before_min_genes - adata.n_obs

    n_genes_dropped_min_cells = 0
    # this filter is not applied in the atlas build--it erases real biology that might materialize at atlas scale
    if cfg.minCellsPerGene > 0:
        before_genes = adata.n_vars
        sc.pp.filter_genes(adata, min_cells=cfg.minCellsPerGene)
        n_genes_dropped_min_cells = before_genes - adata.n_vars

    # scanpy reports pct_counts_* on a 0-100 scale; config ceilings are fractions in (0, 1].
    adata, n_dropped_mito = _drop_cells_by_fraction(adata, "pct_counts_mt", cfg.maxPctMito)
    n_cells_dropped_by_filter["maxPctMito"] = n_dropped_mito

    adata, n_dropped_ribo = _drop_cells_by_fraction(adata, "pct_counts_ribo", cfg.maxPctRibo)
    n_cells_dropped_by_filter["maxPctRibo"] = n_dropped_ribo

    adata, n_dropped_hb = _drop_cells_by_fraction(adata, "pct_counts_hb", cfg.maxPctHb)
    n_cells_dropped_by_filter["maxPctHb"] = n_dropped_hb

    n_cells_after = adata.n_obs
    n_cells_dropped = n_cells_before - n_cells_after
    pct_cells_after = n_cells_after / n_cells_before if n_cells_before > 0 else 0.0

    median_genes = float(adata.obs["n_genes_by_counts"].median()) if n_cells_after > 0 else 0.0
    median_mito = float(adata.obs["pct_counts_mt"].median()) if n_cells_after > 0 else 0.0
    median_ribo = float(adata.obs["pct_counts_ribo"].median()) if n_cells_after > 0 else 0.0
    median_hb = float(adata.obs["pct_counts_hb"].median()) if n_cells_after > 0 else 0.0

    stats = QcStats(
        nCellsBefore=n_cells_before,
        nCellsAfter=n_cells_after,
        nCellsDropped=n_cells_dropped,
        nCellsDroppedByFilter=n_cells_dropped_by_filter,
        nGenesBefore=n_genes_before,
        nGenesAfter=adata.n_vars,
        nGenesDroppedByFilter={"minCellsPerGene": n_genes_dropped_min_cells},
        pctCellsAfter=pct_cells_after,
        postFilterMedians=PostFilterMedians(
            nGenesByCounts=median_genes,
            pctCountsMito=median_mito,
            pctCountsRibo=median_ribo,
            pctCountsHb=median_hb,
        ),
    )

    # mentioned in methods
    if n_cells_after < cfg.minCellsAfterQc:
        raise FileRejected(SkipReason.too_few_cells, qc=stats)
    if pct_cells_after < cfg.minPctCellsAfterQc:
        raise FileRejected(SkipReason.excessive_cell_dropout, qc=stats)

    return adata, stats
