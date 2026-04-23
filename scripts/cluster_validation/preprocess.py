from __future__ import annotations

from dataclasses import dataclass

import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


@dataclass
class PreprocessStats:
    nCellsOriginal: int
    kPrior: int
    nDropped: int
    kFiltered: int


def preprocess(
    adata: sc.AnnData, cfg: ClusterValidationConfig
) -> tuple[sc.AnnData, PreprocessStats]:
    n_cells_original = adata.n_obs
    k_prior = adata.obs["cell_type"].nunique()

    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    adata.var["ribo"] = adata.var["gene_symbols"].str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

    type_counts = adata.obs["cell_type"].value_counts()
    valid_types = type_counts[type_counts >= cfg.minCellsPerType].index
    adata = adata[adata.obs["cell_type"].isin(valid_types)].copy()

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.nTopGenes)

    k_filtered = adata.obs["cell_type"].nunique()
    n_dropped = n_cells_original - adata.n_obs

    return adata, PreprocessStats(
        nCellsOriginal=n_cells_original,
        kPrior=k_prior,
        nDropped=n_dropped,
        kFiltered=k_filtered,
    )
