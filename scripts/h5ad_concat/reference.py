from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import pandas as pd
from anndata._core.merge import gen_reindexer

from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason


@dataclass
class GeneReference:
    ids: list[str]
    var: pd.DataFrame


@dataclass
class AlignStats:
    nGenesMapped: int
    nGenesZeroFilled: int
    nGenesDropped: int
    droppedVarKeys: list[str]


def load_gene_reference(path: Path) -> GeneReference:
    """Load the STAR geneInfo.tab reference and return a GeneReference."""
    table = pd.read_csv(
        path,
        sep="\t",
        skiprows=1,
        header=None,
        names=["ensembl_id", "gene_symbol", "biotype"],
    )
    ids = table["ensembl_id"].astype(str).tolist()
    var = pd.DataFrame(table.set_index("ensembl_id")[["gene_symbol", "biotype"]])
    return GeneReference(ids=ids, var=var)


def align_to_reference(
    adata: ad.AnnData,
    reference: GeneReference,
    *,
    conserve_layers: bool = False,
) -> tuple[ad.AnnData, AlignStats]:
    """Reindex adata onto the canonical reference gene axis and raise FileRejected on zero overlap.

    anndata's Reindexer maps the file's gene columns onto the reference axis in one sparse matmul:
    it reorders shared genes, drops file genes absent from the reference, and zero-fills reference
    genes absent from the file. The file's own var columns are indexed on its gene axis and cannot
    survive the reindex, so the canonical reference annotations, which already line up with the
    reindexed matrix, replace them. reference.var is shared across files, hence the copy.
    AlignStats.droppedVarKeys records the discarded file var columns (per-file QC stats and gene
    annotations absent from the reference) so callers can log what the reindex throws away.

    By default only X is carried over. Set conserve_layers to also reindex every layer of adata
    (for example the STARsolo UniqueAndMult matrices) onto the reference axis and keep them.
    """
    reindexer = gen_reindexer(reference.var.index, adata.var_names)
    n_mapped = len(reindexer.new_pos)
    if n_mapped == 0:
        raise FileRejected(SkipReason.gene_axis_mismatch)

    new_x = reindexer(adata.X, fill_value=0)
    dropped_var_keys = [key for key in adata.var.columns if key not in reference.var.columns]
    aligned = ad.AnnData(X=new_x, obs=pd.DataFrame(adata.obs), var=reference.var.copy())
    if conserve_layers:
        for name, layer in adata.layers.items():
            aligned.layers[name] = reindexer(layer, fill_value=0)
    return aligned, AlignStats(
        nGenesMapped=n_mapped,
        nGenesZeroFilled=len(reference.ids) - n_mapped,
        nGenesDropped=adata.n_vars - n_mapped,
        droppedVarKeys=dropped_var_keys,
    )
