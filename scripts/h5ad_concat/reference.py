from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason


@dataclass
class GeneReference:
    ids: list[str]
    var: pd.DataFrame
    position: dict[str, int]


@dataclass
class AlignStats:
    nGenesMapped: int
    nGenesZeroFilled: int
    nGenesDropped: int


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
    position = {gene_id: idx for idx, gene_id in enumerate(ids)}
    return GeneReference(ids=ids, var=var, position=position)


def align_to_reference(adata: ad.AnnData, reference: GeneReference) -> tuple[ad.AnnData, AlignStats]:
    """Reindex adata to the canonical reference gene axis; raise FileRejected on zero overlap."""
    file_ids = adata.var_names.astype(str).tolist()
    if file_ids == reference.ids:
        aligned = ad.AnnData(
            X=adata.X,
            obs=pd.DataFrame(adata.obs),
            var=reference.var.copy(),
        )
        return aligned, AlignStats(nGenesMapped=len(file_ids), nGenesZeroFilled=0, nGenesDropped=0)

    kept_ids: list[str] = []
    ref_positions: list[int] = []
    for gene_id in file_ids:
        pos = reference.position.get(gene_id)
        if pos is not None:
            kept_ids.append(gene_id)
            ref_positions.append(pos)

    n_mapped = len(kept_ids)
    if n_mapped == 0:
        raise FileRejected(SkipReason.gene_axis_mismatch)

    n_dropped = len(file_ids) - n_mapped
    n_zero_filled = len(reference.ids) - n_mapped

    row = np.arange(n_mapped, dtype=np.int64)
    col = np.asarray(ref_positions, dtype=np.int64)
    selector = sp.coo_matrix(
        (np.ones(n_mapped, dtype=np.float32), (row, col)),
        shape=(n_mapped, len(reference.ids)),
    ).tocsr()

    source_x = adata[:, kept_ids].X
    new_x = source_x @ selector

    aligned = ad.AnnData(X=new_x, obs=pd.DataFrame(adata.obs), var=reference.var.copy())
    return aligned, AlignStats(
        nGenesMapped=n_mapped,
        nGenesZeroFilled=n_zero_filled,
        nGenesDropped=n_dropped,
    )
