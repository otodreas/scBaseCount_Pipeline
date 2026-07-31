import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix


def test_raw_counts_survive_normalize_and_hvg_subset() -> None:
    """Mirror atlas postprocessing: full counts in .raw, HVG-only .X after subset."""
    n_obs, n_vars = 12, 40
    rng = np.random.default_rng(0)
    counts = csr_matrix(rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32))
    adata = sc.AnnData(
        X=counts,
        obs=pd.DataFrame(
            {
                "study_accession": ["A"] * 6 + ["B"] * 6,
            },
            index=pd.Index([f"cell_{i}" for i in range(n_obs)]),
        ),
    )
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    raw_counts = counts.copy()

    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=15, batch_key="study_accession")
    adata = adata[:, adata.var["highly_variable"]].copy()

    assert adata.n_vars == 15
    assert adata.raw is not None
    assert adata.raw.n_vars == n_vars

    full = adata.raw.to_adata()
    assert list(full.var_names) == [f"gene_{i}" for i in range(n_vars)]
    assert (raw_counts != full.X).nnz == 0
