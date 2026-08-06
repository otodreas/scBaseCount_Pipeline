import os
import tempfile
from typing import cast

import h5py
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata.abc import CSRDataset
from anndata.io import write_elem
from cytetype import CyteType, rank_genes_groups_backed
from dotenv import load_dotenv
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.reference import load_gene_reference
from shared.repo import REPO_ROOT
from storage.r2 import upload_to_r2, verify_upload
from storage.transfer import _MD5_METADATA_KEY, _local_md5_b64

load_dotenv()

input_path = REPO_ROOT / "output/atlas/v2/post/production/atlas_v2_post.h5ad"
output_path = REPO_ROOT / "output/atlas/v2/post/production/cytetype/atlas_v2_post_cytetype.h5ad"
r2_key = "atlas/v2/post/production/cytetype/atlas_v2_post_cytetype.h5ad"
group_key = "leiden_atlas"
rank_key = f"rank_genes_{group_key}"

if output_path.exists():
    raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")

output_path.parent.mkdir(parents=True, exist_ok=True)

print(f"Reading input from {input_path}")
atlas = sc.read_h5ad(input_path, backed="r")
if atlas.raw is None:
    raise ValueError("Expected full-gene counts in adata.raw")
print("Loading raw counts, copying X, obs, var")
raw_counts = cast(CSRDataset, atlas.raw.X)
expr = AnnData(
    X=raw_counts.to_memory(),
    obs=cast(pd.DataFrame, atlas.obs).copy(),
    var=atlas.raw.var.copy(),
)
expr.obsm["X_umap"] = atlas.obsm["X_umap"].copy()
expr.layers["counts"] = raw_counts
gene_reference = load_gene_reference(H5adConcatConfig().geneInfoPath)
symbols = gene_reference.var.reindex(expr.var_names)["gene_symbol"]
missing_symbols = symbols.isna()
if bool(missing_symbols.any()):
    missing_ids = expr.var_names[missing_symbols.to_numpy()]
    raise ValueError(
        f"Missing gene symbols for {int(missing_symbols.sum())} Ensembl IDs in "
        f"data/scbasecount/2026-01-12/star_references/Homo_sapiens/hg38_2020/geneInfo.tab; "
        f"examples: {missing_ids[:5].tolist()}"
    )
expr.var["gene_symbol"] = symbols.to_numpy()
print(f"Mapped {expr.n_vars} Ensembl IDs to gene symbols")
sc.pp.normalize_total(expr)
sc.pp.log1p(expr)
print("Normalized and log-transformed expression")

rank_genes_groups_backed(
    expr,
    groupby=group_key,
    use_raw=False,
    key_added=rank_key,
)
print("Ranked genes completed")
with tempfile.TemporaryDirectory(prefix="cytetype-run-", dir=output_path.parent) as temp_dir:
    annotator = CyteType(
        expr,
        group_key=group_key,
        rank_key=rank_key,
        gene_symbols_column="gene_symbol",
        n_top_genes=100,
        auth_token=os.getenv("CYTETYPE_AUTH_TOKEN"),
        vars_h5_path=os.path.join(temp_dir, "vars.h5"),
        obs_duckdb_path=os.path.join(temp_dir, "obs.duckdb"),
    )

    try:
        expr = annotator.run(
            study_context=(
                "Integrated multi-study human lung single-cell RNA-seq atlas "
                "(Homo sapiens). Samples are primarily adult lung parenchyma and "
                "airway-related tissues (including bronchoalveolar lavage, bronchial "
                "epithelium, and related respiratory sites), with a minority of fetal "
                "lung samples. The atlas mixes healthy donors and lung disease "
                "contexts such as idiopathic pulmonary fibrosis, COVID-19 / SARS-CoV-2, "
                "COPD, interstitial lung disease, and lung cancer (including "
                "adenocarcinoma and NSCLC). Data were generated with 10x Genomics "
                "scRNA-seq (mostly 3' GEX, some 5' GEX)."
            ),
            save_query=False,
        )
    finally:
        annotator.cleanup()

result_obs_columns = [
    f"cytetype_annotation_{group_key}",
    f"cytetype_cellOntologyTerm_{group_key}",
    f"cytetype_cellOntologyTermID_{group_key}",
    f"cytetype_cellState_{group_key}",
]

output_obs = atlas.obs.copy()
for column in result_obs_columns:
    output_obs[column] = expr.obs[column]

output_uns = dict(atlas.uns)
output_uns[rank_key] = expr.uns[rank_key]
for key in ("cytetype_results", "cytetype_jobDetails"):
    output_uns[key] = expr.uns[key]

del expr
atlas.file.close()

input_path.replace(output_path)
with h5py.File(output_path, "r+") as output_file:
    write_elem(output_file, "obs", output_obs, dataset_kwargs={"compression": "gzip"})
    write_elem(output_file, "uns", output_uns, dataset_kwargs={"compression": "gzip"})

md5 = _local_md5_b64(output_path)
upload_to_r2(output_path, r2_key, extra_metadata={_MD5_METADATA_KEY: md5})
if not verify_upload(r2_key):
    raise RuntimeError(f"R2 upload verification failed for {r2_key}")
