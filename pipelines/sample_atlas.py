from pathlib import Path

import scanpy as sc

FULL_ADATA_PATH: Path = Path("output/atlas/v1/processed1/atlas_harmony.h5ad")
SAMPLE_ADATA_PATH: Path = Path("output/atlas/v1/processed1/atlas_harmony_sub_100k.h5ad")
N_OBS_SAMPLE: int = 100_000


adata = sc.read_h5ad(FULL_ADATA_PATH, backed="r")
sub = sc.pp.sample(adata, n=N_OBS_SAMPLE, copy=True)
sub.write_h5ad(SAMPLE_ADATA_PATH, compression="gzip")
