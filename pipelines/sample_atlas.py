"""This file is made obsolete by the select_atlas_parameters pipeline, which
calls atlas_postprocessing.sampling.sample_study_proportional.
"""

from pathlib import Path

import scanpy as sc

FULL_ADATA_PATH: Path = Path(__file__).resolve().parents[1] / "output/atlas/v1/processed_1/atlas_harmony.h5ad"
SAMPLE_ADATA_PATH: Path = (
    Path(__file__).resolve().parents[1] / "output/atlas/v1/processed_1/atlas_harmony_sub_100k.h5ad"
)
N_OBS_SAMPLE: int = 100_000


adata = sc.read_h5ad(FULL_ADATA_PATH, backed="r")
print(f"Loaded full adata: {adata.n_obs} observations")
sub = sc.pp.sample(adata, n=N_OBS_SAMPLE, copy=True)
print(f"Subsampled adata: {sub.n_obs} observations, saving to {SAMPLE_ADATA_PATH}")
sub.write_h5ad(SAMPLE_ADATA_PATH, compression="gzip")
