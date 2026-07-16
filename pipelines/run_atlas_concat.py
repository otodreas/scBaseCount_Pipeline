from h5ad_concat import H5adConcatConfig, run_h5ad_concat

cfg = H5adConcatConfig(
    atlasR2Key="lung/atlas_v1.h5ad",
    uploadAtlas=True,
)

run_h5ad_concat(cfg)
