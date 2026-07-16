from h5ad_concat import H5adConcatConfig, run_h5ad_concat

# TODO: confirm contexts on server match local
cfg = H5adConcatConfig(
    atlasR2Key="lung/atlas_v1.h5ad",
)

run_h5ad_concat(cfg)
