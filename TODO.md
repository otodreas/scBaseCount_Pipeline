- [x] Set `verify_md5 = True` in `download_from_r2`[, line 34](scripts/storage/r2.py) to check that downloaded `h5ad` files' md5 hash matches the value stored in metadata on upload
- [ ] Refactor [`run_atlas_harmony.py`](pipelines/run_atlas_harmony.py) to move functionality to [`cluster_validation`](scripts/cluster_validation) where appropriate. there could be an atlas mode
- [ ] Clarify embedding names in [`run_atlas_harmony.py`](pipelines/run_atlas_harmony.py). These are confusing, should be consistent
```py
["X_pca", "X_pca_harmony", "X_umap", "X_umap_uncorrected"]
```