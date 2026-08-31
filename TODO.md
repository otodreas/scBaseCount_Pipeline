# Report wrap-up

- week 1
   - [ ] code review
   - [ ] atlas run
   - [ ] git tag
   - [ ] clean readme
- week 2
   - [ ] report review, update figs
   - [ ] powerpoint presentation


# atlas wrap-up
- [ ] rf merge on atlas, harden leiden sweep methodology (match with small proof of concept runs on individual datasets)
- [ ] make note of missing/deprioritized doublet detection
- [ ] get n cells dropped bc low gene count, n cells dropped bc high mito
- [ ] identify key pathway on atlas to drive home point: see `DailyNotes.md` line 1029
- [ ] mention in the slide deck that rf merge is not included in cluster resolution optimization because the leiden sweep methods tend to overcluster relative to cell_type rather than the latter.

# other
- [x] Set `verify_md5 = True` in `download_from_r2`[, line 34](scripts/storage/r2.py) to check that downloaded `h5ad` files' md5 hash matches the value stored in metadata on upload
- [x] Move atlas Harmony workflow into [`atlas_postprocessing`](scripts/atlas_postprocessing) with twin runners [`select_atlas_parameters.py`](pipelines/select_atlas_parameters.py) and [`run_atlas_postprocessing.py`](pipelines/run_atlas_postprocessing.py)
- [ ] Clarify embedding names in atlas postprocessing. These are confusing, should be consistent
```py
["X_pca", "X_pca_harmony", "X_umap", "X_umap_uncorrected"]
```
- [x] rsync 5k sub atlas to local, use to run atlas pipeline locally
- [ ] **build atlas v2 ([plan](.cursor/plans/include_healthy_atlas_v2_e8ae0bcd.plan.md))**
    - if a single SRX_accession passes (its a lung related sample), then all of its companion datasets that belong to the same study (study accession) should also pass. they can be ruled out later due to poor quality
- [x] Leiden / parameter sweeping for atlas postprocessing (calibration + approved JSON import)
- [ ] ideally, we would do all our ontology pruning BEFORE making the atlas. any pruning that happens during atlas build would involve genuinely fine grained qc stuff. a lot of information like n_obs is present in the metadata pq, along with tissue and disease ontology.
- [ ] Materialize gene symbols upstream in atlas build (`h5ad_concat` / postprocessing). Today `var_names` are Ensembl IDs only with no symbol column in `var`; CyteType failed on [`run_cytetype_on_atlas_lite.py`](pipelines/run_cytetype_on_atlas_lite.py) without remapping. Map from STAR [`geneInfo.tab`](data/scbasecount/2026-01-12/star_references/Homo_sapiens/hg38_2020/geneInfo.tab) during concat so the atlas carries human-readable symbols and downstream tools do not need per-pipeline remaps.

# Clean full-gene atlas postprocessing

## Problem and target

The concatenated atlas enters `run_atlas_postprocessing.py` with full-gene counts in `.X`. `load_and_normalize()` currently copies those counts to `.raw`, normalizes `.X`, then `select_hvgs()` subsets the main object before `scale_and_pca()` overwrites `.X` with scaled HVG expression. The current production artifact is therefore 9,310,348 cells x 2,000 HVGs in `.X` with 36,601-gene counts in `.raw`. CyteType requires full-gene log1p expression in `.X`, so `run_cytetype_on_atlas_lite.py` has to materialize `.raw`, normalize it again, copy embeddings, run annotation, and copy results back.

Rebuild the atlas with this fixed AnnData contract:

- `.X`: sparse, full-gene library-size-normalized `log1p` expression, shape `nCells x nGenes`. This is the direct CyteType and marker-ranking input.
- `.layers["counts"]`: unchanged sparse counts with exactly the same shape, cell order, and gene order as `.X`. Validate that the source `.X` is finite, non-negative, and integer-valued before normalization.
- `.raw`: unset. Counts belong in the explicit `counts` layer, and a duplicate full-gene raw snapshot would be ambiguous and expensive.
- `.var`: the complete gene table plus full-axis HVG statistics and `highly_variable`; exactly `nTopGenes` rows should be selected.
- Scaled HVG expression: use an HVG-only temporary AnnData for scaling and PCA, then discard its dense `.X`. Do not put it in a layer because AnnData layers must have the full `.X` shape. Preserve the relevant scale means/standard deviations on the HVG rows of full `.var`, with missing values for non-HVGs.
- `.obsm["X_pca"]`: PCA scores from the scaled HVG workspace. `.varm["PCs"]` must have `nGenes x nPcsCompute`; copy HVG loadings into their full-gene rows and use zero loadings for genes excluded from PCA. Keep PCA variance metadata in `.uns["pca"]`.
- `.obsm["X_pca_harmony"]`, `.obsm["X_umap"]`, `.obs["leiden_atlas"]`, and the active neighbor graph remain the production integration outputs. Validation may additionally retain `X_umap_uncorrected` and `leiden_uncorrected`.

## Implementation

1. Refactor `scripts/atlas_postprocessing/core.py`.
   - Make `load_and_normalize()` validate the counts input, copy counts to `layers["counts"]`, normalize/log1p the full `.X`, retain `uns["log1p"]`, and never populate `.raw`.
   - Split the current mutating `select_hvgs()` contract into full-axis HVG marking and creation of an HVG-only PCA workspace. Remove the counts layer from temporary workspaces before scaling to avoid unnecessary copies.
   - Make `prepare_pca()` transfer PCA scores, variance, full-axis loadings, and HVG scale statistics back to the canonical full-gene object without replacing its `.X`.
   - Run uncorrected/Harmony neighbors, UMAP, and Leiden on that canonical object after PCA transfer. These stages consume `.obsm`, so their numerical behavior should not depend on whether `.X` remains full-gene.
   - Update every call site in `scripts/atlas_postprocessing/selection.py`: sweeps still operate on normalized full-gene input but create isolated HVG workspaces for each candidate. Keep approved parameter selection numerically equivalent while avoiding repeated copies of the counts layer.

2. Align runner, config, and artifacts.
   - In `scripts/atlas_postprocessing/config.py`, make the production output default agree with the real artifact name `output/atlas/v2/post/production/atlas_v2_post.h5ad`; retain explicit CLI overrides and keep upload opt-in.
   - In `pipelines/run_atlas_postprocessing.py`, add an early input-contract check and log the resolved expression layout and approved parameter source before the expensive run. Keep the approved JSON authoritative for the four calibrated knobs.
   - In `scripts/atlas_postprocessing/core.py::save_atlas()`, report `genes`, `nHighlyVariable`, `xRepresentation`, `countsLayer`, and the relevant shapes instead of inferring gene counts from `.raw` or treating `adata.n_vars` as the HVG count.
   - Leave `AtlasPostprocessingParameters` and calibration JSON focused on algorithm parameters. Do not add a layout version field or create parallel model versions. Enforce the single current layout through validation and explicit manifest fields.
   - Update the atlas postprocessing section of `pipelines/README.md` to describe the new public slot contract and rebuild command.

3. Remove downstream conversion and copy-back paths.
   - Refactor `pipelines/run_cytetype_on_atlas_lite.py` into a callable `main()`. Open the rebuilt atlas, run `rank_genes_groups_backed(..., use_raw=False)` directly against full-gene log1p `.X`, pass the same object to CyteType, and write the annotated output. Remove `.raw.to_adata()`, repeat normalization, embedding copies, result transfer, and `INPUT.unlink()` so the production atlas is never deleted.
   - Update `notebooks/analysis/analyze_atlas_DE.ipynb` to stop rebuilding a full-gene object from `.raw`: use `.X` for full-gene marker ranking and explicitly select `layers["counts"]` for pseudobulk sums. Continue using the existing cluster and embedding keys.
   - Check other consumers that only use `.obs`, `.obsm`, `.obsp`, or `.uns` such as scIB and plotting. They should require no behavior change, but their smoke tests must prove that the full-gene shape does not alter graph outputs.

## Rebuild and migration

- Do not transform the existing HVG-only production file in place. Re-run from `output/atlas/v2/atlas_v2.h5ad` with the already approved `output/atlas/v2/post/parameter_selection/approved_parameters.json`, writing a temporary sibling output first.
- Validate the temporary artifact in backed mode, then replace the production artifact only after all checks pass. Preserve the old artifact until the replacement is verified. Upload remains a separate explicit operator action.
- The accepted artifact must have the source cell/gene identities, exact counts in `layers["counts"]`, full-gene log1p `.X`, the approved HVG count, expected PCA/Harmony/UMAP shapes, and non-null `leiden_atlas`. The approved calibration does not need a new version because the HVG/PCA computation is unchanged.
- Re-run the CyteType atlas entry point on the rebuilt file. Never copy annotations into the obsolete layout.

## Tests and verification

- Replace `tests/pytest/test_atlas_harmony_raw.py` with layout tests proving exact count preservation, independent normalization/log1p equality, full-gene `.X`, no `.raw`, and valid same-shaped layers after an h5ad round trip.
- Extend `tests/pytest/test_atlas_postprocessing.py` to prove that PCA leaves canonical `.X` and counts unchanged; HVG metadata, scale statistics, PCA scores/loadings, Harmony, neighbors, UMAP, and Leiden have valid full-axis shapes; non-HVG loadings are zero; and production versus validation keys remain correct.
- Update `tests/pytest/test_atlas_postprocessing_sampling.py` and parameter-selection tests for the new normalization/HVG helper contracts and manifest fields. Confirm sampling metadata survives.
- Add a small backed-file CyteType smoke test that ranks from `.X` with `use_raw=False` and reaches CyteType validation without constructing a second AnnData or making a network request.
- Run the focused pytest files with `uv run pytest`, then the full pytest suite. Finally inspect the rebuilt production h5ad in backed mode and compare selected PCA/Harmony outputs against a tiny fixture produced by the old computation path within numerical tolerance.

## PR strategy

**Stack:** multi-layer
**Rationale:** The change has real dependencies from the AnnData layout foundation, through production/calibration wiring, to CyteType and analysis consumers, so separate reviewable layers reduce migration risk.
**Shape:**
1. `Define full-gene atlas layout`: `scripts/atlas_postprocessing/core.py`, config/artifact helpers, and core layout tests.
2. `Wire clean atlas production flow`: both atlas pipeline runners, parameter-selection call sites, manifests, runner tests, and `pipelines/README.md`.
3. `Adopt clean atlas in consumers`: `pipelines/run_cytetype_on_atlas_lite.py`, the atlas DE notebook, and downstream smoke tests.