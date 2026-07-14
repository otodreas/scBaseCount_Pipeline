# pytest suite

Tests for the packages under `scripts/`. Run everything from the repo root with `uv run pytest`.

## Test files

| File | Purpose |
|------|---------|
| `test_imports.py` | Import smoke test: every package shipped by the wheel build target imports cleanly |
| `test_h5ad_concat_smoke.py` | In-memory smoke tests for `scripts/h5ad_concat/` that need no R2 download or upload |
| `test_cyteonto_payload.py` | Unit test for `cyteonto.payload.build_payload` label deduplication |
| `test_cluster_validation_regression.py` | Config snapshot test (always on) and golden regression test (opt-in) for `scripts/cluster_validation/` |
| `_clval_capture.py` | Field list and helpers shared by the cluster_validation regression test |

## h5ad_concat smoke tests

`test_h5ad_concat_smoke.py` mirrors the `run_h5ad_concat` data flow from `scripts/h5ad_concat/pipeline.py`. It mocks only the R2/IO boundary (`download_from_r2` writes prebuilt adata to the cache path; `load_contexts_jsonl` returns in-memory contexts) so `prepare_adata` -> `apply_qc_gate` -> `concat_atlas` -> `write_atlas` all run for real. One `_make_adata` helper builds every fixture.

What is checked:

- `run_h5ad_concat`: happy path (two files concatenated, atlas written, unique accession-suffixed barcodes), skip-and-continue when one file is rejected (`cell_type_all_missing`), and `ValueError` when every file is rejected
- `flag_qc_genes`: classifies `mt` / `ribo` / `hb` (case-insensitive, anchored HB regex), prefers `gene_symbols` over `var_names`
- `apply_qc_gate`: filters low-gene cells, optional `maxPctHb` ceiling, `FileRejected(SkipReason.too_few_cells)` when too few cells remain, optional `excessive_cell_dropout` when `maxPctCellsDropped` is set, returns `QcStats`
- `prepare_adata`: download errors map to `SkipReason` (`md5_mismatch` vs `download_failed`) and raw cache paths are deleted

Run it:

```sh
uv run pytest tests/pytest/test_h5ad_concat_smoke.py -q
```

The [`pre-push`](../../.githooks/pre-push) hook runs this file automatically when commits being pushed touch `scripts/h5ad_concat/` or the test itself. No env flag is needed.

## cluster_validation regression tests

Golden regression for `scripts/cluster_validation/`. Compares deterministic `ClusterValidationResult` fields against committed JSON baselines so clustering behavior stays stable when the default weak prior (`weakPriorKey="cell_type"`) is used.

### Layout

| Path | Purpose |
|------|---------|
| `test_cluster_validation_regression.py` | Config snapshot test (always on) and golden regression test (opt-in) |
| `_clval_capture.py` | Field list and helpers shared by the regression test |
| `data/SRX12708356.h5ad` | Committed fixture (~1k cells) |
| `baselines/clval_SRX12708356.json` | Expected deterministic result fields for the fixture |
| `baselines/clval_config_snapshot.json` | Expected default `ClusterValidationConfig.model_dump()` |

### Why numeric regression instead of figure diff

Reference plots under `docs/exploratory_depreciated/cluster_validation_sandbox/.figs/` came from the manual sandbox notebook, not the package. They are useful visually but not a strict oracle.

Figures also stamp a runtime timestamp (`_now()` in `viz.py`), so PNGs are never byte-identical even when clustering is unchanged.

The regression compares deterministic fields from `ClusterValidationResult` that should be stable for the default prior.

### What is checked

24 fields from `_clval_capture.ALL_FIELDS`:

- Scalars: `selectedResolution`, `cumvar`, `nPcs`, `kPrior`, `kFiltered`, `nCellsDropped`, `nCellsFinal`, `nClustersPreMerge`, `nClustersPostMerge`, `clusterKey`, `mergedKey`
- Arrays: `resolutions`, `kArr`, `jaccArr`, `silhouetteArr`, `homogeneityArr`, `completenessArr`, `nmiArr`, `vscoreArr`, `ariArr`, `confMatrix`, `confClasses`
- Maps: `labelMap`, `mergedGroups`

Branch-only metadata (`runTag`, `weakPriorKey`, `adataPath`, etc.) is excluded. Float leaves use `pytest.approx` with `rel=1e-6`, `abs=1e-9`.

The config snapshot test catches drift in default hyperparameters and repo-relative paths without running clustering.

### Environment variables

These are test-only flags. They are not listed in `.env.example` and do not need to live in `.env`.

| Variable | When to set | Effect |
|----------|-------------|--------|
| `RUN_CLVAL_REGRESSION=1` | Local runs or pre-push when you want the full golden test | Enables `test_cluster_validation_matches_baseline` (slow; runs Leiden sweep + RF merge on the fixture) |
| `UPDATE_CLVAL_BASELINES=1` | After an intentional change to config defaults or clustering behavior | Regenerates both `baselines/clval_config_snapshot.json` and `baselines/clval_SRX12708356.json` (the latter only when `RUN_CLVAL_REGRESSION=1` is also set, since that test is opt-in) |

Unset these variables for normal test runs and CI.

### Running

From the repo root:

```sh
# Fast: config snapshot only (same as CI)
uv run pytest tests/pytest/test_cluster_validation_regression.py -q

# Full golden regression
RUN_CLVAL_REGRESSION=1 uv run pytest tests/pytest/test_cluster_validation_regression.py -q
```

The [`pre-push`](../../.githooks/pre-push) hook sets `RUN_CLVAL_REGRESSION=1` automatically when commits being pushed touch `scripts/cluster_validation/`.

Regenerate both baselines after verifying the new output is correct:

```sh
RUN_CLVAL_REGRESSION=1 UPDATE_CLVAL_BASELINES=1 uv run pytest tests/pytest/test_cluster_validation_regression.py -q
git add tests/pytest/baselines/
```

Config snapshot only (no clustering run):

```sh
UPDATE_CLVAL_BASELINES=1 uv run pytest tests/pytest/test_cluster_validation_regression.py::test_cluster_validation_config_matches_snapshot -q
git add tests/pytest/baselines/clval_config_snapshot.json
```

### What this does not check

- Visual parity with sandbox `.figs/` (timestamped plots, different notebook origin)
- Behavior with a non-default weak prior (e.g. CellTypist `predicted_labels`)
- The file-loading entry point `run_cluster_validation()` specifically (the test uses `run_cluster_validation_on_adata()` on the committed h5ad, which runs the same pipeline steps)
