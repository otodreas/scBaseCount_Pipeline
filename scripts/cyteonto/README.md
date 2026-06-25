# cyteonto

Submits cell type annotation labels to the [CyteOnto API](https://cyteonto.nygen.io), polls `/result/` until the CSV is ready, and returns a pandas DataFrame. Designed to be called from `notebooks/pipeline/cyteonto.ipynb` after CyteType has annotated the dataset.

## Usage

```python
from cyteonto import CyteOntoConfig, run_cyteonto, check_pending_runs
from pathlib import Path

cfg = CyteOntoConfig(
    h5adPath=Path("output/cytetype/data/SRX17412841_cytetype_annotated.h5ad"),
)

similarities = run_cyteonto(cfg)
```

Multi-algorithm comparison (for example CellTypist and CyteType against cxg author labels):

```python
cfg = CyteOntoConfig(
    h5adPath=Path("output/celltypist_vs_cxg/data/SRX17412841_cytetype_annotated.h5ad"),
    authorCol="cell_type",
    algorithmCols={
        "celltypist": "predicted_labels",
        "cytetype": "cytetype_annotation_leiden_merged",
    },
)
similarities = run_cyteonto(cfg)
```

`run_cyteonto` returns a `pandas.DataFrame` with one row per unique label combination per algorithm. The CSV is also written to `output/cyteonto/runs/{run_id}.csv` and every step is appended to `logs/cyteonto.log`.

### Interrupting a run

Runs can take a long time. You can safely interrupt the polling loop (e.g. shut down the IDE) without losing the run -- the job continues on the server. On interrupt, a message is written to the log and `run_cyteonto` returns `None`:

```
2026-04-29 14:31:00 INFO polling stopped  runId=run-<uuid>  (run continues on server; call check_pending_runs() to resume)
```

A stub file `output/cyteonto/runs/{run_id}.json` is written immediately after the job is submitted with `completedAt` set to `null`, so it is never lost. When the run finishes the same file is updated in place with a `completedAt` timestamp, and the result CSV lands alongside it at `output/cyteonto/runs/{run_id}.csv`.

### Resuming after a restart

Call `check_pending_runs()` at the top of your session. It scans `output/cyteonto/runs/` for any stub whose `completedAt` is `null`, probes `/result/` for each, fetches any that return HTTP 200, and returns a dict of `{run_id: DataFrame}`:

```python
from cyteonto import check_pending_runs

results = check_pending_runs()
# results is a dict[run_id, pd.DataFrame] for every run that completed
# each run stub is updated in place -- completedAt is stamped on completion
```

## Input conventions

By default the pipeline reads two columns from `adata.obs`:

| Column | Role in payload |
|--------|----------------|
| `cell_type` | `authorLabels` -- the CELLxGENE author reference annotation |
| `cytetype_annotation_leiden_merged` | `algorithms["algo1"]` -- the CyteType annotation |

Override with `authorCol` and `algorithmCols` on `CyteOntoConfig`.

The payload is deduplicated before submission: only unique combinations of `(authorCol, *algorithmCols)` are sent to CyteOnto. Map cytescores back to cells with `attach_cytescores_to_obs`.

## Pipeline steps

| Step | Module | What happens |
|------|--------|--------------|
| Load | `pipeline.py` | Read h5ad with `scanpy` |
| Build payload | `payload.py` | Extract unique label combinations from obs into the API request dict |
| Write payload | `payload.py` | Serialize to `output/cyteonto/payloads/{stem}_annotations.json` |
| Submit | `client.py` | POST `/compare` to the CyteOnto API; receive `run_id` |
| Poll result | `client.py` | GET `/result/{run_id}?format=csv` on `pollIntervalS` cadence; HTTP 409 means still running, HTTP 200 saves CSV to `output/cyteonto/runs/{run_id}.csv` |

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `h5adPath` | required | Path to the annotated h5ad file |
| `authorCol` | `cell_type` | `obs` column used as CyteOnto author labels |
| `algorithmCols` | `{"algo1": "cytetype_annotation_leiden_merged"}` | Map of algorithm key to `obs` column |
| `payloadDir` | `output/cyteonto/payloads` | Directory for the payload JSON |
| `runsDir` | `output/cyteonto/runs` | Directory for run stubs (JSON) and result CSVs |
| `baseUrl` | `https://cyteonto.nygen.io` | CyteOnto service base URL |
| `pollIntervalS` | `10` | Seconds between result polls |
| `pollTimeoutS` | `3600` | Total seconds before a `TimeoutError` is raised |

All default paths are relative to the repo root.

## Output DataFrame columns

| Column | Type | Description |
|--------|------|-------------|
| `run_id` | str | CyteOnto run identifier |
| `algorithm` | str | Algorithm key (`algo1`, `celltypist`, `cytetype`, etc.) |
| `pair_index` | int | Row position in the deduplicated label lists |
| `author_label` | str | Author reference label for this combination |
| `algorithm_label` | str | Algorithm label for this combination |
| `author_ontology_id` | str | Best Cell Ontology match for the author label |
| `author_embedding_similarity` | float | Cosine similarity to that CL term |
| `algorithm_ontology_id` | str | Best Cell Ontology match for the algorithm label |
| `algorithm_embedding_similarity` | float | Cosine similarity to that CL term |
| `cytescore_similarity` | float | Ontology-aware agreement score between the two labels |
| `similarity_method` | str | Scoring method used (`cytescore`, `string_similarity`, etc.) |

Per-cell analysis should join these rows onto `adata.obs` with `attach_cytescores_to_obs`, not treat one CSV row as one cell.

## Logging

Every run appends to `logs/cyteonto.log`:

```
2026-04-24 18:24:01 INFO start  input=output/cytetype/data/SRX17412841_cytetype_annotated.h5ad
2026-04-24 18:24:03 INFO loaded 12847 cells  24280 genes
2026-04-24 18:24:03 INFO payload  author_labels=12847  algorithms=1
2026-04-24 18:24:03 INFO payload written  path=output/cyteonto/payloads/...json
2026-04-24 18:24:04 INFO submitted  runId=run-<uuid>  state=queued
2026-04-24 18:26:01 INFO fetched   runId=run-<uuid>  path=output/cyteonto/runs/run-<uuid>.csv
2026-04-24 18:26:01 INFO completed  runId=run-<uuid>  rows=12847
2026-04-24 18:26:02 INFO done  saved=output/cyteonto/runs/run-<uuid>.csv
```

To manually check whether a run is ready, search for `runId` in `logs/cyteonto.log`, then:

```sh
export CYTEONTO_URL="https://cyteonto.nygen.io"
export RUN_ID="run-<uuid>"
curl -sS "$CYTEONTO_URL/result/$RUN_ID?format=csv" -o /tmp/result.csv -w '%{http_code}\n'
```

HTTP 200 means the CSV is ready; HTTP 409 means still processing.
