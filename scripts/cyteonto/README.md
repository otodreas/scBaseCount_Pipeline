# cyteonto

Submits cell type annotation labels to the [CyteOnto API](https://cyteonto.nygen.io), polls until the run completes, fetches the result as a CSV, and returns a pandas DataFrame. Designed to be called from `notebooks/cytetype.ipynb` after CyteType has annotated the dataset.

## Usage

```python
from cyteonto import CyteOntoConfig, run_cyteonto
from pathlib import Path

cfg = CyteOntoConfig(
    h5adPath=Path("output/cytetype/data/SRX17412841_cytetype_annotated.h5ad"),
)

similarities = run_cyteonto(cfg)
```

`run_cyteonto` returns a `pandas.DataFrame` with one row per `(algorithm, cell)` pair. The CSV is also written to `output/cyteonto/runs/{run_id}.csv` and every step is appended to `logs/cyteonto.log`.

## Input conventions

The pipeline reads two fixed columns from `adata.obs`:

| Column | Role in payload |
|--------|----------------|
| `cell_type` | `authorLabels` -- the STATE reference annotation |
| `cytetype_annotation_leiden_merged` | `algorithms["algo1"]` -- the CyteType annotation |

## Pipeline steps

| Step | Module | What happens |
|------|--------|--------------|
| Load | `pipeline.py` | Read h5ad with `scanpy` |
| Build payload | `payload.py` | Extract the two obs columns into the API request dict |
| Write payload | `payload.py` | Serialize to `output/cyteonto/payloads/{stem}_annotations.json` |
| Submit | `client.py` | POST `/compare` to the CyteOnto API; receive `run_id` |
| Poll | `client.py` | GET `/status/{run_id}` on `pollIntervalS` cadence until `completed` or `failed` |
| Fetch | `client.py` | GET `/result/{run_id}?format=csv`; save to `output/cyteonto/runs/{run_id}.csv` |

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `h5adPath` | required | Path to the annotated h5ad file |
| `payloadDir` | `output/cyteonto/payloads` | Directory for the payload JSON |
| `resultsDir` | `output/cyteonto/runs` | Directory for the fetched CSV |
| `baseUrl` | `https://cyteonto.nygen.io` | CyteOnto service base URL |
| `pollIntervalS` | `10` | Seconds between status polls |
| `pollTimeoutS` | `3600` | Total seconds before a `TimeoutError` is raised |

## Output DataFrame columns

| Column | Type | Description |
|--------|------|-------------|
| `run_id` | str | CyteOnto run identifier |
| `algorithm` | str | Algorithm key (`algo1`) |
| `pair_index` | int | Row position in the label lists |
| `author_label` | str | STATE reference label for this cell |
| `algorithm_label` | str | CyteType label for this cell |
| `author_ontology_id` | str | Best Cell Ontology match for the author label |
| `author_embedding_similarity` | float | Cosine similarity to that CL term |
| `algorithm_ontology_id` | str | Best Cell Ontology match for the algorithm label |
| `algorithm_embedding_similarity` | float | Cosine similarity to that CL term |
| `cytescore_similarity` | float | Ontology-aware agreement score between the two labels |
| `similarity_method` | str | Scoring method used (`cytescore`, `string_similarity`, etc.) |

## Logging

Every run appends to `logs/cyteonto.log`:

```
2026-04-24 18:24:01 INFO start  input=output/cytetype/data/SRX17412841_cytetype_annotated.h5ad
2026-04-24 18:24:03 INFO loaded 12847 cells  24280 genes
2026-04-24 18:24:03 INFO payload  author_labels=12847  algorithms=1
2026-04-24 18:24:03 INFO payload written  path=output/cyteonto/payloads/...json
2026-04-24 18:24:04 INFO submitted  runId=run-<uuid>  state=queued
2026-04-24 18:24:14 INFO status    runId=run-<uuid>  state=running
2026-04-24 18:26:01 INFO status    runId=run-<uuid>  state=completed
2026-04-24 18:26:01 INFO completed  runId=run-<uuid>  rows=12847
2026-04-24 18:26:02 INFO fetched   runId=run-<uuid>  path=output/cyteonto/runs/run-<uuid>.csv
2026-04-24 18:26:02 INFO done  saved=output/cyteonto/runs/run-<uuid>.csv
```
