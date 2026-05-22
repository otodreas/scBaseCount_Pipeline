# cytetype_runner

Runs the CyteType annotation step on a clustered h5ad file and writes the annotated result to disk. Job details (`job_id`, `report_url`, `api_url`) are returned in the result and also embedded in the annotated h5ad under `adata.uns["cytetype_jobDetails"]`.

## Usage

```python
from cytetype_runner import CyteTypeRunnerConfig, run_cytetype

cfg = CyteTypeRunnerConfig(srxAccession="SRX12345678")
result = run_cytetype(cfg, input_path, group_key, study_context)
print(result.outputPath, result.reportUrl)
```

`run_cytetype` returns a `CyteTypeRunResult` with `outputPath` (path to the annotated h5ad), `reportUrl`, `jobId`, and `apiUrl`.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `srxAccession` | required | SRX accession identifier |
| `outputDir` | `output/cytetype/data` | Directory for the annotated h5ad file |

All default paths are relative to the repo root.

## Output

| File | Description |
|------|-------------|
| `output/cytetype/data/{srx}_cytetype_annotated.h5ad` | Annotated h5ad with CyteType labels in `adata.obs` and job details under `adata.uns["cytetype_jobDetails"]` |

## Logging

Steps are appended to `logs/cytetype_runner.log`.
