# cytetype_runner

Runs the CyteType annotation step on a clustered h5ad file and writes the annotated result to disk. Also extracts and persists the job details returned by the CyteType API.

## Usage

```python
from cytetype_runner import CyteTypeRunnerConfig, run_cytetype

cfg = CyteTypeRunnerConfig(srxAccession="SRX12345678")
cytetype_h5ad = run_cytetype(cfg, input_path, group_key, study_context)
```

`run_cytetype` returns the path to the written annotated h5ad file. Job details are written automatically alongside it.

### Extracting job details from an existing h5ad

```python
from cytetype_runner import CyteTypeRunnerConfig, write_job_details

cfg = CyteTypeRunnerConfig(srxAccession="SRX12345678")
cytetype_h5ad = cfg.outputDir / "SRX12345678_cytetype_annotated.h5ad"
write_job_details(cfg, cytetype_h5ad)
```

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `srxAccession` | required | SRX accession identifier |
| `outputDir` | `output/cytetype/data` | Directory for the annotated h5ad file |
| `jobDetailsDir` | `output/cytetype/job_details` | Directory for the job details JSON |

All default paths are relative to the repo root.

## Output

| File | Description |
|------|-------------|
| `output/cytetype/data/{srx}_cytetype_annotated.h5ad` | Annotated h5ad with CyteType labels in `adata.obs` |
| `output/cytetype/job_details/{srx}_cytetype_jobDetails.json` | `{ "<srx>": adata.uns["cytetype_jobDetails"] }` |

## Logging

Steps are appended to `logs/cytetype_runner.log`.
