# h5ad_extractor

Extracts annotation columns from `adata.obs` or `adata.var` in an h5ad file and writes them to a parquet or CSV file.

## Usage

```python
from h5ad_extractor import H5adExtractConfig, extract_annotation_columns
from pathlib import Path

cfg = H5adExtractConfig(
    h5adPath=Path("output/cytetype/data/SRX12345678_cytetype_annotated.h5ad"),
    columnNames=["cell_type", "cytetype_annotation_leiden_merged"],
)

out_path = extract_annotation_columns(cfg)
```

`extract_annotation_columns` returns the path to the written file.

## Config reference

| Field | Default | Description |
|-------|---------|-------------|
| `h5adPath` | required | Path to the h5ad file (absolute or relative to repo root) |
| `columnNames` | required | List of column names to extract (minimum 1) |
| `outputDir` | `output/h5ad_extract` | Directory for the output file when `outputPath` is not set |
| `outputPath` | `None` | Explicit output path; overrides `outputDir` and auto-naming when set |
| `annotationAxis` | `"obs"` | Which axis to read columns from (`"obs"` or `"var"`) |
| `outputFormat` | `"parquet"` | Output format (`"parquet"` or `"csv"`) |

All default paths are relative to the repo root.

## Output

The default output filename is `{h5ad_stem}_{obs|var}_columns.{parquet|csv}` written to `outputDir`.
