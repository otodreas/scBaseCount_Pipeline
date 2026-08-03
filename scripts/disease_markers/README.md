# disease_markers

Per-experiment disease area, disease status, biological control, and atlas eligibility labels derived from `contexts.jsonl` and `atlas.csv`.

## Usage

```python
from pathlib import Path

from disease_markers.labels import build_sample_label_table, coarse_disease_area

label_table = build_sample_label_table(
    Path("output/context/contexts.jsonl"),
    Path("output/atlas/v1/atlas.csv"),
)
```

Join labels onto an AnnData object by experiment accession (for example `SRX_accession` in `obs`). Each label row applies to the full experiment H5AD.

## Outputs

Typical notebook or ad hoc export: `eligibility_labels.csv` with columns `srxAccession`, `studyAccession`, `diseaseRaw`, `diseaseArea`, `diseased`, `isBiologicalControl`, `controlType`, `eligible`, `excludeReason`.

## Config

Labeling logic lives in `labels.py` (`coarse_disease_area`, `build_sample_label_table`). Thresholds for downstream DE (minimum cells per pseudobulk, samples per area, and so on) are set in analysis notebooks.
