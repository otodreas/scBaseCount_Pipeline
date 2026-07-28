# disease_markers

Per-SRX coarse disease area labels and atlas eligibility rules derived from `contexts.jsonl` and `atlas.csv`.

## Usage

```python
from pathlib import Path

from disease_markers.labels import build_sample_label_table, coarse_disease_area

label_table = build_sample_label_table(
    Path("output/context/contexts.jsonl"),
    Path("output/atlas/v1/atlas.csv"),
)
```

Join labels onto an AnnData object by sample accession (for example `SRX_accession` in `obs`).

## Outputs

Typical notebook or ad hoc export: `eligibility_labels.csv` with columns `srxAccession`, `studyAccession`, `diseaseRaw`, `diseaseArea`, `isControl`, `eligible`, `excludeReason`.

## Config

Labeling logic lives in `labels.py` (`coarse_disease_area`, `build_sample_label_table`). Thresholds for downstream DE (minimum cells per pseudobulk, samples per area, and so on) are set in analysis notebooks.
