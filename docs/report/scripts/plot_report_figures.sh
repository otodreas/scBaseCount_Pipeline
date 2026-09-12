#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"

# Committed in git
ATLAS_QC_IN="output/atlas/2026-09-07/atlas_result.json"
ATLAS_QC_OUT="docs/report/tables/table_1.tex"

# Committed in git
SCIB_IN="output/atlas/2026-09-07/post/param_validate/scib/scib_results.csv"
SCIB_OUT="docs/report/figs/batch_benchmark.pdf"

# NOT committed in git -- ensure ATLAS_UMAP_IN exists
ATLAS_UMAP_IN="output/atlas/2026-09-07/post/production/atlas_post.h5ad"
ATLAS_UMAP_OUT="docs/report/figs/atlas_umaps.png"

RESOLUTION_OUT="docs/report/si/resolution_validation.pdf"

missing=0
for path in "$ATLAS_QC_IN" "$SCIB_IN" "$ATLAS_UMAP_IN"; do
    if [[ ! -f "$path" ]]; then
        echo "missing input: $path" >&2
        missing=1
    fi
done
if [[ "$missing" -ne 0 ]]; then
    exit 1
fi

uv run python docs/report/scripts/atlas_qc_table.py -i "$ATLAS_QC_IN" -o "$ATLAS_QC_OUT"
uv run python docs/report/scripts/batch_benchmark.py -i "$SCIB_IN" -o "$SCIB_OUT"
uv run python docs/report/scripts/atlas_umaps.py -i "$ATLAS_UMAP_IN" -o "$ATLAS_UMAP_OUT"
uv run python docs/report/scripts/resolution_validation.py \
    -o "$RESOLUTION_OUT"
