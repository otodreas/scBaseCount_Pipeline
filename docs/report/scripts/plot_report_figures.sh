#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"

ATLAS_QC_IN="output/atlas/2026-09-07/atlas_result.json"
ATLAS_QC_OUT="docs/report/table_1.generated.tex"

SCIB_IN="output/atlas/2026-09-07/post/param_validate/scib/scib_results.csv"
SCIB_OUT="docs/report/figs/batch_benchmark.pdf"

ATLAS_UMAP_IN="output/atlas/2026-09-07/post/production/atlas_post.h5ad"
ATLAS_UMAP_OUT="docs/report/figs/atlas_umaps.png"

RESOLUTION_IN="data/scbasecount/2026-01-12/h5ad/GeneFull/Homo_sapiens/SRX17412841.h5ad"
RESOLUTION_OUT="docs/report/figs/resolution_validation.png"

missing=0
for path in "$ATLAS_QC_IN" "$SCIB_IN" "$ATLAS_UMAP_IN" "$RESOLUTION_IN"; do
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
uv run python docs/report/scripts/resolution_validation.py -i "$RESOLUTION_IN" -o "$RESOLUTION_OUT"
