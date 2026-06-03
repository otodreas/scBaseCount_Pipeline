from __future__ import annotations

import threading
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

UMAP_COLUMNS = [
    "accession",
    "cell_id",
    "umap_x",
    "umap_y",
    "cell_type",
    "cytetype_annotation_leiden_merged",
    "leiden_merged",
    "cytetype_confidence",
    "cytescore_similarity",
]


def open_umap_writer(output_path: Path) -> pq.ParquetWriter:
    """Open a ParquetWriter for streaming per-cell UMAP rows."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    schema = pa.schema(
        [
            ("accession", pa.string()),
            ("cell_id", pa.string()),
            ("umap_x", pa.float64()),
            ("umap_y", pa.float64()),
            ("cell_type", pa.string()),
            ("cytetype_annotation_leiden_merged", pa.string()),
            ("leiden_merged", pa.string()),
            ("cytetype_confidence", pa.string()),
            ("cytescore_similarity", pa.float64()),
        ]
    )
    return pq.ParquetWriter(output_path, schema)


def append_umap_rows(
    writer: pq.ParquetWriter,
    cell_df: pd.DataFrame,
    lock: threading.Lock,
) -> None:
    """Append one accession's per-cell UMAP rows to a shared ParquetWriter."""
    table = pa.Table.from_pandas(cell_df[UMAP_COLUMNS], preserve_index=False)
    with lock:
        writer.write_table(table)
