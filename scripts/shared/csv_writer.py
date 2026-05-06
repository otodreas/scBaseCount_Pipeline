from __future__ import annotations

import csv
from pathlib import Path


def append_csv_row(summary_path: Path, columns: list[str], values: list) -> None:
    write_header = not summary_path.exists()
    with open(summary_path, "a", newline="") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(columns)
        writer.writerow(values)
