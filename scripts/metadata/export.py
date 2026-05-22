from __future__ import annotations

from pathlib import Path

from metadata.config import MetadataConfig
from metadata.filter import FilterResult


def export_datasets(result: FilterResult, cfg: MetadataConfig) -> Path:
    """Write the full lung intersection to outputDir/datasets.csv and return its path."""
    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    datasets_path = cfg.outputDir / "datasets.csv"
    result.lungIntersection[["srx_accession", "file_path", "obs_count"]].to_csv(datasets_path, index=False)
    return datasets_path
