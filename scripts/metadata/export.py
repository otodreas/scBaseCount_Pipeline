from __future__ import annotations

from pathlib import Path

import pandas as pd

from metadata.config import MetadataConfig
from metadata.filter import FilterResult

_QUANTILES = [0.25, 0.33, 0.5, 0.67, 0.75]


def export_datasets(result: FilterResult, cfg: MetadataConfig) -> tuple[Path, Path]:
    cfg.outputDir.mkdir(parents=True, exist_ok=True)

    datasets_path = cfg.outputDir / "datasets.csv"
    df = result.lungIntersection[["srx_accession", "file_path", "obs_count"]]
    df.to_csv(datasets_path, index=False)

    quantiles_path = cfg.outputDir / "quantiles_datasets.csv"
    q_values = df["obs_count"].quantile(_QUANTILES)
    quantiles_df = (
        df[df["obs_count"].isin(q_values)]
        [["srx_accession", "file_path", "obs_count"]]
        .sort_values("obs_count")
        .copy()
    )
    quantiles_df["quantile"] = _QUANTILES
    quantiles_df.to_csv(quantiles_path, index=False)

    return datasets_path, quantiles_path
