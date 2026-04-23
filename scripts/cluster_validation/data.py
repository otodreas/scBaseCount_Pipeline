from __future__ import annotations

import pandas as pd
import scanpy as sc

from cluster_validation.config import ClusterValidationConfig


def load_dataset(cfg: ClusterValidationConfig) -> tuple[sc.AnnData, str, str]:
    catalog = pd.read_csv(cfg.summaryPath).sort_values("quantile").reset_index(drop=True)
    row = _select_row(catalog, cfg)
    srx = str(row["srx_accession"])
    quantile = float(row["quantile"])
    title_suffix = f"{srx} (quantile={quantile})"
    local_path = cfg.localH5adRoot / f"{srx}.h5ad"
    if local_path.exists():
        adata = sc.read(str(local_path))
    else:
        adata = sc.read(str(row["file_path"]))
    adata.obs_names_make_unique()
    return adata, srx, title_suffix


def _select_row(df: pd.DataFrame, cfg: ClusterValidationConfig) -> pd.Series:
    if cfg.srxAccession is not None:
        hit = df.loc[df["srx_accession"] == cfg.srxAccession]
        if hit.empty:
            raise ValueError(
                f"srx_accession {cfg.srxAccession!r} not found in {cfg.summaryPath}; "
                f"choices: {df['srx_accession'].tolist()}"
            )
        return hit.iloc[0]
    idx = int(cfg.datasetIndex)
    if idx < 0 or idx >= len(df):
        raise IndexError(
            f"datasetIndex {idx} out of range for {len(df)} rows in {cfg.summaryPath}"
        )
    return df.iloc[idx]
