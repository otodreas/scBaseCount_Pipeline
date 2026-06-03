from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc
from cyteonto.dedup import dedup_table

from annotation_inspector.confidence import confidence_by_cluster

STATE_COL = "cell_type"
CYTETYPE_COL = "cytetype_annotation_leiden_merged"
CLUSTER_KEY = "leiden_merged"

PAIR_COLUMNS = [
    "accession",
    "cell_type",
    "cytetype_annotation_leiden_merged",
    "leiden_merged",
    "cytetype_confidence",
    "cytescore_similarity",
    "n_cells",
    "report_url",
]


def _report_url(adata: sc.AnnData) -> str:
    job_details = adata.uns.get("cytetype_jobDetails")
    if not isinstance(job_details, dict):
        return ""
    return str(job_details.get("report_url", ""))


def _annotate_obs(
    adata: sc.AnnData,
    srx: str,
    cyteonto_csv_path: Path | None,
) -> pd.DataFrame:
    """Build a per-cell obs frame with confidence and cytescore columns."""
    confidence_map = confidence_by_cluster(adata)
    obs = adata.obs[[STATE_COL, CYTETYPE_COL, CLUSTER_KEY]].copy()
    obs["cytetype_confidence"] = obs[CLUSTER_KEY].astype(str).map(confidence_map)
    obs["pair_label"] = obs[STATE_COL].astype(str) + obs[CYTETYPE_COL].astype(str)

    if cyteonto_csv_path is not None and cyteonto_csv_path.is_file():
        dedup_df = dedup_table(cyteonto_csv_path, srx)
        if dedup_df is not None:
            score_map = dedup_df.set_index("pair_label")["cytescore_similarity"]
            obs["cytescore_similarity"] = obs["pair_label"].map(score_map)
        else:
            obs["cytescore_similarity"] = pd.NA
    else:
        obs["cytescore_similarity"] = pd.NA

    return obs


def build_pair_df(obs: pd.DataFrame, srx: str, report_url: str) -> pd.DataFrame:
    """Aggregate a per-cell obs frame to one row per STATE/CyteType/leiden pair."""
    pair_df = (
        obs.groupby(
            [STATE_COL, CYTETYPE_COL, CLUSTER_KEY, "cytetype_confidence", "cytescore_similarity"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
    )
    pair_df.insert(0, "accession", srx)
    pair_df["report_url"] = report_url
    return pair_df[PAIR_COLUMNS]


def build_cell_df(obs: pd.DataFrame, adata: sc.AnnData, srx: str) -> pd.DataFrame:
    """Build a per-cell UMAP export frame from annotated obs and adata.obsm."""
    umap = adata.obsm["X_umap"]
    cell_df = obs.copy()
    cell_df.insert(0, "cell_id", adata.obs_names.to_numpy())
    cell_df["umap_x"] = umap[:, 0]
    cell_df["umap_y"] = umap[:, 1]
    cell_df.insert(0, "accession", srx)
    return cell_df[
        [
            "accession",
            "cell_id",
            "umap_x",
            "umap_y",
            STATE_COL,
            CYTETYPE_COL,
            CLUSTER_KEY,
            "cytetype_confidence",
            "cytescore_similarity",
        ]
    ]


def inspect_accession(
    srx: str,
    h5ad_path: Path,
    cyteonto_csv_path: Path | None,
    *,
    emit_umap: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Inspect one annotated h5ad; returns pair-level summary and optional per-cell UMAP data."""
    adata = sc.read_h5ad(h5ad_path, backed="r")
    try:
        report_url = _report_url(adata)
        obs = _annotate_obs(adata, srx, cyteonto_csv_path)
        pair_df = build_pair_df(obs, srx, report_url)
        cell_df = build_cell_df(obs, adata, srx) if emit_umap else None
        return pair_df, cell_df
    finally:
        if getattr(adata, "isbacked", False) and adata.file is not None:
            adata.file.close()
