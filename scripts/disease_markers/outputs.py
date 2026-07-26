import json
from pathlib import Path

import pandas as pd
import scanpy as sc
from shared.repo import rel_to_repo

from disease_markers.config import DiseaseMarkersConfig


def write_eligibility_labels(label_table: pd.DataFrame, output_dir: Path) -> Path:
    """Write per-SRX eligibility table to output_dir/eligibility_labels.csv."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "eligibility_labels.csv"
    label_table.to_csv(path, index=False)
    return path


def write_area_cluster_counts(adata: sc.AnnData, cfg: DiseaseMarkersConfig, output_dir: Path) -> Path:
    """Write sample and study counts per cluster and disease area."""
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped = (
        adata.obs.groupby([cfg.clusterKey, "diseaseArea"], observed=True)
        .agg(
            nCells=("diseaseArea", "size"),
            nSamples=(cfg.sampleKey, "nunique"),
            nStudies=(cfg.studyKey, "nunique"),
        )
        .reset_index()
    )
    path = output_dir / "area_cluster_counts.csv"
    grouped.to_csv(path, index=False)
    return path


def write_marker_tables(
    markers: dict[str, dict[str, pd.DataFrame]],
    output_dir: Path,
) -> list[Path]:
    """Write DE hit tables under output_dir/markers/."""
    markers_dir = output_dir / "markers"
    markers_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for cluster, area_tables in markers.items():
        for area, table in area_tables.items():
            safe_area = area.replace("/", "_")
            path = markers_dir / f"{cluster}__{safe_area}.csv"
            table.to_csv(path, index=False)
            paths.append(path)
    return paths


def write_summary_json(
    cfg: DiseaseMarkersConfig,
    label_table: pd.DataFrame,
    adata: sc.AnnData,
    marker_paths: list[Path],
    output_dir: Path,
) -> Path:
    """Write a JSON run summary with config paths and high-level counts."""
    summary = {
        "inputAtlas": rel_to_repo(cfg.inputAtlasH5ad),
        "harmonyAtlas": rel_to_repo(cfg.harmonyAtlasH5ad),
        "outputDir": rel_to_repo(output_dir),
        "nCellsAnnotated": int(adata.n_obs),
        "nGenes": int(adata.n_vars),
        "nEligibleSamples": int(label_table["eligible"].sum()),
        "nMarkerTables": len(marker_paths),
        "clusterKey": cfg.clusterKey,
        "minCellsPerProfile": cfg.minCellsPerProfile,
        "minSamplesPerArea": cfg.minSamplesPerArea,
        "minStudiesPerArea": cfg.minStudiesPerArea,
    }
    path = output_dir / "summary.json"
    path.write_text(json.dumps(summary, indent=2))
    return path
