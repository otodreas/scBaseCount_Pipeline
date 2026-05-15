from __future__ import annotations

from pathlib import Path

import scanpy as sc
from shared.repo import REPO_ROOT

from h5ad_extractor.config import H5adExtractConfig


def _resolved_repo_path(p: Path) -> Path:
    if p.is_absolute():
        return p
    return (REPO_ROOT / p).resolve()


def _h5ad_stem(p: Path) -> str:
    return p.stem


def _resolve_output_path(cfg: H5adExtractConfig, h5ad_path: Path) -> Path:
    if cfg.outputPath is not None:
        return _resolved_repo_path(cfg.outputPath)
    ext = "parquet" if cfg.outputFormat == "parquet" else "csv"
    name = f"{_h5ad_stem(h5ad_path)}_{cfg.annotationAxis}_columns.{ext}"
    out_dir = _resolved_repo_path(cfg.outputDir)
    return out_dir / name


def extract_annotation_columns(cfg: H5adExtractConfig) -> Path:
    h5ad_path = _resolved_repo_path(cfg.h5adPath)
    out_path = _resolve_output_path(cfg, h5ad_path)
    adata = sc.read(str(h5ad_path), backed="r")
    try:
        ann = adata.obs if cfg.annotationAxis == "obs" else adata.var
        available = set(ann.columns)
        missing = [c for c in cfg.columnNames if c not in available]
        if missing:
            raise ValueError(
                f"Missing {cfg.annotationAxis} column(s) {missing!r} in {h5ad_path}. Available: {sorted(available)!r}"
            )
        df = ann[cfg.columnNames].copy()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Note that when brought back into pandas, some columns data types may be coerced differently between parquet vs csv
        if cfg.outputFormat == "parquet":
            df.to_parquet(out_path, index=False)
        else:
            df.to_csv(out_path, index=False)
        return out_path
    finally:
        if getattr(adata, "isbacked", False) and adata.file is not None:
            adata.file.close()
