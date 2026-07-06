from __future__ import annotations

from pathlib import Path

import orjson
import scanpy as sc

_DEFAULT_AUTHOR_COL = "cell_type"
_DEFAULT_ALGORITHM_COLS = {"algo1": "cytetype_annotation_leiden_merged"}


def build_payload(
    adata: sc.AnnData,
    author_col: str = _DEFAULT_AUTHOR_COL,
    algorithm_cols: dict[str, str] | None = None,
) -> dict:
    """Build a deduplicated CyteOnto payload from unique label combinations in obs."""
    algo_cols = algorithm_cols if algorithm_cols is not None else _DEFAULT_ALGORITHM_COLS
    obs_cols = [author_col, *algo_cols.values()]
    missing = [col for col in obs_cols if col not in adata.obs]
    if missing:
        raise ValueError(f"column(s) not found in adata.obs: {missing}")

    unique = adata.obs[obs_cols].astype(str).drop_duplicates()
    return {
        "authorLabels": unique[author_col].tolist(),
        "algorithms": {name: unique[col].tolist() for name, col in algo_cols.items()},
    }


def write_payload(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(orjson.dumps(payload))
