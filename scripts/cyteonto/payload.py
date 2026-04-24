from __future__ import annotations

import json
from pathlib import Path

import scanpy as sc

_AUTHOR_COL = "cell_type"
_ALGO_COL = "cytetype_annotation_leiden_merged"


def build_payload(adata: sc.AnnData) -> dict:
    return {
        "authorLabels": adata.obs[_AUTHOR_COL].tolist(),
        "algorithms": {
            "algo1": adata.obs[_ALGO_COL].tolist(),
        },
    }


def write_payload(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f)
