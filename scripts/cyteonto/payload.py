from __future__ import annotations

from pathlib import Path

import orjson
import scanpy as sc

_AUTHOR_COL = "cell_type"
_ALGO_COL = "cytetype_annotation_leiden_merged"


def build_payload(adata: sc.AnnData) -> dict:
    return {
        "authorLabels": adata.obs[_AUTHOR_COL].to_list(),
        "algorithms": {
            "algo1": adata.obs[_ALGO_COL].to_list(),
        },
    }


def write_payload(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(orjson.dumps(payload))  # , option=orjson.OPT_SERIALIZE_NUMPY))
