from __future__ import annotations

import pandas as pd
import scanpy as sc
from cyteonto.payload import build_payload


def test_build_payload_deduplicates_label_combinations() -> None:
    obs = pd.DataFrame(
        {
            "cell_type": ["A", "A", "B", "B"],
            "predicted_labels": ["x", "x", "y", "z"],
            "cytetype_annotation_leiden_merged": ["p", "p", "q", "q"],
        }
    )
    adata = sc.AnnData(obs=obs)

    payload = build_payload(
        adata,
        author_col="cell_type",
        algorithm_cols={
            "celltypist": "predicted_labels",
            "cytetype": "cytetype_annotation_leiden_merged",
        },
    )

    assert payload["authorLabels"] == ["A", "B", "B"]
    assert payload["algorithms"]["celltypist"] == ["x", "y", "z"]
    assert payload["algorithms"]["cytetype"] == ["p", "q", "q"]
    assert len(payload["authorLabels"]) == 3
