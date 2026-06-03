from __future__ import annotations

import json

import scanpy as sc


def confidence_by_cluster(adata: sc.AnnData) -> dict[str, str]:
    """Map leiden cluster ids to CyteType review confidence labels from adata.uns."""
    payload = adata.uns["cytetype_results"]["result"]
    cytetype_result = json.loads(payload) if isinstance(payload, str) else payload
    return {
        str(cluster_id): entry["latest"]["review"]["confidence"]
        for cluster_id, entry in cytetype_result["raw_annotations"].items()
    }
