from __future__ import annotations

import json

import anndata as ad


def get_link(adata: ad.AnnData) -> str | None:
    """Return the CyteType report URL from adata.uns, or None when missing."""
    try:
        return adata.uns["cytetype_jobDetails"]["report_url"]
    except KeyError:
        print("No report URL found in adata.uns['cytetype_jobDetails']")
        return None


def confidence_by_cluster(adata: ad.AnnData) -> dict[str, str]:
    """Map leiden cluster ids to CyteType review confidence labels from adata.uns."""
    payload = adata.uns["cytetype_results"]["result"]
    cytetype_result = json.loads(payload) if isinstance(payload, str) else payload
    return {
        str(cluster_id): entry["latest"]["review"]["confidence"]
        for cluster_id, entry in cytetype_result["raw_annotations"].items()
    }


def add_confidence_column(adata: ad.AnnData, cluster_key: str = "leiden_merged") -> ad.AnnData | None:
    """Add CyteType confidence column to adata object."""
    try:
        confidence_map = confidence_by_cluster(adata)
        adata.obs["cytetype_confidence"] = adata.obs[cluster_key].map(lambda x: confidence_map[str(x)])
        return adata
    except KeyError:
        print("No CyteType results found in adata.uns['cytetype_results']")
        return None
