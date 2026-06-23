"""Shared capture helpers for cluster_validation golden regression tests."""

from __future__ import annotations

from typing import Any

from cluster_validation.models import ClusterValidationResult

_NUMERIC_SCALARS = ("selectedResolution", "cumvar")
_INT_SCALARS = (
    "nPcs",
    "kPrior",
    "kFiltered",
    "nCellsDropped",
    "nCellsFinal",
    "nClustersPreMerge",
    "nClustersPostMerge",
)
_STR_SCALARS = ("clusterKey", "mergedKey")
_LIST_FIELDS = (
    "resolutions",
    "kArr",
    "jaccArr",
    "silhouetteArr",
    "homogeneityArr",
    "completenessArr",
    "nmiArr",
    "vscoreArr",
    "ariArr",
    "confMatrix",
    "confClasses",
)
_MAP_FIELDS = ("labelMap", "mergedGroups")

ALL_FIELDS = _NUMERIC_SCALARS + _INT_SCALARS + _STR_SCALARS + _LIST_FIELDS + _MAP_FIELDS


def capture_fields(result: ClusterValidationResult) -> dict[str, Any]:
    """Return deterministic ClusterValidationResult fields for baseline comparison."""
    dump = result.model_dump(mode="json")
    return {field: dump.get(field) for field in ALL_FIELDS}
