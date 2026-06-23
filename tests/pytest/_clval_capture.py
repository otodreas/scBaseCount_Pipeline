"""Shared capture helpers for cluster_validation golden regression tests."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from cluster_validation.config import ClusterValidationConfig
from cluster_validation.models import ClusterValidationResult
from shared.repo import rel_to_repo

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


def config_snapshot(cfg: ClusterValidationConfig) -> dict[str, Any]:
    """Return the config's model_dump with Path fields normalized to repo-relative strings."""
    dump = cfg.model_dump(mode="python")
    # Make sure all Path fields are relative to the repo root, not absolute
    return {k: (rel_to_repo(v) if isinstance(v, Path) else v) for k, v in dump.items()}
