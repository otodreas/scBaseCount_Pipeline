"""Smoke test that every package shipped by the wheel build target imports cleanly.

Catches missing __init__.py re-exports, broken cross-package imports, accidentally
renamed symbols, and the wheel-target list in pyproject.toml drifting from reality.
"""

from __future__ import annotations

import importlib

import pytest

PACKAGES = [
    "shared",
    "study_context",
    "cluster_validation",
    "metadata",
    "cyteonto",
    "h5ad_extractor",
    "storage",
    "cytetype_runner",
    "annotation_inspector",
]


@pytest.mark.parametrize("pkg", PACKAGES)
def test_package_imports(pkg: str) -> None:
    importlib.import_module(pkg)
