"""Smoke tests for h5ad_concat that run entirely on in-memory AnnData (no R2 I/O).

Covers concat_atlas / write_atlas and the adata-level prepare helpers
(cell_type_all_missing, fill_cell_type, validate_single_accession).
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import SkipReason
from h5ad_concat.prepare import cell_type_all_missing, fill_cell_type, validate_single_accession

_LOG = logging.getLogger("h5ad_concat_smoke")


def _make_adata(accession: str, cell_types: list[str | None], var_names: list[str]) -> ad.AnnData:
    """Build a tiny AnnData with SRX_accession and cell_type obs columns and shared barcodes."""
    n_obs = len(cell_types)
    x = np.arange(n_obs * len(var_names), dtype=np.float32).reshape(n_obs, len(var_names))
    obs = pd.DataFrame(
        {"SRX_accession": [accession] * n_obs, "cell_type": cell_types},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=list(var_names))
    return ad.AnnData(X=x, obs=obs, var=var)


def test_concat_atlas_stacks_and_suffixes_barcodes() -> None:
    cfg = H5adConcatConfig()
    genes = ["g1", "g2", "g3"]
    a1 = _make_adata("SRX1", ["T cell", "B cell"], genes)
    a2 = _make_adata("SRX2", ["NK cell", "T cell", "B cell"], genes)

    atlas = concat_atlas([a1, a2], ["SRX1", "SRX2"], cfg, _LOG)

    assert atlas.n_obs == a1.n_obs + a2.n_obs
    assert atlas.obs_names.is_unique
    assert list(atlas.obs_names) == [
        "cell0_SRX1",
        "cell1_SRX1",
        "cell0_SRX2",
        "cell1_SRX2",
        "cell2_SRX2",
    ]
    assert set(atlas.obs["SRX_accession"]) == {"SRX1", "SRX2"}


def test_concat_atlas_inner_join_intersects_vars() -> None:
    cfg = H5adConcatConfig(join="inner")
    a1 = _make_adata("SRX1", ["T cell"], ["g1", "g2", "g3"])
    a2 = _make_adata("SRX2", ["B cell"], ["g2", "g3", "g4"])

    atlas = concat_atlas([a1, a2], ["SRX1", "SRX2"], cfg, _LOG)

    assert list(atlas.var_names) == ["g2", "g3"]


def test_concat_atlas_outer_join_unions_vars() -> None:
    cfg = H5adConcatConfig(join="outer")
    a1 = _make_adata("SRX1", ["T cell"], ["g1", "g2"])
    a2 = _make_adata("SRX2", ["B cell"], ["g2", "g3"])

    atlas = concat_atlas([a1, a2], ["SRX1", "SRX2"], cfg, _LOG)

    assert set(atlas.var_names) == {"g1", "g2", "g3"}


def test_write_atlas_roundtrip(tmp_path) -> None:
    out = tmp_path / "atlas.h5ad"
    cfg = H5adConcatConfig(outputPath=out)
    a1 = _make_adata("SRX1", ["T cell", "B cell"], ["g1", "g2"])
    a2 = _make_adata("SRX2", ["NK cell"], ["g1", "g2"])
    atlas = concat_atlas([a1, a2], ["SRX1", "SRX2"], cfg, _LOG)

    path = write_atlas(atlas, cfg, _LOG)

    assert path == out
    assert out.exists()
    reloaded = ad.read_h5ad(out)
    assert reloaded.n_obs == atlas.n_obs
    assert list(reloaded.var_names) == list(atlas.var_names)
    assert reloaded.obs_names.is_unique


def test_cell_type_all_missing_when_column_absent() -> None:
    adata = _make_adata("SRX1", ["T cell"], ["g1"])
    del adata.obs["cell_type"]
    assert cell_type_all_missing(adata, "cell_type") is True


def test_cell_type_all_missing_when_all_blank_or_nan() -> None:
    adata = _make_adata("SRX1", [None, "", "  "], ["g1"])
    assert cell_type_all_missing(adata, "cell_type") is True


def test_cell_type_all_missing_false_when_some_present() -> None:
    adata = _make_adata("SRX1", [None, "T cell"], ["g1"])
    assert cell_type_all_missing(adata, "cell_type") is False


def test_fill_cell_type_replaces_blanks_and_keeps_labels() -> None:
    cfg = H5adConcatConfig()
    adata = _make_adata("SRX1", ["T cell", None, ""], ["g1"])

    fill_cell_type(adata, cfg)

    assert list(adata.obs["cell_type"]) == ["T cell", cfg.missingLabel, cfg.missingLabel]


def test_validate_single_accession_passes_on_match() -> None:
    cfg = H5adConcatConfig()
    adata = _make_adata("SRX1", ["T cell", "B cell"], ["g1"])
    validate_single_accession(adata, "SRX1", cfg)


def test_validate_single_accession_rejects_multiple_values() -> None:
    cfg = H5adConcatConfig()
    adata = _make_adata("SRX1", ["T cell", "B cell"], ["g1"])
    adata.obs["SRX_accession"] = ["SRX1", "SRX2"]

    with pytest.raises(FileRejected) as excinfo:
        validate_single_accession(adata, "SRX1", cfg)

    assert excinfo.value.reason is SkipReason.accession_mismatch
    assert isinstance(excinfo.value.__cause__, ValueError)


def test_validate_single_accession_rejects_wrong_value() -> None:
    cfg = H5adConcatConfig()
    adata = _make_adata("SRX9", ["T cell"], ["g1"])

    with pytest.raises(FileRejected) as excinfo:
        validate_single_accession(adata, "SRX1", cfg)

    assert excinfo.value.reason is SkipReason.accession_mismatch
