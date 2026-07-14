"""Smoke tests for h5ad_concat that run entirely on in-memory AnnData (no real R2 I/O).

Covers concat_atlas / write_atlas, the adata-level prepare helpers
(cell_type_all_missing, fill_cell_type, validate_single_accession), and the
prepare_adata download-failure handling (with download_from_r2 mocked out).
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from botocore.exceptions import BotoCoreError, ClientError
from h5ad_concat import prepare
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import SkipReason
from h5ad_concat.prepare import (
    cell_type_all_missing,
    fill_cell_type,
    prepare_adata,
    validate_single_accession,
)

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


def _client_error() -> ClientError:
    """Build a minimal botocore ClientError for a failed GetObject call."""
    return ClientError({"Error": {"Code": "500", "Message": "boom"}}, "GetObject")


def _run_prepare_with_failing_download(
    monkeypatch: pytest.MonkeyPatch, tmp_path, download_exc: Exception
) -> tuple[FileRejected, list]:
    """Run prepare_adata with download_from_r2 raising download_exc; return (raised, deleted paths)."""
    cfg = H5adConcatConfig(cacheDir=tmp_path)
    deleted: list = []

    def _raise_download(*_args, **_kwargs) -> None:
        raise download_exc

    monkeypatch.setattr(prepare, "resolve_batch_key", lambda accession, contexts: "STUDY1")
    monkeypatch.setattr(prepare, "download_from_r2", _raise_download)
    monkeypatch.setattr(prepare, "safe_delete", lambda path, log: deleted.append(path))

    with pytest.raises(FileRejected) as excinfo:
        prepare_adata("prefix/SRX1.h5ad", "SRX1", cfg, {}, _LOG)

    return excinfo.value, deleted


def test_prepare_adata_md5_mismatch_rejected(monkeypatch, tmp_path) -> None:
    exc = ValueError("Download MD5 mismatch for r2://bucket/key: local=a stored=b")
    rejected, deleted = _run_prepare_with_failing_download(monkeypatch, tmp_path, exc)

    assert rejected.reason is SkipReason.md5_mismatch
    assert rejected.__cause__ is exc
    raw_path = tmp_path / "raw" / "SRX1.h5ad"
    assert raw_path in deleted


def test_prepare_adata_other_value_error_rejected_as_download_failed(monkeypatch, tmp_path) -> None:
    exc = ValueError("unexpected value error")
    rejected, deleted = _run_prepare_with_failing_download(monkeypatch, tmp_path, exc)

    assert rejected.reason is SkipReason.download_failed
    assert rejected.__cause__ is exc
    assert (tmp_path / "raw" / "SRX1.h5ad") in deleted


def test_prepare_adata_client_error_rejected_as_download_failed(monkeypatch, tmp_path) -> None:
    exc = _client_error()
    rejected, deleted = _run_prepare_with_failing_download(monkeypatch, tmp_path, exc)

    assert rejected.reason is SkipReason.download_failed
    assert rejected.__cause__ is exc
    assert (tmp_path / "raw" / "SRX1.h5ad") in deleted


def test_prepare_adata_botocore_error_rejected_as_download_failed(monkeypatch, tmp_path) -> None:
    exc = BotoCoreError()
    rejected, deleted = _run_prepare_with_failing_download(monkeypatch, tmp_path, exc)

    assert rejected.reason is SkipReason.download_failed
    assert rejected.__cause__ is exc
    assert (tmp_path / "raw" / "SRX1.h5ad") in deleted


def test_prepare_adata_os_error_rejected_as_download_failed(monkeypatch, tmp_path) -> None:
    exc = OSError("disk full")
    rejected, deleted = _run_prepare_with_failing_download(monkeypatch, tmp_path, exc)

    assert rejected.reason is SkipReason.download_failed
    assert rejected.__cause__ is exc
    assert (tmp_path / "raw" / "SRX1.h5ad") in deleted
