from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from h5ad_concat.checkpoint import MANIFEST_KEY, load_checkpoint, write_checkpoint
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import fold_atlas
from h5ad_concat.models import ConcatManifest, ManifestEntry, SkipReason
from h5ad_concat.pipeline import run_h5ad_concat

_GENES = ["G1", "G2", "G3"]


def _make_adata(n_obs: int, study: str, cell_type: str = "T cell") -> ad.AnnData:
    x = sp.csr_matrix(np.ones((n_obs, len(_GENES)), dtype=np.float32))
    return ad.AnnData(
        X=x,
        obs={
            "cell_type": [cell_type] * n_obs,
            "study_accession": [study] * n_obs,
        },
        var={"gene_ids": _GENES},
    )


def _manifest_with(entries: list[ManifestEntry]) -> ConcatManifest:
    return ConcatManifest(join="inner", batchKey="study_accession", entries=entries)


def test_write_and_load_checkpoint_roundtrip(tmp_path: Path) -> None:
    atlas = _make_adata(3, "PRJONE")
    manifest = _manifest_with(
        [
            ManifestEntry(
                r2Key="key/a.h5ad",
                accession="a",
                concatenated=True,
                study="PRJONE",
            )
        ]
    )
    cfg = H5adConcatConfig(
        r2Keys=["key/a.h5ad"],
        outputPath=tmp_path / "atlas.h5ad",
        resume=True,
    )
    log = logging.getLogger("test_h5ad_concat_checkpoint")

    write_checkpoint(atlas, manifest, cfg, log)
    loaded_atlas, loaded_manifest = load_checkpoint(cfg, log)

    assert loaded_atlas is not None
    assert loaded_atlas.n_obs == 3
    assert loaded_manifest.processedKeys() == {"key/a.h5ad"}
    assert loaded_atlas.uns[MANIFEST_KEY] == manifest.model_dump_json()
    assert not (tmp_path / "atlas.tmp.h5ad").exists()


def test_load_checkpoint_rejects_join_mismatch(tmp_path: Path) -> None:
    atlas = _make_adata(2, "PRJONE")
    manifest = _manifest_with([])
    cfg = H5adConcatConfig(
        r2Keys=["key/a.h5ad"],
        outputPath=tmp_path / "atlas.h5ad",
        join="inner",
    )
    write_checkpoint(atlas, manifest, cfg, logging.getLogger("test_h5ad_concat_checkpoint"))

    cfg_outer = H5adConcatConfig(
        r2Keys=["key/a.h5ad"],
        outputPath=tmp_path / "atlas.h5ad",
        join="outer",
        resume=True,
    )
    with pytest.raises(ValueError, match="Cannot resume"):
        load_checkpoint(cfg_outer, logging.getLogger("test_h5ad_concat_checkpoint"))


def test_fold_atlas_matches_single_concat(tmp_path: Path) -> None:
    cfg = H5adConcatConfig(r2Keys=["k1", "k2", "k3"])
    adatas = [_make_adata(2, "PRJONE"), _make_adata(3, "PRJTWO"), _make_adata(1, "PRJTHREE")]

    folded = fold_atlas(None, [adatas[0]], cfg)
    folded = fold_atlas(folded, [adatas[1]], cfg)
    folded = fold_atlas(folded, [adatas[2]], cfg)

    direct = ad.concat(adatas, axis="obs", join="inner")
    assert folded.n_obs == direct.n_obs
    assert folded.n_vars == direct.n_vars


def test_run_h5ad_concat_resumes_after_simulated_crash(tmp_path: Path) -> None:
    keys = ["key/a.h5ad", "key/b.h5ad", "key/c.h5ad"]
    cfg = H5adConcatConfig(
        r2Keys=keys,
        outputPath=tmp_path / "atlas.h5ad",
        cacheDir=tmp_path / "cache",
        checkpointEvery=2,
        resume=True,
        verifyMd5=False,
    )
    prepared = {
        "key/a.h5ad": (_make_adata(2, "PRJONE"), "PRJONE"),
        "key/b.h5ad": (_make_adata(3, "PRJTWO"), "PRJTWO"),
        "key/c.h5ad": (_make_adata(1, "PRJTHREE"), "PRJTHREE"),
    }
    call_log: list[str] = []

    def fake_prepare(r2_key: str, *_args, **_kwargs):
        call_log.append(r2_key)
        return prepared[r2_key]

    with patch("h5ad_concat.pipeline.prepare_adata", side_effect=fake_prepare):
        with patch("h5ad_concat.pipeline.load_contexts_jsonl", return_value={}):
            first_result = run_h5ad_concat(cfg)

    assert first_result.nObs == 6
    assert first_result.nFilesConcatenated == 3
    assert call_log == keys

    call_log.clear()
    with patch("h5ad_concat.pipeline.prepare_adata", side_effect=fake_prepare):
        with patch("h5ad_concat.pipeline.load_contexts_jsonl", return_value={}):
            resumed_result = run_h5ad_concat(cfg)

    assert resumed_result.nObs == 6
    assert resumed_result.nFilesConcatenated == 3
    assert call_log == []


def test_run_h5ad_concat_records_skipped_in_manifest(tmp_path: Path) -> None:
    keys = ["key/good.h5ad", "key/bad.h5ad"]
    cfg = H5adConcatConfig(
        r2Keys=keys,
        outputPath=tmp_path / "atlas.h5ad",
        cacheDir=tmp_path / "cache",
        checkpointEvery=1,
        resume=True,
        verifyMd5=False,
    )

    def fake_prepare(r2_key: str, *_args, **_kwargs):
        if r2_key.endswith("bad.h5ad"):
            raise FileRejected(SkipReason.cell_type_all_missing)
        return _make_adata(2, "PRJONE"), "PRJONE"

    with patch("h5ad_concat.pipeline.prepare_adata", side_effect=fake_prepare):
        with patch("h5ad_concat.pipeline.load_contexts_jsonl", return_value={}):
            result = run_h5ad_concat(cfg)

    assert result.nFilesConcatenated == 1
    assert len(result.skipped) == 1
    assert result.skipped[0].reason == SkipReason.cell_type_all_missing

    _, manifest = load_checkpoint(cfg, logging.getLogger("test_h5ad_concat_checkpoint"))
    assert manifest.processedKeys() == set(keys)
