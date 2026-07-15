"""Smoke tests for h5ad_concat aligned with run_h5ad_concat (pipeline.py).

Mocks only the R2/IO boundary (download_from_r2, load_contexts_jsonl) so
prepare_adata -> apply_qc_gate -> concat_atlas -> write_atlas run for real.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from botocore.exceptions import BotoCoreError, ClientError
from h5ad_concat import pipeline, prepare
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.models import SkipReason
from h5ad_concat.prepare import accession_from_r2_key, prepare_adata
from h5ad_concat.qc import QcStats, apply_qc_gate, flag_qc_genes
from h5ad_concat.reference import GeneReference
from study_context.models import ExperimentContext, StudyContext

ad.settings.allow_write_nullable_strings = True

_LOG = logging.getLogger("h5ad_concat_smoke")


def _cfg(**overrides) -> H5adConcatConfig:
    """Build a minimal H5adConcatConfig for tests."""
    defaults = {
        "r2Keys": ["k"],
        "datasetsPath": None,
        "minCellsAfterQc": 1,
        "minPctCellsAfterQc": 0.0,
    }
    defaults.update(overrides)
    return H5adConcatConfig(**defaults)


def _make_adata(
    accession: str = "SRX1",
    cell_types: Sequence[str | None] = ("T cell",),
    gene_names: Sequence[str] = ("GAPDH",),
    counts: np.ndarray | None = None,
    *,
    use_gene_symbols: bool = False,
    pad_to_genes: int = 0,
) -> ad.AnnData:
    """Build one in-memory AnnData for h5ad_concat tests."""
    names = list(gene_names)
    explicit_counts = counts is not None
    if counts is None:
        n_obs = len(cell_types)
    else:
        counts = np.asarray(counts, dtype=np.float32)
        n_obs = counts.shape[0]
        if len(cell_types) != n_obs:
            cell_types = ["T cell"] * n_obs
    if pad_to_genes > len(names):
        n_pad = pad_to_genes - len(names)
        names.extend(f"PAD{i}" for i in range(n_pad))
        if explicit_counts:
            counts = np.hstack([counts, np.zeros((n_obs, n_pad), dtype=np.float32)])
    if counts is None:
        counts = np.ones((n_obs, len(names)), dtype=np.float32)
    obs = pd.DataFrame(
        {"SRX_accession": [accession] * n_obs, "cell_type": list(cell_types)},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(len(names))])
    if use_gene_symbols:
        var["gene_symbols"] = names
    else:
        var.index = pd.Index(names)
    return ad.AnnData(X=counts, obs=obs, var=var)


def _contexts_for_keys(r2_keys: Sequence[str]) -> dict[str, ExperimentContext]:
    """Return study contexts keyed by accession for the given R2 keys."""
    return {
        accession_from_r2_key(key): ExperimentContext(
            accession=accession_from_r2_key(key),
            study=StudyContext(studyAccession=f"STUDY_{accession_from_r2_key(key)}"),
        )
        for key in r2_keys
    }


def _reference_for_adatas(adatas_by_key: dict[str, ad.AnnData]) -> GeneReference:
    """Build a minimal GeneReference covering all var_names in the test adatas."""
    ids = sorted({str(gene_id) for adata in adatas_by_key.values() for gene_id in adata.var_names})
    var = pd.DataFrame(
        {"gene_symbol": ids, "biotype": ["protein_coding"] * len(ids)},
        index=ids,
    )
    position = {gene_id: idx for idx, gene_id in enumerate(ids)}
    return GeneReference(ids=ids, var=var, position=position)


def _run_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
    adatas_by_key: dict[str, ad.AnnData],
    cfg: H5adConcatConfig,
):
    """Run run_h5ad_concat with download_from_r2 and load_contexts_jsonl mocked."""

    def fake_download(r2_key: str, raw_path, **_kwargs) -> None:
        adatas_by_key[r2_key].write_h5ad(raw_path)

    monkeypatch.setattr(prepare, "download_from_r2", fake_download)
    monkeypatch.setattr(pipeline, "load_contexts_jsonl", lambda _path: _contexts_for_keys(adatas_by_key))
    monkeypatch.setattr(
        pipeline,
        "load_gene_reference",
        lambda _path: _reference_for_adatas(adatas_by_key),
    )
    return pipeline.run_h5ad_concat(cfg)


def _qc_ready_adata(
    gene_names: Sequence[str],
    counts: np.ndarray,
    *,
    use_gene_symbols: bool = False,
) -> ad.AnnData:
    """Build an adata with explicit counts and enough genes for scanpy QC metrics."""
    return _make_adata(
        gene_names=gene_names,
        counts=counts,
        use_gene_symbols=use_gene_symbols,
        pad_to_genes=500,
    )


# --- run_h5ad_concat (end-to-end) ---


def test_run_h5ad_concat_happy_path(monkeypatch, tmp_path) -> None:
    key1, key2 = "prefix/SRX1.h5ad", "prefix/SRX2.h5ad"
    a1 = _make_adata("SRX1", ["T cell", "B cell"], pad_to_genes=500)
    a2 = _make_adata("SRX2", ["NK cell", "T cell", "B cell"], pad_to_genes=500)
    out = tmp_path / "atlas.h5ad"
    cfg = _cfg(
        r2Keys=[key1, key2],
        outputPath=out,
        cacheDir=tmp_path,
        minGenesPerCell=10,
        minCellsPerGene=0,
        maxPctMito=1.0,
    )

    result = _run_pipeline(monkeypatch, tmp_path, {key1: a1, key2: a2}, cfg)

    assert result.nFilesConcatenated == 2
    assert result.nObs == a1.n_obs + a2.n_obs
    assert result.skipped == []
    assert result.studiesSeen == ["STUDY_SRX1", "STUDY_SRX2"]
    assert out.exists()
    reloaded = ad.read_h5ad(out)
    assert reloaded.n_obs == result.nObs
    assert reloaded.obs_names.is_unique
    assert set(reloaded.obs_names.str.split("_").str[-1]) == {"SRX1", "SRX2"}


def test_run_h5ad_concat_skips_rejected_and_continues(monkeypatch, tmp_path) -> None:
    good_key, bad_key = "prefix/SRX1.h5ad", "prefix/SRX2.h5ad"
    good = _make_adata("SRX1", ["T cell"], pad_to_genes=500)
    bad = _make_adata("SRX2", [None, ""], pad_to_genes=500)
    out = tmp_path / "atlas.h5ad"
    cfg = _cfg(
        r2Keys=[good_key, bad_key],
        outputPath=out,
        cacheDir=tmp_path,
        minGenesPerCell=10,
        minCellsPerGene=0,
        maxPctMito=1.0,
    )

    result = _run_pipeline(monkeypatch, tmp_path, {good_key: good, bad_key: bad}, cfg)

    assert result.nFilesConcatenated == 1
    assert result.nObs == good.n_obs
    assert len(result.skipped) == 1
    assert result.skipped[0].accession == "SRX2"
    assert result.skipped[0].reason is SkipReason.cell_type_all_missing
    assert out.exists()


def test_run_h5ad_concat_raises_when_all_rejected(monkeypatch, tmp_path) -> None:
    key = "prefix/SRX1.h5ad"
    rejected = _make_adata("SRX1", [None, ""], pad_to_genes=500)
    cfg = _cfg(
        r2Keys=[key],
        outputPath=tmp_path / "atlas.h5ad",
        cacheDir=tmp_path,
        minGenesPerCell=10,
        minCellsPerGene=0,
        maxPctMito=1.0,
    )

    with pytest.raises(ValueError, match="No files passed validation"):
        _run_pipeline(monkeypatch, tmp_path, {key: rejected}, cfg)


# --- qc.py (unit) ---


def test_flag_qc_genes_classifies_and_anchors() -> None:
    genes = ["GAPDH", "MT-ND1", "mt-nd2", "RPS18", "RPL10", "HBA1", "HBB", "HBP1", "HBEGF"]
    adata = _make_adata(gene_names=genes, pad_to_genes=0)

    flag_qc_genes(adata)

    assert list(adata.var["mt"]) == [False, True, True, False, False, False, False, False, False]
    assert list(adata.var["ribo"]) == [False, False, False, True, True, False, False, False, False]
    assert list(adata.var["hb"]) == [False, False, False, False, False, True, True, False, False]


def test_flag_qc_genes_prefers_gene_symbols_column() -> None:
    genes = ["MT-ND1", "GAPDH"]
    adata = _make_adata(gene_names=genes, use_gene_symbols=True, pad_to_genes=0)

    flag_qc_genes(adata)

    assert list(adata.var_names) == ["g0", "g1"]
    assert list(adata.var["mt"]) == [True, False]


def test_apply_qc_gate_filters_and_reports() -> None:
    genes = ["G1", "G2", "G3"]
    counts = np.array([[1, 0, 0], [1, 1, 1]], dtype=np.float32)
    adata = _qc_ready_adata(genes, counts)
    cfg = _cfg(minGenesPerCell=2, minCellsPerGene=0, maxPctMito=1.0)

    filtered, stats = apply_qc_gate(adata, cfg)

    assert filtered.n_obs == 1
    assert list(filtered.obs_names) == ["cell1"]
    assert isinstance(stats, QcStats)
    assert stats.nCellsBefore == 2
    assert stats.nCellsAfter == 1
    assert stats.nGenesBefore == 500
    assert stats.nGenesAfter == 500
    assert stats.pctCellsAfter == 0.5
    assert stats.medianGenesPerCell == 3.0
    assert stats.medianPctMito == 0.0
    assert stats.medianPctRibo == 0.0
    assert stats.medianPctHb == 0.0


def test_apply_qc_gate_max_pct_hb_opt_in() -> None:
    genes = ["GAPDH", "HBA1"]
    counts = np.array([[20, 80], [50, 50]], dtype=np.float32)
    adata = _qc_ready_adata(genes, counts)
    cfg_no_hb = _cfg(minGenesPerCell=1, minCellsPerGene=0, maxPctMito=1.0, maxPctHb=1.0)

    filtered_no_hb, _ = apply_qc_gate(adata, cfg_no_hb)
    assert filtered_no_hb.n_obs == 2

    cfg_with_hb = _cfg(minGenesPerCell=1, minCellsPerGene=0, maxPctMito=1.0, maxPctHb=0.6)
    filtered_with_hb, stats = apply_qc_gate(adata, cfg_with_hb)

    assert filtered_with_hb.n_obs == 1
    assert list(filtered_with_hb.obs_names) == ["cell1"]
    assert stats.nCellsAfter == 1


def test_apply_qc_gate_rejects_when_no_cells_remain() -> None:
    genes = ["GAPDH", "MT-ND1"]
    counts = np.array([[10, 90], [10, 90]], dtype=np.float32)
    adata = _qc_ready_adata(genes, counts)
    cfg = _cfg(minGenesPerCell=1, minCellsPerGene=0, maxPctMito=0.2, minCellsAfterQc=1)

    with pytest.raises(FileRejected) as excinfo:
        apply_qc_gate(adata, cfg)

    assert excinfo.value.reason is SkipReason.too_few_cells


def _counts_by_genes_per_cell(genes_per_cell: Sequence[int], *, n_genes: int = 500) -> np.ndarray:
    """Build a count matrix where row i has genes_per_cell[i] genes detected."""
    counts = np.zeros((len(genes_per_cell), n_genes), dtype=np.float32)
    for row, n_detected in enumerate(genes_per_cell):
        counts[row, :n_detected] = 1.0
    return counts


def test_apply_qc_gate_rejects_excessive_dropout() -> None:
    genes = [f"g{i}" for i in range(500)]
    counts = _counts_by_genes_per_cell([250, 250, 10, 10, 10])
    adata = _make_adata(gene_names=genes, counts=counts, pad_to_genes=0)
    cfg = _cfg(
        minGenesPerCell=200,
        minCellsPerGene=0,
        maxPctMito=1.0,
        minCellsAfterQc=1,
        minPctCellsAfterQc=0.5,
    )

    with pytest.raises(FileRejected) as excinfo:
        apply_qc_gate(adata, cfg)

    assert excinfo.value.reason is SkipReason.excessive_cell_dropout


def test_apply_qc_gate_dropout_gate_off_by_default() -> None:
    genes = [f"g{i}" for i in range(500)]
    counts = _counts_by_genes_per_cell([250, 250, 10, 10, 10])
    adata = _make_adata(gene_names=genes, counts=counts, pad_to_genes=0)
    cfg = _cfg(
        minGenesPerCell=200,
        minCellsPerGene=0,
        maxPctMito=1.0,
        minCellsAfterQc=1,
        minPctCellsAfterQc=0.0,
    )

    filtered, stats = apply_qc_gate(adata, cfg)

    assert filtered.n_obs == 2
    assert stats.pctCellsAfter == 0.4


# --- prepare_adata download failures (unit) ---


def _client_error() -> ClientError:
    return ClientError({"Error": {"Code": "500", "Message": "boom"}}, "GetObject")


@pytest.mark.parametrize(
    ("download_exc", "expected_reason"),
    [
        (ValueError("Download MD5 mismatch for r2://bucket/key: local=a stored=b"), SkipReason.md5_mismatch),
        (ValueError("unexpected value error"), SkipReason.download_failed),
        (_client_error(), SkipReason.download_failed),
        (BotoCoreError(), SkipReason.download_failed),
        (OSError("disk full"), SkipReason.download_failed),
    ],
)
def test_prepare_adata_maps_download_errors_to_skip_reason(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
    download_exc: Exception,
    expected_reason: SkipReason,
) -> None:
    cfg = H5adConcatConfig(cacheDir=tmp_path)
    deleted: list = []
    contexts = _contexts_for_keys(["prefix/SRX1.h5ad"])
    reference = _reference_for_adatas({"prefix/SRX1.h5ad": _make_adata()})

    def _raise_download(*_args, **_kwargs) -> None:
        raise download_exc

    monkeypatch.setattr(prepare, "download_from_r2", _raise_download)
    monkeypatch.setattr(prepare, "safe_delete", lambda path, log: deleted.append(path))

    with pytest.raises(FileRejected) as excinfo:
        prepare_adata("prefix/SRX1.h5ad", "SRX1", cfg, contexts, reference, _LOG)

    assert excinfo.value.reason is expected_reason
    assert excinfo.value.__cause__ is download_exc
    assert tmp_path / "raw" / "SRX1.h5ad" in deleted
