import argparse
import importlib.util
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_MIGRATE_PATH = _REPO_ROOT / "pipelines" / "migrate_gcs_to_r2.py"

dotenv_module = MagicMock()
dotenv_module.load_dotenv = MagicMock()
sys.modules["dotenv"] = dotenv_module

_SPEC = importlib.util.spec_from_file_location("migrate_gcs_to_r2", _MIGRATE_PATH)
assert _SPEC and _SPEC.loader
migrate = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(migrate)


def _source_csv(path: Path, rows: list[dict[str, str]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def _baseline_csv(path: Path, accessions: list[str]) -> None:
    pd.DataFrame(
        {
            "srx_accession": accessions,
            "file_path": [f"gs://bucket/old/{accession}.h5ad" for accession in accessions],
        }
    ).to_csv(path, index=False)


def test_select_migration_rows_preserves_order_and_excludes_baseline() -> None:
    source = pd.DataFrame(
        {
            "srx_accession": ["SRX1", "SRX2", "SRX3", "SRX4"],
            "file_path": [
                "gs://bucket/a/SRX1.h5ad",
                "gs://bucket/a/SRX2.h5ad",
                "gs://bucket/a/SRX3.h5ad",
                "gs://bucket/a/SRX4.h5ad",
            ],
        }
    )
    baseline = pd.DataFrame({"srx_accession": ["SRX2", "SRX4"], "file_path": ["x", "y"]})

    selected = migrate.select_migration_rows(source, baseline)

    assert selected["srx_accession"].tolist() == ["SRX1", "SRX3"]


def test_validate_datasets_csv_rejects_missing_columns(tmp_path: Path) -> None:
    csv_path = tmp_path / "datasets.csv"
    pd.DataFrame({"srx_accession": ["SRX1"]}).to_csv(csv_path, index=False)

    with pytest.raises(ValueError, match="missing required columns"):
        migrate.validate_datasets_csv(pd.read_csv(csv_path), csv_path, require_unique_accessions=True)


def test_validate_datasets_csv_rejects_blank_values(tmp_path: Path) -> None:
    csv_path = tmp_path / "datasets.csv"
    pd.DataFrame({"srx_accession": ["SRX1", ""], "file_path": ["gs://bucket/a.h5ad", "gs://bucket/b.h5ad"]}).to_csv(
        csv_path, index=False
    )

    with pytest.raises(ValueError, match="blank srx_accession"):
        migrate.validate_datasets_csv(pd.read_csv(csv_path), csv_path, require_unique_accessions=True)


def test_validate_datasets_csv_rejects_duplicate_accessions(tmp_path: Path) -> None:
    csv_path = tmp_path / "datasets.csv"
    pd.DataFrame(
        {
            "srx_accession": ["SRX1", "SRX1"],
            "file_path": ["gs://bucket/a.h5ad", "gs://bucket/b.h5ad"],
        }
    ).to_csv(csv_path, index=False)

    with pytest.raises(ValueError, match="duplicate accession"):
        migrate.validate_datasets_csv(pd.read_csv(csv_path), csv_path, require_unique_accessions=True)


@pytest.mark.parametrize("baseline_args", [[], ["--baseline"]], ids=["omitted", "blank"])
def test_optional_baseline_excludes_nothing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, baseline_args: list[str]
) -> None:
    source_path = tmp_path / "source.csv"
    _source_csv(source_path, [{"srx_accession": "SRX1", "file_path": "gs://bucket/a/SRX1.h5ad"}])

    monkeypatch.setattr(
        sys,
        "argv",
        ["migrate_gcs_to_r2.py", "--datasets", str(source_path), *baseline_args, "--dry-run"],
    )
    monkeypatch.setattr(migrate, "gcs_uri_to_r2_raw_key", lambda _gs_uri: "raw/key")
    monkeypatch.setattr(migrate, "gcs_blob_md5", lambda _gs_uri: "md5-ok")
    monkeypatch.setattr(migrate, "r2_raw_matches_gcs", lambda _key, _md5: False)
    monkeypatch.setattr(migrate, "RUN_OUTPUT_DIR", tmp_path / "migration")

    args = migrate._parse_args()
    exit_code = migrate.run_migration(args, run_timestamp="test_run")

    summary = pd.read_csv(tmp_path / "migration" / "test_run" / "run.csv")
    assert args.baseline is None
    assert summary["status"].tolist() == ["dry_run"]
    assert exit_code == 0


def test_run_migration_accepts_blank_baseline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    source_path = tmp_path / "source.csv"
    baseline_path = tmp_path / "baseline.csv"
    _source_csv(source_path, [{"srx_accession": "SRX1", "file_path": "gs://bucket/a/SRX1.h5ad"}])
    baseline_path.touch()

    monkeypatch.setattr(migrate, "gcs_uri_to_r2_raw_key", lambda _gs_uri: "raw/key")
    monkeypatch.setattr(migrate, "gcs_blob_md5", lambda _gs_uri: "md5-ok")
    monkeypatch.setattr(migrate, "r2_raw_matches_gcs", lambda _key, _md5: False)
    monkeypatch.setattr(migrate, "RUN_OUTPUT_DIR", tmp_path / "migration")

    args = argparse.Namespace(datasets=source_path, baseline=baseline_path, dry_run=True)
    exit_code = migrate.run_migration(args, run_timestamp="test_run")

    summary = pd.read_csv(tmp_path / "migration" / "test_run" / "run.csv")
    assert summary["status"].tolist() == ["dry_run"]
    assert exit_code == 0


def test_run_migration_continues_after_precheck_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    source_path = tmp_path / "source.csv"
    baseline_path = tmp_path / "baseline.csv"
    _source_csv(
        source_path,
        [
            {"srx_accession": "SRX1", "file_path": "gs://bucket/a/SRX1.h5ad"},
            {"srx_accession": "SRX2", "file_path": "gs://bucket/a/SRX2.h5ad"},
        ],
    )
    _baseline_csv(baseline_path, [])

    calls: list[str] = []

    def fake_gcs_blob_md5(gs_uri: str) -> str:
        calls.append(gs_uri)
        if gs_uri.endswith("SRX1.h5ad"):
            raise RuntimeError("GCS unavailable")
        return "md5-ok"

    monkeypatch.setattr(migrate, "gcs_uri_to_r2_raw_key", lambda gs_uri: gs_uri.replace("gs://bucket/", "raw/"))
    monkeypatch.setattr(migrate, "gcs_blob_md5", fake_gcs_blob_md5)
    monkeypatch.setattr(migrate, "r2_raw_matches_gcs", lambda _key, _md5: False)
    monkeypatch.setattr(migrate, "RUN_OUTPUT_DIR", tmp_path / "migration")

    args = argparse.Namespace(datasets=source_path, baseline=baseline_path, dry_run=True)
    exit_code = migrate.run_migration(args, run_timestamp="test_run")

    summary = pd.read_csv(tmp_path / "migration" / "test_run" / "run.csv")
    assert calls == ["gs://bucket/a/SRX1.h5ad", "gs://bucket/a/SRX2.h5ad"]
    assert summary["status"].tolist() == ["failed", "dry_run"]
    assert exit_code == 1


def test_run_migration_returns_zero_on_success(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    source_path = tmp_path / "source.csv"
    baseline_path = tmp_path / "baseline.csv"
    _source_csv(source_path, [{"srx_accession": "SRX1", "file_path": "gs://bucket/a/SRX1.h5ad"}])
    _baseline_csv(baseline_path, [])

    monkeypatch.setattr(migrate, "gcs_uri_to_r2_raw_key", lambda gs_uri: "raw/key")
    monkeypatch.setattr(migrate, "gcs_blob_md5", lambda _gs_uri: "md5-ok")
    monkeypatch.setattr(migrate, "r2_raw_matches_gcs", lambda _key, _md5: True)
    monkeypatch.setattr(migrate, "RUN_OUTPUT_DIR", tmp_path / "migration")

    args = argparse.Namespace(datasets=source_path, baseline=baseline_path, dry_run=False)
    exit_code = migrate.run_migration(args, run_timestamp="test_run")

    summary = pd.read_csv(tmp_path / "migration" / "test_run" / "run.csv")
    assert summary["status"].tolist() == ["skipped"]
    assert exit_code == 0


def test_run_migration_delta_against_repo_metadata() -> None:
    source = pd.read_csv(_REPO_ROOT / "output/metadata/datasets_v2.csv")
    baseline = pd.read_csv(_REPO_ROOT / "output/metadata/datasets.csv")

    selected = migrate.select_migration_rows(source, baseline)

    assert len(selected) == 1048
    assert selected["srx_accession"].is_unique
