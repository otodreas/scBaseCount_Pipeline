from pathlib import Path
from unittest.mock import MagicMock, patch

import scanpy as sc
from atlas_postprocessing.scib import run_scib_benchmark


def test_run_scib_benchmark_writes_artifacts(tmp_path: Path) -> None:
    adata = sc.AnnData()
    out_dir = tmp_path / "scib"

    fake_df = MagicMock()
    bm = MagicMock()
    bm.get_results.return_value = fake_df

    with patch("atlas_postprocessing.scib.Benchmarker", return_value=bm) as benchmarker:
        csv_path = run_scib_benchmark(
            adata,
            outDir=out_dir,
            batchKey="study_accession",
            labelKey="cell_type",
            force=True,
        )

    benchmarker.assert_called_once()
    kwargs = benchmarker.call_args.kwargs
    assert kwargs["batch_key"] == "study_accession"
    assert kwargs["label_key"] == "cell_type"
    assert kwargs["embedding_obsm_keys"] == ["X_pca", "X_pca_harmony"]
    assert kwargs["pre_integrated_embedding_obsm_key"] == "X_pca"
    bm.benchmark.assert_called_once()
    fake_df.to_csv.assert_called_once_with(csv_path)
    bm.plot_results_table.assert_called_once_with(show=False, save_dir=str(out_dir))


def test_run_scib_benchmark_skips_when_artifacts_exist(tmp_path: Path) -> None:
    out_dir = tmp_path / "scib"
    out_dir.mkdir()
    (out_dir / "scib_results.csv").write_text("ok")
    (out_dir / "scib_results.svg").write_text("<svg/>")

    with patch("atlas_postprocessing.scib.Benchmarker") as benchmarker:
        path = run_scib_benchmark(sc.AnnData(), outDir=out_dir, force=False)

    assert path == out_dir / "scib_results.csv"
    benchmarker.assert_not_called()
