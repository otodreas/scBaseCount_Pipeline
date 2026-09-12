"""
This script calls the same functions as the single-SRX cluster validation
notebook (REPO_ROOT/notebooks/utility/single_srx_cluster_validation.ipynb) to
generate the figures for Supporting information S1 in the report. It loads
a reproducible five-dataset sample from the report cohort that is tracked in git,
runs the cluster validation, and saves the figures to a PDF.
"""

import argparse
from pathlib import Path

import matplotlib
import pandas as pd
from dotenv import load_dotenv

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData
from cluster_validation import ClusterValidationConfig, ClusterValidationResult, run_cluster_validation_on_adata
from cluster_validation.viz import (
    plot_composition_bars,
    plot_metrics,
    plot_pca_cumvar,
    plot_resolution_sweep,
    plot_umap_merged,
)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from shared.repo import REPO_ROOT
from storage import download_from_r2, gcs_local_path, gcs_uri_to_r2_raw_key

_DEFAULT_DATASETS = REPO_ROOT / "output" / "metadata" / "datasets_v2.csv"
_DEFAULT_DATA_ROOT = REPO_ROOT / "data"
_SAMPLE_SIZE = 5
_SAMPLE_SEED = 42
_REQUIRED_COLUMNS = ("srx_accession", "file_path")
_CELL_TYPE_KEY = "cell_type"
_RESOLUTIONS = [i / 10 for i in range(1, 20)]
_DPI = 200


def _pdf_path(value: str) -> Path:
    path = Path(value)
    if path.suffix.lower() != ".pdf":
        raise argparse.ArgumentTypeError("output path must end in .pdf")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the 25-page single-SRX cluster-validation PDF.")
    parser.add_argument(
        "--datasets",
        type=Path,
        default=_DEFAULT_DATASETS,
        help="dataset catalog sampled for validation",
    )
    parser.add_argument("--data-root", type=Path, default=_DEFAULT_DATA_ROOT, help="local h5ad root")
    parser.add_argument("-o", "--output", type=_pdf_path, required=True, help="output PDF path")
    args = parser.parse_args()
    load_dotenv(REPO_ROOT / ".env")
    inputs = prepare_inputs(args.datasets, args.data_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output) as pdf:
        for accession, path in inputs:
            adata = load_srx(path, accession)
            adata, result = run_validation(adata, accession)
            for plot_name, fig in validation_figures(adata, result):
                pdf.savefig(fig, dpi=_DPI, bbox_inches="tight")
                plt.close(fig)
                print(f"Added {accession}: {plot_name}")
            print(f"Completed {accession} (selected resolution {result.selectedResolution:g})")
    print(f"Wrote {args.output}")


def sample_datasets(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"dataset catalog not found: {path}")
    datasets = pd.read_csv(path)
    missing_columns = [column for column in _REQUIRED_COLUMNS if column not in datasets.columns]
    if missing_columns:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing_columns)}")
    if len(datasets) < _SAMPLE_SIZE:
        raise ValueError(f"{path}: expected at least {_SAMPLE_SIZE} datasets, found {len(datasets)}")
    if datasets[list(_REQUIRED_COLUMNS)].isna().to_numpy().any():
        raise ValueError(f"{path}: required columns contain missing values")
    return datasets.sample(n=_SAMPLE_SIZE, random_state=_SAMPLE_SEED).reset_index(drop=True)


def prepare_inputs(datasets_path: Path, data_root: Path) -> list[tuple[str, Path]]:
    sampled = sample_datasets(datasets_path)
    accessions = sampled["srx_accession"].astype(str).tolist()
    duplicates = sorted({accession for accession in accessions if accessions.count(accession) > 1})
    if duplicates:
        raise ValueError(f"sampled duplicate accessions: {', '.join(duplicates)}")

    inputs: list[tuple[str, Path]] = []
    for (_, row), accession in zip(sampled.iterrows(), accessions, strict=True):
        gs_uri = str(row["file_path"])
        local_path = gcs_local_path(gs_uri, data_root)
        if local_path.stem != accession:
            raise ValueError(f"{datasets_path}: accession {accession!r} does not match file path {gs_uri!r}")
        if local_path.is_file():
            print(f"Using local h5ad for {accession}: {local_path}")
        else:
            r2_key = gcs_uri_to_r2_raw_key(gs_uri)
            print(f"Downloading {accession} from R2: {r2_key}")
            try:
                download_from_r2(r2_key, local_path, verify_md5=True)
            except Exception:
                local_path.unlink(missing_ok=True)
                raise
        inputs.append((accession, local_path))
    return inputs


def load_srx(path: Path, expected_accession: str) -> AnnData:
    if not path.is_file():
        raise FileNotFoundError(f"h5ad not found: {path}")
    if path.stem != expected_accession:
        raise ValueError(f"expected filename stem {expected_accession!r}, got {path.stem!r}")
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    if _CELL_TYPE_KEY not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing {_CELL_TYPE_KEY!r}")
    if "gene_symbols" not in adata.var.columns:
        raise ValueError("adata.var is missing 'gene_symbols'")
    if "SRX_accession" in adata.obs.columns:
        found = sorted(adata.obs["SRX_accession"].astype(str).unique().tolist())
        if found != [expected_accession]:
            raise ValueError(f"SRX_accession values {found} do not match {expected_accession}")
    return adata


def run_validation(
    adata: AnnData,
    accession: str,
) -> tuple[AnnData, ClusterValidationResult]:
    cfg = ClusterValidationConfig(srxAccession=accession, resolutions=_RESOLUTIONS)
    return run_cluster_validation_on_adata(
        adata,
        cfg,
        srx=accession,
        write_outputs=False,
        plot=False,
    )


def validation_figures(
    adata: AnnData,
    result: ClusterValidationResult,
) -> tuple[tuple[str, Figure], ...]:
    return (
        ("PCA cumulative variance", plot_pca_cumvar(adata, result)),
        ("resolution sweep", plot_resolution_sweep(result)),
        ("selected and RF-merged UMAPs", plot_umap_merged(adata, result)),
        ("cluster composition", plot_composition_bars(adata, result)),
        ("cluster similarity metrics", plot_metrics(result)),
    )


if __name__ == "__main__":
    main()
