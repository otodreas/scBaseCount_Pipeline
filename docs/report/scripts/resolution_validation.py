"""
This script calls the same functions as the single-SRX cluster validation
notebook (REPO_ROOT/notebooks/utility/single_srx_cluster_validation.ipynb) to
generate the figures for Supporting information S1 in the report. It loads
the five validation SRXs, runs the cluster validation, and saves the figures
to a PDF.
"""

import argparse
from pathlib import Path

import matplotlib

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

_EXPECTED_ACCESSIONS = (
    "SRX22996378",
    "SRX12366723",
    "SRX17412841",
    "SRX24313469",
    "SRX13198730",
)
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
        "-i",
        "--input",
        type=Path,
        nargs=len(_EXPECTED_ACCESSIONS),
        required=True,
        metavar="H5AD",
        help="raw h5ad paths for the five validation SRXs",
    )
    parser.add_argument("-o", "--output", type=_pdf_path, required=True, help="output PDF path")
    args = parser.parse_args()
    inputs = _index_inputs(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output) as pdf:
        for accession in _EXPECTED_ACCESSIONS:
            adata = load_srx(inputs[accession], accession)
            adata, result = run_validation(adata, accession)
            for plot_name, fig in validation_figures(adata, result):
                pdf.savefig(fig, dpi=_DPI, bbox_inches="tight")
                plt.close(fig)
                print(f"Added {accession}: {plot_name}")
            print(f"Completed {accession} (selected resolution {result.selectedResolution:g})")
    print(f"Wrote {args.output}")


def _index_inputs(paths: list[Path]) -> dict[str, Path]:
    missing_paths = [str(path) for path in paths if not path.is_file()]
    if missing_paths:
        raise FileNotFoundError(f"h5ad files not found: {', '.join(missing_paths)}")

    accessions = [path.stem for path in paths]
    duplicates = sorted({accession for accession in accessions if accessions.count(accession) > 1})
    if duplicates:
        raise ValueError(f"duplicate input accessions: {', '.join(duplicates)}")

    expected: set[str] = set(_EXPECTED_ACCESSIONS)
    found = set(accessions)
    if found != expected:
        missing = sorted(expected - found)
        unexpected = sorted(found - expected)
        raise ValueError(f"input accessions do not match validation set; missing={missing}, unexpected={unexpected}")
    return dict(zip(accessions, paths, strict=True))


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
