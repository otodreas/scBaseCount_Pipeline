import argparse
import datetime
from pathlib import Path

from disease_markers.config import DiseaseMarkersConfig
from disease_markers.pipeline import log, run_disease_markers
from shared.logger import log_run_separator

_DEFAULT_CFG = DiseaseMarkersConfig()


def _parse_args() -> argparse.Namespace:
    d = _DEFAULT_CFG
    parser = argparse.ArgumentParser(
        description="Transfer Harmony clusters, annotate disease areas, and run pseudobulk one-vs-rest DE.",
    )
    parser.add_argument(
        "--input-atlas", type=Path, default=d.inputAtlasH5ad, metavar="PATH", help="Full-gene atlas h5ad"
    )
    parser.add_argument(
        "--harmony-atlas",
        type=Path,
        default=d.harmonyAtlasH5ad,
        metavar="PATH",
        help="Harmony atlas h5ad with cluster labels",
    )
    parser.add_argument("--contexts", type=Path, default=d.contextsPath, metavar="PATH", help="contexts.jsonl path")
    parser.add_argument("--atlas-csv", type=Path, default=d.atlasCsvPath, metavar="PATH", help="atlas.csv manifest")
    parser.add_argument("--output-dir", type=Path, default=d.outputDir, metavar="PATH", help="Output directory")
    parser.add_argument(
        "--cluster-key", type=str, default=d.clusterKey, metavar="COL", help="Leiden cluster obs column"
    )
    parser.add_argument(
        "--min-cells-per-profile",
        type=int,
        default=d.minCellsPerProfile,
        metavar="N",
        help="Minimum cells per pseudobulk profile",
    )
    parser.add_argument(
        "--min-samples-per-area",
        type=int,
        default=d.minSamplesPerArea,
        metavar="N",
        help="Minimum samples in the one group",
    )
    parser.add_argument(
        "--min-studies-per-area",
        type=int,
        default=d.minStudiesPerArea,
        metavar="N",
        help="Minimum studies in the one group",
    )
    parser.add_argument("--labels-only", action="store_true", help="Write eligibility labels and exit before DE")
    parser.add_argument(
        "--no-write-transferred",
        action="store_true",
        help="Skip writing atlas_with_clusters.h5ad after cluster transfer",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    cfg = _DEFAULT_CFG.model_copy(
        update={
            "inputAtlasH5ad": args.input_atlas,
            "harmonyAtlasH5ad": args.harmony_atlas,
            "contextsPath": args.contexts,
            "atlasCsvPath": args.atlas_csv,
            "outputDir": args.output_dir,
            "clusterKey": args.cluster_key,
            "minCellsPerProfile": args.min_cells_per_profile,
            "minSamplesPerArea": args.min_samples_per_area,
            "minStudiesPerArea": args.min_studies_per_area,
            "writeTransferredAtlas": not args.no_write_transferred,
        }
    )

    log_run_separator(log)
    started = datetime.datetime.now()
    run_disease_markers(cfg, runDe=not args.labels_only, labelsOnly=args.labels_only)
    print(f"Disease markers finished in {datetime.datetime.now() - started}")


if __name__ == "__main__":
    main()
