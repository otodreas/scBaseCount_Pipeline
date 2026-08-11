# %% Exploratory entry point for sampled atlas discovery
"""Sample-oriented wrapper around the production atlas DE modules.

Use pipelines/run_atlas_de_analysis.py for full-atlas runs.
This script keeps a local exploratory path pointing at the 100k sample.
"""

from disease_markers.analysis import analyze_atlas
from disease_markers.config import AtlasDeAnalysisConfig
from shared.repo import REPO_ROOT

ATLAS_PATH = REPO_ROOT / "output/atlas/v2/post/production/atlas_v2_post_sample.h5ad"
OUTPUT_DIR = REPO_ROOT / "output/atlas/v2/agent_analysis/"


def main() -> None:
    cfg = AtlasDeAnalysisConfig(
        atlasPath=ATLAS_PATH,
        outputDir=OUTPUT_DIR,
        primaryBudget=20,
        extendedBudget=60,
    )
    print(f"Running exploratory atlas DE on {cfg.atlasPath}")
    print(f"Writing outputs under {cfg.outputDir}")
    summary = analyze_atlas(cfg, reuseCheckpoint=True)
    print(summary)
    print(
        "Done. Review noteworthy_gene_shortlist.csv and noteworthy_gene_extended.csv. "
        "Hand-authored literature proposals are intentionally not generated here."
    )


if __name__ == "__main__":
    main()
