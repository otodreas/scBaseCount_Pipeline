from dotenv import load_dotenv
from h5ad_concat import H5adConcatConfig, run_h5ad_concat
from shared.repo import REPO_ROOT

load_dotenv(REPO_ROOT / ".env")

# Fields not defined here are left as their default values.
# See the definition of H5adConcatConfig in scripts/h5ad_concat/config.py for details.
# There, you will find qc settings.
cfg = H5adConcatConfig(
    datasetsPath=REPO_ROOT / "output/metadata/datasets_v2.csv",
    atlasR2Key="atlas/2026-08-12/atlas.h5ad",
    uploadAtlas=True,
    minPctCellsAfterQc=0.5,
    maxPctRibo=0.5,
    outputPath=REPO_ROOT / "output/atlas/2026-08-12/atlas.h5ad",
)

print("Atlas concat pipeline started")
run_h5ad_concat(cfg)
print("Atlas concat pipeline completed")
