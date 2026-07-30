from dotenv import load_dotenv
from h5ad_concat import H5adConcatConfig, run_h5ad_concat
from shared.repo import REPO_ROOT

load_dotenv(REPO_ROOT / ".env")

cfg = H5adConcatConfig(
    datasetsPath=REPO_ROOT / "output/metadata/datasets_v2.csv",
    atlasR2Key="lung/atlas_v2.h5ad",
    uploadAtlas=True,
    minPctCellsAfterQc=0.5,
    outputPath=REPO_ROOT / "output/atlas/v2/atlas_v2.h5ad",
)

print("Atlas concat pipeline started")
run_h5ad_concat(cfg)
print("Atlas concat pipeline completed")
