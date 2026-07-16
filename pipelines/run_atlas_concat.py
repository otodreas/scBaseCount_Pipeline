from dotenv import load_dotenv
from h5ad_concat import H5adConcatConfig, run_h5ad_concat
from shared.repo import REPO_ROOT

load_dotenv(REPO_ROOT / ".env")

cfg = H5adConcatConfig(
    datasetsPath=REPO_ROOT / "tests/datasets_sample20.csv",
    # atlasR2Key="lung/atlas_v1.h5ad",
    uploadAtlas=False,
    outputPath=REPO_ROOT / "output" / "atlas" / "data" / "atlas_sample20.h5ad",
)

run_h5ad_concat(cfg)
