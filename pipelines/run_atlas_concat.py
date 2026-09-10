import argparse
from pathlib import Path

from dotenv import load_dotenv
from h5ad_concat import H5adConcatConfig, run_h5ad_concat
from shared.repo import REPO_ROOT

load_dotenv(REPO_ROOT / ".env")

parser = argparse.ArgumentParser()
parser.add_argument("--datasets", help="Path to the datasets CSV file", type=Path, required=True)
parser.add_argument("--output", help="Path to the output atlas h5ad file", type=Path, required=True)
args = parser.parse_args()
DATASETS_FILE = args.datasets
OUTPUT_FILE = args.output

if DATASETS_FILE.suffix != ".csv":
    raise ValueError("Datasets file must have a .csv suffix")

if not DATASETS_FILE.exists():
    raise ValueError("Datasets file does not exist.")

if OUTPUT_FILE.suffix != ".h5ad":
    raise ValueError("Output file must have a .h5ad suffix")

if OUTPUT_FILE.exists():
    raise ValueError("File already exists and will be overwritten. Aborting.")

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# Fields not defined here are left as their default values.
# See the definition of H5adConcatConfig in scripts/h5ad_concat/config.py for details.
# There, you will find qc settings.
cfg = H5adConcatConfig(
    datasetsPath=DATASETS_FILE,
    atlasR2Key=OUTPUT_FILE,
    uploadAtlas=True,
    minPctCellsAfterQc=0.5,
    maxPctRibo=0.5,
    outputPath=OUTPUT_FILE,
)

print("Atlas concat pipeline started")
run_h5ad_concat(cfg)
print("Atlas concat pipeline completed")
