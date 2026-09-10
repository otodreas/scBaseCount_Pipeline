import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import batched
from pathlib import Path
from time import sleep

import pandas as pd
from metadata.regexes import LUNG_TISSUE_RE
from shared.repo import REPO_ROOT
from study_context import fetch_study_accession

parser = argparse.ArgumentParser()
parser.add_argument("--output", help="Path to the output CSV file", type=Path, required=True)
args = parser.parse_args()
OUTPUT_FILE = args.output

if OUTPUT_FILE.suffix != ".csv":
    raise ValueError("Output file must have a .csv suffix")

if OUTPUT_FILE.exists():
    raise ValueError("File already exists and will be overwritten. Aborting.")

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)


# # Debug version of fetch_study_accession
# import numpy as np
# rng = np.random.default_rng()
# def fetch_study_accession(accession: str) -> str:
#     return rng.choice(np.arange(1000))

# OUTPUT_FILE = REPO_ROOT / "output/metadata/datasets_v2.csv"

ENA_MAX_REQUESTS_PER_SECOND = 50
MIN_OBS_COUNT_PER_STUDY = 1_000
EXCLUDED_STUDY_ACCESSIONS = (
    # Short ena context listed below
    "PRJEB51634",  # The goal of this project is to perform a systematic comparison of immune cell lineages across human tissues. To this end, we collected up to 16 tissues from twelve adult deceased organ donors, isolated immune cells and profiled them using single-cell RNA sequencing and VDJ sequencing generating a dataset of around 360,000 cells.
    "PRJNA1005589",  # Through analysis of gd T cells in mucosal and lymphoid tissues across the human lifespan ...
    "PRJNA1179423",  # Here, using single nucleus (sn)RNA-seq on temporal lobe tissue
    "PRJNA1188170",  # unclear
    "PRJNA1215450",  # Here, we took a systems approach to comprehensively profile RNA and surface protein expression of over 1.25 million immune cells isolated from blood, lymphoid organs, and mucosal tissues
    "PRJNA657844",  # Pancreatic Cancer Biopsies
    "PRJNA902813",  # Overall design: We sorted CD3+ T cells from four patients with lung cancer
)

metadata_pq: pd.DataFrame = pd.read_parquet(
    REPO_ROOT
    / "data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_sample_metadata.parquet"
)
metadata_pq = metadata_pq.set_index("srx_accession")

# Drop accessions that are not SRA/ENA
metadata_pq = metadata_pq.loc[~metadata_pq.index.str.startswith("NRX"), :]


# Check for lung tissue
metadata_pq["is_lung"] = metadata_pq["tissue"].str.contains(LUNG_TISSUE_RE, na=False)

# Get study accession (ENA portal: up to 50 read_experiment calls per second)
study_accessions: list[str | None] = []
with ThreadPoolExecutor(max_workers=ENA_MAX_REQUESTS_PER_SECOND) as pool:
    batches = list(batched(metadata_pq.index.tolist(), ENA_MAX_REQUESTS_PER_SECOND))
    for batch_number, batch in enumerate(batches):
        print(f"Getting {ENA_MAX_REQUESTS_PER_SECOND} study accessions for batch {batch_number + 1} of {len(batches)}")
        study_accessions.extend(pool.map(fetch_study_accession, batch))
        if batch_number + 1 < len(batches):
            sleep(1)

metadata_pq["study_accession"] = study_accessions

# Keep lung accessions, excluding known broad multi-organ studies
metadata_pq = metadata_pq.loc[
    metadata_pq["is_lung"] & ~metadata_pq["study_accession"].isin(EXCLUDED_STUDY_ACCESSIONS)
].dropna(subset=["study_accession"])
# Drop studies with less than MIN_OBS_COUNT_PER_STUDY observations
metadata_pq = metadata_pq.loc[
    metadata_pq.groupby("study_accession")["obs_count"].transform("sum") > MIN_OBS_COUNT_PER_STUDY,
    :,
]

metadata_pq.to_csv(OUTPUT_FILE)
print(f"Saved datasets with {len(metadata_pq.loc[metadata_pq['is_lung']])} lung samples to {OUTPUT_FILE}")
