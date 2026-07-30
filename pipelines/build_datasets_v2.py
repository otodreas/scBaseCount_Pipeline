from concurrent.futures import ThreadPoolExecutor
from itertools import batched
from time import sleep

import pandas as pd
from metadata.regexes import LUNG_TISSUE_RE
from shared.repo import REPO_ROOT
from study_context import fetch_study_accession

# # Debug version of fetch_study_accession
# import numpy as np
# rng = np.random.default_rng()
# def fetch_study_accession(accession: str) -> str:
#     return rng.choice(np.arange(1000))

ENA_MAX_REQUESTS_PER_SECOND = 50
MIN_OBS_COUNT_PER_STUDY = 1_000
EXCLUDED_STUDY_ACCESSIONS = (
    "PRJEB51634",
    "PRJNA1005589",
    "PRJNA1179423",
    "PRJNA1188170",
    "PRJNA1215450",
    "PRJNA657844",
    "PRJNA902813",
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


metadata_pq.to_csv(REPO_ROOT / "output/metadata/datasets_v2.csv")
print(
    f"Saved datasets with {len(metadata_pq.loc[metadata_pq['is_lung']])} lung samples to {REPO_ROOT / 'output/metadata/datasets_v2.csv'}"
)
