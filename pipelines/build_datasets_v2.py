from concurrent.futures import ThreadPoolExecutor
from itertools import batched
from time import sleep

import pandas as pd
from metadata.regexes import LUNG_TISSUE_RE
from shared.repo import REPO_ROOT
from study_context import fetch_study_accession

# # Debug version of fetch_study_accession
# def fetch_study_accession(accession: str) -> str:
#     return accession

ENA_MAX_REQUESTS_PER_SECOND = 50

metadata_pq = pd.read_parquet(
    REPO_ROOT
    / "data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_sample_metadata.parquet"
)
metadata_pq = metadata_pq[["srx_accession", "disease", "tissue", "czi_collection_id", "czi_collection_name"]].set_index(
    "srx_accession"
)

# Check for lung tissue
metadata_pq["is_lung"] = metadata_pq["tissue"].str.contains(LUNG_TISSUE_RE, na=False)

# Get study accession (ENA portal: up to 50 read_experiment calls per second)
ignore_accessions = metadata_pq.index.to_series().str.startswith("NRX")
fetch_accessions = ignore_accessions.index[~ignore_accessions]
study_accessions: dict[str, str | None] = {}
with ThreadPoolExecutor(max_workers=ENA_MAX_REQUESTS_PER_SECOND) as pool:
    batches = list(batched(fetch_accessions, ENA_MAX_REQUESTS_PER_SECOND))
    for batch_number, batch in enumerate(batches):
        print(f"Getting {ENA_MAX_REQUESTS_PER_SECOND} study accessions for batch {batch_number + 1} of {len(batches)}")
        study_accessions.update(
            zip(
                batch,
                pool.map(fetch_study_accession, batch),
                strict=True,
            )
        )
        if batch_number + 1 < len(batches):
            sleep(1)

metadata_pq["study_accession"] = metadata_pq.index.to_series().map(study_accessions)

lung_studies = metadata_pq.loc[metadata_pq["is_lung"], "study_accession"].dropna().unique()
metadata_pq["in_atlas"] = metadata_pq["is_lung"] | metadata_pq["study_accession"].isin(lung_studies)

metadata_pq.to_csv(REPO_ROOT / "output/metadata/datasets_v2.csv")
