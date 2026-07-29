from concurrent.futures import ThreadPoolExecutor
from itertools import batched
from time import sleep

import pandas as pd
from metadata.regexes import LUNG_TISSUE_RE
from shared.repo import REPO_ROOT
from study_context import fetch_study_accession

ENA_MAX_REQUESTS_PER_SECOND = 50

metadata_pq = pd.read_parquet(
    REPO_ROOT
    / "data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens/scbasecount_2026-01-12_metadata_GeneFull_Homo_sapiens_sample_metadata.parquet"
)
metadata_pq = metadata_pq[["srx_accession", "disease", "tissue"]].set_index("srx_accession")

# Check for lung tissue
metadata_pq["isLung"] = metadata_pq["tissue"].str.contains(LUNG_TISSUE_RE, na=False)

# Get study accession (ENA portal: up to 50 read_experiment calls per second)
study_accessions: list[str | None] = []
accession_list = metadata_pq.index.tolist()
with ThreadPoolExecutor(max_workers=ENA_MAX_REQUESTS_PER_SECOND) as pool:
    batches = list(batched(accession_list, ENA_MAX_REQUESTS_PER_SECOND))
    for batch_index, batch in enumerate(batches):
        print(f"Getting {ENA_MAX_REQUESTS_PER_SECOND} study accessions for batch {batch_index + 1} of {len(batches)}")
        study_accessions.extend(pool.map(fetch_study_accession, batch))
        if batch_index + 1 < len(batches):
            sleep(1)

metadata_pq["study_accession"] = study_accessions

metadata_pq.to_csv(REPO_ROOT / "output/metadata/datasets_v2.csv")
