import json
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from ena_context import ExperimentContext, pipeline_for_accession_list

load_dotenv()

output_path = Path("output/contexts.jsonl")

READ_FROM_CACHE = False

if READ_FROM_CACHE and output_path.exists():
    with open(output_path) as f:
        contexts = [ExperimentContext.model_validate_json(line) for line in f]
else:
    SAMPLE_SIZE = False
    data = pd.read_csv("metadata_analysis/v2_lung/datasets.csv")
    if SAMPLE_SIZE:
        data = data.sample(n=SAMPLE_SIZE, random_state=42)
    accessions = data["srx_accession"].tolist()
    contexts = pipeline_for_accession_list(accessions)
    with open(output_path, "w") as f:
        for ctx in contexts:
            f.write(ctx.model_dump_json() + "\n")


# descriptions = [ctx.study.studyDescription for ctx in contexts]

# data = pd.DataFrame({"srx_accession": accessions, "study_description": descriptions})
# data.to_csv("output/datasets_with_study_description.csv", index=False)

# missing_descriptions = []
# for ctx in contexts:
#     if ctx.study.studyDescription is None:
#         missing_descriptions.append(ctx.accession)
# if len(missing_descriptions) > 0:
#     print(f"Missing descriptions for {len(missing_descriptions)} accessions:")
#     for accession in missing_descriptions:
#         print(accession)
# else:
#     print("All accessions have descriptions")