import pandas as pd
from ena_context import pipeline_for_accession, pipeline_for_accession_list

SAMPLE_SIZE = 0

data = pd.read_csv("metadata_analysis/v2_lung/datasets.csv")
if SAMPLE_SIZE:
    data = data.sample(n=SAMPLE_SIZE, random_state=42)

accessions = data["srx_accession"].tolist()

contexts = pipeline_for_accession_list(accessions)
descriptions = [ctx.study.studyDescription for ctx in contexts]
data = pd.DataFrame({"srx_accession": accessions, "study_description": descriptions})
data.to_csv("output/datasets_with_study_description.csv", index=False)

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