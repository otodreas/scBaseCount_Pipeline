# Bioinformatics Master's Project Plan

## General Information

Name of student: **Oliver E. Todreas**
Contact email: **[ol6437ek-s@student.lu.se](mailto:ol6437ek-s@student.lu.se)**
Name of supervisor: **Parashar Dhapola**
Name of department: **Nygen Analytics AB**
Contact email: **[parashar@nygen.io](mailto:parashar@nygen.io)**
Course code: **BINP39**
Start and end date (expected): **March 23 2026–August 28 2026**

## Abstract

Disease-focused single-cell atlases will be built from SC Base Count, a large standardized scRNA-seq resource that currently lacks systematic clustering and annotation. Disease priorities will be informed by a structured survey of CELLxGENE Census metadata. The work implements a reproducible pipeline (scVI, scanpy, CyteType; AnnData; GitHub) to curate 2–3 disease areas, process datasets individually over an extended per-dataset phase, integrate batches, and benchmark scGPT embeddings against scVI-, PCA/Harmony-, and STATE (Arc Institute) annotation workflows. Evaluation combines concordance metrics (e.g. ARI, NMI) with documented edge cases for manuscript-ready comparison.

## Project Plan

### Introduction

Aggregating many scRNA-seq studies requires consistent preprocessing, batch integration, and interpretable cell typing [1,2]. scBaseCount standardizes data at scale but does not yet provide systematic clustering and annotation, limiting cross-cohort disease comparison. Nygen’s CyteType delivers evidence-traced, ontology-linked annotations suitable for auditing at scale. How variational integration (scVI) [3], linear reductions with batch correction (scanpy/Harmony) [1], and foundation embeddings (scGPT) [4] interact with CyteType on this resource is an open, testable question [5]. Candidate disease areas will be evaluated using CELLxGENE Census metadata (tissues, dataset counts, cells, available fields) before joint selection of 2–3 priorities aligned with the company’s therapeutic focus.

**References.** [1] Wolf, Angerer & Theis (2018) *Genome Biology*. [2] Regev *et al.* (2017) *Nature*. [3] Lopez *et al.* (2018) *Nature Methods*. [4] Cui *et al.* (2024) *Nature Methods*. [5] Heumos *et al.* (2023) *Nature Biotechnology* (best-practice single-cell analysis context).

### Hypothesis and Specific Aims

**H1:** For at least one prioritized disease area, scVI integration will yield CyteType label assignments that differ from per-dataset processing in a pattern measurable by concordance metrics (e.g. ARI, NMI); whether integration improves or harms annotation consistency will be assessed per dataset and subset.

**H2:** CyteType on scGPT-derived structures will not be identical to scVI- or PCA/Harmony-based workflows on the same cells; concordance and discordance versus STATE (Arc Institute) will be quantified and interpreted for biological plausibility.

**Aims:** (A) Survey CELLxGENE Census metadata; assess 5–6 candidate disease areas; jointly select 2–3 priority areas and subset SC Base Count accordingly. (B) Per dataset: scVI, multi-resolution clustering, subclustering as needed, CyteType, documented quality. (C) Integrate per disease area; compare integrated versus individual CyteType outputs. (D) Benchmark scGPT + CyteType against the above and STATE in a structured four-way comparison.

### Methods

**Data:** SC Base Count disease subsets; AnnData; parameterized pipelines with logging and version control. Phase 1 uses CELLxGENE Census metadata to shortlist disease contexts with sufficient datasets and cells.

**Procedure:** Phases follow the extended host specification—**Phase 1 (weeks 1–3)** disease atlas curation and selection; **Phase 2 (weeks 4–10)** per-dataset scVI, clustering/subclustering, manual cluster QC, CyteType, edge-case logs; **Phase 3 (weeks 11–16)** scVI integration per disease, CyteType on integrated objects, comparison to per-dataset annotations, manuscript-style notes; **Phase 4 (weeks 17–20)** scGPT embeddings and CyteType versus scVI-based (Phases 2–3), scanpy PCA/Harmony + CyteType, and STATE (Arc Institute). Tools: Python, scanpy, scvi-tools, scGPT, CyteType SDK.

**Evaluation:** Hypotheses are tested with (i) quantitative agreement between partitions or label vectors where cells are comparable; (ii) qualitative review of discordant clusters. **Support H1:** sustained cross-batch coherence of labels after integration without loss of disease-relevant separation. **Refute H1:** integration collapses biologically distinct states or degrades metrics without compensating batch alignment. **Support H2:** non-random, interpretable differences across embedding sources at matched resolution. **Refute H2:** near-identical outputs across scGPT, scVI, and PCA/Harmony at chosen clustering depth. The pipeline serves these biological comparisons; success is defined by evaluable metrics and documented interpretation, not by tool deployment alone.

### Time plan


| Weeks | Phase   | Activities                                                                                                                                                                                             |
| ----- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1–3   | Phase 1 | Literature; **CELLxGENE Census** metadata survey; structured assessment of 5–6 candidate disease areas (tissues, dataset counts, cells, metadata); joint selection of 2–3 priorities                   |
| 4–10  | Phase 2 | Per-dataset scVI training; multi-resolution clustering and subclustering where warranted; manual cluster inspection; CyteType; document annotation quality, edge cases, and CyteType improvement areas |
| 11–16 | Phase 3 | scVI batch integration per disease area; CyteType on integrated data; compare with individual-dataset annotations; evaluate where integration helps or hurts; manuscript-ready documentation           |
| 17–20 | Phase 4 | scGPT embeddings on constructed atlases; CyteType on scGPT-derived clusters; systematic comparison to scVI-, PCA/Harmony-, and STATE-based workflows; structured comparison for manuscript             |


**Course deliverables:** Reserve time in the final project weeks for synthesis, final report drafting, and GitHub + README polish; submit the **final report and repository one week before the seminar presentation** (per course guidelines). Weekly progress reviews with the supervisor continue throughout; remote or in-office work; optional Nygen Lund office days; cloud, compute, and LLM resources as provided by the host.