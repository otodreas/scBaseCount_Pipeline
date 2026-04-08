# Project Plan Submission v2

The file `Project_Plan_v2.pdf` is the project plan.

## Changes from submission v1

### Requested changes

1. The entire project hinges on: compare methods and see if they differ. There is no way to evaluate if the results are correct. I don't see any suggestion on how sides are right in case of conflict.
2. There is a circular evaluation problem, beause CyteType depends on the embedding structure. You are testing embedding &rarr; annotation pipeline, But evaluating using annotation outputs influenced by the same embedding.
3. Hypothesis is weakly falsifiable. You hypothesize that methods will differ in measurable ways. This is obvious: they operate completely differently. What are you testing here? That the methods are different or that they can measure the same data in the same way.
4. Overall, the hypothesis is:
  - Not biologically meaningful
  - Not risky
  - Not strongly falsifiable
  Even the “refutation” case (methods agree) is unlikely.
5. There is no formal statistical framework; there certainly is no biological validation layer (it's not mandatory, but it is added to the concerns of no validation).
6. The company needs to provide a statement that they agree that the outcome of this project will be published on GitHub.
7. The proposed project is too ambitious.

### Changes

1. `Project_Plan_v2.pdf` now treats expert-curated author cell-type labels from CELLxGENE-aligned data as the primary quantitative reference. Agreement is defined per cell type (not only method-vs-method), and conflicting pipeline outputs are scored against those labels; STATE provides an additional embedding-independent check (Abstract, Introduction, Aim (D), Methods — Evaluation).
2. The plan states that CyteType uses cluster-level differential expression, marker genes, and LLM interpretation with tissue context, not a direct readout of embedding coordinates, and frames the scientific question as whether embeddings shift cluster boundaries enough to change markers and biological calls. STATE is used as an external reference that does not depend on the embedding under test (Abstract, Introduction, Phase 4).
3. The hypothesis was replaced: annotation agreement with author labels is expected to be robust to embedding choice for well-characterized populations but may diverge for rare or transitional states. Either strong cross-method agreement with references everywhere or systematic divergence for ambiguous populations counts as an informative outcome. Support/refute criteria in Methods — Evaluation are aligned with that framing.
4. The comparison is anchored in clinically relevant resolution (e.g. T cell exhaustion, regulatory populations, polarized macrophages), and the plan states that a null result (such as foundation embeddings not improving disease-relevant separation over standard PCA) is scientifically meaningful (Abstract). The falsifiable hypothesis and evaluation bullets tie results to author labels and marker programs, not only to pairwise method difference.
5. Methods — Evaluation specifies adjusted Rand index, normalized mutual information, and per-cell-type F1 against author labels (with ontology or naming harmonization where types are merged or split), plus literature-defined gene-set scores for expected programs (exhaustion, cytotoxicity, polarization, etc.) on annotated populations, without wet-lab validation.
6. The file `project_specification_v2.pdf` (which was submitted in the first round) is signed by the CEO of Nygen and explicitly grants permission for the project to be published on GitHub (see §3.3). The work is also hosted in [this repo](https://github.com/otodreas/scBaseCount_Pipeline). The same permission is summarized in `Project_Plan_v2.pdf` under General Information (**Publication**).
7. Aim (A), the Introduction, the Abstract, and the Phase 1 row in the time plan now commit to **two** priority disease areas, with a **third only if schedule permits**, reducing default scope while keeping an escape hatch if timelines allow.