# disease_markers

Per-experiment disease area, disease status, biological control, and atlas eligibility labels derived from ontology metadata (`disease_ontology_term_id`, `tissue_ontology_term_id`), with narrow free-text and study-consensus fallbacks. The package also provides atlas-cluster concordance, candidate screening, and study-aware validation helpers used by analysis notebooks.

## Usage

```python
from pathlib import Path

from disease_markers import (
    annotate_obs_with_labels,
    build_sample_label_table,
    classify_cluster_candidates,
    cluster_support_table,
    concordance_summary,
    same_study_contrast_support,
)
from metadata.config import MetadataConfig

label_table = build_sample_label_table(
    Path("output/context/contexts_v2.jsonl"),
    Path("output/atlas/v2/atlas_v2.csv"),
    MetadataConfig().sampleParquetPath,
)
obs = annotate_obs_with_labels(adata.obs, label_table)
support = cluster_support_table(obs)
contrast = same_study_contrast_support(obs)
candidates = classify_cluster_candidates(support, contrast)
concordance = concordance_summary(
    obs,
    clusterKey="leiden_atlas",
    labelKey="cell_type",
    studyKey="study_accession",
)
```

Join labels onto an AnnData object by experiment accession (for example `SRX_accession` in `obs`). Each label row applies to the full experiment H5AD.

## Modules

| Module | Role |
|--------|------|
| `labels.py` | Build per-SRX disease and eligibility labels from MONDO/UBERON caches |
| `status.py` | Infer diseased / control status and comparator arms |
| `concordance.py` | Compare integrated clusters with preserved source labels |
| `candidates.py` | Screen clusters for shared and disease-associated support |
| `validation.py` | Same-study abundance and pseudobulk DE helpers |
| `config.py` | `AtlasDeAnalysisConfig` for full-atlas discovery runs |
| `aggregation.py` | Copy-conscious sparse sample x Leiden pseudobulk aggregation |
| `specificity.py` | Tau specificity, home-cluster, and study-direction helpers |
| `ranking.py` | Adaptive evidence scoring and review-budget shortlists |
| `plots.py` | Review figures driven by the shortlist |
| `analysis.py` | Analyze-stage orchestration from a pseudobulk checkpoint |
| `memory.py` | RAM preflight and process telemetry |

## Outputs

Production discovery artifacts under `output/atlas/v2/analysis/production/`:

- `checkpoints/pseudobulk.h5ad` and `checkpoints/aggregate_fingerprint.json`
- `noteworthy_gene_shortlist.csv`, `noteworthy_gene_extended.csv`, `candidate_thresholds.csv`
- `de_hits.csv`, `de_results.parquet`, `restricted_genes.csv`, `unexpected_expression_candidates.csv`
- `figures/` with score distributions, heatmaps, volcanoes, and study evidence panels

Exploratory sample runs can still write under `output/atlas/v2/agent_analysis/` via `notebooks/analysis/analyze_atlas_DE.py`.

## Config

Labeling logic lives in `labels.py` (`build_sample_label_table`) and uses `ontology_lookup.OntologyLookupConfig` for release-pinned MONDO/UBERON caches under `data/ontologies/`. Full-atlas discovery defaults live in `AtlasDeAnalysisConfig` (`primaryBudget=20`, `extendedBudget=60`, `memoryReserveBytes` sized for a ~2 TiB server). Candidate thresholds and DE filters in `candidates.py` / `validation.py` remain overridable by callers.
