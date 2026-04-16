from __future__ import annotations

from ena_context.fetch import fetch_experiment_context
from ena_context.models import BiologicalContext, ExperimentContext, StudyContext, TechnicalContext


def pipeline_for_accession_list(accessions: list[str]) -> list[ExperimentContext]:
    return [fetch_experiment_context(acc) for acc in accessions]


__all__ = [
    "pipeline_for_accession_list",
    "ExperimentContext",
    "StudyContext",
    "TechnicalContext",
    "BiologicalContext",
]
