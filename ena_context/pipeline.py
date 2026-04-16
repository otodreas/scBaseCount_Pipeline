from __future__ import annotations

from ena_context.fetch import fetch_experiment_context
from ena_context.models import ExperimentContext


def pipeline_for_accession(accession: str) -> ExperimentContext:
    return fetch_experiment_context(accession)


def pipeline_for_accession_list(accessions: list[str]) -> list[ExperimentContext]:
    return [fetch_experiment_context(acc) for acc in accessions]
