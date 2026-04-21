from __future__ import annotations

import logging

from study_context.fetch import fetch_experiment_context
from study_context.models import BiologicalContext, ExperimentContext, StudyContext, TechnicalContext

_log = logging.getLogger(__name__)


def pipeline_for_accession_list(accessions: list[str]) -> list[ExperimentContext]:
    _log.info("New context acquisition run started (%d accession(s))", len(accessions))
    return [fetch_experiment_context(acc) for acc in accessions]


__all__ = [
    "pipeline_for_accession_list",
    "ExperimentContext",
    "StudyContext",
    "TechnicalContext",
    "BiologicalContext",
]
