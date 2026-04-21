from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor

from study_context.fetch import fetch_experiment_context
from study_context.models import BiologicalContext, ExperimentContext, StudyContext, TechnicalContext

_log = logging.getLogger(__name__)


def pipeline_for_accession_list(
    accessions: list[str], *, max_workers: int = 8
) -> list[ExperimentContext]:
    _log.info("New context acquisition run started (%d accession(s))", len(accessions))
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        return list(pool.map(fetch_experiment_context, accessions))


__all__ = [
    "pipeline_for_accession_list",
    "ExperimentContext",
    "StudyContext",
    "TechnicalContext",
    "BiologicalContext",
]
