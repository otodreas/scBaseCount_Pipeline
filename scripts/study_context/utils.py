from __future__ import annotations

from study_context.models import ExperimentContext


def experiment_context_summary(ctx: ExperimentContext) -> str:
    parts = []
    if ctx.biological.tissueType:
        parts.append(f"Tissue: {ctx.biological.tissueType}")
    if ctx.technical.libraryStrategy:
        parts.append(f"Library strategy: {ctx.technical.libraryStrategy}")
    if ctx.technical.libraryConstructionProtocol:
        parts.append(f"Library prep: {ctx.technical.libraryConstructionProtocol}")
    if ctx.study and ctx.study.studyDescription:
        parts.append(f"Description: {ctx.study.studyDescription}")
    if ctx.study and ctx.study.pubmedAbstract:
        parts.append(f"Abstract: {ctx.study.pubmedAbstract}")

    return ". ".join(parts) if parts else "No descriptive fields found."
