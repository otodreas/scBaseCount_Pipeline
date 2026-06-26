from __future__ import annotations

from pathlib import Path

from shared.repo import REPO_ROOT

from study_context.models import ExperimentContext

CONTEXTS_JSONL_PATH = REPO_ROOT / "output" / "context" / "contexts.jsonl"


def load_contexts_jsonl(path: Path = CONTEXTS_JSONL_PATH) -> dict[str, ExperimentContext]:
    """Load contexts.jsonl into a dict keyed by accession, raising if the file is missing."""
    if not path.is_file():
        raise FileNotFoundError(f"study context file not found at {path}")
    contexts: dict[str, ExperimentContext] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            ctx = ExperimentContext.model_validate_json(line)
            contexts[ctx.accession] = ctx
    return contexts


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
