from __future__ import annotations

import os
from pathlib import Path

import scanpy as sc
from cytetype import CyteType, rank_genes_groups_backed
from shared.logger import configure_file_logger

from cytetype_runner.config import N_TOP_GENES, CyteTypeRunnerConfig, CyteTypeRunResult

_log = configure_file_logger("cytetype_runner.log", __name__)


def require_api_key() -> str:
    """Return CYTETYPE_API_KEY from the environment, raising if missing or empty."""
    api_key = os.environ.get("CYTETYPE_API_KEY", "").strip()
    if not api_key:
        raise RuntimeError(
            "CYTETYPE_API_KEY is not set or empty. "
            "Add it to your .env file or export it in your shell before running."
        )
    return api_key


def run_cytetype(
    cfg: CyteTypeRunnerConfig,
    input_path: Path,
    group_key: str,
    study_context: str,
    metadata: dict[str, str] | None = None,
) -> CyteTypeRunResult:
    output_path = cfg.outputDir / f"{cfg.srxAccession}_cytetype_annotated.h5ad"
    adata = sc.read_h5ad(input_path, backed="r")
    rank_genes_groups_backed(adata, groupby=group_key, use_raw=False, key_added=f"rank_genes_{group_key}")
    annotator = CyteType(adata, group_key, rank_key=f"rank_genes_{group_key}", n_top_genes=N_TOP_GENES)
    run_kwargs: dict = {"study_context": study_context, "auth_token": require_api_key()}
    if metadata:
        run_kwargs["metadata"] = metadata
    adata = annotator.run(**run_kwargs)

    job_details = adata.uns.get("cytetype_jobDetails")
    if not isinstance(job_details, dict):
        job_details = {}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    _log.info("%s: annotation written to %s", cfg.srxAccession, output_path)

    return CyteTypeRunResult(
        outputPath=output_path,
        reportUrl=str(job_details.get("report_url", "")),
        jobId=str(job_details.get("job_id", "")),
        apiUrl=str(job_details.get("api_url", "")),
    )
