from __future__ import annotations

import os
from pathlib import Path

import orjson
import scanpy as sc
from cytetype import CyteType, rank_genes_groups_backed
from shared.logger import configure_file_logger

from cytetype_runner.config import N_TOP_GENES, CyteTypeRunnerConfig

_log = configure_file_logger("cytetype_runner.log", __name__)


def write_job_details(cfg: CyteTypeRunnerConfig, h5ad_path: Path) -> Path:
    adata = sc.read_h5ad(h5ad_path, backed="r")
    job_details_path = cfg.jobDetailsDir / f"{cfg.srxAccession}_cytetype_jobDetails.json"
    job_details_path.parent.mkdir(parents=True, exist_ok=True)
    job_details_path.write_bytes(orjson.dumps({cfg.srxAccession: adata.uns.get("cytetype_jobDetails")}))
    _log.info("%s: job details written to %s", cfg.srxAccession, job_details_path)
    return job_details_path


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
) -> Path:
    output_path = cfg.outputDir / f"{cfg.srxAccession}_cytetype_annotated.h5ad"
    adata = sc.read_h5ad(input_path, backed="r")
    rank_genes_groups_backed(adata, groupby=group_key, use_raw=False, key_added=f"rank_genes_{group_key}")
    annotator = CyteType(adata, group_key, rank_key=f"rank_genes_{group_key}", n_top_genes=N_TOP_GENES)
    run_kwargs: dict = {"study_context": study_context, "auth_token": require_api_key()}
    if metadata:
        run_kwargs["metadata"] = metadata
    adata = annotator.run(**run_kwargs)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    _log.info("%s: annotation written to %s", cfg.srxAccession, output_path)

    write_job_details(cfg, output_path)

    return output_path
