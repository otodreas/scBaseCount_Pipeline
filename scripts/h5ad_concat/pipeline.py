from __future__ import annotations

from pathlib import Path

from shared.files import safe_delete
from shared.logger import configure_file_logger
from study_context.utils import load_contexts_jsonl

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_prepared, read_h5ad_shape
from h5ad_concat.models import H5adConcatResult, SkippedFile
from h5ad_concat.prepare import accession_from_r2_key, prepare_accession

_log = configure_file_logger("h5ad_concat.log", __name__)


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    # TODO(datasets-csv): resolve cfg.r2Keys from output/metadata/datasets.csv accessions
    # mapped to R2 raw keys (see pipelines/run_clustering_pipeline.py).
    contexts = load_contexts_jsonl(cfg.contextsPath)
    skipped: list[SkippedFile] = []
    prepared_paths: list[tuple[Path, str]] = []

    # TODO(stream-pipeline): interleave prepare and merge so peak staging is one batch of prepared
    # files, not all passing files. After mergeBatchSize files pass validation, fold them into the
    # partial atlas via concat_prepared, delete that batch, and continue. Move orchestration here or
    # into a new merge entry point that accepts an r2-key iterator instead of a full prepared list.

    for r2_key in cfg.r2Keys:
        accession = accession_from_r2_key(r2_key)
        try:
            path, study_accession = prepare_accession(r2_key, cfg, contexts, _log)
            prepared_paths.append((path, study_accession))
        except FileRejected as exc:
            _log.warning("%s: skipped (%s)", accession, exc.reason.value)
            skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))

    if not prepared_paths:
        msg = "No files passed validation; nothing to concatenate"
        raise ValueError(msg)

    output_path = concat_prepared(prepared_paths, cfg, _log)

    # TODO(upload-atlas): upload output_path to R2 via upload_to_r2 with _local_md5_b64 metadata.

    n_obs, n_vars = read_h5ad_shape(output_path)
    studies_seen = sorted({study for _, study in prepared_paths})

    for path in cfg.cacheDir.glob("**/*"):
        if path.is_file() and path != output_path:
            safe_delete(path, _log)

    return H5adConcatResult(
        outputPath=output_path,
        nObs=n_obs,
        nVars=n_vars,
        nFilesConcatenated=len(prepared_paths),
        studiesSeen=studies_seen,
        skipped=skipped,
    )
