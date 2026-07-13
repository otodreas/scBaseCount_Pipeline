from __future__ import annotations

from shared.logger import configure_file_logger
from study_context.utils import load_contexts_jsonl

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import H5adConcatResult, SkippedFile
from h5ad_concat.prepare import accession_from_r2_key, prepare_adata

_log = configure_file_logger("h5ad_concat.log", __name__)


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    # TODO(datasets-csv): resolve cfg.r2Keys from output/metadata/datasets.csv accessions
    # mapped to R2 raw keys (see pipelines/run_clustering_pipeline.py).
    contexts = load_contexts_jsonl(cfg.contextsPath)
    skipped: list[SkippedFile] = []
    adatas = []
    studies: list[str] = []

    for r2_key in cfg.r2Keys:
        accession = accession_from_r2_key(r2_key)
        try:
            adata, study_accession = prepare_adata(r2_key, accession, cfg, contexts, _log)
            adatas.append(adata)
            studies.append(study_accession)
        except FileRejected as exc:
            detail = f": {exc.__cause__}" if exc.__cause__ else ""
            # Warn on skipped files, continue pipeline
            _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
            skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))

    if not adatas:
        msg = "No files passed validation; nothing to concatenate"
        raise ValueError(msg)

    try:
        atlas = concat_atlas(adatas, cfg, _log)
    except Exception:
        _log.exception("concat failed")
        raise

    try:
        output_path = write_atlas(atlas, cfg, _log)
    except Exception:
        _log.exception("write failed")
        raise

    # TODO(upload-atlas): upload output_path to R2 via upload_to_r2 with _local_md5_b64 metadata.

    return H5adConcatResult(
        outputPath=output_path,
        nObs=atlas.n_obs,
        nVars=atlas.n_vars,
        nFilesConcatenated=len(adatas),
        studiesSeen=sorted(set(studies)),
        skipped=skipped,
    )
