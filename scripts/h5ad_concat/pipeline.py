from __future__ import annotations

import anndata as ad
from shared.logger import configure_file_logger
from study_context.utils import load_contexts_jsonl

from h5ad_concat.checkpoint import load_checkpoint
from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import fold_atlas, write_atlas
from h5ad_concat.models import ConcatManifest, H5adConcatResult, ManifestEntry
from h5ad_concat.prepare import accession_from_r2_key, prepare_adata

_log = configure_file_logger("h5ad_concat.log", __name__)


def _maybe_checkpoint(
    atlas: ad.AnnData | None,
    pending: list[ad.AnnData],
    manifest: ConcatManifest,
    cfg: H5adConcatConfig,
) -> ad.AnnData | None:
    """Fold pending into atlas and write when an atlas exists."""
    if pending:
        atlas = fold_atlas(atlas, pending, cfg)
        pending.clear()
    if atlas is not None:
        write_atlas(atlas, manifest, cfg, _log)
    return atlas


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    # TODO(datasets-csv): resolve cfg.r2Keys from output/metadata/datasets.csv accessions
    # mapped to R2 raw keys (see pipelines/run_clustering_pipeline.py).
    contexts = load_contexts_jsonl(cfg.contextsPath)
    atlas, manifest = load_checkpoint(cfg, _log)
    processed = manifest.processedKeys()
    pending: list[ad.AnnData] = []
    since_checkpoint = 0
    new_this_run = 0

    for r2_key in cfg.r2Keys:
        if r2_key in processed:
            continue

        accession = accession_from_r2_key(r2_key)
        try:
            adata, study_accession = prepare_adata(r2_key, cfg, contexts, _log)
            pending.append(adata)
            manifest.entries.append(
                ManifestEntry(
                    r2Key=r2_key,
                    accession=accession,
                    concatenated=True,
                    study=study_accession,
                )
            )
        except FileRejected as exc:
            _log.warning("%s: skipped (%s)", accession, exc.reason.value)
            manifest.entries.append(
                ManifestEntry(
                    r2Key=r2_key,
                    accession=accession,
                    concatenated=False,
                    reason=exc.reason,
                )
            )

        new_this_run += 1
        since_checkpoint += 1
        if cfg.checkpointEvery > 0 and since_checkpoint >= cfg.checkpointEvery:
            atlas = _maybe_checkpoint(atlas, pending, manifest, cfg)
            since_checkpoint = 0

    if atlas is None and not pending and not any(entry.concatenated for entry in manifest.entries):
        msg = "No files passed validation; nothing to concatenate"
        raise ValueError(msg)

    if pending:
        atlas = fold_atlas(atlas, pending, cfg)
        pending.clear()

    if new_this_run > 0:
        if atlas is None:
            msg = "No files passed validation; nothing to concatenate"
            raise ValueError(msg)
        write_atlas(atlas, manifest, cfg, _log)

    if atlas is None:
        msg = "No files passed validation; nothing to concatenate"
        raise ValueError(msg)

    # TODO(upload-atlas): upload output_path to R2 via upload_to_r2 with _local_md5_b64 metadata.

    return H5adConcatResult(
        outputPath=cfg.outputPath,
        nObs=atlas.n_obs,
        nVars=atlas.n_vars,
        nFilesConcatenated=manifest.concatenatedCount(),
        studiesSeen=manifest.studiesSeen(),
        skipped=manifest.skippedFiles(),
    )
