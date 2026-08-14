from pathlib import Path

import pandas as pd
from shared.logger import configure_file_logger
from shared.repo import rel_to_repo
from storage import gcs_uri_to_r2_raw_key

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import (
    CELL_FILTER_ORDER,
    FileRecord,
    H5adConcatResult,
    QcSummary,
    SkippedFile,
)
from h5ad_concat.outputs import (
    append_file_record,
    ensure_atlas_targets_absent,
    file_log_path,
    finalize_outputs,
    init_file_log,
    write_config_manifest,
)
from h5ad_concat.prepare import prepare_adata
from h5ad_concat.reference import load_gene_reference

_log = configure_file_logger("h5ad_concat.log", __name__)


def load_concat_inputs(datasets_path: Path) -> list[tuple[str, str, str]]:
    """Load R2 key, experiment accession, and study accession from a datasets CSV."""
    if not datasets_path.is_file():
        raise FileNotFoundError(f"datasets file not found at {datasets_path}")

    datasets = pd.read_csv(datasets_path, dtype="string")
    required_columns = ["file_path", "srx_accession", "study_accession"]
    missing_columns = [column for column in required_columns if column not in datasets.columns]
    if missing_columns:
        raise ValueError(f"{datasets_path}: missing required columns: {', '.join(missing_columns)}")

    for column in required_columns:
        datasets[column] = datasets[column].str.strip()
        if bool(datasets[column].isna().any()) or bool(datasets[column].eq("").any()):
            raise ValueError(f"{datasets_path}: blank values in {column}")

    inputs = [
        (gcs_uri_to_r2_raw_key(gs_uri), accession, study_accession)
        for gs_uri, accession, study_accession in datasets[required_columns].itertuples(index=False, name=None)
    ]
    _log.info("Loaded %d input(s) from %s", len(inputs), rel_to_repo(datasets_path))
    return inputs


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    ensure_atlas_targets_absent(cfg)

    if cfg.datasetsPath is None:
        raise ValueError("datasetsPath is required")

    inputs = load_concat_inputs(cfg.datasetsPath)
    _log.info("Starting h5ad_concat run: %d input(s)", len(inputs))

    log_path = file_log_path(cfg.outputPath)
    init_file_log(log_path)
    config_path = write_config_manifest(cfg, _log)
    reference = load_gene_reference(cfg.geneInfoPath)
    skipped: list[SkippedFile] = []
    file_records: list[FileRecord] = []
    qc_summary = QcSummary()
    adatas = []
    studies: list[str] = []
    accessions: list[str] = []

    try:
        for r2_key, accession, study_accession in inputs:
            try:
                adata, study_accession, qc_stats = prepare_adata(
                    r2_key, accession, study_accession, cfg, reference, _log
                )
                adatas.append(adata)
                studies.append(study_accession)
                accessions.append(accession)
                record = FileRecord(
                    accession=accession,
                    studyAccession=study_accession,
                    r2Key=r2_key,
                    status="success",
                    skipReason=None,
                    qc=qc_stats,
                )
                if qc_stats is not None:
                    qc_summary.concatenatedFiles.add(qc_stats)
                    qc_summary.allQcProcessedFiles.add(qc_stats)
                append_file_record(log_path, record)
                file_records.append(record)
            except FileRejected as exc:
                detail = f": {exc.__cause__}" if exc.__cause__ else ""
                _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
                skipped.append(
                    SkippedFile(
                        r2Key=r2_key,
                        accession=accession,
                        reason=exc.reason,
                        studyAccession=study_accession,
                        qc=exc.qc,
                    )
                )
                record = FileRecord(
                    accession=accession,
                    studyAccession=study_accession,
                    r2Key=r2_key,
                    status="skip",
                    skipReason=exc.reason,
                    qc=exc.qc,
                )
                if exc.qc is not None:
                    qc_summary.allQcProcessedFiles.add(exc.qc)
                    dropped_by_filter = ", ".join(
                        f"{name}={count}" for name, count in exc.qc.nCellsDroppedByFilter.items()
                    )
                    _log.info(
                        "%s: QC before rejection kept %d/%d cells; dropped by filter: %s",
                        accession,
                        exc.qc.nCellsAfter,
                        exc.qc.nCellsBefore,
                        dropped_by_filter,
                    )
                append_file_record(log_path, record)
                file_records.append(record)

        if not adatas:
            msg = "No files passed validation; nothing to concatenate"
            raise ValueError(msg)

        try:
            atlas = concat_atlas(adatas, accessions, cfg, _log)
        except Exception:
            _log.exception("concat failed")
            raise

        try:
            output_path = write_atlas(atlas, cfg, _log)
        except Exception:
            _log.exception("write failed")
            raise

        result = H5adConcatResult(
            outputPath=output_path,
            nObs=atlas.n_obs,
            nVars=atlas.n_vars,
            nFilesConcatenated=len(adatas),
            nFilesSkipped=len(skipped),
            studiesSeen=sorted(set(studies)),
            skipped=skipped,
            cellFilterOrder=list(CELL_FILTER_ORDER),
            qcSummary=qc_summary,
            files=file_records,
            fileLogPath=log_path,
            configPath=config_path,
            atlasR2Key=cfg.atlasR2Key if cfg.uploadAtlas else None,
            conserveLayers=cfg.conserveLayers,
        )

        finalize_outputs(cfg, output_path, log_path, result, _log)

        _log.info(
            "h5ad_concat run complete: %d concatenated, %d skipped; cell drops by filter: %s",
            len(adatas),
            len(skipped),
            qc_summary.concatenatedFiles.nCellsDroppedByFilter,
        )
        return result
    except KeyboardInterrupt:
        _log.warning("h5ad_concat run interrupted (KeyboardInterrupt)")
        raise
