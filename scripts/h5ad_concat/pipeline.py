from pathlib import Path

import pandas as pd
from shared.logger import configure_file_logger
from storage import gcs_uri_to_r2_raw_key

from h5ad_concat.config import H5adConcatConfig
from h5ad_concat.exceptions import FileRejected
from h5ad_concat.merge import concat_atlas, write_atlas
from h5ad_concat.models import H5adConcatResult, SkippedFile
from h5ad_concat.outputs import (
    append_status_row,
    finalize_outputs,
    init_status_csv,
    status_csv_path,
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
    _log.info("Loaded %d input(s) from %s", len(inputs), datasets_path)
    return inputs


def run_h5ad_concat(cfg: H5adConcatConfig) -> H5adConcatResult:
    """Download, validate, and concatenate h5ad files from R2 into a local atlas."""
    inputs = load_concat_inputs(cfg.datasetsPath)
    _log.info("Starting h5ad_concat run: %d input(s)", len(inputs))

    csv_path = status_csv_path(cfg.outputPath)
    init_status_csv(csv_path)
    config_path = write_config_manifest(cfg, _log)
    reference = load_gene_reference(cfg.geneInfoPath)
    skipped: list[SkippedFile] = []
    adatas = []
    studies: list[str] = []
    accessions: list[str] = []

    try:
        for r2_key, accession, study_accession in inputs:
            try:
                adata, study_accession = prepare_adata(r2_key, accession, study_accession, cfg, reference, _log)
                adatas.append(adata)
                studies.append(study_accession)
                accessions.append(accession)
                append_status_row(csv_path, accession, r2_key, "success", "", study_accession)
            except FileRejected as exc:
                detail = f": {exc.__cause__}" if exc.__cause__ else ""
                _log.warning("%s: skipped (%s)%s", accession, exc.reason.value, detail)
                skipped.append(SkippedFile(r2Key=r2_key, accession=accession, reason=exc.reason))
                append_status_row(csv_path, accession, r2_key, "skip", exc.reason.value, "")

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
            studiesSeen=sorted(set(studies)),
            skipped=skipped,
            statusCsvPath=csv_path,
            configPath=config_path,
            atlasR2Key=cfg.atlasR2Key if cfg.uploadAtlas else None,
            conserveLayers=cfg.conserveLayers,
        )

        finalize_outputs(cfg, output_path, csv_path, result, _log)

        _log.info(
            "h5ad_concat run complete: %d concatenated, %d skipped",
            len(adatas),
            len(skipped),
        )
        return result
    except KeyboardInterrupt:
        _log.warning("h5ad_concat run interrupted (KeyboardInterrupt)")
        raise
