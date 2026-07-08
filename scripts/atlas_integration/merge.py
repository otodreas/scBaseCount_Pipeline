import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from study_context.models import ExperimentContext
from study_context.utils import load_contexts_jsonl

from atlas_integration.config import AtlasIntegrationConfig
from atlas_integration.models import MergeStats

_log = logging.getLogger(__name__)

_UNKNOWN_CELL_TYPE_VALUES = {"", "nan", "none", "na", "unknown"}


def read_datasets_csv(path: Path) -> pd.DataFrame:
    """Load the datasets catalog CSV and return rows with srx_accession and file_path."""
    df = pd.read_csv(path)
    required = {"srx_accession", "file_path"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return df


def build_accession_study_map(
    accessions: list[str],
    contexts_path: Path,
) -> dict[str, str]:
    """Map each accession to studyAccession, falling back to the accession when study context is missing."""
    contexts = load_contexts_jsonl(contexts_path)
    mapping: dict[str, str] = {}
    for accession in accessions:
        mapping[accession] = study_for_accession(accession, contexts)
    return mapping


def study_for_accession(accession: str, contexts: dict[str, ExperimentContext]) -> str:
    """Return studyAccession for an accession, or the accession itself when no study block is present."""
    ctx = contexts.get(accession)
    if ctx is None:
        _log.warning("%s: missing from contexts.jsonl; using accession as batch key", accession)
        return accession
    if ctx.study is None or not ctx.study.studyAccession:
        _log.warning("%s: study context missing studyAccession; using accession as batch key", accession)
        return accession
    return ctx.study.studyAccession


def normalize_cell_type_labels(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> None:
    """Fill missing or blank cell_type values with cfg.missingLabel in adata.obs."""
    if cfg.cellTypeKey not in adata.obs:
        adata.obs[cfg.cellTypeKey] = cfg.missingLabel
        return

    raw = adata.obs[cfg.cellTypeKey]
    labels = pd.Series(raw, index=raw.index, dtype="object").map(
        lambda value: None if pd.isna(value) else str(value).strip()
    )
    normalized = labels.fillna("")
    missing = labels.isna() | normalized.eq("") | normalized.str.lower().isin(_UNKNOWN_CELL_TYPE_VALUES)
    adata.obs[cfg.cellTypeKey] = labels.where(~missing, cfg.missingLabel)


def qc_accession(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> sc.AnnData:
    """Apply per-accession QC filters aligned with cluster_validation thresholds."""
    if "gene_symbols" not in adata.var:
        raise ValueError("adata.var is missing required column 'gene_symbols'")

    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    adata.var["ribo"] = adata.var["gene_symbols"].str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

    sc.pp.filter_cells(adata, min_genes=cfg.minGenesPerCell)
    sc.pp.filter_genes(adata, min_cells=cfg.minCellsPerGene)
    adata = adata[adata.obs["pct_counts_mt"] < cfg.maxPctCountsMt].copy()
    return adata


def load_accession_h5ad(path: Path, accession: str, gs_uri: str | None = None) -> sc.AnnData:
    """Read a local or remote h5ad for one accession and make obs names unique."""
    if path.exists():
        adata = sc.read(str(path))
    elif gs_uri is not None:
        adata = sc.read(gs_uri)
    else:
        raise FileNotFoundError(f"No h5ad found for {accession} at {path}")

    adata.obs_names_make_unique()
    return adata


def prepare_accession_adata(
    adata: sc.AnnData,
    accession: str,
    study: str,
    cfg: AtlasIntegrationConfig,
) -> sc.AnnData | None:
    """QC one accession, attach metadata columns, and optionally subsample cells per study."""
    n_cells_original = adata.n_obs
    adata = qc_accession(adata, cfg)
    if adata.n_obs < cfg.minCellsTotal:
        _log.warning(
            "%s: skipped after QC (%d cells < minCellsTotal=%d)",
            accession,
            adata.n_obs,
            cfg.minCellsTotal,
        )
        return None

    adata.obs[cfg.accessionKey] = accession
    adata.obs[cfg.batchKey] = study
    normalize_cell_type_labels(adata, cfg)

    if cfg.subsampleN is not None and adata.n_obs > cfg.subsampleN:
        sc.pp.subsample(adata, n_obs=cfg.subsampleN, random_state=0)

    adata.layers["counts"] = adata.X.copy()
    _log.info(
        "%s: kept %d/%d cells after QC (study=%s)",
        accession,
        adata.n_obs,
        n_cells_original,
        study,
    )
    return adata


def concat_accession_adatas(adatas: list[sc.AnnData]) -> sc.AnnData:
    """Concatenate per-accession AnnData objects on the intersection of shared genes."""
    if not adatas:
        raise ValueError("concat_accession_adatas requires at least one AnnData object")
    merged = ad.concat(adatas, join="inner", merge="same", label="atlas_batch", index_unique="-")
    merged.obs_names_make_unique()
    return merged


def build_merged_adata(cfg: AtlasIntegrationConfig) -> tuple[sc.AnnData, MergeStats]:
    """Load, QC, and concatenate all accessions listed in cfg.datasetsCsvPath."""
    datasets = read_datasets_csv(cfg.datasetsCsvPath)
    accessions = datasets["srx_accession"].astype(str).tolist()
    study_map = build_accession_study_map(accessions, cfg.contextsPath)

    merged_parts: list[sc.AnnData] = []
    skipped: list[str] = []

    for _, row in datasets.iterrows():
        accession = str(row["srx_accession"])
        local_path = cfg.localH5adRoot / f"{accession}.h5ad"
        try:
            adata = load_accession_h5ad(local_path, accession, str(row["file_path"]))
        except FileNotFoundError:
            _log.warning("%s: h5ad not found locally or remotely", accession)
            skipped.append(accession)
            continue

        prepared = prepare_accession_adata(adata, accession, study_map[accession], cfg)
        if prepared is None:
            skipped.append(accession)
            continue
        merged_parts.append(prepared)

    if not merged_parts:
        raise ValueError("No accessions passed QC; atlas merge produced zero cells")

    merged = concat_accession_adatas(merged_parts)
    stats = MergeStats(
        nAccessionsRequested=len(accessions),
        nAccessionsMerged=len(merged_parts),
        nAccessionsSkipped=len(skipped),
        nCellsFinal=merged.n_obs,
        nGenesFinal=merged.n_vars,
        nStudies=int(merged.obs[cfg.batchKey].nunique()),
        skippedAccessions=skipped,
    )
    return merged, stats
