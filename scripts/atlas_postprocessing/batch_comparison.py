import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from shared.repo import REPO_ROOT, rel_to_repo

from atlas_postprocessing.artifacts import write_json
from atlas_postprocessing.config import AtlasPostprocessingConfig
from atlas_postprocessing.core import run_harmony_on_pcs, timed
from atlas_postprocessing.scib import run_scib_benchmark

log = logging.getLogger(__name__)

STUDY_KEY = "study_accession"
SRX_KEY = "SRX_accession"
TECH_KEY = "tech_10x"
STUDY_TECH_KEY = "batch_study_tech"
EXISTING_HARMONY_KEY = "X_pca_harmony"
PCA_KEY = "X_pca"

EXPECTED_BASELINE = {
    "nTopGenes": 2000,
    "nPcs": 50,
    "nNeighbors": 15,
    "resolution": 0.8,
}


def slugify(value: str) -> str:
    """Filesystem-safe slug for embedding and evaluation labels."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "unnamed"


def harmony_obsm_key(batchKey: str) -> str:
    return f"X_pca_harmony__{slugify(batchKey)}"


def load_run_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Run JSON not found: {rel_to_repo(path)}")
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Run JSON must be an object: {rel_to_repo(path)}")
    return payload


def validate_baseline_manifests(
    *,
    subsetRun: dict[str, Any],
    productionRun: dict[str, Any],
    expected: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Require subset and production manifests to share the comparison baseline."""
    resolved_expected = EXPECTED_BASELINE if expected is None else expected
    fields = ("nTopGenes", "nPcs", "nNeighbors", "resolution")
    subset_values = {field: subsetRun.get(field) for field in fields}
    production_values = {field: productionRun.get(field) for field in fields}
    mismatches: list[str] = []
    for field in fields:
        if subset_values[field] != production_values[field]:
            mismatches.append(f"{field}: subset={subset_values[field]!r} production={production_values[field]!r}")
        if subset_values[field] != resolved_expected[field]:
            mismatches.append(f"{field}: subset={subset_values[field]!r} expected={resolved_expected[field]!r}")
    if mismatches:
        raise ValueError(
            "Subset and production baselines disagree with the fixed comparison baseline: " + "; ".join(mismatches)
        )
    return {
        "expected": dict(resolved_expected),
        "subset": subset_values,
        "production": production_values,
    }


def attach_tech_10x(
    adata: sc.AnnData,
    datasetsPath: Path,
    *,
    accessionKey: str = SRX_KEY,
    studyKey: str = STUDY_KEY,
    techKey: str = TECH_KEY,
) -> dict[str, Any]:
    """Join ``tech_10x`` onto ``adata.obs`` by experiment accession with strict validation."""
    if accessionKey not in adata.obs:
        raise ValueError(f"adata.obs is missing accession key {accessionKey!r}")
    if studyKey not in adata.obs:
        raise ValueError(f"adata.obs is missing study key {studyKey!r}")
    if not datasetsPath.is_file():
        raise FileNotFoundError(f"datasets CSV not found: {rel_to_repo(datasetsPath)}")

    catalog = pd.read_csv(
        datasetsPath,
        usecols=["srx_accession", "study_accession", techKey],
        dtype="string",
    )
    catalog = catalog.rename(
        columns={
            "srx_accession": accessionKey,
            "study_accession": f"{studyKey}_catalog",
        }
    )
    if catalog[accessionKey].isna().any() or (catalog[accessionKey].str.strip() == "").any():
        raise ValueError(f"{rel_to_repo(datasetsPath)} has blank srx_accession values")
    duplicated = catalog.loc[catalog[accessionKey].duplicated(), accessionKey].tolist()
    if duplicated:
        raise ValueError(f"{rel_to_repo(datasetsPath)} has duplicate srx_accession values: {duplicated[:5]}")

    obs = adata.obs[[accessionKey, studyKey]].copy()
    obs[accessionKey] = obs[accessionKey].astype("string")
    obs[studyKey] = obs[studyKey].astype("string")
    joined = obs.merge(catalog, on=accessionKey, how="left", validate="many_to_one")

    missing_catalog = int(joined[f"{studyKey}_catalog"].isna().sum())
    if missing_catalog:
        missing_ids = sorted(joined.loc[joined[f"{studyKey}_catalog"].isna(), accessionKey].unique().tolist())
        raise ValueError(f"{missing_catalog} cells lack a catalog row for {accessionKey}; examples: {missing_ids[:5]}")
    missing_tech = int(joined[techKey].isna().sum()) + int((joined[techKey].str.strip() == "").sum())
    if missing_tech:
        blank_ids = sorted(
            joined.loc[joined[techKey].isna() | (joined[techKey].str.strip() == ""), accessionKey].unique().tolist()
        )
        raise ValueError(f"{missing_tech} cells have blank {techKey}; examples: {blank_ids[:5]}")

    study_mismatches = int((joined[studyKey] != joined[f"{studyKey}_catalog"]).sum())
    if study_mismatches:
        bad = joined.loc[
            joined[studyKey] != joined[f"{studyKey}_catalog"], [accessionKey, studyKey, f"{studyKey}_catalog"]
        ]
        examples = bad.drop_duplicates().head(5).to_dict(orient="records")
        raise ValueError(f"{study_mismatches} cells have study_accession mismatches; examples: {examples}")

    adata.obs[techKey] = joined[techKey].to_numpy()
    n_experiments = int(obs[accessionKey].nunique())
    n_study_tech = int(joined[[studyKey, techKey]].drop_duplicates().shape[0])
    audit = {
        "datasetsPath": rel_to_repo(datasetsPath),
        "nCells": int(adata.n_obs),
        "nExperiments": n_experiments,
        "nStudyTechBatches": n_study_tech,
        "techValues": sorted(joined[techKey].astype(str).unique().tolist()),
        "missingCatalog": 0,
        "missingTech": 0,
        "studyMismatches": 0,
    }
    log.info(
        "Attached %s for %s experiments (%s study×tech batches)",
        techKey,
        n_experiments,
        n_study_tech,
    )
    return audit


def make_study_tech_batch_key(
    adata: sc.AnnData,
    *,
    studyKey: str = STUDY_KEY,
    techKey: str = TECH_KEY,
    outCol: str = STUDY_TECH_KEY,
    separator: str = "|",
) -> dict[str, Any]:
    """Create a collision-safe study×technology batch column."""
    if studyKey not in adata.obs:
        raise ValueError(f"adata.obs is missing study key {studyKey!r}")
    if techKey not in adata.obs:
        raise ValueError(f"adata.obs is missing tech key {techKey!r}")
    if adata.obs[studyKey].astype(str).str.contains(re.escape(separator), regex=True).any():
        raise ValueError(f"{studyKey} values contain separator {separator!r}")
    if adata.obs[techKey].astype(str).str.contains(re.escape(separator), regex=True).any():
        raise ValueError(f"{techKey} values contain separator {separator!r}")

    composed = adata.obs[studyKey].astype(str) + separator + adata.obs[techKey].astype(str)
    adata.obs[outCol] = composed.astype("string")
    n_batches = int(adata.obs[outCol].nunique())
    log.info("Created %s with %s unique batches", outCol, n_batches)
    return {"column": outCol, "nBatches": n_batches, "separator": separator}


def ensure_shared_pca(adata: sc.AnnData, *, nPcs: int) -> None:
    """Confirm the frozen subset already carries a shared PCA embedding."""
    if PCA_KEY not in adata.obsm:
        raise ValueError(f"adata.obsm[{PCA_KEY!r}] is required for fixed-subset comparison")
    n_computed = int(adata.obsm[PCA_KEY].shape[1])
    if n_computed < nPcs:
        raise ValueError(f"Shared PCA has {n_computed} components; need at least {nPcs}")
    if EXISTING_HARMONY_KEY not in adata.obsm:
        raise ValueError(f"adata.obsm[{EXISTING_HARMONY_KEY!r}] is required to reuse study Harmony")


def preserve_study_harmony_embedding(
    adata: sc.AnnData,
    *,
    studyKey: str = STUDY_KEY,
    sourceKey: str = EXISTING_HARMONY_KEY,
) -> str:
    """Copy the existing study Harmony embedding under an explicit keyed name."""
    target = harmony_obsm_key(studyKey)
    adata.obsm[target] = np.asarray(adata.obsm[sourceKey]).copy()
    log.info("Preserved %s as %s", sourceKey, target)
    return target


def run_harmony_variants(
    adata: sc.AnnData,
    cfg: AtlasPostprocessingConfig,
    *,
    batchKeys: list[str],
    nPcs: int,
    skipExisting: bool = True,
) -> dict[str, str]:
    """Run Harmony for each batch key into distinct ``obsm`` slots."""
    embedding_map: dict[str, str] = {}
    for batch_key in batchKeys:
        obsm_key = harmony_obsm_key(batch_key)
        embedding_map[batch_key] = obsm_key
        if skipExisting and obsm_key in adata.obsm:
            log.info("Reuse existing Harmony embedding %s for batch key %s", obsm_key, batch_key)
            continue
        corrected = timed(
            f"Harmony batch={batch_key}",
            lambda key=batch_key: run_harmony_on_pcs(adata, cfg, nPcs=nPcs, batchKey=key),
            logger=log,
        )
        adata.obsm[obsm_key] = corrected
    return embedding_map


def batch_cardinalities(adata: sc.AnnData, batchKeys: list[str]) -> dict[str, int]:
    return {key: int(adata.obs[key].nunique()) for key in batchKeys}


def run_scib_evaluation_grid(
    adata: sc.AnnData,
    *,
    outRoot: Path,
    embeddingKeys: list[str],
    evalBatchKeys: list[str],
    labelKey: str,
    nJobs: int,
    force: bool,
) -> list[dict[str, Any]]:
    """Run one scIB benchmark per evaluation batch key across all embeddings."""
    outRoot.mkdir(parents=True, exist_ok=True)
    runs: list[dict[str, Any]] = []
    for eval_key in evalBatchKeys:
        run_dir = outRoot / f"eval_batch={slugify(eval_key)}"
        timed(
            f"scIB eval_batch={eval_key}",
            lambda run_dir=run_dir, eval_key=eval_key: run_scib_benchmark(
                adata,
                outDir=run_dir,
                batchKey=eval_key,
                labelKey=labelKey,
                embeddingKeys=embeddingKeys,
                preIntegratedKey=PCA_KEY,
                nJobs=nJobs,
                force=force,
            ),
            logger=log,
        )
        runs.append(
            {
                "evalBatchKey": eval_key,
                "csv": rel_to_repo(run_dir / "scib_results.csv"),
                "svg": rel_to_repo(run_dir / "scib_results.svg"),
                "embeddingKeys": list(embeddingKeys),
            }
        )
    return runs


def aggregate_scib_results(
    runs: list[dict[str, Any]],
    *,
    embeddingMap: dict[str, str],
    outCsv: Path,
) -> Path:
    """Combine per-evaluation scIB tables into one long-form CSV."""
    reverse_map = {obsm_key: batch_key for batch_key, obsm_key in embeddingMap.items()}
    frames: list[pd.DataFrame] = []
    for run in runs:
        csv_path = Path(run["csv"])
        if not csv_path.is_absolute():
            csv_path = REPO_ROOT / csv_path
        frame = pd.read_csv(csv_path, index_col=0)
        frame = frame.reset_index(names="embedding")
        frame.insert(0, "evalBatchKey", run["evalBatchKey"])
        frame.insert(
            1,
            "harmonyBatchKey",
            frame["embedding"].map(lambda key: reverse_map.get(str(key), "uncorrected" if key == PCA_KEY else None)),
        )
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    outCsv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(outCsv, index=False)
    log.info("Wrote aggregated scIB matrix %s", rel_to_repo(outCsv))
    return outCsv


def write_comparison_summary(path: Path, payload: dict[str, Any]) -> Path:
    return write_json(path, payload)
