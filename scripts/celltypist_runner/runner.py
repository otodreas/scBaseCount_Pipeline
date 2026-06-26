from __future__ import annotations

import scanpy as sc
from celltypist import annotate, models
from shared.logger import configure_file_logger

from celltypist_runner.config import CellTypistRunnerConfig

_log = configure_file_logger("celltypist_runner.log", __name__)


def annotate_celltypist(adata: sc.AnnData, cfg: CellTypistRunnerConfig) -> sc.AnnData:
    """Run CellTypist on a copy of adata and write predicted labels onto the input obs."""
    if cfg.geneSymbolCol not in adata.var:
        raise ValueError(f"gene symbol column {cfg.geneSymbolCol!r} not found in adata.var")

    if cfg.downloadIfMissing:
        models.download_models()

    model = models.Model.load(model=cfg.modelName)

    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=cfg.targetSum)
    sc.pp.log1p(work.X)

    # CellTypist expects gene symbols as var_names; keep Ensembl IDs before swapping.
    if "ensembl_id" not in work.var:
        work.var["ensembl_id"] = work.var_names
    work.var_names = work.var[cfg.geneSymbolCol].astype(str).tolist()
    work.var_names_make_unique()

    predictions = annotate(work, model, majority_voting=False)
    predicted = predictions.to_adata().obs["predicted_labels"].astype(str)
    adata.obs[cfg.predictedLabelKey] = predicted.reindex(index=adata.obs_names).values

    _log.info(
        "annotated %d cells with model=%s  n_labels=%d",
        adata.n_obs,
        cfg.modelName,
        adata.obs[cfg.predictedLabelKey].nunique(),
    )
    return adata
