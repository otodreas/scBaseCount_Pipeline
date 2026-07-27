import pandas as pd
import scanpy as sc

from disease_markers.config import DiseaseMarkersConfig
from disease_markers.labels import sample_labels_by_srx


def annotate_atlas(
    adata: sc.AnnData,
    label_table: pd.DataFrame,
    cfg: DiseaseMarkersConfig,
) -> sc.AnnData:
    """Map diseaseArea/isControl onto cells and keep only eligible samples."""
    labels = sample_labels_by_srx(label_table)
    sample_key = cfg.sampleKey

    disease_areas: list[str] = []
    is_control: list[bool] = []
    eligible_mask: list[bool] = []
    for srx in adata.obs[sample_key].astype(str):
        row = labels.get(srx)
        if row is None:
            disease_areas.append("Other")
            is_control.append(False)
            eligible_mask.append(False)
            continue
        disease_areas.append(row.diseaseArea)
        is_control.append(row.isControl)
        eligible_mask.append(row.eligible)

    out = adata.copy()
    out.obs["diseaseArea"] = pd.Series(disease_areas, index=out.obs_names, dtype="category")
    out.obs["isControl"] = pd.Series(is_control, index=out.obs_names)
    keep = pd.Series(eligible_mask, index=out.obs_names)
    return out[keep.to_numpy()].copy()
