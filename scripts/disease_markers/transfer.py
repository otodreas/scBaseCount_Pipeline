from pathlib import Path

import scanpy as sc


def transfer_leiden_clusters(
    adata: sc.AnnData,
    harmonyH5ad: Path,
    clusterKey: str = "leiden_atlas",
) -> sc.AnnData:
    """Attach cluster labels from the Harmony h5ad onto adata by obs_names without changing X."""
    harmony = sc.read_h5ad(harmonyH5ad, backed="r")
    if clusterKey not in harmony.obs.columns:
        msg = f"{clusterKey!r} missing from Harmony atlas obs"
        raise KeyError(msg)

    labels = harmony.obs[clusterKey]
    if not adata.obs_names.equals(harmony.obs_names):
        labels = labels.reindex(adata.obs_names)
    missing = int(labels.isna().sum())
    if missing:
        msg = f"cluster join left {missing} cells without a {clusterKey} label"
        raise ValueError(msg)

    out = adata.copy()
    out.obs[clusterKey] = labels.astype(str).astype("category")
    return out


def load_full_atlas_transfer_clusters(
    inputAtlasH5ad: Path,
    harmonyH5ad: Path,
    clusterKey: str = "leiden_atlas",
) -> sc.AnnData:
    """Read the full-gene atlas and attach Harmony cluster labels."""
    adata = sc.read_h5ad(inputAtlasH5ad)
    return transfer_leiden_clusters(adata, harmonyH5ad, clusterKey=clusterKey)
