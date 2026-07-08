import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

from atlas_integration.config import AtlasIntegrationConfig
from atlas_integration.models import BatchMixingMetrics, ClusterConservationMetrics


def compute_batch_mixing(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> BatchMixingMetrics:
    """Estimate batch mixing as the mean fraction of same-study neighbors in PCA vs Harmony space."""
    sample_idx = _sample_indices(adata.n_obs, cfg.metricsSampleSize)
    batch_labels = adata.obs[cfg.batchKey].astype(str).to_numpy()[sample_idx]

    uncorrected = _mean_same_label_neighbor_fraction(adata.obsm[cfg.pcaKey][sample_idx], batch_labels)
    corrected = _mean_same_label_neighbor_fraction(adata.obsm[cfg.harmonyKey][sample_idx], batch_labels)
    return BatchMixingMetrics(
        meanSameStudyNeighborFractionUncorrected=uncorrected,
        meanSameStudyNeighborFractionCorrected=corrected,
    )


def compute_cluster_conservation(adata: sc.AnnData, cfg: AtlasIntegrationConfig) -> ClusterConservationMetrics:
    """Compare uncorrected and corrected Leiden partitions and silhouette scores."""
    uncorrected_labels = adata.obs[cfg.leidenKeyUncorrected].astype(str)
    corrected_labels = adata.obs[cfg.leidenKeyAtlas].astype(str)
    sample_idx = _sample_indices(adata.n_obs, cfg.metricsSampleSize)

    uncorrected_sample = uncorrected_labels.iloc[sample_idx]
    corrected_sample = corrected_labels.iloc[sample_idx]

    return ClusterConservationMetrics(
        ariUncorrectedVsCorrected=float(adjusted_rand_score(uncorrected_sample, corrected_sample)),
        nmiUncorrectedVsCorrected=float(normalized_mutual_info_score(uncorrected_sample, corrected_sample)),
        silhouetteUncorrectedEmbedding=_safe_silhouette(adata.obsm[cfg.pcaKey][sample_idx], uncorrected_sample),
        silhouetteCorrectedEmbedding=_safe_silhouette(adata.obsm[cfg.harmonyKey][sample_idx], corrected_sample),
    )


def _sample_indices(n_obs: int, sample_size: int) -> np.ndarray:
    if n_obs <= sample_size:
        return np.arange(n_obs)
    rng = np.random.default_rng(0)
    return np.sort(rng.choice(n_obs, size=sample_size, replace=False))


def _mean_same_label_neighbor_fraction(embedding: np.ndarray, labels: np.ndarray, n_neighbors: int = 15) -> float:
    n_neighbors = min(n_neighbors, len(labels) - 1)
    if n_neighbors < 1:
        return 1.0

    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric="euclidean")
    nn.fit(embedding)
    neighbor_idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]

    same_label_fractions: list[float] = []
    for row_idx, neighbors in enumerate(neighbor_idx):
        neighbor_labels = labels[neighbors]
        same_label_fractions.append(float(np.mean(neighbor_labels == labels[row_idx])))
    return float(np.mean(same_label_fractions))


def _safe_silhouette(embedding: np.ndarray, labels: pd.Series) -> float:
    if labels.nunique() <= 1:
        return float("nan")
    return float(silhouette_score(embedding, labels))
