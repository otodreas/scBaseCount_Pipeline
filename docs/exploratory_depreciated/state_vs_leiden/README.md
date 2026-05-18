# STATE vs leiden clustering assessment

Exploratory analysis on `STATE x leiden_cluster` matrices produced by [`pipelines/cluster_stats.py`](../../../pipelines/cluster_stats.py), run `cluster_stats/clustered_20260509`. Source matrices live in [`output/cluster_stats/clustered_20260509/cluster_stats.json`](../../../output/cluster_stats/clustered_20260509/cluster_stats.json).

Leiden clusters from the clustering pipeline are passed to CyteType for annotation. STATE labels are treated as weak priors during clustering, informing the clustering resolution selection process. The goal here is to flag cell types whose STATE labels overlap poorly with leiden clusters, since those are the ones most likely to disagree with CyteType.

Analysis script: [state_vs_leiden.py](state_vs_leiden.py).

## Figures

<!-- Lognormal fit of leiden cluster sizes pooled across accessions, plus the KS goodness-of-fit statistic.

<img src=".figs/cluster_size_lognormal_fit.png" width="600"/> -->

KDE of STATE-mixing NSE within each cluster, split by log2(cluster size) quartile, with a one-way ANOVA across quartiles.

<img src=".figs/nse_leiden_kde_by_size_quartile.png" width="500"/>

NSE against log2(group size) for STATEs (left) and leiden clusters (right) with linear and LOWESS fits and Spearman rho.

<img src=".figs/nse_vs_size_hexbin.png" width="700"/>

Distribution of `nse_STATE` conditioned on the NSE tertile of the leiden cluster it appears in.

<!-- <img src=".figs/nse_state_kde_by_nse_leiden_bin.png" width="500"/>

Per-STATE share of cells and clusters that fall in low, mid, and high NSE leiden tertiles. STATEs sorted by share of cells in the low tertile. -->

<!-- <img src=".figs/state_behavior_across_nse_tertiles.png" width="700"/>

Per-STATE distribution of `nse_STATE` across accessions, sorted by median (top) and variance (bottom), colored by cell category. -->

<img src=".figs/nse_state_by_state_boxplots.png" width="600"/>
