# STATE vs Leiden disagreement

Exploratory analysis on `STATE x leiden_cluster` matrices produced by [`pipelines/cluster_stats.py`](../../pipelines/cluster_stats.py), run `cluster_stats/clustered_20260509`. Source matrices live in [`output/cluster_stats/clustered_20260509/cluster_stats.json`](../../output/cluster_stats/clustered_20260509/cluster_stats.json).

Leiden clusters from the clustering pipeline are passed to CyteType for annotation. STATE labels are treated as weak priors, informing which clustering resolution is selected. This report showcases the agreement between STATE clusters and STATE-informed leiden clusters on a cluster-by-cluster basis.

The number of STATE names that appear across the data make gleaning any memorable insights difficult, since a reader is unlikely to remember the distribution of STATE names according to a given metric. The figures presented in this report are intended to serve as a proof of concept that normalized Shannon entropy on the STATE x Leiden matrix can highlight the extent to which STATE labels agree or disagree with Leiden cluster labels.

Analysis script: [state_vs_leiden.py](state_vs_leiden.py).

## Leiden cluster size does not inform entropy of STATE labels within Leiden clusters

![KDE of STATE-mixing NSE by cluster-size quartile](.figs/nse_leiden_kde_by_size_quartile.png)

***Figure 1.*** *KDE of within-cluster STATE-label NSE (`nse_leiden`), split by log2(Leiden cluster size) quartile.*

The distribution of NSE is consistent across Leiden cluster quartiles, suggesting that the number of cells in a Leiden cluster does not influence the entropy of STATE labels in the cluster.

## Two complementary views on STATE-Leiden agreement

The next two figures are both organized one row per STATE, but they pivot on different entropies and answer different questions. STATE labels represented in fewer than 7 datasets are not shown.

![Per-STATE distribution of nse_STATE across accessions](.figs/nse_state_by_state_boxplots.png)

***Figure 2.*** *Per-STATE distribution of `nse_STATE` (entropy across Leiden clusters within a single accession), one observation per accession, sorted by median (top) and variance (bottom) and colored by cell category.*

A STATE with low median NSE and low variance, like neutrophils or basophils, consistently concentrates in only a few Leiden clusters across accessions. That makes the STATE easy for a clustering-based annotator to "find", provided the cluster it ends up in is also dominated by it.

At the other extreme, some STATE labels sit at consistently high NSE in every accession: their cells are spread broadly across Leiden clusters in every dataset where they appear. This is the pattern expected of STATE labels that bundle several transcriptionally distinct populations under one umbrella, or whose expression signature overlaps with that of co-occurring STATEs.

The normalization in NSE applies to each individual measurement, not to the distribution of medians across STATEs. In principle, every STATE could have a median near zero (every label consistently concentrated in a few Leiden clusters) or every STATE could have a median near one (every label consistently spread across many clusters). What we see instead is a continuum: some STATE labels are reliably well-isolated by Leiden, others are reliably mixed in with other STATEs, and many sit somewhere in between. The most natural reading is that STATE labels in this collection are not all at the same biological granularity: some name a tight, transcriptionally coherent population that Leiden picks out, others act as broad bins that Leiden splits across several clusters.

The variance panel layers a second axis on top of this: low variance means a STATE behaves the same way in every accession (uniformly concentrated or uniformly spread), high variance means its placement in the clustering shifts with the dataset, which is consistent with accession-specific composition, tissue context, or batch effects influencing where its cells land.

![Per-STATE cell and cluster shares across nse_leiden quartiles](.figs/state_behavior_across_nse_quartiles.png)

***Figure 3.*** *Per-STATE share of cells (left) and Leiden clusters (right) across `nse_leiden` quartiles, ordered by ascending median `nse_STATE` to match ***Figure 2*** (top).*

Reading top to bottom, STATEs progress from concentrating within few Leiden clusters per accession to spreading across many. A STATE label near the top with most of its cells in q1 sits in Leiden clusters dominated by a single STATE label, so downstream analyses that only see the cluster-level expression effectively have access to STATE label signal. A STATE whose cells fall mostly in q3 or q4 sits in mixed clusters, meaning that downstream analyses using Leiden clusters may be "blind" to STATE labels because the cluster's expression is driven by whichever STATE dominates that mix.

There is no clear vertical gradient in **Figure 3**: a STATE's median `nse_STATE` does not predict the `nse_leiden` quartiles of the Leiden clusters its cells occupy. Concentrating within few Leiden clusters per accession (low `nse_STATE`) does not imply that those clusters are themselves dominated by a single STATE (low `nse_leiden`), and the two entropies should be read as complementary rather than redundant.

Read together: **Figure 2** highlights whether a STATE label *lands in few Leiden clusters*, while **Figure 3** highlights whether *those clusters are clean*.

## Conclusion

STATE labels are not ground truth; they are weak priors used to guide clustering resolution. For this reason, disagreement between STATE labels and Leiden clusters is not a cause for concern. However, it is noteworthy that there are differences in agreement across STATE labels (**Figure 2**, **3**), suggesting that certain STATE labels systemically exhibit less agreement to Leiden clusters than others.

The three figures together provide that characterization:

- **Figure 1** establishes that within-cluster STATE entropy (`nse_leiden`) is largely independent of Leiden cluster size, with a modest tendency for the largest clusters to be more STATE-pure. Cluster size on its own is therefore a poor proxy for STATE-Leiden agreement.
- **Figure 2** ranks STATE labels by how concentrated they are across Leiden clusters within an accession (`nse_STATE`). STATEs at the top of the median panel (e.g. neutrophils, basophils) consistently land in a few clusters across accessions; STATE labels further down the panel consistently spread across many Leiden clusters. NSE is normalized per measurement, so a continuum of medians is not a foregone conclusion: the medians could all be near zero or all near one. The fact that they instead span the full range suggests STATE labels here are not at uniform biological granularity, with some naming a tight cell type that Leiden isolates and others acting as broad bins that Leiden splits across many clusters.
- **Figure 3**, ordered to match the top panel of **Figure 2, shows where the cells of each STATE land along the `nse_leiden` axis. The absence of a vertical gradient confirms that being concentrated (low `nse_STATE`) does not imply being in clean clusters (low `nse_leiden`); the two entropies have to be checked together.

note caveats like  put leave interpretation on the x axis

hand draw what we could have seen, and if so, what could we expect to see? x axis of KDE, mixing independent of state (no trend on fig 2 x)

