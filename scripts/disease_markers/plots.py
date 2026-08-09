"""Review-oriented plots for atlas noteworthy-gene discovery."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from disease_markers.candidates import study_balanced_weights
from disease_markers.concordance import as_string_series
from disease_markers.validation import same_study_case_control_profiles


def plot_heatmap(matrix: pd.DataFrame, path: Path, title: str, *, cmap: str = "viridis") -> None:
    if matrix.empty:
        return
    fig_w = max(4.0, 0.35 * matrix.shape[1] + 2.0)
    fig_h = max(3.5, 0.28 * matrix.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap=cmap)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=90, fontsize=7)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str), fontsize=7)
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_unexpected_dotplot(frame: pd.DataFrame, path: Path) -> None:
    if frame.empty:
        return
    plot = frame.copy()
    plot["label"] = plot["geneSymbol"].fillna(plot["gene"]).astype(str) + " | c" + plot["cluster"].astype(str)
    fig, ax = plt.subplots(figsize=(7.5, max(3.0, 0.28 * len(plot) + 1.0)))
    sizes = 40 + 200 * plot["detectionDelta"].abs().clip(upper=1.0)
    colors = np.where(plot["log2FoldChange"] >= 0, "C0", "C3")
    ax.scatter(plot["log2FoldChange"], range(len(plot)), s=sizes, c=colors, alpha=0.85)
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels(plot["label"], fontsize=7)
    ax.axvline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("disease vs control log2FC")
    ax.set_title("Unexpected expression candidates")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_score_distributions(
    candidates: pd.DataFrame,
    thresholds: pd.DataFrame,
    path: Path,
) -> None:
    if candidates.empty:
        return
    classes = sorted(candidates["proposalClass"].dropna().astype(str).unique())
    n = max(len(classes), 1)
    fig, axes = plt.subplots(n, 1, figsize=(7.0, max(2.5, 1.8 * n)), sharex=False)
    if n == 1:
        axes = [axes]
    threshold_map = (
        thresholds.set_index("proposalClass")
        if not thresholds.empty and "proposalClass" in thresholds.columns
        else None
    )
    for ax, proposal_class in zip(axes, classes, strict=True):
        scores = candidates.loc[candidates["proposalClass"].astype(str).eq(proposal_class), "evidenceScore"]
        ax.hist(scores.astype(float), bins=20, color="0.55", edgecolor="white")
        if threshold_map is not None and proposal_class in threshold_map.index:
            primary_cut = threshold_map.loc[proposal_class, "primaryCutoff"]
            extended_cut = threshold_map.loc[proposal_class, "extendedCutoff"]
            if np.isfinite(primary_cut):
                ax.axvline(float(primary_cut), color="C0", ls="--", label="primary cutoff")
            if np.isfinite(extended_cut):
                ax.axvline(float(extended_cut), color="C1", ls=":", label="extended cutoff")
            ax.legend(fontsize=7, frameon=False)
        ax.set_title(proposal_class, fontsize=9)
        ax.set_ylabel("count")
    axes[-1].set_xlabel("evidence score")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(
    frame: pd.DataFrame,
    path: Path,
    title: str,
    labelGenes: list[str],
    *,
    padjThreshold: float = 0.05,
    lfcThreshold: float = 1.0,
) -> None:
    plot = frame.copy()
    plot["log2FoldChange"] = pd.to_numeric(plot["log2FoldChange"], errors="coerce")
    plot["padj"] = pd.to_numeric(plot["padj"], errors="coerce")
    log2FoldChange = plot["log2FoldChange"].to_numpy(dtype=float)
    adjustedPValue = plot["padj"].to_numpy(dtype=float)
    finite = np.isfinite(log2FoldChange) & np.isfinite(adjustedPValue) & (adjustedPValue >= 0)
    plot = plot.loc[finite].copy()
    if plot.empty:
        return

    log2FoldChange = log2FoldChange[finite]
    adjustedPValue = adjustedPValue[finite]
    positivePValues = adjustedPValue[adjustedPValue > 0]
    pValueFloor = (
        max(float(positivePValues.min()) / 10, np.finfo(float).tiny) if positivePValues.size else np.finfo(float).tiny
    )
    minusLog10Padj = -np.log10(np.clip(adjustedPValue, pValueFloor, None))
    significant = (adjustedPValue < padjThreshold) & (np.abs(log2FoldChange) >= lfcThreshold)
    gain = significant & (log2FoldChange > 0)
    depletion = significant & (log2FoldChange < 0)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.scatter(
        log2FoldChange[~significant],
        minusLog10Padj[~significant],
        s=10,
        color="0.78",
        alpha=0.55,
        linewidths=0,
        label="Below thresholds",
        rasterized=True,
    )
    ax.scatter(
        log2FoldChange[depletion],
        minusLog10Padj[depletion],
        s=13,
        color="#C44E52",
        alpha=0.75,
        linewidths=0,
        label="Significant depletion",
        rasterized=True,
    )
    ax.scatter(
        log2FoldChange[gain],
        minusLog10Padj[gain],
        s=13,
        color="#4C72B0",
        alpha=0.75,
        linewidths=0,
        label="Significant gain",
        rasterized=True,
    )
    ax.axvline(-lfcThreshold, color="0.35", lw=0.8, ls="--")
    ax.axvline(lfcThreshold, color="0.35", lw=0.8, ls="--")
    ax.axhline(-np.log10(padjThreshold), color="0.35", lw=0.8, ls="--")

    geneSymbols = plot["geneSymbol"].fillna("").astype(str).to_numpy()
    labelPoints: list[tuple[str, float, float]] = []
    for geneSymbol in labelGenes:
        matches = np.flatnonzero((geneSymbols == geneSymbol) & significant)
        if not matches.size:
            continue
        pointIndex = int(matches[0])
        labelPoints.append((geneSymbol, float(log2FoldChange[pointIndex]), float(minusLog10Padj[pointIndex])))

    xRange = max(float(np.ptp(log2FoldChange)), 1.0)
    yRange = max(float(np.ptp(minusLog10Padj)), 1.0)
    for side in (-1, 1):
        previousTextY = -np.inf
        sidePoints = sorted(
            (point for point in labelPoints if np.sign(point[1]) == side),
            key=lambda point: point[2],
        )
        for geneSymbol, x, y in sidePoints:
            textY = max(y + 0.02 * yRange, previousTextY + 0.05 * yRange)
            textX = x + side * 0.02 * xRange
            ax.annotate(
                geneSymbol,
                xy=(x, y),
                xytext=(textX, textY),
                textcoords="data",
                ha="left" if side > 0 else "right",
                va="bottom",
                fontsize=7,
                arrowprops={"arrowstyle": "-", "color": "0.45", "lw": 0.5},
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.8},
                clip_on=False,
            )
            previousTextY = textY

    ax.set_xlabel("disease vs control log2FC")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title, fontsize=10)
    ax.margins(x=0.08, y=0.15)
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def auto_label_genes(
    contrastResults: pd.DataFrame,
    *,
    padj: float,
    lfc: float,
    maxLabels: int = 8,
    preferred: list[str] | None = None,
) -> list[str]:
    preferred = preferred or []
    frame = contrastResults.copy()
    frame["padj"] = pd.to_numeric(frame["padj"], errors="coerce")
    frame["log2FoldChange"] = pd.to_numeric(frame["log2FoldChange"], errors="coerce")
    significant = frame[frame["padj"].notna() & (frame["padj"] <= padj) & (frame["log2FoldChange"].abs() >= lfc)].copy()
    labels: list[str] = []
    gene_symbol_set = set(significant["geneSymbol"].fillna("").astype(str))
    for gene in preferred:
        if gene in gene_symbol_set and gene not in labels:
            labels.append(gene)
    ranked = significant.assign(absLfc=lambda df: df["log2FoldChange"].abs()).sort_values(
        ["padj", "absLfc"], ascending=[True, False]
    )
    for gene in ranked["geneSymbol"].fillna(ranked["gene"]).astype(str):
        if gene and gene not in labels:
            labels.append(gene)
        if len(labels) >= maxLabels:
            break
    return labels


def plot_study_evidence_panel(
    pdata: AnnData,
    *,
    gene: str,
    area: str,
    cluster: str,
    path: Path,
    title: str,
    clusterKey: str = "leiden_atlas",
    studyKey: str = "study_accession",
) -> None:
    selected = same_study_case_control_profiles(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=clusterKey,
        studyKey=studyKey,
    )
    if not bool(selected.any()) or gene not in pdata.var_names:
        return
    sub = pdata[selected.to_numpy(), [gene]].copy()
    props = np.asarray(sub.layers["psbulk_props"]).ravel()
    obs = pd.DataFrame(sub.obs).copy()
    obs["detection"] = props
    obs["arm"] = np.where(obs["diseased"].astype("boolean").eq(True).fillna(False), "case", "control")
    obs[studyKey] = as_string_series(pd.Series(obs[studyKey]))
    summary = obs.groupby([studyKey, "arm"], observed=True)["detection"].mean().reset_index(name="meanDetection")
    if summary.empty:
        return
    studies = sorted(summary[studyKey].unique())
    fig, ax = plt.subplots(figsize=(max(4.5, 0.55 * len(studies) + 2.0), 3.5))
    width = 0.35
    x = np.arange(len(studies))
    for offset, arm, color in ((-width / 2, "control", "0.65"), (width / 2, "case", "C0")):
        values = [
            float(summary.loc[(summary[studyKey] == study) & (summary["arm"] == arm), "meanDetection"].mean())
            if ((summary[studyKey] == study) & (summary["arm"] == arm)).any()
            else np.nan
            for study in studies
        ]
        ax.bar(x + offset, values, width=width, label=arm, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(studies, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("mean detection")
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_review_figures(
    *,
    figDir: Path,
    pdata: AnnData,
    shortlist: pd.DataFrame,
    extended: pd.DataFrame,
    thresholds: pd.DataFrame,
    restricted: pd.DataFrame,
    deHits: pd.DataFrame,
    deResults: pd.DataFrame,
    unexpected: pd.DataFrame,
    clusterKey: str,
    studyKey: str,
    padj: float,
    lfc: float,
    maxVolcanoPlots: int,
    maxEvidencePanels: int,
) -> None:
    figDir.mkdir(parents=True, exist_ok=True)
    review = pd.concat([shortlist, extended], ignore_index=True, sort=False)
    plot_score_distributions(review, thresholds, figDir / "score_distributions.png")

    if not restricted.empty:
        top_restricted = restricted.sort_values(["tau", "detectionDifference"], ascending=[False, False]).head(30)
        gene_list = top_restricted["gene"].astype(str).tolist()
        cluster_order = sorted(
            top_restricted["topCluster"].astype(str).unique(),
            key=lambda x: int(x) if x.isdigit() else x,
        )
        props = np.asarray(pdata[:, gene_list].layers["psbulk_props"], dtype=float)
        labels = top_restricted["geneSymbol"].fillna(top_restricted["gene"]).astype(str).tolist()
        det_mat = pd.DataFrame(index=pd.Index(labels), columns=pd.Index(cluster_order), dtype=float)
        weights = study_balanced_weights(pd.Series(pdata.obs[studyKey])).to_numpy(dtype=np.float64)
        for j, _gene in enumerate(gene_list):
            for cluster in cluster_order:
                mask = as_string_series(pd.Series(pdata.obs[clusterKey])).to_numpy() == cluster
                if not mask.any():
                    det_mat.iloc[j, det_mat.columns.get_loc(cluster)] = np.nan
                    continue
                w = weights[mask]
                det_mat.iloc[j, det_mat.columns.get_loc(cluster)] = float(np.average(props[mask, j], weights=w))
        plot_heatmap(
            det_mat,
            figDir / "heatmap_restricted_genes.png",
            "Restricted gene detection across home clusters",
            cmap="magma",
        )

    if not deHits.empty:
        strong = deHits.sort_values("padj").groupby(["gene", "diseaseArea"], observed=True).head(1)
        strong = strong.sort_values("padj").head(40)
        labels = strong["geneSymbol"].fillna(strong["gene"]).astype(str) + " | c" + strong["cluster"].astype(str)
        pivot = (
            pd.DataFrame(
                {
                    "label": labels.to_numpy(),
                    "diseaseArea": strong["diseaseArea"].astype(str).to_numpy(),
                    "log2FoldChange": strong["log2FoldChange"].astype(float).to_numpy(),
                }
            )
            .pivot_table(index="label", columns="diseaseArea", values="log2FoldChange", aggfunc="first")
            .fillna(0.0)
        )
        plot_heatmap(pivot, figDir / "heatmap_disease_effects.png", "Disease vs control log2FC", cmap="coolwarm")

    if not unexpected.empty:
        plot_unexpected_dotplot(unexpected.head(40), figDir / "dotplot_unexpected_expression.png")

    if not shortlist.empty and not deResults.empty and maxVolcanoPlots > 0:
        contrast_keys = (
            shortlist.dropna(subset=["cluster", "diseaseArea"])
            .loc[:, ["cluster", "diseaseArea"]]
            .astype(str)
            .drop_duplicates()
            .head(maxVolcanoPlots)
        )
        for row in contrast_keys.to_dict(orient="records"):
            cluster = str(row["cluster"])
            disease_area = str(row["diseaseArea"])
            mask = deResults["cluster"].astype(str).eq(cluster) & deResults["diseaseArea"].astype(str).eq(disease_area)
            frame = deResults.loc[mask].copy()
            preferred = (
                shortlist.loc[
                    shortlist["cluster"].astype(str).eq(cluster)
                    & shortlist["diseaseArea"].astype(str).eq(disease_area),
                    "geneSymbol",
                ]
                .fillna("")
                .astype(str)
                .tolist()
            )
            labels = auto_label_genes(frame, padj=padj, lfc=lfc, preferred=preferred)
            safe_area = "".join(ch if ch.isalnum() else "_" for ch in disease_area).strip("_").lower()
            plot_volcano(
                frame,
                figDir / f"volcano_{safe_area}_cluster_{cluster}.png",
                f"{disease_area} cluster {cluster}",
                labels,
                padjThreshold=padj,
                lfcThreshold=lfc,
            )

    if not shortlist.empty and maxEvidencePanels > 0:
        panel_dir = figDir / "evidence_panels"
        panel_dir.mkdir(parents=True, exist_ok=True)
        panels = shortlist.dropna(subset=["gene", "cluster", "diseaseArea"]).head(maxEvidencePanels)
    for row in panels.to_dict(orient="records"):
        symbol = str(row.get("geneSymbol") or row["gene"])
        safe_area = "".join(ch if ch.isalnum() else "_" for ch in str(row["diseaseArea"])).strip("_").lower()
        plot_study_evidence_panel(
            pdata,
            gene=str(row["gene"]),
            area=str(row["diseaseArea"]),
            cluster=str(row["cluster"]),
            path=panel_dir / f"{symbol}_c{row['cluster']}_{safe_area}.png",
            title=f"{symbol} | cluster {row['cluster']} | {row['diseaseArea']}",
            clusterKey=clusterKey,
            studyKey=studyKey,
        )
