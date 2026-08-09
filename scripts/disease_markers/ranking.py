"""Adaptive noteworthy-gene ranking with review-budget diversity constraints."""

from __future__ import annotations

import numpy as np
import pandas as pd

EVIDENCE_CLASSES = (
    "clusterRestricted",
    "replicatedDiseaseGain",
    "replicatedDiseaseDepletion",
    "sharedDiseaseProgram",
    "oppositeDiseaseEffect",
    "unexpectedExpression",
)


def empirical_percentile_scores(values: pd.Series) -> pd.Series:
    """Rank-based percentile in [0, 1]. Ties share the average rank. All-NaN -> 0."""
    numeric = pd.Series(pd.to_numeric(values, errors="coerce"), index=values.index, dtype="float64")
    if int(numeric.notna().sum()) == 0:
        return pd.Series(np.zeros(len(values), dtype=float), index=values.index)
    ranks = numeric.rank(method="average", na_option="keep")
    valid = ranks.notna()
    if int(valid.sum()) <= 1:
        out = pd.Series(np.zeros(len(values), dtype=float), index=values.index)
        out.loc[valid] = 1.0
        return out
    scaled = (ranks - 1.0) / (float(valid.sum()) - 1.0)
    return scaled.fillna(0.0)


def _score_frame(frame: pd.DataFrame, columns: list[str], weights: dict[str, float]) -> pd.Series:
    score = pd.Series(0.0, index=frame.index, dtype=float)
    weight_sum = 0.0
    for column in columns:
        if column not in frame.columns:
            continue
        weight = float(weights.get(column, 1.0))
        score = score + weight * empirical_percentile_scores(frame[column])
        weight_sum += weight
    if weight_sum <= 0:
        return score
    return score / weight_sum


def _ensure_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column not in out.columns:
            out[column] = pd.NA
    return out


def build_evidence_pools(
    *,
    restricted: pd.DataFrame,
    deHits: pd.DataFrame,
    sharedGenes: pd.DataFrame,
    geneClass: pd.DataFrame,
    unexpected: pd.DataFrame,
    padj: float,
    lfc: float,
    minDetectionDelta: float,
    minTau: float,
    minTargetDetection: float,
    maxBackgroundDetection: float,
    minStudiesForSpecificity: int,
) -> dict[str, pd.DataFrame]:
    pools: dict[str, pd.DataFrame] = {name: pd.DataFrame() for name in EVIDENCE_CLASSES}

    if not restricted.empty:
        keep = restricted[
            (restricted["tau"] >= minTau)
            & (restricted["meanDetectionTop"] >= minTargetDetection)
            & (restricted["maxDetectionBackground"] <= maxBackgroundDetection)
            & (restricted["nStudiesAgreeTop"] >= minStudiesForSpecificity)
            & (restricted["nStudiesScored"] >= minStudiesForSpecificity)
        ].copy()
        if not keep.empty:
            keep["proposalClass"] = "clusterRestricted"
            keep["cluster"] = keep["topCluster"].astype(str)
            keep["diseaseArea"] = pd.NA
            keep["evidenceScore"] = _score_frame(
                keep,
                ["tau", "detectionDifference", "nStudiesAgreeTop", "meanDetectionTop"],
                {
                    "tau": 1.2,
                    "detectionDifference": 1.0,
                    "nStudiesAgreeTop": 1.0,
                    "meanDetectionTop": 0.6,
                },
            )
            pools["clusterRestricted"] = keep

    if not deHits.empty:
        gain = deHits[
            (deHits["log2FoldChange"] >= lfc)
            & (deHits["padj"] <= padj)
            & (deHits["detectionDelta"] >= minDetectionDelta)
        ].copy()
        if not gain.empty:
            gain = gain.sort_values("padj").groupby(["gene", "cluster", "diseaseArea"], observed=True).head(1)
            gain["proposalClass"] = "replicatedDiseaseGain"
            gain["negLog10Padj"] = -np.log10(np.clip(gain["padj"].astype(float), 1e-300, None))
            gain["absLog2FoldChange"] = gain["log2FoldChange"].abs()
            gain["absDetectionDelta"] = gain["detectionDelta"].abs()
            gain["evidenceScore"] = _score_frame(
                gain,
                ["negLog10Padj", "absLog2FoldChange", "absDetectionDelta", "nStudies"],
                {
                    "negLog10Padj": 1.2,
                    "absLog2FoldChange": 1.0,
                    "absDetectionDelta": 1.0,
                    "nStudies": 0.8,
                },
            )
            pools["replicatedDiseaseGain"] = gain

        depletion = deHits[
            (deHits["log2FoldChange"] <= -lfc)
            & (deHits["padj"] <= padj)
            & (deHits["detectionDelta"] <= -minDetectionDelta)
        ].copy()
        if not depletion.empty:
            depletion = depletion.sort_values("padj").groupby(["gene", "cluster", "diseaseArea"], observed=True).head(1)
            depletion["proposalClass"] = "replicatedDiseaseDepletion"
            depletion["negLog10Padj"] = -np.log10(np.clip(depletion["padj"].astype(float), 1e-300, None))
            depletion["absLog2FoldChange"] = depletion["log2FoldChange"].abs()
            depletion["absDetectionDelta"] = depletion["detectionDelta"].abs()
            depletion["evidenceScore"] = _score_frame(
                depletion,
                ["negLog10Padj", "absLog2FoldChange", "absDetectionDelta", "nStudies"],
                {
                    "negLog10Padj": 1.2,
                    "absLog2FoldChange": 1.0,
                    "absDetectionDelta": 1.0,
                    "nStudies": 0.8,
                },
            )
            pools["replicatedDiseaseDepletion"] = depletion

    if not sharedGenes.empty:
        shared = sharedGenes.copy()
        shared["proposalClass"] = "sharedDiseaseProgram"
        shared["cluster"] = pd.NA
        shared["diseaseArea"] = pd.NA
        shared["absMeanLog2FoldChange"] = shared["meanLog2FoldChange"].abs()
        shared["evidenceScore"] = _score_frame(
            shared,
            ["nDiseaseAreas", "absMeanLog2FoldChange", "nClusters"],
            {"nDiseaseAreas": 1.4, "absMeanLog2FoldChange": 1.0, "nClusters": 0.6},
        )
        pools["sharedDiseaseProgram"] = shared

    if not geneClass.empty:
        opposite = geneClass[geneClass["evidenceClass"].eq("oppositeDirection")].copy()
        if not opposite.empty:
            opposite["proposalClass"] = "oppositeDiseaseEffect"
            opposite["diseaseArea"] = opposite.get("diseaseAreas", pd.NA)
            opposite["negLog10Padj"] = -np.log10(np.clip(opposite["minPadj"].astype(float), 1e-300, None))
            opposite["absMeanLog2FoldChange"] = opposite["meanLog2FoldChange"].abs()
            opposite["evidenceScore"] = _score_frame(
                opposite,
                ["negLog10Padj", "absMeanLog2FoldChange", "nDiseaseAreasHit"],
                {"negLog10Padj": 1.0, "absMeanLog2FoldChange": 1.0, "nDiseaseAreasHit": 1.2},
            )
            pools["oppositeDiseaseEffect"] = opposite

    if not unexpected.empty:
        unexp = unexpected[unexpected["detectionDelta"].abs() >= minDetectionDelta].copy()
        if not unexp.empty:
            unexp["proposalClass"] = "unexpectedExpression"
            unexp["negLog10Padj"] = -np.log10(np.clip(unexp["padj"].astype(float), 1e-300, None))
            unexp["absLog2FoldChange"] = unexp["log2FoldChange"].abs()
            unexp["absDetectionDelta"] = unexp["detectionDelta"].abs()
            unexp["resolvedBonus"] = (
                unexp.get("interpretationStatus", pd.Series("unresolved", index=unexp.index))
                .astype(str)
                .eq("resolved")
                .astype(float)
            )
            unexp["evidenceScore"] = _score_frame(
                unexp,
                ["negLog10Padj", "absLog2FoldChange", "absDetectionDelta", "resolvedBonus"],
                {
                    "negLog10Padj": 1.0,
                    "absLog2FoldChange": 0.8,
                    "absDetectionDelta": 1.2,
                    "resolvedBonus": 0.5,
                },
            )
            pools["unexpectedExpression"] = unexp

    return pools


def _event_key(row: pd.Series) -> str:
    gene = str(row.get("geneSymbol") or row.get("gene") or "")
    cluster = str(row.get("cluster") or "")
    area = str(row.get("diseaseArea") or "")
    return f"{gene}|{cluster}|{area}"


def merge_duplicate_evidence(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return candidates.copy()
    frame = candidates.copy()
    frame["eventKey"] = frame.apply(_event_key, axis=1)
    rows: list[dict[str, object]] = []
    for _, group in frame.groupby("eventKey", sort=False):
        ordered = group.sort_values("evidenceScore", ascending=False)
        primary = ordered.iloc[0].to_dict()
        classes = sorted({str(value) for value in ordered["proposalClass"].dropna().astype(str)})
        primary["proposalClass"] = classes[0]
        primary["allEvidenceClasses"] = ",".join(classes)
        primary["nEvidenceClasses"] = len(classes)
        primary["evidenceScore"] = float(ordered["evidenceScore"].max())
        rows.append(primary)
    out = pd.DataFrame(rows)
    return out.drop(columns=["eventKey"], errors="ignore")


def _passes_diversity(
    row: pd.Series,
    selected: pd.DataFrame,
    *,
    maxPerGene: int,
    maxPerCluster: int,
    maxPerDiseaseArea: int,
    maxPerClass: int,
) -> bool:
    if selected.empty:
        return True
    gene = str(row.get("geneSymbol") or row.get("gene") or "")
    cluster = str(row.get("cluster") or "")
    area = str(row.get("diseaseArea") or "")
    proposal_class = str(row.get("proposalClass") or "")
    gene_count = int((selected["geneSymbol"].fillna(selected["gene"]).astype(str) == gene).sum()) if gene else 0
    if gene and gene_count >= maxPerGene:
        return False
    if cluster and cluster != "nan" and int((selected["cluster"].astype(str) == cluster).sum()) >= maxPerCluster:
        return False
    if (
        area
        and area not in {"", "nan", "<NA>"}
        and int((selected["diseaseArea"].astype(str) == area).sum()) >= maxPerDiseaseArea
    ):
        return False
    if proposal_class and int((selected["proposalClass"].astype(str) == proposal_class).sum()) >= maxPerClass:
        return False
    return True


def select_review_queue(
    pools: dict[str, pd.DataFrame],
    *,
    primaryBudget: int = 20,
    extendedBudget: int = 60,
    maxPerClassPrimary: int = 6,
    maxPerClassExtended: int = 15,
    maxPerGene: int = 2,
    maxPerCluster: int = 4,
    maxPerDiseaseArea: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return primary shortlist, extended queue, and per-class cutoff table.

    Cutoffs are derived from score order statistics needed to fill the budget.
    Rows that fail validity floors are never added merely to fill a quota.
    """
    frames = [frame for frame in pools.values() if frame is not None and not frame.empty]
    if not frames:
        empty = pd.DataFrame()
        return empty, empty, empty

    combined = merge_duplicate_evidence(pd.concat(frames, ignore_index=True, sort=False))
    combined = _ensure_columns(
        combined,
        [
            "gene",
            "geneSymbol",
            "cluster",
            "diseaseArea",
            "proposalClass",
            "evidenceScore",
            "allEvidenceClasses",
            "interpretationStatus",
            "interpretedCellType",
        ],
    )
    combined = combined.sort_values(["evidenceScore", "proposalClass"], ascending=[False, True]).reset_index(drop=True)

    threshold_rows: list[dict[str, object]] = []
    for proposal_class, group in combined.groupby("proposalClass", observed=True):
        scores = group["evidenceScore"].sort_values(ascending=False).to_numpy(dtype=float)
        primary_cut = float(scores[min(len(scores), maxPerClassPrimary) - 1]) if len(scores) else float("nan")
        extended_cut = float(scores[min(len(scores), maxPerClassExtended) - 1]) if len(scores) else float("nan")
        threshold_rows.append(
            {
                "proposalClass": str(proposal_class),
                "nCandidates": int(len(scores)),
                "primaryCutoff": primary_cut,
                "extendedCutoff": extended_cut,
                "scoreMedian": float(np.median(scores)) if len(scores) else float("nan"),
                "scoreP90": float(np.quantile(scores, 0.9)) if len(scores) else float("nan"),
            }
        )
    thresholds = pd.DataFrame(threshold_rows)

    primary_rows: list[pd.Series] = []
    extended_rows: list[pd.Series] = []
    primary = pd.DataFrame(columns=combined.columns)
    extended = pd.DataFrame(columns=combined.columns)

    for _, row in combined.iterrows():
        if len(primary_rows) < primaryBudget and _passes_diversity(
            row,
            primary,
            maxPerGene=maxPerGene,
            maxPerCluster=maxPerCluster,
            maxPerDiseaseArea=maxPerDiseaseArea,
            maxPerClass=maxPerClassPrimary,
        ):
            primary_rows.append(row)
            primary = pd.DataFrame(primary_rows)
            continue
        if len(primary_rows) + len(extended_rows) < extendedBudget and _passes_diversity(
            row,
            pd.concat([primary, extended], ignore_index=True) if len(extended_rows) else primary,
            maxPerGene=maxPerGene,
            maxPerCluster=maxPerCluster,
            maxPerDiseaseArea=maxPerDiseaseArea,
            maxPerClass=maxPerClassExtended,
        ):
            extended_rows.append(row)
            extended = pd.DataFrame(extended_rows)

    if primary_rows:
        primary = pd.DataFrame(primary_rows).copy()
        primary["reviewTier"] = "primary"
        primary["reviewRank"] = np.arange(1, len(primary) + 1)
        primary["reviewPriority"] = np.where(primary["reviewRank"] <= max(5, primaryBudget // 4), "high", "medium")
    else:
        primary = pd.DataFrame()

    if extended_rows:
        extended = pd.DataFrame(extended_rows).copy()
        extended["reviewTier"] = "extended"
        start = len(primary) + 1
        extended["reviewRank"] = np.arange(start, start + len(extended))
        extended["reviewPriority"] = "medium"
    else:
        extended = pd.DataFrame()

    return primary, extended, thresholds
