"""Analyze-stage orchestration for atlas disease DE and noteworthy ranking."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from h5ad_concat.reference import load_gene_reference

from disease_markers.candidates import same_study_contrast_support
from disease_markers.config import AtlasDeAnalysisConfig
from disease_markers.labels import OTHER_AREA
from disease_markers.memory import format_bytes, log_memory, snapshot_memory
from disease_markers.plots import write_review_figures
from disease_markers.ranking import build_evidence_pools, select_review_queue
from disease_markers.specificity import (
    attach_gene_symbols,
    cluster_interpretation_table,
    control_home_cluster,
    gene_specificity_table,
    study_direction_agreement,
)
from disease_markers.validation import (
    disease_vs_control_deseq2,
    filter_two_sided_de,
    same_study_case_control_profiles,
    shared_direction_genes,
)

log = logging.getLogger(__name__)


def _contrast_key(cluster: str, area: str) -> str:
    safe_area = "".join(ch if ch.isalnum() else "_" for ch in area).strip("_")
    return f"cluster_{cluster}__{safe_area}"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n")


def _load_or_run_contrast(
    pdata: AnnData,
    *,
    cluster: str,
    area: str,
    cfg: AtlasDeAnalysisConfig,
    fingerprintSha: str,
) -> pd.DataFrame:
    cfg.deCheckpointsDir.mkdir(parents=True, exist_ok=True)
    key = _contrast_key(cluster, area)
    result_path = cfg.deCheckpointsDir / f"{key}.parquet"
    meta_path = cfg.deCheckpointsDir / f"{key}.json"
    if result_path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if (
            meta.get("fingerprintSha") == fingerprintSha
            and meta.get("cluster") == cluster
            and meta.get("diseaseArea") == area
        ):
            log.info("Reusing DE checkpoint %s", result_path.name)
            return pd.read_parquet(result_path)

    full = disease_vs_control_deseq2(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=cfg.clusterKey,
        studyKey=cfg.studyKey,
    )
    if full.empty:
        full.to_parquet(result_path, index=False)
    else:
        full.to_parquet(result_path, index=False)
    _write_json(
        meta_path,
        {
            "cluster": cluster,
            "diseaseArea": area,
            "fingerprintSha": fingerprintSha,
            "nRows": int(len(full)),
            "ran": not full.empty,
        },
    )
    return full


def _attach_detection_metrics(
    full: pd.DataFrame,
    pdata: AnnData,
    *,
    area: str,
    cluster: str,
    cfg: AtlasDeAnalysisConfig,
) -> pd.DataFrame:
    if full.empty:
        return full
    gene_idx = {gene: i for i, gene in enumerate(pdata.var_names)}
    selected = same_study_case_control_profiles(
        pdata,
        area=area,
        cluster=cluster,
        clusterKey=cfg.clusterKey,
        studyKey=cfg.studyKey,
    )
    sub = pdata[selected.to_numpy()]
    diseased = sub.obs["diseased"].astype("boolean").eq(True).fillna(False).to_numpy()
    props = np.asarray(sub.layers["psbulk_props"], dtype=float)
    case_det = props[diseased].mean(axis=0) if diseased.any() else np.full(props.shape[1], np.nan)
    control_det = props[~diseased].mean(axis=0) if (~diseased).any() else np.full(props.shape[1], np.nan)
    out = full.copy()
    positions = [gene_idx[str(gene)] for gene in out["gene"].astype(str)]
    out["caseDetection"] = case_det[positions]
    out["controlDetection"] = control_det[positions]
    out["detectionDelta"] = out["caseDetection"] - out["controlDetection"]
    return out


def _build_unexpected(
    de_hits: pd.DataFrame,
    pdata: AnnData,
    cluster_interp: pd.DataFrame,
    cfg: AtlasDeAnalysisConfig,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if de_hits.empty:
        return pd.DataFrame()
    home_cache: dict[str, tuple[str | None, float]] = {}
    for row in de_hits.to_dict(orient="records"):
        detection_delta = float(row["detectionDelta"])
        if not np.isfinite(detection_delta) or abs(detection_delta) < cfg.minDetectionDelta:
            continue
        gene = str(row["gene"])
        if gene not in home_cache:
            home_cache[gene] = control_home_cluster(
                pdata,
                gene=gene,
                studyKey=cfg.studyKey,
                clusterKey=cfg.clusterKey,
            )
        home_cluster, home_detection = home_cache[gene]
        if home_cluster is None or home_cluster == str(row["cluster"]):
            continue
        home_meta = cluster_interp.loc[cluster_interp["cluster"].astype(str) == home_cluster]
        recipient_meta = cluster_interp.loc[cluster_interp["cluster"].astype(str) == str(row["cluster"])]
        rows.append(
            {
                "gene": gene,
                "geneSymbol": row.get("geneSymbol"),
                "biotype": row.get("biotype"),
                "cluster": str(row["cluster"]),
                "diseaseArea": str(row["diseaseArea"]),
                "log2FoldChange": float(row["log2FoldChange"]),
                "padj": float(row["padj"]),
                "caseDetection": float(row["caseDetection"]),
                "controlDetection": float(row["controlDetection"]),
                "detectionDelta": detection_delta,
                "nStudies": int(row["nStudies"]),
                "interpretedCellType": row.get("interpretedCellType"),
                "interpretationStatus": row.get("interpretationStatus"),
                "topSourceLabel": row.get("topSourceLabel"),
                "topSourceLabelFraction": row.get("topSourceLabelFraction"),
                "homeCluster": home_cluster,
                "homeClusterControlDetection": home_detection,
                "homeInterpretedCellType": None if home_meta.empty else home_meta["interpretedCellType"].iloc[0],
                "homeTopLabel": None if home_meta.empty else home_meta["topLabel"].iloc[0],
                "homeTopLabelFraction": None if home_meta.empty else float(home_meta["topLabelFraction"].iloc[0]),
                "recipientLabelEntropy": None
                if recipient_meta.empty
                else float(recipient_meta["labelEntropy"].iloc[0]),
                "direction": "gain" if float(row["log2FoldChange"]) > 0 else "loss",
            }
        )
    unexpected = pd.DataFrame(rows)
    if unexpected.empty:
        return unexpected
    return unexpected.sort_values(["direction", "detectionDelta", "padj"], ascending=[True, False, True]).reset_index(
        drop=True
    )


def _annotate_study_agreement(
    shortlist: pd.DataFrame,
    pdata: AnnData,
    cfg: AtlasDeAnalysisConfig,
) -> pd.DataFrame:
    if shortlist.empty:
        return shortlist
    rows: list[dict[str, object]] = []
    for payload in shortlist.to_dict(orient="records"):
        area = payload.get("diseaseArea")
        cluster = payload.get("cluster")
        gene = payload.get("gene")
        if pd.isna(area) or pd.isna(cluster) or pd.isna(gene) or not gene:
            payload["studyDirectionAgreement"] = pd.NA
            payload["nStudiesAgree"] = pd.NA
            payload["nStudiesScoredAgreement"] = pd.NA
        else:
            agreement = study_direction_agreement(
                pdata,
                gene=str(gene),
                area=str(area),
                cluster=str(cluster),
                clusterKey=cfg.clusterKey,
                studyKey=cfg.studyKey,
            )
            payload["studyDirectionAgreement"] = agreement["studyDirectionAgreement"]
            payload["nStudiesAgree"] = agreement["nStudiesAgree"]
            payload["nStudiesScoredAgreement"] = agreement["nStudiesScored"]
        rows.append(payload)
    return pd.DataFrame(rows)


def analyze_from_pseudobulk(
    cfg: AtlasDeAnalysisConfig,
    pdata: AnnData,
    *,
    clusterObs: pd.DataFrame | None = None,
    fingerprintSha: str = "",
) -> dict[str, Any]:
    """Run specificity, DE, ranking, and review figures from a pseudobulk object."""
    cfg.outputDir.mkdir(parents=True, exist_ok=True)
    cfg.figuresDir.mkdir(parents=True, exist_ok=True)
    log_memory("analyze-start", logger=log)

    gene_reference = load_gene_reference(cfg.geneInfoPath)
    gene_map = gene_reference.var.reset_index(names="ensembl_id").rename(
        columns={"gene_symbol": "geneSymbol", "biotype": "biotype"}
    )
    pdata = pdata.copy()
    gene_annotations = gene_map.set_index("ensembl_id").reindex(pdata.var_names)
    pdata.var["geneSymbol"] = gene_annotations["geneSymbol"].to_numpy()
    pdata.var["biotype"] = gene_annotations["biotype"].to_numpy()

    cluster_interp_path = cfg.outputDir / "cluster_interpretation.csv"
    if cluster_interp_path.exists():
        cluster_interp = pd.read_csv(cluster_interp_path)
        cluster_interp["cluster"] = cluster_interp["cluster"].astype(str)
    else:
        obs_for_interp = clusterObs if clusterObs is not None else pdata.obs
        cluster_interp = cluster_interpretation_table(
            pd.DataFrame(obs_for_interp),
            clusterKey=cfg.clusterKey,
            labelKey=cfg.labelKey,
            ontologyKey=cfg.ontologyKey,
            studyKey=cfg.studyKey,
            sampleKey=cfg.sampleKey,
            highPurity=cfg.highPurity,
            resolvedMinStudies=cfg.resolvedMinStudies,
        )
        cluster_interp.to_csv(cluster_interp_path, index=False)

    log.info("Computing tau specificity across Leiden clusters...")
    specificity = gene_specificity_table(
        pdata,
        clusterKey=cfg.clusterKey,
        studyKey=cfg.studyKey,
        minStudies=cfg.minStudiesForSpecificity,
        minProfilesForGene=cfg.minProfilesForGene,
        minTotalCounts=cfg.minTotalCountsForGene,
        geneChunkSize=cfg.geneChunkSize,
    )
    specificity = attach_gene_symbols(specificity, gene_map)
    specificity = specificity.merge(
        cluster_interp[
            [
                "cluster",
                "interpretedCellType",
                "interpretationStatus",
                "topLabel",
                "topLabelFraction",
                "labelEntropy",
            ]
        ].rename(columns={"cluster": "topCluster"}),
        on="topCluster",
        how="left",
    )
    specificity.to_parquet(cfg.outputDir / "gene_specificity_tau.parquet", index=False)

    restricted = specificity[
        (specificity["tau"] >= cfg.minTau)
        & (specificity["meanDetectionTop"] >= cfg.minTargetDetection)
        & (specificity["maxDetectionBackground"] <= cfg.maxBackgroundDetection)
        & (specificity["nStudiesAgreeTop"] >= cfg.minStudiesForSpecificity)
        & (specificity["nStudiesScored"] >= cfg.minStudiesForSpecificity)
    ].copy()
    restricted = restricted.sort_values(["tau", "detectionDifference"], ascending=[False, False])
    restricted.to_csv(cfg.outputDir / "restricted_genes.csv", index=False)

    contrast_path = cfg.outputDir / "same_study_contrast_support.csv"
    if contrast_path.exists():
        contrast_support = pd.read_csv(contrast_path)
        contrast_support["cluster"] = contrast_support["cluster"].astype(str)
        if "interpretedCellType" not in contrast_support.columns:
            contrast_support = contrast_support.merge(
                cluster_interp[
                    ["cluster", "interpretedCellType", "interpretationStatus", "topLabel", "topLabelFraction"]
                ],
                on="cluster",
                how="left",
            )
    else:
        contrast_support = same_study_contrast_support(
            pd.DataFrame(pdata.obs),
            clusterKey=cfg.clusterKey,
            sampleKey=cfg.sampleKey,
            studyKey=cfg.studyKey,
            minCellsPerProfile=1,
        )
        contrast_support = contrast_support.merge(
            cluster_interp[["cluster", "interpretedCellType", "interpretationStatus", "topLabel", "topLabelFraction"]],
            on="cluster",
            how="left",
        )
        contrast_support.to_csv(contrast_path, index=False)

    eligible_contrasts = contrast_support[
        contrast_support["eligibleForContrast"].astype(bool)
        & (contrast_support["nOverlapStudies"] >= cfg.minOverlapStudies)
        & contrast_support["diseaseArea"].astype(str).ne(OTHER_AREA)
    ].copy()
    log.info("Eligible disease contrasts: %s", len(eligible_contrasts))

    matrix = np.asarray(pdata.X)
    gene_profile_support = (matrix > 0).sum(axis=0)
    gene_total_counts = matrix.sum(axis=0)
    keep_var = (gene_profile_support >= cfg.minProfilesForGene) & (gene_total_counts >= cfg.minTotalCountsForGene)
    pdata_de = pdata[:, keep_var].copy()
    log.info("Genes retained for DE: %s / %s", f"{pdata_de.n_vars:,}", f"{pdata.n_vars:,}")

    summary_rows: list[dict[str, object]] = []
    result_frames: list[pd.DataFrame] = []
    hit_frames: list[pd.DataFrame] = []
    for row in eligible_contrasts.to_dict(orient="records"):
        cluster = str(row["cluster"])
        area = str(row["diseaseArea"])
        log.info("DESeq2 cluster=%s area=%s studies=%s", cluster, area, row["nOverlapStudies"])
        full = _load_or_run_contrast(
            pdata_de,
            cluster=cluster,
            area=area,
            cfg=cfg,
            fingerprintSha=fingerprintSha,
        )
        if full.empty:
            summary_rows.append(
                {
                    "cluster": cluster,
                    "diseaseArea": area,
                    "interpretedCellType": row.get("interpretedCellType"),
                    "interpretationStatus": row.get("interpretationStatus"),
                    "nOverlapStudies": int(row["nOverlapStudies"]),
                    "nCaseProfiles": int(row["nCaseProfiles"]),
                    "nControlProfiles": int(row["nControlProfiles"]),
                    "nHits": 0,
                    "ran": False,
                }
            )
            continue
        full = _attach_detection_metrics(full, pdata_de, area=area, cluster=cluster, cfg=cfg)
        full["interpretedCellType"] = row.get("interpretedCellType")
        full["interpretationStatus"] = row.get("interpretationStatus")
        full["topSourceLabel"] = row.get("topLabel")
        full["topSourceLabelFraction"] = row.get("topLabelFraction")
        full = attach_gene_symbols(full, gene_map)
        hits = filter_two_sided_de(full, padj=cfg.padj, lfc=cfg.lfc)
        summary_rows.append(
            {
                "cluster": cluster,
                "diseaseArea": area,
                "interpretedCellType": row.get("interpretedCellType"),
                "interpretationStatus": row.get("interpretationStatus"),
                "nOverlapStudies": int(row["nOverlapStudies"]),
                "nCaseProfiles": int(row["nCaseProfiles"]),
                "nControlProfiles": int(row["nControlProfiles"]),
                "nHits": int(len(hits)),
                "ran": True,
            }
        )
        result_frames.append(full)
        if not hits.empty:
            hit_frames.append(hits)

    de_summary = pd.DataFrame(summary_rows)
    de_results = pd.concat(result_frames, ignore_index=True) if result_frames else pd.DataFrame()
    de_hits = pd.concat(hit_frames, ignore_index=True) if hit_frames else pd.DataFrame()
    shared_genes = shared_direction_genes(de_hits, minDiseaseAreas=2) if not de_hits.empty else pd.DataFrame()
    if not shared_genes.empty:
        shared_genes = attach_gene_symbols(shared_genes, gene_map)

    if de_hits.empty:
        gene_class = pd.DataFrame()
    else:
        class_rows: list[dict[str, object]] = []
        for key_values, group in de_hits.groupby(["gene", "cluster"], observed=True):
            gene, cluster = tuple(key_values)  # type: ignore[misc]
            directions = pd.Series(np.sign(group["log2FoldChange"].astype(float)).astype(int))
            areas = group["diseaseArea"].astype(str).tolist()
            if directions.nunique() > 1:
                evidence_class = "oppositeDirection"
            elif len(areas) >= 2:
                evidence_class = "sharedSameDirection"
            else:
                evidence_class = "areaSelective"
            class_rows.append(
                {
                    "gene": str(gene),
                    "cluster": str(cluster),
                    "evidenceClass": evidence_class,
                    "nDiseaseAreasHit": int(group["diseaseArea"].nunique()),
                    "diseaseAreas": ",".join(sorted(set(areas))),
                    "meanLog2FoldChange": float(group["log2FoldChange"].mean()),
                    "minPadj": float(group["padj"].min()),
                    "interpretedCellType": group["interpretedCellType"].iloc[0],
                    "geneSymbol": group["geneSymbol"].iloc[0] if "geneSymbol" in group.columns else None,
                }
            )
        gene_class = pd.DataFrame(class_rows)

    de_summary.to_csv(cfg.outputDir / "de_summary.csv", index=False)
    de_hits.to_csv(cfg.outputDir / "de_hits.csv", index=False)
    if not de_results.empty:
        de_results.to_parquet(cfg.outputDir / "de_results.parquet", index=False)
    shared_genes.to_csv(cfg.outputDir / "shared_direction_genes.csv", index=False)
    gene_class.to_csv(cfg.outputDir / "gene_evidence_classes.csv", index=False)

    unexpected = _build_unexpected(de_hits, pdata_de, cluster_interp, cfg)
    unexpected.to_csv(cfg.outputDir / "unexpected_expression_candidates.csv", index=False)

    pools = build_evidence_pools(
        restricted=restricted,
        deHits=de_hits,
        sharedGenes=shared_genes,
        geneClass=gene_class,
        unexpected=unexpected,
        padj=cfg.padj,
        lfc=cfg.lfc,
        minDetectionDelta=cfg.minDetectionDelta,
        minTau=cfg.minTau,
        minTargetDetection=cfg.minTargetDetection,
        maxBackgroundDetection=cfg.maxBackgroundDetection,
        minStudiesForSpecificity=cfg.minStudiesForSpecificity,
    )
    shortlist, extended, thresholds = select_review_queue(
        pools,
        primaryBudget=cfg.primaryBudget,
        extendedBudget=cfg.extendedBudget,
        maxPerClassPrimary=cfg.maxPerClassPrimary,
        maxPerClassExtended=cfg.maxPerClassExtended,
        maxPerGene=cfg.maxPerGene,
        maxPerCluster=cfg.maxPerCluster,
        maxPerDiseaseArea=cfg.maxPerDiseaseArea,
    )
    shortlist = _annotate_study_agreement(shortlist, pdata_de, cfg)
    extended = _annotate_study_agreement(extended, pdata_de, cfg)

    shortlist.to_csv(cfg.outputDir / "noteworthy_gene_shortlist.csv", index=False)
    extended.to_csv(cfg.outputDir / "noteworthy_gene_extended.csv", index=False)
    thresholds.to_csv(cfg.outputDir / "candidate_thresholds.csv", index=False)
    pd.concat([shortlist, extended], ignore_index=True, sort=False).to_csv(
        cfg.outputDir / "noteworthy_gene_candidates.csv",
        index=False,
    )

    write_review_figures(
        figDir=cfg.figuresDir,
        pdata=pdata_de,
        shortlist=shortlist,
        extended=extended,
        thresholds=thresholds,
        restricted=restricted,
        deHits=de_hits,
        deResults=de_results,
        unexpected=unexpected,
        clusterKey=cfg.clusterKey,
        studyKey=cfg.studyKey,
        padj=cfg.padj,
        lfc=cfg.lfc,
        maxVolcanoPlots=cfg.maxVolcanoPlots,
        maxEvidencePanels=cfg.maxEvidencePanels,
    )

    snap = snapshot_memory()
    summary = {
        "atlasPath": str(cfg.atlasPath),
        "outputDir": str(cfg.outputDir),
        "pseudobulkProfiles": int(pdata.n_obs),
        "genesScoredForSpecificity": int(len(specificity)),
        "restrictedGenes": int(len(restricted)),
        "eligibleContrasts": int(len(eligible_contrasts)),
        "deContrastsRan": int(de_summary["ran"].sum()) if not de_summary.empty else 0,
        "deHits": int(len(de_hits)),
        "sharedDirectionGenes": int(len(shared_genes)),
        "unexpectedCandidates": int(len(unexpected)),
        "primaryCandidates": int(len(shortlist)),
        "extendedCandidates": int(len(extended)),
        "peakRssBytes": snap.peakRssBytes,
        "peakRssGiB": format_bytes(snap.peakRssBytes),
        "memoryReserveBytes": cfg.memoryReserveBytes,
        "thresholds": {
            "minCellsPerProfile": cfg.minCellsPerProfile,
            "minOverlapStudies": cfg.minOverlapStudies,
            "padj": cfg.padj,
            "lfc": cfg.lfc,
            "minTau": cfg.minTau,
            "minTargetDetection": cfg.minTargetDetection,
            "maxBackgroundDetection": cfg.maxBackgroundDetection,
            "minDetectionDelta": cfg.minDetectionDelta,
            "highPurity": cfg.highPurity,
            "primaryBudget": cfg.primaryBudget,
            "extendedBudget": cfg.extendedBudget,
        },
        "note": (
            "Automated discovery shortlist for manual review. "
            "Expression evidence alone does not establish novelty or mechanism. "
            "Source cell_type labels interpret clusters and do not define pseudobulk groups."
        ),
    }
    _write_json(cfg.outputDir / "analysis_summary.json", summary)
    log.info("Wrote shortlist=%s extended=%s", len(shortlist), len(extended))
    log_memory("analyze-end", logger=log)
    return summary


def analyze_atlas(cfg: AtlasDeAnalysisConfig, *, reuseCheckpoint: bool = True) -> dict[str, Any]:
    from disease_markers.aggregation import aggregate_atlas, build_aggregate_fingerprint

    fingerprint = build_aggregate_fingerprint(cfg)
    pdata = aggregate_atlas(cfg, reuseCheckpoint=reuseCheckpoint)
    return analyze_from_pseudobulk(
        cfg,
        pdata,
        clusterObs=pd.DataFrame(pdata.obs),
        fingerprintSha=str(fingerprint["sha256"]),
    )
