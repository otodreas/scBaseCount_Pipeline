from __future__ import annotations

import csv
import logging
import os
import sys
from pathlib import Path
from urllib.parse import urlparse

import boto3
import scanpy as sc
from cytetype import CyteType, rank_genes_groups_backed
from dotenv import load_dotenv
from google.cloud import storage as gcs

from cluster_validation import ClusterValidationConfig, run_cluster_validation
from shared.logger import configure_file_logger
from shared.repo import REPO_ROOT
from study_context import ExperimentContext, experiment_context_summary

load_dotenv()

DATASETS_CSV = REPO_ROOT / "output" / "metadata" / "datasets.csv"
CONTEXTS_JSONL = REPO_ROOT / "output" / "context" / "contexts.jsonl"
CYTETYPE_OUTPUT_DIR = REPO_ROOT / "output" / "cytetype" / "data"

R2_PREFIX = "cytetype"

GROUP_KEY_FALLBACK = "leiden_merged"
N_TOP_GENES = 100

log = configure_file_logger("annotation_pipeline.log", __name__)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def _r2_client() -> "boto3.client":
    return boto3.client(
        "s3",
        endpoint_url=os.environ["ENDPOINT_URL"],
        aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"],
    )


def fetch_uploaded_r2_keys() -> set[str]:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    keys: set[str] = set()
    paginator = client.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket):
        for obj in page.get("Contents") or []:
            keys.add(obj["Key"])
    log.info("R2 bucket contains %d existing object(s)", len(keys))
    return keys


def read_datasets() -> list[dict[str, str]]:
    with open(DATASETS_CSV, newline="") as fh:
        return list(csv.DictReader(fh))


def load_contexts() -> dict[str, ExperimentContext]:
    if not CONTEXTS_JSONL.exists():
        log.warning("contexts.jsonl not found at %s; study context will be empty for all samples", CONTEXTS_JSONL)
        return {}
    contexts: dict[str, ExperimentContext] = {}
    with open(CONTEXTS_JSONL) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ctx = ExperimentContext.model_validate_json(line)
            contexts[ctx.accession] = ctx
    return contexts


def download_from_gcs(gs_uri: str, dest: Path) -> None:
    parsed = urlparse(gs_uri)
    bucket_name = parsed.netloc
    blob_name = parsed.path.lstrip("/")
    dest.parent.mkdir(parents=True, exist_ok=True)
    client = gcs.Client(project=os.environ["GCP_PROJECT"])
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    log.info("Downloading %s -> %s", gs_uri, dest)
    blob.download_to_filename(str(dest))


def run_cytetype_step(
    input_path: Path,
    group_key: str,
    study_context: str,
    output_path: Path,
) -> None:
    adata = sc.read_h5ad(input_path, backed="r")
    rank_genes_groups_backed(adata, groupby=group_key, use_raw=False, key_added=f"rank_genes_{group_key}")
    annotator = CyteType(adata, group_key, rank_key=f"rank_genes_{group_key}", n_top_genes=N_TOP_GENES)
    adata = annotator.run(study_context=study_context, auth_token=os.environ["CYTETYPE_API_KEY"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)


def upload_to_r2(local_path: Path, r2_key: str) -> None:
    bucket = os.environ["BUCKET"]
    log.info("Uploading %s -> r2://%s/%s", local_path.name, bucket, r2_key)
    _r2_client().upload_file(str(local_path), bucket, r2_key)


def _safe_delete(path: Path) -> None:
    try:
        if path.exists():
            path.unlink()
            log.debug("Deleted %s", path)
    except OSError as exc:
        log.warning("Could not delete %s: %s", path, exc)


def process_accession(
    srx: str,
    gs_uri: str,
    contexts: dict[str, ExperimentContext],
) -> None:
    cfg = ClusterValidationConfig(srxAccession=srx)
    raw_h5ad = cfg.localH5adRoot / f"{srx}.h5ad"
    cytetype_h5ad = CYTETYPE_OUTPUT_DIR / f"{srx}_cytetype_annotated.h5ad"
    r2_key = f"{R2_PREFIX}/{srx}_cytetype_annotated.h5ad"

    try:
        download_from_gcs(gs_uri, raw_h5ad)

        _, result = run_cluster_validation(cfg)
        _safe_delete(raw_h5ad)

        ctx = contexts.get(srx)
        study_context = experiment_context_summary(ctx) if ctx else ""
        if not ctx:
            log.warning("%s: no study context found in contexts.jsonl; proceeding with empty context", srx)

        group_key = result.clusterKey or GROUP_KEY_FALLBACK
        run_cytetype_step(result.adataPath, group_key, study_context, cytetype_h5ad)
        _safe_delete(result.adataPath)

        upload_to_r2(cytetype_h5ad, r2_key)
        _safe_delete(cytetype_h5ad)

        log.info("%s: done", srx)

    except Exception:
        log.exception("%s: pipeline failed, skipping", srx)
        _safe_delete(raw_h5ad)
        _safe_delete(cfg.outputDir / f"{srx}_clustered.h5ad")
        _safe_delete(cytetype_h5ad)


def main() -> None:
    datasets = read_datasets()
    log.info("Loaded %d accession(s) from datasets.csv", len(datasets))

    uploaded = fetch_uploaded_r2_keys()
    contexts = load_contexts()

    skipped = 0
    for row in datasets:
        srx = row["srx_accession"]
        gs_uri = row["file_path"]
        r2_key = f"{R2_PREFIX}/{srx}_cytetype_annotated.h5ad"

        if r2_key in uploaded:
            log.info("%s: already uploaded, skipping", srx)
            skipped += 1
            continue

        log.info("%s: starting", srx)
        process_accession(srx, gs_uri, contexts)

    log.info("Pipeline complete. %d skipped, %d processed.", skipped, len(datasets) - skipped)


if __name__ == "__main__":
    main()
