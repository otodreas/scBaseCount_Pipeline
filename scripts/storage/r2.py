from __future__ import annotations

import os
from pathlib import Path

import boto3
from botocore.exceptions import ClientError
from shared.logger import configure_file_logger

_log = configure_file_logger("r2.log", __name__)


def _r2_client() -> boto3.client:
    return boto3.client(
        "s3",
        endpoint_url=os.environ["ENDPOINT_URL"],
        aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"],
    )


def r2_key_exists(r2_key: str) -> bool:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    try:
        client.head_object(Bucket=bucket, Key=r2_key)
        return True
    except ClientError as exc:
        if exc.response["Error"]["Code"] in ("404", "NoSuchKey"):
            return False
        raise


def download_from_r2(r2_key: str, local_path: Path) -> None:
    bucket = os.environ["BUCKET"]
    local_path.parent.mkdir(parents=True, exist_ok=True)
    _log.info("Downloading r2://%s/%s -> %s", bucket, r2_key, local_path)
    _r2_client().download_file(bucket, r2_key, str(local_path))


def fetch_uploaded_r2_keys(prefix: str | None = None) -> set[str]:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    keys: set[str] = set()
    paginator = client.get_paginator("list_objects_v2")
    paginate_kwargs: dict = {"Bucket": bucket}
    if prefix is not None:
        paginate_kwargs["Prefix"] = prefix
    for page in paginator.paginate(**paginate_kwargs):
        for obj in page.get("Contents") or []:
            keys.add(obj["Key"])
    if prefix is None:
        _log.info("R2 bucket contains %d existing object(s)", len(keys))
    else:
        _log.info("R2 prefix %r contains %d existing object(s)", prefix, len(keys))
    return keys


def r2_object_md5(r2_key: str) -> str | None:
    from storage.transfer import _MD5_METADATA_KEY

    bucket = os.environ["BUCKET"]
    try:
        resp = _r2_client().head_object(Bucket=bucket, Key=r2_key)
        return resp.get("Metadata", {}).get(_MD5_METADATA_KEY)
    except ClientError as exc:
        if exc.response["Error"]["Code"] in ("404", "NoSuchKey"):
            return None
        raise


def upload_to_r2(local_path: Path, r2_key: str, extra_metadata: dict[str, str] | None = None) -> None:
    bucket = os.environ["BUCKET"]
    _log.info("Uploading %s -> r2://%s/%s", local_path.name, bucket, r2_key)
    kwargs: dict = {"Filename": str(local_path), "Bucket": bucket, "Key": r2_key}
    if extra_metadata:
        kwargs["ExtraArgs"] = {"Metadata": extra_metadata}
    _r2_client().upload_file(**kwargs)


def verify_upload(r2_key: str) -> bool:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    try:
        client.head_object(Bucket=bucket, Key=r2_key)
        _log.info("Upload verified: r2://%s/%s", bucket, r2_key)
        return True
    except ClientError as exc:
        if exc.response["Error"]["Code"] in ("404", "NoSuchKey"):
            _log.warning("Upload verification failed: r2://%s/%s not found", bucket, r2_key)
            return False
        raise
