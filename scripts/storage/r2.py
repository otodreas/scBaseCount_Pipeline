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


def download_from_r2(r2_key: str, local_path: Path, verify_md5: bool = False) -> None:
    bucket = os.environ["BUCKET"]
    local_path.parent.mkdir(parents=True, exist_ok=True)  # local path is a file, not the location of the file
    _log.info("Downloading r2://%s/%s -> %s", bucket, r2_key, local_path)
    _r2_client().download_file(bucket, r2_key, str(local_path))
    if verify_md5:
        _verify_downloaded_md5(r2_key, local_path)


def _verify_downloaded_md5(r2_key: str, local_path: Path) -> None:
    from storage.transfer import _local_md5_b64

    stored = r2_object_md5(r2_key)
    if stored is None:
        _log.warning("No stored MD5 for r2://%s/%s, skipping download verification", os.environ["BUCKET"], r2_key)
        return
    local = _local_md5_b64(local_path)
    if local != stored:
        raise ValueError(
            f"Download MD5 mismatch for r2://{os.environ['BUCKET']}/{r2_key}: local={local} stored={stored}"
        )
    _log.info("Download MD5 verified: r2://%s/%s", os.environ["BUCKET"], r2_key)


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


def delete_from_r2(r2_key: str) -> None:
    bucket = os.environ["BUCKET"]
    client = _r2_client()

    try:
        client.head_object(Bucket=bucket, Key=r2_key)
    except ClientError as exc:
        if exc.response["Error"]["Code"] in ("404", "NoSuchKey"):
            raise FileNotFoundError(f"R2 object not found: r2://{bucket}/{r2_key}") from None
        raise
    _log.info("Deleting r2://%s/%s", bucket, r2_key)
    client.delete_object(Bucket=bucket, Key=r2_key)


def delete_r2_prefix(prefix: str) -> list[str]:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    keys = sorted(fetch_uploaded_r2_keys(prefix))
    if not keys:
        _log.info("No objects under prefix %r to delete", prefix)
        return []
    failed: set[str] = set()
    for start in range(0, len(keys), 1000):
        batch = keys[start : start + 1000]
        resp = client.delete_objects(
            Bucket=bucket,
            Delete={"Objects": [{"Key": key} for key in batch], "Quiet": True},
        )
        for err in resp.get("Errors", []):
            _log.error("Failed to delete r2://%s/%s: %s", bucket, err["Key"], err["Message"])
            failed.add(err["Key"])
    deleted = [key for key in keys if key not in failed]
    _log.info("Deleted %d object(s) under prefix %r", len(deleted), prefix)
    return deleted
