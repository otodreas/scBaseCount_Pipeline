from __future__ import annotations

import os
from pathlib import Path

import boto3

from shared.logger import configure_file_logger

_log = configure_file_logger("r2.log", __name__)


def _r2_client() -> boto3.client:
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
    _log.info("R2 bucket contains %d existing object(s)", len(keys))
    return keys


def upload_to_r2(local_path: Path, r2_key: str) -> None:
    bucket = os.environ["BUCKET"]
    _log.info("Uploading %s -> r2://%s/%s", local_path.name, bucket, r2_key)
    _r2_client().upload_file(str(local_path), bucket, r2_key)


def verify_upload(r2_key: str) -> bool:
    bucket = os.environ["BUCKET"]
    client = _r2_client()
    try:
        client.head_object(Bucket=bucket, Key=r2_key)
        _log.info("Upload verified: r2://%s/%s", bucket, r2_key)
        return True
    except client.exceptions.ClientError as exc:
        if exc.response["Error"]["Code"] == "404":
            _log.warning("Upload verification failed: r2://%s/%s not found", bucket, r2_key)
            return False
        raise
