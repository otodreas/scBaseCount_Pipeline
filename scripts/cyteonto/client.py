from __future__ import annotations

import logging
import time
from pathlib import Path

import httpx
import orjson
import pandas as pd


def _check_response(resp: httpx.Response, context: str) -> None:
    if resp.is_success:
        return
    body = resp.text[:300].strip() if resp.text else "(empty body)"
    raise RuntimeError(
        f"CyteOnto API error [{context}]: HTTP {resp.status_code}  {body}"
    )


def get_status(run_id: str, base_url: str) -> dict:
    resp = httpx.get(f"{base_url}/status/{run_id}", timeout=30)
    _check_response(resp, f"status/{run_id}")
    return resp.json()


def submit(payload: dict, base_url: str, log: logging.Logger) -> str:
    resp = httpx.post(
        f"{base_url}/compare",
        content=orjson.dumps(payload),
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    _check_response(resp, "submit")
    data = resp.json()
    log.info("submitted  runId=%s  state=%s", data["runId"], data["state"])
    return data["runId"]


def poll(
    run_id: str,
    base_url: str,
    interval_s: int,
    timeout_s: int,
    log: logging.Logger,
) -> dict:
    deadline = time.time() + timeout_s
    last_state = None
    while time.time() < deadline:
        resp = httpx.get(f"{base_url}/status/{run_id}", timeout=30)
        if resp.status_code >= 500:
            body = resp.text[:300].strip() if resp.text else "(empty body)"
            log.warning("poll: HTTP %d  retrying  %s", resp.status_code, body)
            print(f"[cyteonto] server error {resp.status_code} while polling -- retrying")
            time.sleep(interval_s)
            continue
        _check_response(resp, f"poll/{run_id}")
        status = resp.json()
        if status["state"] != last_state:
            log.info("status    runId=%s  state=%s", run_id, status["state"])
            print(f"[cyteonto] {run_id}  ->  {status['state']}")
            last_state = status["state"]
        if status["state"] in ("completed", "failed"):
            return status
        time.sleep(interval_s)
    raise TimeoutError(f"run {run_id} did not finish within {timeout_s}s")


def fetch_result(
    run_id: str,
    base_url: str,
    out_path: Path,
    log: logging.Logger,
) -> pd.DataFrame:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    resp = httpx.get(
        f"{base_url}/result/{run_id}",
        params={"format": "csv"},
        timeout=60,
    )
    _check_response(resp, f"fetch/{run_id}")
    out_path.write_bytes(resp.content)
    log.info("fetched   runId=%s  path=%s", run_id, out_path)
    return pd.read_csv(out_path)
