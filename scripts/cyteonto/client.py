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
    raise RuntimeError(f"CyteOnto API error [{context}]: HTTP {resp.status_code}  {body}")


def _response_detail(resp: httpx.Response, max_len: int = 200) -> str:
    if not resp.text:
        return "(empty body)"
    return resp.text[:max_len].strip()


def submit(payload: dict, base_url: str, log: logging.Logger) -> str:
    """POST /compare and return the run id."""
    resp = httpx.post(
        f"{base_url}/compare",
        content=orjson.dumps(payload),  # , option=orjson.OPT_SERIALIZE_NUMPY),
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    _check_response(resp, "submit")
    data = resp.json()
    log.info("submitted  runId=%s  state=%s", data["runId"], data["state"])
    return data["runId"]


def poll_result(
    run_id: str,
    base_url: str,
    out_path: Path,
    interval_s: int,
    timeout_s: int,
    log: logging.Logger,
) -> pd.DataFrame:
    """Poll /result until 200 with CSV body; 409 means still processing."""
    deadline = time.time() + timeout_s
    last_pending_detail: str | None = None

    while time.time() < deadline:
        try:
            resp = httpx.get(
                f"{base_url}/result/{run_id}",
                params={"format": "csv"},
                timeout=60,
            )
        except httpx.RequestError as exc:
            log.warning("poll_result: request error for %s  retrying  %s", run_id, exc)
            print("[cyteonto] request error while polling -- retrying")
            time.sleep(interval_s)
            continue

        if resp.status_code == 200:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(resp.content)
            log.info("fetched   runId=%s  path=%s", run_id, out_path)
            return pd.read_csv(out_path)

        if resp.status_code == 409:
            detail = _response_detail(resp)
            if detail != last_pending_detail:
                log.info("poll_result  runId=%s  still running  %s", run_id, detail)
                print(f"[cyteonto] {run_id}  ->  still running")
                last_pending_detail = detail
        else:
            detail = _response_detail(resp)
            log.warning(
                "poll_result: HTTP %d for %s  retrying  %s",
                resp.status_code,
                run_id,
                detail,
            )
            print(f"[cyteonto] HTTP {resp.status_code} while polling -- retrying")

        time.sleep(interval_s)

    raise TimeoutError(f"run {run_id} did not finish within {timeout_s}s")
