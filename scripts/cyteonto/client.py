from __future__ import annotations

import time
from pathlib import Path

import httpx
import pandas as pd


def submit(payload: dict, base_url: str) -> str:
    resp = httpx.post(f"{base_url}/compare", json=payload, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    print(f"[cyteonto] submitted  runId={data['runId']}  state={data['state']}")
    return data["runId"]


def poll(run_id: str, base_url: str, interval_s: int, timeout_s: int) -> dict:
    deadline = time.time() + timeout_s
    last_state = None
    while time.time() < deadline:
        resp = httpx.get(f"{base_url}/status/{run_id}", timeout=30)
        resp.raise_for_status()
        status = resp.json()
        if status["state"] != last_state:
            print(f"[cyteonto] status    {status['state']}")
            last_state = status["state"]
        if status["state"] in ("completed", "failed"):
            return status
        time.sleep(interval_s)
    raise TimeoutError(f"run {run_id} did not finish within {timeout_s}s")


def fetch_result(run_id: str, base_url: str, out_path: Path) -> pd.DataFrame:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    resp = httpx.get(
        f"{base_url}/result/{run_id}",
        params={"format": "csv"},
        timeout=60,
    )
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    print(f"[cyteonto] saved     {out_path}")
    return pd.read_csv(out_path)
