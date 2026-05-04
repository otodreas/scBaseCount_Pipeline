from __future__ import annotations

import datetime
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

from cyteonto.client import fetch_result, get_status, poll, submit
from cyteonto.config import CyteOntoConfig, _REPO_ROOT
from cyteonto.payload import build_payload, write_payload
from shared.logger import configure_file_logger

_log = configure_file_logger("cyteonto.log", __name__)

_DEFAULT_BASE_URL = "https://cyteonto.nygen.io"
_DEFAULT_RUNS_DIR = _REPO_ROOT / "output" / "cyteonto" / "runs"


# --- run stub helpers (one JSON per run in runs/, completedAt null until done) ---

def _stub_path(run_id: str, runs_dir: Path) -> Path:
    return runs_dir / f"{run_id}.json"


def _write_run_stub(run_id: str, h5ad_stem: str, payload_path: Path, runs_dir: Path) -> None:
    runs_dir.mkdir(parents=True, exist_ok=True)
    stub = {
        "runId": run_id,
        "h5adStem": h5ad_stem,
        "payloadPath": str(payload_path.relative_to(_REPO_ROOT)),
        "submittedAt": datetime.datetime.now().isoformat(timespec="seconds"),
        "completedAt": None,
    }
    _stub_path(run_id, runs_dir).write_text(json.dumps(stub, indent=2))


def _mark_completed(run_id: str, runs_dir: Path) -> None:
    path = _stub_path(run_id, runs_dir)
    if not path.exists():
        return
    stub = json.loads(path.read_text())
    stub["completedAt"] = datetime.datetime.now().isoformat(timespec="seconds")
    path.write_text(json.dumps(stub, indent=2))


def _load_pending(runs_dir: Path) -> list[dict]:
    if not runs_dir.exists():
        return []
    stubs = [json.loads(p.read_text()) for p in sorted(runs_dir.glob("*.json"))]
    return [s for s in stubs if s.get("completedAt") is None]


# --- public API ---

def run_cyteonto(cfg: CyteOntoConfig) -> pd.DataFrame | None:
    _log.info("start  input=%s", cfg.h5adPath)
    print(f"[cyteonto] input     {cfg.h5adPath.name}")

    adata = sc.read_h5ad(cfg.h5adPath)
    _log.info("loaded %d cells  %d genes", adata.n_obs, adata.n_vars)
    print(f"[cyteonto] loaded    {adata.n_obs} cells  {adata.n_vars} genes")

    payload = build_payload(adata)
    _log.info(
        "payload  author_labels=%d  algorithms=%d",
        len(payload["authorLabels"]),
        len(payload["algorithms"]),
    )

    payload_path = cfg.payloadDir / f"{cfg.h5adPath.stem}_annotations.json"
    write_payload(payload, payload_path)
    _log.info("payload written  path=%s", payload_path)

    run_id = submit(payload, cfg.baseUrl, _log)
    print(f"[cyteonto] submitted  {run_id}")
    _write_run_stub(run_id, cfg.h5adPath.stem, payload_path, cfg.runsDir)
    _log.info("run stub written  runId=%s  dir=%s", run_id, cfg.runsDir)

    try:
        status = poll(run_id, cfg.baseUrl, cfg.pollIntervalS, cfg.pollTimeoutS, _log)
    except KeyboardInterrupt:
        _log.info(
            "polling stopped  runId=%s  (run continues on server; call check_pending_runs() to resume)",
            run_id,
        )
        print(f"[cyteonto] polling stopped  {run_id}  (call check_pending_runs() to resume)")
        return None

    if status["state"] == "failed":
        _mark_completed(run_id, cfg.runsDir)
        _log.error("failed  runId=%s  error=%s", run_id, status.get("error"))
        raise RuntimeError(f"CyteOnto run failed: {status.get('error')}")

    _log.info("completed  runId=%s  rows=%s", run_id, status.get("numRows"))

    out_path = cfg.runsDir / f"{run_id}.csv"
    df = fetch_result(run_id, cfg.baseUrl, out_path, _log)
    _mark_completed(run_id, cfg.runsDir)
    _log.info("done  saved=%s", out_path)
    print(f"[cyteonto] done      {run_id}  rows={status.get('numRows')}  saved={out_path.name}")
    return df


def check_pending_runs(
    base_url: str = _DEFAULT_BASE_URL,
    runs_dir: Path = _DEFAULT_RUNS_DIR,
) -> dict[str, pd.DataFrame]:
    runs = _load_pending(runs_dir)
    if not runs:
        _log.info("check_pending_runs: no pending runs")
        print("[cyteonto] no pending runs")
        return {}

    _log.info("check_pending_runs: checking %d run(s)", len(runs))
    print(f"[cyteonto] checking {len(runs)} pending run(s)")

    completed: dict[str, pd.DataFrame] = {}

    for entry in runs:
        run_id = entry["runId"]
        h5ad_stem = entry.get("h5adStem", "unknown")
        status = get_status(run_id, base_url)
        state = status["state"]
        _log.info("check  runId=%s  h5adStem=%s  state=%s", run_id, h5ad_stem, state)
        print(f"[cyteonto] {run_id}  ({h5ad_stem})  ->  {state}")

        if state == "completed":
            out_path = runs_dir / f"{run_id}.csv"
            completed[run_id] = fetch_result(run_id, base_url, out_path, _log)
            _mark_completed(run_id, runs_dir)
        elif state == "failed":
            _log.error(
                "failed  runId=%s  h5adStem=%s  error=%s",
                run_id,
                h5ad_stem,
                status.get("error"),
            )
            _mark_completed(run_id, runs_dir)

    still_pending = len(_load_pending(runs_dir))
    _log.info(
        "check_pending_runs: completed=%d  still_pending=%d",
        len(completed),
        still_pending,
    )
    print(f"[cyteonto] done  completed={len(completed)}  still_pending={still_pending}")
    return completed
