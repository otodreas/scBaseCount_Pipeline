from __future__ import annotations

import datetime
import json
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc

from cyteonto.client import fetch_result, get_status, poll, submit
from cyteonto.config import CyteOntoConfig, _REPO_ROOT
from cyteonto.payload import build_payload, write_payload

_LOG_PATH = Path(__file__).parents[2] / "logs" / "cyteonto.log"
_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(_LOG_PATH, mode="a", encoding="utf-8"),
    ],
)

_log = logging.getLogger(__name__)

_DEFAULT_BASE_URL = "https://cyteonto.nygen.io"
_DEFAULT_COMPLETED_DIR = _REPO_ROOT / "output" / "cyteonto" / "runs" / "completed"
_DEFAULT_PENDING_DIR = _REPO_ROOT / "output" / "cyteonto" / "runs" / "pending"


# --- pending-run helpers (one JSON stub per run in runs/pending/) ---

def _pending_stub_path(run_id: str, pending_dir: Path) -> Path:
    return pending_dir / f"{run_id}.json"


def _add_pending_run(run_id: str, h5ad_stem: str, pending_dir: Path) -> None:
    pending_dir.mkdir(parents=True, exist_ok=True)
    stub = {
        "runId": run_id,
        "h5adStem": h5ad_stem,
        "submittedAt": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    _pending_stub_path(run_id, pending_dir).write_text(json.dumps(stub, indent=2))


def _remove_pending_run(run_id: str, pending_dir: Path) -> None:
    stub = _pending_stub_path(run_id, pending_dir)
    if stub.exists():
        stub.unlink()


def _load_pending(pending_dir: Path) -> list[dict]:
    if not pending_dir.exists():
        return []
    return [json.loads(p.read_text()) for p in sorted(pending_dir.glob("*.json"))]


# --- public API ---

def run_cyteonto(cfg: CyteOntoConfig) -> pd.DataFrame | None:
    _log.info("start  input=%s", cfg.h5adPath)

    adata = sc.read_h5ad(cfg.h5adPath)
    _log.info("loaded %d cells  %d genes", adata.n_obs, adata.n_vars)

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
    _add_pending_run(run_id, cfg.h5adPath.stem, cfg.pendingDir)
    _log.info("pending stub written  runId=%s  dir=%s", run_id, cfg.pendingDir)

    try:
        status = poll(run_id, cfg.baseUrl, cfg.pollIntervalS, cfg.pollTimeoutS, _log)
    except KeyboardInterrupt:
        _log.info(
            "polling stopped  runId=%s  (run continues on server; call check_pending_runs() to resume)",
            run_id,
        )
        return None

    if status["state"] == "failed":
        _remove_pending_run(run_id, cfg.pendingDir)
        _log.error("failed  runId=%s  error=%s", run_id, status.get("error"))
        raise RuntimeError(f"CyteOnto run failed: {status.get('error')}")

    _log.info("completed  runId=%s  rows=%s", run_id, status.get("numRows"))

    out_path = cfg.completedDir / f"{run_id}.csv"
    df = fetch_result(run_id, cfg.baseUrl, out_path, _log)
    _remove_pending_run(run_id, cfg.pendingDir)
    _log.info("done  saved=%s", out_path)
    return df


def check_pending_runs(
    base_url: str = _DEFAULT_BASE_URL,
    completed_dir: Path = _DEFAULT_COMPLETED_DIR,
    pending_dir: Path = _DEFAULT_PENDING_DIR,
) -> dict[str, pd.DataFrame]:
    runs = _load_pending(pending_dir)
    if not runs:
        _log.info("check_pending_runs: no pending runs")
        return {}

    _log.info("check_pending_runs: checking %d run(s)", len(runs))

    completed: dict[str, pd.DataFrame] = {}

    for entry in runs:
        run_id = entry["runId"]
        h5ad_stem = entry.get("h5adStem", "unknown")
        status = get_status(run_id, base_url)
        state = status["state"]
        _log.info("check  runId=%s  h5adStem=%s  state=%s", run_id, h5ad_stem, state)

        if state == "completed":
            out_path = completed_dir / f"{run_id}.csv"
            completed[run_id] = fetch_result(run_id, base_url, out_path, _log)
            _remove_pending_run(run_id, pending_dir)
        elif state == "failed":
            _log.error(
                "failed  runId=%s  h5adStem=%s  error=%s",
                run_id,
                h5ad_stem,
                status.get("error"),
            )
            _remove_pending_run(run_id, pending_dir)

    still_pending = len(_load_pending(pending_dir))
    _log.info(
        "check_pending_runs: completed=%d  still_pending=%d",
        len(completed),
        still_pending,
    )
    return completed
