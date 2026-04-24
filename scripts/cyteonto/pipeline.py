from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import scanpy as sc

from cyteonto.client import fetch_result, poll, submit
from cyteonto.config import CyteOntoConfig
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


def run_cyteonto(cfg: CyteOntoConfig) -> pd.DataFrame:
    _log.info("start  input=%s", cfg.h5adPath)

    adata = sc.read_h5ad(cfg.h5adPath)
    _log.info("loaded %d cells  %d genes", adata.n_obs, adata.n_vars)

    payload = build_payload(adata)
    n_author = len(payload["authorLabels"])
    n_algos = len(payload["algorithms"])
    _log.info("payload  author_labels=%d  algorithms=%d", n_author, n_algos)

    payload_path = cfg.payloadDir / f"{cfg.h5adPath.stem}_annotations.json"
    write_payload(payload, payload_path)
    _log.info("payload written  path=%s", payload_path)

    run_id = submit(payload, cfg.baseUrl, _log)
    status = poll(run_id, cfg.baseUrl, cfg.pollIntervalS, cfg.pollTimeoutS, _log)

    if status["state"] == "failed":
        _log.error("failed  runId=%s  error=%s", run_id, status.get("error"))
        raise RuntimeError(f"CyteOnto run failed: {status.get('error')}")

    num_rows = status.get("numRows")
    _log.info("completed  runId=%s  rows=%s", run_id, num_rows)

    out_path = cfg.resultsDir / f"{run_id}.csv"
    df = fetch_result(run_id, cfg.baseUrl, out_path, _log)
    _log.info("done  saved=%s", out_path)
    return df
