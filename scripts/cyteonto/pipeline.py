from __future__ import annotations

import scanpy as sc
import pandas as pd

from cyteonto.client import fetch_result, poll, submit
from cyteonto.config import CyteOntoConfig
from cyteonto.payload import build_payload, write_payload


def run_cyteonto(cfg: CyteOntoConfig) -> pd.DataFrame:
    adata = sc.read_h5ad(cfg.h5adPath)

    payload = build_payload(adata)

    payload_path = cfg.payloadDir / f"{cfg.h5adPath.stem}_annotations.json"
    write_payload(payload, payload_path)
    print(f"[cyteonto] payload   {payload_path}")

    run_id = submit(payload, cfg.baseUrl)
    status = poll(run_id, cfg.baseUrl, cfg.pollIntervalS, cfg.pollTimeoutS)

    if status["state"] == "failed":
        raise RuntimeError(f"CyteOnto run failed: {status.get('error')}")

    out_path = cfg.resultsDir / f"{run_id}.csv"
    return fetch_result(run_id, cfg.baseUrl, out_path)
