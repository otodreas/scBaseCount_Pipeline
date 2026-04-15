from __future__ import annotations

import subprocess
import time


def curl_get(url: str, *, retries: int = 3, timeout_s: int = 30) -> str:
    last_stderr = ""
    for attempt in range(retries + 1):
        result = subprocess.run(
            ["curl", "-sf", "--max-time", str(timeout_s), url],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return result.stdout
        last_stderr = result.stderr.strip()
        if attempt < retries:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"curl failed for {url!r} (exit {result.returncode}): {last_stderr}")
