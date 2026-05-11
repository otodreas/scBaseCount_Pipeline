from __future__ import annotations

import logging
from pathlib import Path


def safe_delete(path: Path, logger: logging.Logger) -> None:
    try:
        if path.exists():
            path.unlink()
            logger.debug("Deleted %s", path)
    except OSError as exc:
        logger.warning("Could not delete %s: %s", path, exc)
