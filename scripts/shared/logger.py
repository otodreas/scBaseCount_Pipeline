from __future__ import annotations

import logging

from shared.repo import REPO_ROOT


def configure_file_logger(log_filename: str, logger_name: str) -> logging.Logger:
    log_path = REPO_ROOT / "logs" / log_filename
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_path, mode="a", encoding="utf-8"),
        ],
    )
    return logging.getLogger(logger_name)
