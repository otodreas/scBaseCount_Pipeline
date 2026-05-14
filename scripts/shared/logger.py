from __future__ import annotations

import logging
import sys

from shared.repo import REPO_ROOT


def configure_file_logger(log_filename: str, logger_name: str) -> logging.Logger:
    log_path = REPO_ROOT / "logs" / log_filename
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_path.resolve()) for h in logger.handlers):
        handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(handler)

    return logger


def add_stdout_handler() -> None:
    root = logging.getLogger()
    if not any(isinstance(h, logging.StreamHandler) and getattr(h, "stream", None) is sys.stdout for h in root.handlers):
        root.addHandler(logging.StreamHandler(sys.stdout))


def log_run_separator(logger: logging.Logger) -> None:
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write("\n")
