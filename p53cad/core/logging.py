import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional, Union


DEFAULT_LOG_FILE = Path("logs/p53cad_workflow.log")


def setup_logging(level=logging.INFO, log_file: Optional[Union[Path, str]] = None) -> Path:
    """
    Sets up logging for the p53cad package.
    Always logs to console and a persistent workflow file.
    """
    env_log_file = os.getenv("P53CAD_LOG_FILE", "").strip()
    if env_log_file:
        resolved_log_file = Path(env_log_file).expanduser()
    elif log_file is not None:
        resolved_log_file = Path(log_file).expanduser()
    else:
        resolved_log_file = DEFAULT_LOG_FILE

    resolved_log_file.parent.mkdir(parents=True, exist_ok=True)

    handlers = [
        logging.StreamHandler(sys.stdout),
        RotatingFileHandler(
            resolved_log_file,
            maxBytes=10 * 1024 * 1024,
            backupCount=5,
            encoding="utf-8",
        ),
    ]

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )
    logging.captureWarnings(True)
    return resolved_log_file

def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)
