from __future__ import annotations

import logging
from typing import Optional

_LOGGING_CONFIGURED = False


def setup_logging(level: str = "INFO") -> logging.Logger:
    """Configure root logging once and return the root logger."""
    global _LOGGING_CONFIGURED
    if not _LOGGING_CONFIGURED:
        logging.basicConfig(
            level=level,
            format="%(asctime)s %(levelname)s %(name)s %(message)s",
        )
        _LOGGING_CONFIGURED = True
    return logging.getLogger()


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a named logger, defaulting to the root logger."""
    return logging.getLogger(name)
