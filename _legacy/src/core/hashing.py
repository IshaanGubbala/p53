from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return SHA256 hash for a file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text_sha256(text: str) -> str:
    """Return SHA256 hash for a text string."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def json_sha256(obj: Any) -> str:
    """Return SHA256 hash for JSON-serializable content."""
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def dict_sha256(data: dict[str, Any]) -> str:
    """Return SHA256 hash for a dict using stable JSON encoding."""
    return json_sha256(data)
