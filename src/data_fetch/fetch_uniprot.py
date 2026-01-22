from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

USER_AGENT = "p53-stabilimut/0.1"


def _download(url: str, out_path: Path) -> None:
    if out_path.exists():
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            payload = response.read()
    except URLError as exc:
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    out_path.write_bytes(payload)


def download_uniprot_entry(uniprot_id: str, out_path: Path) -> None:
    """Download the UniProt text entry for a protein."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    _download(url, out_path)


def extract_sequence(entry_text: str) -> str:
    """Extract the amino acid sequence from a UniProt text entry."""
    lines = entry_text.splitlines()
    seq_started = False
    seq_parts: list[str] = []

    for line in lines:
        if line.startswith("SQ"):
            seq_started = True
            continue
        if seq_started:
            if line.startswith("//"):
                break
            seq_parts.append("".join(ch for ch in line if ch.isalpha()))

    return "".join(seq_parts)


def _parse_feature_position(token: str) -> int | None:
    token = token.strip()
    if not token:
        return None
    token = token.replace("<", "").replace(">", "").replace("?", "")
    return int(token) if token.isdigit() else None


def extract_features(entry_text: str) -> dict[str, Any]:
    """Extract feature annotations from a UniProt text entry."""
    features: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for line in entry_text.splitlines():
        if not line.startswith("FT"):
            continue

        key = line[5:13].strip()
        start_txt = line[14:20].strip()
        end_txt = line[21:27].strip()
        desc = line[34:].strip()

        if key:
            if current:
                features.append(current)
            current = {
                "type": key,
                "start": _parse_feature_position(start_txt),
                "end": _parse_feature_position(end_txt),
                "description": desc,
            }
        elif current is not None:
            if desc:
                current["description"] = (current.get("description", "") + " " + desc).strip()

    if current:
        features.append(current)

    return {"features": features}


def save_fasta(seq: str, out_fasta: Path) -> None:
    """Write a FASTA file with 60-character line wrapping."""
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    header = f">{out_fasta.stem}"
    with out_fasta.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for i in range(0, len(seq), 60):
            handle.write(seq[i : i + 60] + "\n")


def save_annotations(features: dict[str, Any], out_json: Path) -> None:
    """Write UniProt feature annotations to JSON."""
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as handle:
        json.dump(features, handle, indent=2, sort_keys=True)


def extract_uniprot_id(entry_text: str) -> str | None:
    """Try to capture the UniProt accession from the entry text."""
    match = re.search(r"^AC\s+([A-Z0-9]+);", entry_text, re.MULTILINE)
    return match.group(1) if match else None
