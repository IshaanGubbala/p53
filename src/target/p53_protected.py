from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable


def get_zinc_binding_residues(p53_cfg: dict) -> set[int]:
    residues = p53_cfg.get("zinc_binding_residues", [])
    return {int(pos) for pos in residues}


def _feature_positions(feature: dict) -> set[int]:
    start = feature.get("start")
    end = feature.get("end")
    if start is None or end is None:
        return set()
    return set(range(int(start), int(end) + 1))


def get_dna_contact_residues(annotations_path: Path) -> set[int]:
    if not annotations_path.exists():
        return set()

    payload = json.loads(annotations_path.read_text(encoding="utf-8"))
    features = payload.get("features", [])
    residues: set[int] = set()

    for feature in features:
        ftype = (feature.get("type") or "").upper()
        desc = (feature.get("description") or "").lower()
        if ftype in {"DNA_BIND", "DNA_BINDING"}:
            residues.update(_feature_positions(feature))
        elif "dna-binding" in desc or "dna binding" in desc:
            residues.update(_feature_positions(feature))

    return residues


def get_hotspot_residues(hotspots_path: Path) -> set[int]:
    if not hotspots_path.exists():
        return set()
    payload = json.loads(hotspots_path.read_text(encoding="utf-8"))
    hotspots = payload.get("hotspots", [])
    residues: set[int] = set()
    for entry in hotspots:
        pos = entry.get("pos")
        if pos is None:
            continue
        residues.add(int(pos))
    return residues


def combine_protected_sets(*sets: Iterable[int]) -> set[int]:
    combined: set[int] = set()
    for values in sets:
        combined.update(int(pos) for pos in values)
    return combined
