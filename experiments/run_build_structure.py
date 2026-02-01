from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.hashing import dict_sha256, file_sha256
from src.core.logging import get_logger
from src.structure.clean_structure import (
    extract_domain_pdb,
    remove_altloc,
    standardize_residue_numbering,
    validate_structure,
)
from src.structure.residue_geometry import compute_distance_matrix
from src.structure.sasa_burial import classify_burial, compute_sasa


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
    }


def _find_latest_pdb(raw_dir: Path) -> Path:
    candidates = list(raw_dir.glob("*.pdb"))
    if not candidates:
        raise RuntimeError(f"No PDB files found in {raw_dir}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _structure_basename(uniprot_id: str, start: int, end: int) -> str:
    return f"{uniprot_id}_core_{start}_{end}"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_distance_map(dist_map: dict[tuple[int, int], float], out_path: Path) -> None:
    rows = []
    for (pos_i, pos_j), dist in dist_map.items():
        if pos_i >= pos_j:
            continue
        rows.append({"pos_i": pos_i, "pos_j": pos_j, "distance": dist})

    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))

    p53_cfg = configs.get("p53", {})
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    core_domain = p53_cfg.get("core_domain", {})
    start = int(core_domain.get("start", 1))
    end = int(core_domain.get("end", 9999))

    structure_cfg = p53_cfg.get("structure", {})
    chain_id = structure_cfg.get("chain_id")
    atom_name = structure_cfg.get("atom_name", "CA")
    sasa_thresholds = structure_cfg.get("sasa_thresholds")

    raw_alphafold_dir = paths["raw"] / "alphafold"
    base_pdb = _find_latest_pdb(raw_alphafold_dir)

    out_dir = paths["interim"] / "structure"
    base_name = _structure_basename(uniprot_id, start, end)
    domain_pdb = out_dir / f"{base_name}.pdb"
    summary_path = out_dir / f"{base_name}_summary.json"
    sasa_path = out_dir / f"{base_name}_sasa.json"
    burial_path = out_dir / f"{base_name}_burial.json"
    dist_path = out_dir / f"{base_name}_distances.parquet"
    signature_path = out_dir / f"{base_name}_signature.json"

    recompute = bool(structure_cfg.get("recompute", False))
    signature = {
        "pdb_sha256": file_sha256(base_pdb),
        "domain_start": start,
        "domain_end": end,
        "chain_id": chain_id,
        "atom_name": atom_name,
        "sasa_thresholds": sasa_thresholds,
    }
    signature_hash = dict_sha256(signature)
    if not recompute and signature_path.exists():
        existing = _read_json(signature_path)
        if (
            existing.get("signature_sha256") == signature_hash
            and domain_pdb.exists()
            and summary_path.exists()
            and sasa_path.exists()
            and burial_path.exists()
            and dist_path.exists()
        ):
            logger.info("Structure outputs up to date; skipping rebuild")
            return 0

    logger.info("Extracting core domain PDB")
    extract_domain_pdb(base_pdb, start, end, domain_pdb)
    remove_altloc(domain_pdb)
    standardize_residue_numbering(domain_pdb)

    summary = validate_structure(domain_pdb)
    _write_json(summary_path, summary)

    logger.info("Computing SASA and burial")
    sasa_map = compute_sasa(domain_pdb, chain_id=chain_id)
    burial_map = classify_burial(sasa_map, thresholds=sasa_thresholds)

    _write_json(sasa_path, {str(k): v for k, v in sasa_map.items()})
    _write_json(burial_path, {str(k): v for k, v in burial_map.items()})

    logger.info("Computing distance map")
    dist_map = compute_distance_matrix(domain_pdb, chain_id=chain_id, atom_name=atom_name)
    _write_distance_map(dist_map, dist_path)

    _write_json(signature_path, {"signature_sha256": signature_hash, "inputs": signature})

    logger.info("Structure build complete")
    return 0
