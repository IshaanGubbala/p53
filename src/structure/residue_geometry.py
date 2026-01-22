from __future__ import annotations

import math
from pathlib import Path

from Bio.PDB import PDBParser


def _select_chain(structure, chain_id: str | None):
    models = list(structure)
    if not models:
        raise RuntimeError("No models found in PDB")
    model = models[0]

    if chain_id:
        return model[chain_id]

    best_chain = None
    best_count = -1
    for chain in model:
        residues = [res for res in chain if res.id[0] == " "]
        if len(residues) > best_count:
            best_chain = chain
            best_count = len(residues)

    if best_chain is None:
        raise RuntimeError("No suitable chain found in PDB")
    return best_chain


def _residue_coord(residue, atom_name: str) -> tuple[float, float, float] | None:
    if atom_name in residue:
        return tuple(residue[atom_name].coord)
    if atom_name == "CA" and "CB" in residue:
        return tuple(residue["CB"].coord)
    return None


def compute_distance_matrix(
    pdb: Path,
    chain_id: str | None = None,
    atom_name: str = "CA",
) -> dict[tuple[int, int], float]:
    """Compute pairwise distances between residue atoms."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(pdb))
    chain = _select_chain(structure, chain_id)

    coords: dict[int, tuple[float, float, float]] = {}
    for residue in chain:
        if residue.id[0] != " ":
            continue
        resseq = residue.id[1]
        coord = _residue_coord(residue, atom_name)
        if coord is None:
            continue
        coords[int(resseq)] = coord

    positions = sorted(coords.keys())
    dist_map: dict[tuple[int, int], float] = {}
    for i, pos_i in enumerate(positions):
        xi, yi, zi = coords[pos_i]
        for pos_j in positions[i + 1 :]:
            xj, yj, zj = coords[pos_j]
            dist = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            dist_map[(pos_i, pos_j)] = dist
            dist_map[(pos_j, pos_i)] = dist
    return dist_map


def min_distance_to_set(
    pos: int,
    protected: set[int],
    dist_map: dict[tuple[int, int], float],
) -> float:
    distances = []
    for prot in protected:
        if pos == prot:
            return 0.0
        dist = dist_map.get((pos, prot))
        if dist is None:
            dist = dist_map.get((prot, pos))
        if dist is not None:
            distances.append(dist)
    return min(distances) if distances else float("inf")


def risk_from_distance(distance: float, cutoff: float) -> float:
    if distance >= cutoff:
        return 0.0
    if distance <= 0.0:
        return 1.0
    return max(0.0, 1.0 - distance / cutoff)
