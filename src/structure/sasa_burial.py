from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBParser, ShrakeRupley


def _select_chain(structure, chain_id: str | None):
    models = list(structure)
    if not models:
        raise RuntimeError("No models found in PDB")
    model = models[0]

    if chain_id:
        chain = model[chain_id]
        return chain

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


def compute_sasa(pdb: Path, chain_id: str | None = None) -> dict[int, float]:
    """Compute per-residue solvent accessible surface area (Å^2)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(pdb))

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    chain = _select_chain(structure, chain_id)
    sasa_map: dict[int, float] = {}
    for residue in chain:
        if residue.id[0] != " ":
            continue
        resseq = residue.id[1]
        sasa = getattr(residue, "sasa", None)
        if sasa is None:
            continue
        sasa_map[int(resseq)] = float(sasa)
    return sasa_map


def classify_burial(
    sasa_map: dict[int, float],
    thresholds: dict[str, float] | None = None,
) -> dict[int, str]:
    """Classify residues into buried/partial/exposed using SASA thresholds."""
    if thresholds is None:
        thresholds = {"buried": 10.0, "exposed": 30.0}

    buried_cutoff = thresholds.get("buried", 10.0)
    exposed_cutoff = thresholds.get("exposed", 30.0)

    burial_map: dict[int, str] = {}
    for pos, sasa in sasa_map.items():
        if sasa <= buried_cutoff:
            burial_map[pos] = "buried"
        elif sasa >= exposed_cutoff:
            burial_map[pos] = "exposed"
        else:
            burial_map[pos] = "partial"
    return burial_map
