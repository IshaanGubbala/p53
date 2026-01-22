from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Iterable

from Bio.PDB import PDBIO, PDBParser, Select


def extract_domain_pdb(in_pdb: Path, start: int, end: int, out_pdb: Path) -> None:
    """Extract a residue range from a PDB into a new file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(in_pdb))

    class DomainSelect(Select):
        def accept_residue(self, residue) -> int:
            hetfield, resseq, _ = residue.id
            if hetfield.strip() and residue.resname != "ZN":
                return 0
            if residue.resname == "HOH":
                return 0
            if isinstance(resseq, int) and start <= resseq <= end:
                return 1
            return 0

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), select=DomainSelect())


def remove_altloc(pdb_path: Path) -> None:
    """Remove alternate location atoms, keeping blank or 'A'."""
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()
    cleaned: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
            altloc = line[16]
            if altloc not in {" ", "A"}:
                continue
            if altloc == "A":
                line = line[:16] + " " + line[17:]
        cleaned.append(line)

    _replace_file(pdb_path, cleaned)


def standardize_residue_numbering(pdb_path: Path) -> None:
    """Strip insertion codes while preserving residue numbering."""
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()
    updated: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
            if line[26] != " ":
                line = line[:26] + " " + line[27:]
        updated.append(line)

    _replace_file(pdb_path, updated)


def validate_structure(pdb_path: Path) -> dict[str, object]:
    """Summarize chains, residue counts, and gaps from a PDB file."""
    chain_residues: dict[str, set[int]] = {}
    atom_count = 0

    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_count += 1
        if len(line) < 26:
            continue
        chain_id = line[21].strip() or "?"
        resseq_text = line[22:26].strip()
        if not resseq_text.isdigit():
            continue
        resseq = int(resseq_text)
        chain_residues.setdefault(chain_id, set()).add(resseq)

    chains = sorted(chain_residues.keys())
    residues = {resseq for values in chain_residues.values() for resseq in values}
    residue_numbers = sorted(residues)

    missing_total = 0
    for chain_id, values in chain_residues.items():
        ordered = sorted(values)
        if not ordered:
            continue
        prev = ordered[0]
        for resseq in ordered[1:]:
            gap = resseq - prev
            if gap > 1:
                missing_total += gap - 1
            prev = resseq

    return {
        "chains": chains,
        "residues": len(residues),
        "min_residue": min(residue_numbers) if residue_numbers else None,
        "max_residue": max(residue_numbers) if residue_numbers else None,
        "missing_residues": missing_total,
        "atom_count": atom_count,
    }


def _replace_file(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        delete=False,
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
    ) as tmp:
        for line in lines:
            tmp.write(line + "\n")
        tmp_path = Path(tmp.name)

    try:
        tmp_path.replace(path)
    except OSError:
        path.write_text(tmp_path.read_text(encoding="utf-8"), encoding="utf-8")
        tmp_path.unlink(missing_ok=True)
