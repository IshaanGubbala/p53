from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

_MUT_RE = re.compile(r"^([A-Za-z*])(\d+)([A-Za-z*])$")


def _parse_mutation_position(token: str) -> int:
    match = _MUT_RE.match(token.strip())
    if not match:
        raise ValueError(f"Invalid mutation token: {token}")
    return int(match.group(2))


def _parse_ca_coords(pdb_path: Path) -> tuple[np.ndarray, np.ndarray]:
    coords: list[list[float]] = []
    residues: list[int] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        if line[12:16].strip() != "CA":
            continue
        resseq_text = line[22:26].strip()
        if not resseq_text.isdigit():
            continue
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        coords.append([x, y, z])
        residues.append(int(resseq_text))
    if not coords:
        raise ValueError(f"No CA atoms found in PDB: {pdb_path}")
    return np.array(coords), np.array(residues)


def _set_equal_aspect(ax, coords: np.ndarray) -> None:
    span = coords.max(axis=0) - coords.min(axis=0)
    max_range = float(span.max()) / 2.0
    mid = (coords.max(axis=0) + coords.min(axis=0)) / 2.0
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)


def render_rotation_gif(
    pdb_path: Path,
    target_mutation: str,
    rescue_mutations: Iterable[str],
    out_path: Path,
    title: str | None = None,
    frames: int = 72,
    fps: int = 18,
    dpi: int = 140,
    elev: float = 22.0,
    azim_offset: float = 0.0,
) -> None:
    coords, residues = _parse_ca_coords(pdb_path)
    target_pos = _parse_mutation_position(target_mutation)
    rescue_positions = [_parse_mutation_position(mut) for mut in rescue_mutations if mut]

    fig = plt.figure(figsize=(6.6, 6.6))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")

    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#bdbdbd", linewidth=1.0, alpha=0.6, zorder=1)
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color="#d9d9d9", s=12, alpha=0.5, zorder=2)

    target_mask = residues == target_pos
    if target_mask.any():
        target_coords = coords[target_mask]
        ax.scatter(
            target_coords[:, 0],
            target_coords[:, 1],
            target_coords[:, 2],
            color="#e63946",
            s=160,
            edgecolor="#1d3557",
            linewidth=0.6,
            alpha=0.95,
            zorder=4,
        )

    if rescue_positions:
        rescue_mask = np.isin(residues, rescue_positions)
        rescue_coords = coords[rescue_mask]
        if rescue_coords.size:
            ax.scatter(
                rescue_coords[:, 0],
                rescue_coords[:, 1],
                rescue_coords[:, 2],
                color="#2a9d8f",
                s=120,
                edgecolor="#1d3557",
                linewidth=0.5,
                alpha=0.9,
                zorder=3,
            )

    if title:
        ax.set_title(title, fontsize=10, weight="bold", pad=12)

    ax.set_axis_off()
    _set_equal_aspect(ax, coords)

    def _update(frame_idx: int):
        angle = azim_offset + (360.0 * frame_idx / max(frames, 1))
        ax.view_init(elev=elev, azim=angle)
        return []

    anim = animation.FuncAnimation(fig, _update, frames=frames, interval=1000 / max(fps, 1), blit=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = animation.PillowWriter(fps=fps)
    anim.save(out_path, writer=writer, dpi=dpi)
    plt.close(fig)
