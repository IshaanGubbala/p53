from __future__ import annotations

from typing import Iterable


def get_core_domain(p53_cfg: dict) -> tuple[int, int]:
    core = p53_cfg.get("core_domain", {})
    start = int(core.get("start", 1))
    end = int(core.get("end", 9999))
    if start > end:
        raise ValueError("core_domain start must be <= end")
    return start, end


def get_designable_region(
    p53_cfg: dict,
    protected: set[int],
    residues: Iterable[int],
) -> list[int]:
    start, end = get_core_domain(p53_cfg)
    positions = [pos for pos in residues if start <= pos <= end and pos not in protected]
    return sorted(positions)
