from __future__ import annotations

from typing import Iterable

from src.core.hashing import json_sha256


def format_mutation(pos: int, ref: str, alt: str) -> str:
    return f"{ref.upper()}{int(pos)}{alt.upper()}"


def canonicalize_set(muts: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({mut.strip() for mut in muts if mut}))


def mutation_set_id(muts: tuple[str, ...]) -> str:
    return json_sha256({"mutations": muts})
