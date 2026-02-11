#!/usr/bin/env python3
"""Non-destructive workspace organizer for generated artifacts.

Default mode is dry-run. Use --apply to perform moves.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from typing import Iterable, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
LEGACY_DIR = ROOT / "_legacy"

CAMPAIGNS_DIR = ROOT / "data" / "campaigns"
EXPORTS_DIR = ROOT / "data" / "processed" / "exports"
ARCHIVE_BASE = ROOT / "data" / "processed" / "archive"

ROOT_ALLOWLIST = {
    "results_raw.csv",
    "generative_design_candidates.csv",
    "candidates.csv",
    "candidates_analyzed.csv",
    "top30.csv",
    "top30.parquet",
    "campaign_summary.md",
}

PROCESSED_EXPORT_KEYWORDS = (
    "top30",
    "shortlist",
    "presentation",
    "report",
)

PROCESSED_MOVE_SUFFIXES = {".csv", ".json", ".parquet", ".md", ".txt"}


def _safe_destination(dst: Path) -> Path:
    if not dst.exists():
        return dst
    stem = dst.stem
    suffix = dst.suffix
    parent = dst.parent
    idx = 1
    while True:
        candidate = parent / f"{stem}_{idx}{suffix}"
        if not candidate.exists():
            return candidate
        idx += 1


def _in_legacy(path: Path) -> bool:
    try:
        path.resolve().relative_to(LEGACY_DIR.resolve())
        return True
    except Exception:
        return False


def _collect_root_moves(archive_dir: Path) -> List[Tuple[Path, Path]]:
    moves: List[Tuple[Path, Path]] = []
    for name in sorted(ROOT_ALLOWLIST):
        src = ROOT / name
        if src.exists() and src.is_file() and not _in_legacy(src):
            dst = _safe_destination(archive_dir / src.name)
            moves.append((src, dst))
    return moves


def _collect_processed_moves(archive_dir: Path) -> List[Tuple[Path, Path]]:
    moves: List[Tuple[Path, Path]] = []
    processed_dir = ROOT / "data" / "processed"
    if not processed_dir.exists():
        return moves

    blocked_dirs = {
        (processed_dir / "exports").resolve(),
        (processed_dir / "archive").resolve(),
    }

    for src in sorted(processed_dir.rglob("*")):
        if not src.is_file():
            continue
        if _in_legacy(src):
            continue
        src_resolved = src.resolve()
        if any(str(src_resolved).startswith(str(blocked)) for blocked in blocked_dirs):
            continue
        if src.suffix.lower() not in PROCESSED_MOVE_SUFFIXES:
            continue

        lower_name = src.name.lower()
        if any(k in lower_name for k in PROCESSED_EXPORT_KEYWORDS):
            dst = _safe_destination(EXPORTS_DIR / src.name)
        else:
            dst = _safe_destination(archive_dir / src.name)
        moves.append((src, dst))
    return moves


def _print_moves(moves: Iterable[Tuple[Path, Path]]) -> None:
    for src, dst in moves:
        print(f"MOVE {src.relative_to(ROOT)} -> {dst.relative_to(ROOT)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Organize generated campaign/result files.")
    parser.add_argument("--apply", action="store_true", help="Perform moves (default is dry-run).")
    args = parser.parse_args()

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    archive_dir = ARCHIVE_BASE / timestamp

    # Ensure directory layout exists regardless of dry-run.
    CAMPAIGNS_DIR.mkdir(parents=True, exist_ok=True)
    EXPORTS_DIR.mkdir(parents=True, exist_ok=True)
    archive_dir.mkdir(parents=True, exist_ok=True)

    moves = _collect_root_moves(archive_dir)
    moves.extend(_collect_processed_moves(archive_dir))

    if not moves:
        print("No generated files matched cleanup rules. Workspace already organized.")
        return

    print(f"Planned moves: {len(moves)}")
    _print_moves(moves)

    if not args.apply:
        print("Dry run complete. Re-run with --apply to execute moves.")
        return

    for src, dst in moves:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))

    print(f"Applied {len(moves)} moves.")


if __name__ == "__main__":
    main()
