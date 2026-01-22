from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Callable

from src.core.config import load_configs, normalize_paths
from src.core.logging import setup_logging


CONFIG_NAMES = ["p53", "paths", "scoring", "optimizer", "report"]
RUNNERS = {
    "fetch": ("experiments.run_fetch", "run"),
    "build-dataset": ("experiments.run_build_dataset", "run"),
    "build-structure": ("experiments.run_build_structure", "run"),
    "score-variants": ("experiments.run_score_variants", "run"),
    "design": ("experiments.run_design_rescues", "run"),
    "report": ("visualization.create_visualizations", "run"),
    "export-benchmark": ("experiments.run_export_benchmark_repro", "run"),
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="p53 StabiliMut pipeline")
    parser.add_argument(
        "--config-dir",
        default=str(repo_root() / "configs"),
        help="Path to configuration directory",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (e.g., INFO, DEBUG)",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch = subparsers.add_parser("fetch", help="Download raw datasets")
    fetch.add_argument("--all", action="store_true", help="Fetch all sources")
    fetch.add_argument("--uniprot", action="store_true", help="Fetch UniProt entry")
    fetch.add_argument("--alphafold", action="store_true", help="Fetch AlphaFold PDB")
    fetch.add_argument("--clinvar", action="store_true", help="Fetch ClinVar release")
    fetch.add_argument("--nci", action="store_true", help="Fetch NCI TP53 tables")

    subparsers.add_parser("build-dataset", help="Parse and map variants")

    subparsers.add_parser("build-structure", help="Prepare structure geometry outputs")

    score = subparsers.add_parser("score-variants", help="Score variants with EvoEF2")
    score.add_argument("--engine", choices=["evoef2"], help="Scoring engine override")
    score.add_argument("--parallel", type=int, help="Number of parallel workers")

    design = subparsers.add_parser("design", help="Design rescue mutation sets")
    design.add_argument("--targets", nargs="+", help="Target destabilizing mutations (e.g., R175H)")
    design.add_argument("--max-muts", type=int, help="Maximum number of rescue mutations")
    design.add_argument("--beam-width", type=int, help="Beam width for candidate expansion")
    design.add_argument("--depth", type=int, help="Search depth for rescue sets")
    design.add_argument("--parallel", type=int, help="Number of parallel workers")

    subparsers.add_parser("report", help="Generate figures and summary tables")

    subparsers.add_parser("export-benchmark", help="Export benchmark reproducibility artifacts")

    return parser


def _load_runner(module_name: str, func_name: str) -> Callable:
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name == module_name:
            raise RuntimeError(f"Runner module missing: {module_name}") from exc
        raise

    runner = getattr(module, func_name, None)
    if runner is None:
        raise RuntimeError(f"Runner function missing: {module_name}.{func_name}")
    return runner


def main(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    setup_logging(args.log_level)

    config_dir = Path(args.config_dir)
    configs = load_configs(config_dir, CONFIG_NAMES)
    configs["paths"] = normalize_paths(configs.get("paths", {}), repo_root())

    module_name, func_name = RUNNERS[args.command]
    try:
        runner = _load_runner(module_name, func_name)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    result = runner(args, configs)
    return int(result or 0)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
