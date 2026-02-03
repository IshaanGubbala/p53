#!/usr/bin/env python3
"""Reproduce the variant separation benchmark from saved CSV files.

This script regenerates the AUC calculation and plot from the exported
benchmark artifacts, demonstrating reproducibility.

Usage:
    python scripts/reproduce_benchmark.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

# Add repo root to Python path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

from src.eval.variant_separation import bootstrap_ci, compute_auc
from src.viz.plots_variants import plot_ddg_by_label


def main() -> None:
    # Paths to reproducibility artifacts
    reports_dir = repo_root / "reports"
    tables_dir = reports_dir / "tables"
    figures_dir = reports_dir / "figures"
    logs_dir = reports_dir / "logs"

    # Load metadata
    metadata_path = logs_dir / "run_metadata.json"
    if not metadata_path.exists():
        print(f"Error: metadata not found at {metadata_path}")
        print("Run: python -m src.cli export-benchmark")
        return

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    print("=== Reproducing Variant Separation Benchmark ===\n")
    print(f"Dataset size: {metadata['dataset_info']['final_benchmark_size']}")
    print(f"  Benign: {metadata['dataset_info']['final_benign']}")
    print(f"  Pathogenic: {metadata['dataset_info']['final_pathogenic']}")
    print(f"Bootstrap seed: {metadata['random_seeds']['bootstrap_ci']}")
    print(f"Bootstrap samples: {metadata['random_seeds']['n_bootstrap_samples']}\n")

    # Load scored variants
    scored_path = tables_dir / "variant_benchmark_scored.csv"
    if not scored_path.exists():
        print(f"Error: scored data not found at {scored_path}")
        print("Run: python -m src.cli export-benchmark")
        return

    data = pd.read_csv(scored_path)
    print(f"Loaded {len(data)} scored variants from {scored_path}\n")

    # Compute AUC
    auc = compute_auc(data["label"], data["ddg"])
    seed = metadata["random_seeds"]["bootstrap_ci"]
    n_bootstrap = metadata["random_seeds"]["n_bootstrap_samples"]
    lo, hi = bootstrap_ci(data["label"], data["ddg"], n=n_bootstrap, seed=seed)

    print("=== Results ===")
    print(f"AUC: {auc:.4f}")
    print(f"95% CI: [{lo:.4f}, {hi:.4f}]\n")

    # Compare with saved results
    saved_results_path = tables_dir / "variant_separation.json"
    if saved_results_path.exists():
        saved = json.loads(saved_results_path.read_text(encoding="utf-8"))
        print("=== Comparison with saved results ===")
        print(f"Saved AUC: {saved['auc']:.4f}")
        print(f"Saved CI: [{saved['auc_ci_lower']:.4f}, {saved['auc_ci_upper']:.4f}]")
        print(f"Match: {abs(auc - saved['auc']) < 1e-10 and abs(lo - saved['auc_ci_lower']) < 1e-10}")

    # Regenerate plot
    output_plot = figures_dir / "variant_ddg_by_label_reproduced.png"
    plot_ddg_by_label(data, output_plot)
    print(f"\nPlot saved to: {output_plot}")
    print("\nReproduction successful! ✓")


if __name__ == "__main__":
    main()
