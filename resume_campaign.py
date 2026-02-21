#!/usr/bin/env python3
"""Resume a specific campaign run with multicore processing enabled."""

import argparse
import os
import sys
import time

# Enable multicore CPU processing for PyTorch/NumPy/MKL
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = str(os.cpu_count() or 8)
os.environ["MKL_NUM_THREADS"] = str(os.cpu_count() or 8)
os.environ["NUMEXPR_NUM_THREADS"] = str(os.cpu_count() or 8)
os.environ["OPENBLAS_NUM_THREADS"] = str(os.cpu_count() or 8)

# Ensure project root is on path
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)

from p53cad.core.runtime import bootstrap_runtime
from p53cad.core.logging import setup_logging

bootstrap_runtime()
setup_logging()


def main():
    parser = argparse.ArgumentParser(description="Resume a specific p53 campaign run")
    parser.add_argument("--run-id", required=True, help="Run ID to resume (e.g., campaign_20260216_090733)")
    args = parser.parse_args()

    from p53cad.engine.campaign import CampaignRunner

    # Enable PyTorch multicore processing
    import torch
    num_threads = os.cpu_count() or 8
    torch.set_num_threads(num_threads)
    torch.set_num_interop_threads(num_threads)

    print(f"\n{'='*70}")
    print(f"  p53 Rescue Campaign — RESUMING RUN")
    print(f"{'='*70}")
    print(f"  Run ID: {args.run_id}")
    print(f"  PyTorch threads: {num_threads} (multicore enabled)")
    print(f"  Resume: True")
    print()

    t0 = time.time()
    runner = CampaignRunner()

    # Report model info
    if runner.embedder is not None:
        print(f"  ESM-2 model: {runner.embedder.model_name}")
    if runner.oracle is not None:
        arch = "attention_pooling" if hasattr(runner.oracle.model, "attn") else "legacy_mlp"
        print(f"  Oracle arch: {arch}")
    print()

    # Resume the specific run with multicore enabled
    result = runner.run(
        run_id=args.run_id,
        budget="medium",
        seed=42,
        include_pairs=True,
        shortlist_n=30,
        with_clinical=True,
    )

    elapsed = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  CAMPAIGN COMPLETE — {elapsed/60:.1f} min")
    print(f"{'='*70}")
    print(f"  Run ID: {result['run_id']}")
    print(f"  Candidates: {result['n_candidates']}")
    print(f"  Shortlist: {result['n_shortlist']}")
    print(f"  Run dir: {result['run_dir']}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
