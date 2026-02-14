#!/usr/bin/env python3
"""Full campaign with all v2 improvements enabled.

Runs the complete p53 rescue candidate campaign using:
- 650M ESM-2 model (Fix A)
- Attention-pooling oracle if trained (Fix B)
- Double-mutant DMS epistasis data (Fix D)
- Conditional DMS rescue scoring (Fix E)
- Pareto multi-objective ranking (Fix F)
- DBD structural confidence in loop (Fix G)
- Autoregressive trial sampling (Fix I)
- Full post-campaign analysis with ESMFold + MD scripts

Usage:
    KMP_DUPLICATE_LIB_OK=TRUE python scripts/run_full_campaign.py [--budget BUDGET]
"""

import argparse
import json
import os
import sys
import time

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

# Ensure project root is on path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from p53cad.core.runtime import bootstrap_runtime

bootstrap_runtime()


def main():
    parser = argparse.ArgumentParser(description="Run full p53 rescue campaign")
    parser.add_argument("--budget", default="medium", choices=["fast", "medium", "high"],
                        help="Compute budget (fast=~5min, medium=~80min, high=~4hr)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--shortlist-n", type=int, default=30, help="Number of shortlist candidates")
    parser.add_argument("--no-pairs", action="store_true", help="Skip multi-hotspot combos")
    parser.add_argument("--hotspots", nargs="+", default=None,
                        help="Specific hotspots (default: all BIG8)")
    args = parser.parse_args()

    from p53cad.engine.campaign import CampaignRunner

    print(f"\n{'='*70}")
    print(f"  p53 Rescue Campaign — Full v2 Pipeline")
    print(f"{'='*70}")
    print(f"  Budget: {args.budget}")
    print(f"  Seed: {args.seed}")
    print(f"  Include pairs: {not args.no_pairs}")
    print(f"  Shortlist: {args.shortlist_n}")
    if args.hotspots:
        print(f"  Hotspots: {args.hotspots}")
    print()

    t0 = time.time()
    runner = CampaignRunner()

    # Report model info
    if runner.embedder is not None:
        print(f"  ESM-2 model: {runner.embedder.model_name}")
        print(f"  Hidden dim: {runner.embedder.hidden_size}")
    if runner.oracle is not None:
        print(f"  Oracle input_dim: {runner.oracle.input_dim}")
        arch = "attention_pooling" if hasattr(runner.oracle.model, "attn") else "legacy_mlp"
        print(f"  Oracle arch: {arch}")
    print(f"  Pairwise DMS entries: {len(runner._pairwise_dms)}")
    print(f"  Contact map entries: {len(runner._wt_contacts)}")
    print()

    run_kwargs = dict(
        budget=args.budget,
        seed=args.seed,
        include_pairs=not args.no_pairs,
        shortlist_n=args.shortlist_n,
        with_clinical=True,
    )
    if args.hotspots:
        run_kwargs["hotspots"] = args.hotspots

    result = runner.run(**run_kwargs)
    elapsed = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  CAMPAIGN COMPLETE — {elapsed/60:.1f} min")
    print(f"{'='*70}")
    print(f"  Run ID: {result['run_id']}")
    print(f"  Scenarios: {result['n_scenarios']}")
    print(f"  Candidates: {result['n_candidates']}")
    print(f"  Shortlist: {result['n_shortlist']}")
    print(f"  Run dir: {result['run_dir']}")

    # ── Detailed analysis ──
    import numpy as np
    import pandas as pd

    cand = pd.read_parquet(os.path.join(result["run_dir"], "candidates.parquet"))
    deep = cand[cand["pass_name"] == "deep"].copy() if "pass_name" in cand.columns else cand

    print(f"\n  Deep candidates: {len(deep)}")

    # Target retention
    retains = 0
    for _, row in deep.iterrows():
        targets = json.loads(str(row.get("targets_json", "[]")))
        muts = json.loads(str(row.get("mutations_json", "[]")))
        if set(targets).issubset(set(muts)):
            retains += 1
    print(f"  Target retention: {retains}/{len(deep)} ({100*retains/max(len(deep),1):.1f}%)")

    # Pareto ranking
    if "pareto_rank" in deep.columns:
        rank1 = (deep["pareto_rank"] == 1).sum()
        print(f"  Pareto rank-1 candidates: {rank1}")
        max_rank = int(deep["pareto_rank"].replace([np.inf], np.nan).dropna().max())
        print(f"  Pareto fronts: {max_rank}")

    # DMS quality
    if "rescue_dms_mean" in deep.columns:
        valid_dms = deep[deep["rescue_dms_mean"].notna() & (deep["rescue_dms_mean"] != 0)]
        print(f"\n  Candidates with DMS data: {len(valid_dms)}")
        if len(valid_dms):
            print(f"  Mean rescue DMS Z: {valid_dms['rescue_dms_mean'].mean():.3f}")
            func = (valid_dms["rescue_dms_mean"] < 0).sum()
            print(f"  Functional rescues (Z<0): {func}/{len(valid_dms)} ({100*func/len(valid_dms):.0f}%)")

    # Evidence scores
    if "score" in deep.columns:
        print(f"\n  Score stats:")
        print(f"    Mean: {deep['score'].mean():.4f}")
        print(f"    Max:  {deep['score'].max():.4f}")
        print(f"    Min:  {deep['score'].min():.4f}")

    # Shortlist analysis
    from p53cad.results.schema import select_presentation_shortlist

    top = select_presentation_shortlist(deep, top_n=args.shortlist_n)
    print(f"\n  Shortlist: {len(top)}")

    if len(top):
        # Top candidates table
        print(f"\n  {'#':<4} {'Target':<20} {'Rescue':<25} {'Oracle':<8} {'DMS-Z':<8} {'Pareto':<8} {'Delivery':<16}")
        print(f"  {'-'*4} {'-'*20} {'-'*25} {'-'*8} {'-'*8} {'-'*8} {'-'*16}")
        for _, row in top.head(20).iterrows():
            targets = json.loads(str(row.get("targets_json", "[]")))
            muts = json.loads(str(row.get("mutations_json", "[]")))
            rescue = [m for m in muts if m not in targets]
            dms_z = row.get("rescue_dms_mean", None)
            dms_str = f"{dms_z:+.2f}" if dms_z and dms_z != 0 else "N/A"
            pr = row.get("pareto_rank", "?")
            pr_str = f"{int(pr)}" if pr and pr != np.inf else "?"
            print(
                f"  {int(row.get('presentation_rank', 0)):<4} "
                f"{row['target_label']:<20} "
                f"{'+'.join(rescue):<25} "
                f"{float(row['score']):<8.3f} "
                f"{dms_str:<8} "
                f"{pr_str:<8} "
                f"{row['delivery_method']:<16}"
            )

        # Delivery distribution
        delivery_counts = top["delivery_method"].value_counts()
        print(f"\n  Delivery distribution:")
        for d, c in delivery_counts.items():
            print(f"    {d}: {c}")

        # Target diversity
        unique_targets = top["target_label"].nunique()
        print(f"\n  Unique target combos: {unique_targets}")

    print(f"\n{'='*70}")
    print(f"  Artifacts: {result['run_dir']}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
