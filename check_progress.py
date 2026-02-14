#!/usr/bin/env python
"""Live campaign progress monitor — refreshes every 15 seconds."""

import json
import os
import time
import sys
from collections import Counter, defaultdict
from datetime import datetime

CAMPAIGN_DIR = "data/campaigns"
REFRESH = 15  # seconds


def load_jsonl(path):
    items = []
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    items.append(json.loads(line))
    return items


def fmt_score(s):
    return f"{s:+.4f}" if s is not None else "N/A"


def run():
    while True:
        os.system("cls" if os.name == "nt" else "clear")
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # --- Find latest campaign ------------------------------------------------
        if not os.path.isdir(CAMPAIGN_DIR):
            print(f"[{now}]  Waiting for campaign to start...")
            print(f"         (no {CAMPAIGN_DIR}/ directory found)")
            time.sleep(REFRESH)
            continue

        campaigns = sorted(
            [d for d in os.listdir(CAMPAIGN_DIR) if d.startswith("campaign_")]
        )
        if not campaigns:
            print(f"[{now}]  No campaigns found yet.")
            time.sleep(REFRESH)
            continue

        latest = campaigns[-1]
        manifest_path = f"{CAMPAIGN_DIR}/{latest}/manifest.json"
        cp = f"{CAMPAIGN_DIR}/{latest}/checkpoints"

        # --- Manifest ------------------------------------------------------------
        with open(manifest_path) as f:
            m = json.load(f)

        status = m.get("status", "unknown")
        budget = m.get("config", {}).get("budget", "?")
        seed = m.get("config", {}).get("seed", "?")
        pass_a_steps = m.get("config", {}).get("pass_a", {}).get("steps", "?")
        pass_a_restarts = m.get("config", {}).get("pass_a", {}).get("restarts", "?")
        pass_b_steps = m.get("config", {}).get("pass_b", {}).get("steps", "?")
        pass_b_restarts = m.get("config", {}).get("pass_b", {}).get("restarts", "?")
        shortlist = m.get("config", {}).get("shortlist", "?")

        # --- Header --------------------------------------------------------------
        print("=" * 70)
        print("  p53 RESCUE CAMPAIGN — LIVE MONITOR")
        print("=" * 70)
        print(f"  Time        : {now}")
        print(f"  Campaign    : {latest}")
        print(f"  Status      : {status.upper()}")
        print(f"  Budget      : {budget}")
        print(f"  Seed        : {seed}")
        print(f"  Pass A      : {pass_a_steps} steps x {pass_a_restarts} restart(s)")
        print(f"  Pass B      : {pass_b_steps} steps x {pass_b_restarts} restart(s)")
        print(f"  Shortlist   : top {shortlist} candidates")
        print("-" * 70)

        # --- Scenarios -----------------------------------------------------------
        scenarios = load_jsonl(f"{cp}/scenario_status.jsonl")
        n_done = len(scenarios)
        n_total = 108
        pct = 100 * n_done / n_total if n_total else 0

        bar_width = 40
        filled = int(bar_width * pct / 100)
        bar = "#" * filled + "-" * (bar_width - filled)

        print(f"\n  PASS A SCREENING")
        print(f"  Progress    : [{bar}] {pct:.1f}%")
        print(f"  Scenarios   : {n_done} / {n_total}")

        if scenarios:
            print(f"  Latest done : {scenarios[-1].get('scenario_id', '?')}")

        # --- Per-hotspot breakdown -----------------------------------------------
        scenarios_by_target = defaultdict(list)
        for s in scenarios:
            sid = s.get("scenario_id", "")
            parts = sid.split("__")
            target = parts[0] if parts else "?"
            delivery = parts[1] if len(parts) > 1 else "?"
            scenarios_by_target[target].append(delivery)

        all_hotspots = [
            "R175H", "R248Q", "R248W", "R273H",
            "R273C", "G245S", "Y220C", "V157F",
        ]

        print(f"\n  HOTSPOT BREAKDOWN")
        print(f"  {'Hotspot':<12} {'gene':>8} {'mrna':>8} {'protein':>8} {'pairs':>8}")
        print(f"  {'-'*12} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
        for h in all_hotspots:
            deliveries = scenarios_by_target.get(h, [])
            gene = "done" if "gene_therapy" in deliveries else "..."
            mrna = "done" if "mrna_therapy" in deliveries else "..."
            prot = "done" if "protein_therapy" in deliveries else "..."
            pair_count = sum(1 for d in deliveries if "+" in h or d not in
                            ["gene_therapy", "mrna_therapy", "protein_therapy"])
            print(f"  {h:<12} {gene:>8} {mrna:>8} {prot:>8}")

        # Count pair scenarios
        pair_scenarios = [s for s in scenarios if "+" in s.get("scenario_id", "")]
        if pair_scenarios:
            print(f"\n  Pair scenarios completed: {len(pair_scenarios)}")
            for ps in pair_scenarios[-5:]:
                print(f"    {ps.get('scenario_id', '?')}")

        # --- Candidates ----------------------------------------------------------
        candidates = load_jsonl(f"{cp}/candidates.jsonl")
        if candidates:
            scores = [c.get("score", 0) for c in candidates]
            print(f"\n  CANDIDATES")
            print(f"  Total       : {len(candidates)}")
            print(f"  Best score  : {fmt_score(max(scores))}")
            print(f"  Worst score : {fmt_score(min(scores))}")
            print(f"  Mean score  : {fmt_score(sum(scores) / len(scores))}")
            print(f"  Median      : {fmt_score(sorted(scores)[len(scores) // 2])}")

            # Per-target stats
            target_scores = defaultdict(list)
            for c in candidates:
                target_scores[c.get("target_label", "?")].append(c.get("score", 0))

            print(f"\n  PER-HOTSPOT SCORES")
            print(f"  {'Hotspot':<12} {'Count':>6} {'Best':>10} {'Mean':>10}")
            print(f"  {'-'*12} {'-'*6} {'-'*10} {'-'*10}")
            for t in sorted(target_scores.keys()):
                sc = target_scores[t]
                print(f"  {t:<12} {len(sc):>6} {fmt_score(max(sc)):>10} {fmt_score(sum(sc)/len(sc)):>10}")

            # Top 10 candidates
            top = sorted(candidates, key=lambda c: c.get("score", 0), reverse=True)[:10]
            print(f"\n  TOP 10 CANDIDATES")
            print(f"  {'Rank':<6} {'Target':<12} {'Delivery':<16} {'Score':>10} {'Mutations':>10}")
            print(f"  {'-'*6} {'-'*12} {'-'*16} {'-'*10} {'-'*10}")
            for i, c in enumerate(top, 1):
                target = c.get("target_label", "?")
                delivery = c.get("delivery", c.get("scenario_id", "?").split("__")[-1] if "__" in c.get("scenario_id", "") else "?")
                score = c.get("score", 0)
                n_muts = c.get("n_mutations", c.get("num_mutations", "?"))
                print(f"  {i:<6} {target:<12} {delivery:<16} {fmt_score(score):>10} {str(n_muts):>10}")
        else:
            print(f"\n  CANDIDATES: none yet")

        # --- Trajectories --------------------------------------------------------
        trajectories = load_jsonl(f"{cp}/trajectories.jsonl")
        if trajectories:
            print(f"\n  OPTIMIZATION TRAJECTORIES")
            print(f"  Total entries: {len(trajectories)}")

        # --- Timing estimate -----------------------------------------------------
        if n_done >= 2:
            # Try to estimate from manifest timestamps
            first_scenario_id = scenarios[0].get("scenario_id", "")
            elapsed_since_start = None

            # Use file modification time as rough proxy
            status_file = f"{cp}/scenario_status.jsonl"
            manifest_mtime = os.path.getmtime(manifest_path)
            status_mtime = os.path.getmtime(status_file)
            elapsed_min = (status_mtime - manifest_mtime) / 60

            if elapsed_min > 0:
                rate = elapsed_min / n_done
                remaining = (n_total - n_done) * rate
                print(f"\n  TIMING")
                print(f"  Elapsed     : {elapsed_min:.1f} min")
                print(f"  Rate        : {rate:.1f} min/scenario")
                print(f"  ETA (Pass A): {remaining:.0f} min ({remaining/60:.1f} hrs)")

        # --- Pass B check --------------------------------------------------------
        pass_b_file = f"{cp}/pass_b_results.jsonl"
        if os.path.exists(pass_b_file):
            pass_b = load_jsonl(pass_b_file)
            print(f"\n  PASS B DEEP REFINEMENT")
            print(f"  Results: {len(pass_b)}")
            if pass_b:
                b_scores = [r.get("score", 0) for r in pass_b]
                print(f"  Best Pass B score: {fmt_score(max(b_scores))}")

        # --- Final results -------------------------------------------------------
        results_file = f"{CAMPAIGN_DIR}/{latest}/results.parquet"
        results_csv = f"{CAMPAIGN_DIR}/{latest}/results.csv"
        if os.path.exists(results_file) or os.path.exists(results_csv):
            print(f"\n  *** FINAL RESULTS AVAILABLE ***")
            if os.path.exists(results_file):
                print(f"  Parquet: {results_file}")
            if os.path.exists(results_csv):
                print(f"  CSV: {results_csv}")

        if status == "completed":
            print(f"\n  *** CAMPAIGN COMPLETED ***")
            print(f"\n  Refreshing stopped.")
            break

        print(f"\n{'=' * 70}")
        print(f"  Refreshing every {REFRESH}s ... Press Ctrl+C to stop.")
        print(f"{'=' * 70}")

        time.sleep(REFRESH)


if __name__ == "__main__":
    try:
        run()
    except KeyboardInterrupt:
        print("\n\nMonitor stopped.")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
