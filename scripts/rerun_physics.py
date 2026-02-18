#!/usr/bin/env python3
"""
Re-run physics validation for a completed campaign.

Usage (in openmm-cuda conda env):
    python scripts/rerun_physics.py --run-id campaign_20260217_060021
    python scripts/rerun_physics.py  # auto-picks most recent campaign
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import torch BEFORE any p53cad/OpenMM imports — OpenMM loads Intel OpenMP DLLs
# which conflict with fbgemm.dll if torch is imported afterwards (WinError 127).
import torch  # noqa: E402

from p53cad.core.logging import setup_logging
setup_logging()

from p53cad.core.runtime import bootstrap_runtime
bootstrap_runtime()


def find_latest_run(campaigns_dir: Path) -> str:
    """Return the run_id of the most recently modified completed campaign."""
    runs = sorted(
        [d for d in campaigns_dir.iterdir() if d.is_dir() and (d / "top30.parquet").exists()],
        key=lambda d: d.stat().st_mtime,
        reverse=True,
    )
    if not runs:
        raise RuntimeError(f"No completed campaigns found in {campaigns_dir}")
    return runs[0].name


def main():
    parser = argparse.ArgumentParser(description="Re-run physics validation on a completed campaign")
    parser.add_argument("--run-id", default=None,
                        help="Campaign run ID (default: most recent completed campaign)")
    parser.add_argument("--top-n", type=int, default=None,
                        help="Validate only the top N candidates (default: all 30)")
    args = parser.parse_args()

    campaigns_dir = project_root / "data" / "campaigns"

    run_id = args.run_id or find_latest_run(campaigns_dir)
    run_dir = campaigns_dir / run_id

    if not run_dir.exists():
        print(f"ERROR: Campaign directory not found: {run_dir}")
        sys.exit(1)

    top30_path = run_dir / "top30.parquet"
    if not top30_path.exists():
        print(f"ERROR: top30.parquet not found in {run_dir}")
        print("       Campaign must be fully complete before running physics validation.")
        sys.exit(1)

    import pandas as pd
    top_df = pd.read_parquet(top30_path)
    print(f"\n{'='*65}")
    print(f"  Physics Validation Re-run")
    print(f"{'='*65}")
    print(f"  Campaign : {run_id}")
    print(f"  Candidates: {len(top_df)}")
    top_n = args.top_n or len(top_df)
    print(f"  Validating: top {top_n}")

    # Device selection
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"  Device   : {device}")
    if device == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        print(f"  TF32     : enabled")
    print()

    # Build cancer sequences for DDG reference
    from p53cad.data.dms import P53_WT, apply_mutation
    cancer_sequences = {}
    if "target_label" in top_df.columns:
        for tl in top_df["target_label"].unique():
            target_mut = str(tl).split("+")[0].strip()
            cancer_seq = apply_mutation(P53_WT, target_mut)
            if cancer_seq:
                cancer_sequences[str(tl)] = cancer_seq
    print(f"  Cancer reference sequences: {len(cancer_sequences)}")

    from p53cad.engine.physics_validation import PhysicsValidationPipeline

    pipeline = PhysicsValidationPipeline(
        device=device,
        cache_dir=run_dir / "esmfold_cache",
    )

    t0 = time.time()
    print(f"\n  Starting validation — ESMFold loads once (~1-2 min first time)...\n")

    report = pipeline.validate_campaign(
        run_id=run_id,
        top_df=top_df,
        output_dir=run_dir,
        wt_sequence=P53_WT,
        cancer_sequences=cancer_sequences,
        tier1_top_n=top_n,
        tier2_top_n=0,      # skip MD (expensive)
        skip_md=True,
        skip_esmfold=False,
        skip_energy=False,
        skip_dna=False,
    )

    elapsed = time.time() - t0

    # Overwrite the previous (empty) physics_validation.json
    val_path = run_dir / "physics_validation.json"
    val_path.write_text(json.dumps(report.to_dict(), indent=2), encoding="utf-8")

    print(f"\n{'='*65}")
    print(f"  Physics Validation Complete — {elapsed/60:.1f} min")
    print(f"{'='*65}")
    print(f"  Candidates validated : {report.n_candidates}")
    print(f"  Results saved to     : {val_path}")
    print()

    # Print per-candidate summary (report.candidates is a list of dicts)
    candidates = report.candidates
    candidates_sorted = sorted(candidates, key=lambda r: r.get("overall_physics_score", 0), reverse=True)
    print(f"  {'Rank':<5} {'Target':<22} {'pLDDT':<7} {'Energy':<14} {'Score':<7} {'Verdict'}")
    print(f"  {'─'*4} {'─'*22} {'─'*7} {'─'*14} {'─'*7} {'─'*12}")
    for r in candidates_sorted:
        esm = r.get("esmfold") or {}
        eng = r.get("energy") or {}
        plddt_str = f"{esm['mean_plddt']:.1f}" if esm.get("mean_plddt") is not None else "N/A"
        energy_str = f"{eng['potential_energy_kcal']:.0f} kcal" if eng.get("potential_energy_kcal") is not None else "N/A"
        print(f"  {r.get('candidate_rank', '?'):<5} {r.get('target_label', '?'):<22} "
              f"{plddt_str:<7} {energy_str:<14} "
              f"{r.get('overall_physics_score', 0):<7.1f} {r.get('verdict', '?')}")

    print(f"\n  View results: python -m streamlit run scripts/visualize_results.py")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
