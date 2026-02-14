#!/usr/bin/env python3
"""Validate the p53-proteoMgCAD pipeline against experimentally proven rescue mutations.

Tests whether the pipeline's scoring (oracle + PLL) correctly identifies known
second-site suppressor mutations as beneficial.  This is the ground-truth
validation that links computational predictions to wet-lab evidence.

Literature sources:
    - Nikolova et al. 2000 (EMBO J)
    - Baroni et al. 2004 (PNAS)
    - Otsuka et al. 2007 (JBC)
    - Suad et al. 2009 (Cancer Research)

Metric design notes:
    - PLL is compared at the SAME positions (cancer site only) to measure
      whether the rescue mutation improves ESM-2's confidence at the damage
      site.  Comparing at different position counts is apples-to-oranges.
    - DMS measures single-mutation effects on WT, so known rescues being LoF
      on WT background is EXPECTED — they're compensatory, not independently
      beneficial.  We test for this epistatic signature (DMS Z > 0 = confirms
      context-dependent rescue).

Usage:
    KMP_DUPLICATE_LIB_OK=TRUE python scripts/validate_literature_rescues.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Ensure project root is on the path
_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np

from p53cad.core.runtime import bootstrap_runtime

bootstrap_runtime()

import torch
from p53cad.data.dms import P53_WT, apply_mutation, parse_single_mutation
from p53cad.engine.latent import ManifoldEmbedder
from p53cad.engine.oracle import FunctionalOracle, compute_masked_marginal_pll

# ─── Known literature rescues ────────────────────────────────────────────
KNOWN_RESCUES: Dict[str, List[Dict[str, str]]] = {
    "R175H": [
        {"rescue": "N268D", "ref": "Nikolova 2000"},
        {"rescue": "H178Y", "ref": "Baroni 2004"},
        {"rescue": "N239Y", "ref": "Nikolova 2000"},
        {"rescue": "T284R", "ref": "Otsuka 2007"},
    ],
    "R248Q": [
        {"rescue": "T123A", "ref": "Baroni 2004"},
        {"rescue": "H178Y", "ref": "Baroni 2004"},
    ],
    "R249S": [
        {"rescue": "N239Y", "ref": "Nikolova 2000"},
        {"rescue": "H168R", "ref": "Baroni 2004"},
    ],
    "R273H": [
        {"rescue": "T284R", "ref": "Otsuka 2007"},
        {"rescue": "S240R", "ref": "Otsuka 2007"},
    ],
    "G245S": [
        {"rescue": "R249M", "ref": "Baroni 2004"},
        {"rescue": "T123A", "ref": "Baroni 2004"},
    ],
}


def _apply_mutations(seq: str, mutations: List[str]) -> Optional[str]:
    """Apply a series of mutations to a sequence."""
    for m in mutations:
        result = apply_mutation(seq, m)
        if result is None:
            return None
        seq = result
    return seq


def _mutated_positions(mutations: List[str]) -> List[int]:
    """Extract 1-indexed positions from mutation strings."""
    positions = []
    for m in mutations:
        parsed = parse_single_mutation(m)
        if parsed is not None:
            positions.append(parsed[1])
    return positions


def score_sequence(
    embedder: ManifoldEmbedder,
    oracle: FunctionalOracle,
    sequence: str,
) -> float:
    """Get the oracle score for a sequence."""
    from p53cad.engine.oracle import AttentionPoolingNet

    emb = embedder.get_embeddings(sequence).detach()
    with torch.no_grad():
        z, _, _ = embedder.latent_forward_ascent(emb)
        if isinstance(oracle.model, AttentionPoolingNet):
            z_oracle = z
            if z_oracle.shape[-1] != oracle.input_dim:
                z_oracle = z_oracle[:, :, :oracle.input_dim]
            return float(oracle.model(z_oracle).squeeze(-1).item())
        pooled = z.mean(dim=1)
        if pooled.shape[-1] != oracle.input_dim:
            pooled = pooled[:, :oracle.input_dim]
        return float(oracle.model(pooled).squeeze(-1).item())


def load_dms_lookup() -> Dict[Tuple[int, str], float]:
    """Load raw DMS Z-scores into a (position, amino_acid) → Z-score dict."""
    import pandas as pd

    dms_path = _root / "data" / "raw" / "p53_DMS_Giacomelli_2018.csv"
    if not dms_path.exists():
        print(f"  [WARN] DMS CSV not found at {dms_path}")
        return {}
    raw = pd.read_csv(dms_path)
    score_col = "A549_p53WT_Nutlin-3_Z-score"
    if score_col not in raw.columns:
        return {}
    lookup: Dict[Tuple[int, str], float] = {}
    for _, row in raw.iterrows():
        pos = int(row.get("Position", 0))
        aa = str(row.get("AA_variant", ""))
        z = row.get(score_col)
        if not np.isnan(z) and 1 <= pos <= len(P53_WT) and len(aa) == 1:
            lookup[(pos, aa.upper())] = float(z)
    return lookup


def main() -> None:
    print("=" * 70)
    print("  p53-proteoMgCAD  |  Literature Rescue Validation")
    print("=" * 70)
    print()

    # Load models
    print("Loading ESM-2 and oracle...")
    embedder = ManifoldEmbedder()
    oracle = FunctionalOracle(model_path=_root / "data" / "models" / "functional_oracle.pt")
    dms_lookup = load_dms_lookup()
    print(f"  DMS lookup: {len(dms_lookup)} entries")
    print()

    # Score wild-type baseline
    wt_score = score_sequence(embedder, oracle, P53_WT)
    # PLL at a few hotspot positions for reference
    hotspot_positions = [175, 248, 249, 273, 245]
    wt_pll_hotspots = compute_masked_marginal_pll(embedder, P53_WT, hotspot_positions)
    print(f"WT baseline:  oracle = {wt_score:.4f}  PLL(hotspots) = {wt_pll_hotspots:.2f}")
    print()

    print("Metrics:")
    print("  Oracle Δ    : score(cancer+rescue) - score(cancer)  [+ = rescue helps]")
    print("  PLL Δ       : PLL at cancer site in rescue seq vs cancer seq  [+ = rescue helps]")
    print("  Specificity : oracle Δ on cancer - oracle Δ on WT  [+ = context-dependent]")
    print("  DMS epistasis: rescue Z > 0 on WT = LoF alone, confirming compensatory role")
    print()

    total_tests = 0
    passed_oracle = 0
    passed_pll = 0
    passed_specificity = 0
    passed_epistatic_dms = 0

    results_table: List[Dict] = []

    for cancer_mut, rescues in KNOWN_RESCUES.items():
        cancer_seq = apply_mutation(P53_WT, cancer_mut)
        if cancer_seq is None:
            print(f"[SKIP] Could not apply {cancer_mut}")
            continue

        cancer_score = score_sequence(embedder, oracle, cancer_seq)
        cancer_positions = _mutated_positions([cancer_mut])
        # PLL at the cancer mutation site(s) only
        cancer_pll = compute_masked_marginal_pll(embedder, cancer_seq, cancer_positions)

        print(f"─── {cancer_mut} ───")
        print(f"  Cancer-only:  oracle = {cancer_score:.4f}  PLL(pos {cancer_positions}) = {cancer_pll:.4f}")

        for entry in rescues:
            rescue_label = entry["rescue"]
            ref = entry["ref"]

            # Build cancer + rescue sequence
            rescue_seq = _apply_mutations(P53_WT, [cancer_mut, rescue_label])
            if rescue_seq is None:
                print(f"  [SKIP] Could not apply {rescue_label}")
                continue

            total_tests += 1

            # 1. Oracle delta
            rescue_score = score_sequence(embedder, oracle, rescue_seq)
            oracle_delta = rescue_score - cancer_score
            oracle_pass = oracle_delta > 0

            # 2. PLL delta at the SAME positions (cancer site only)
            #    Does the rescue mutation improve ESM-2's confidence at the
            #    cancer damage site?  This is the fair comparison.
            rescue_pll = compute_masked_marginal_pll(
                embedder, rescue_seq, cancer_positions
            )
            pll_delta = rescue_pll - cancer_pll
            pll_pass = pll_delta > 0

            # Also compute PLL at the rescue site for context
            rescue_positions = _mutated_positions([rescue_label])
            rescue_site_pll = compute_masked_marginal_pll(
                embedder, rescue_seq, rescue_positions
            )

            # 3. Rescue specificity: rescue helps cancer > rescue helps WT
            wt_plus_rescue = apply_mutation(P53_WT, rescue_label)
            if wt_plus_rescue is not None:
                wt_rescue_score = score_sequence(embedder, oracle, wt_plus_rescue)
                wt_rescue_delta = wt_rescue_score - wt_score
                specificity = oracle_delta - wt_rescue_delta
                specificity_pass = specificity > 0
            else:
                wt_rescue_delta = 0.0
                specificity = 0.0
                specificity_pass = False

            # 4. DMS epistatic signature: rescue mutation LoF on WT = compensatory
            #    Known rescues are expected to have POSITIVE Z (LoF on WT)
            #    because they're beneficial only in the cancer context.
            parsed = parse_single_mutation(rescue_label)
            dms_z = None
            dms_epistatic = False
            if parsed is not None:
                _, pos, var_aa = parsed
                dms_z = dms_lookup.get((pos, var_aa))
                if dms_z is not None:
                    dms_epistatic = dms_z > 0  # LoF on WT = compensatory

            # Tally
            if oracle_pass:
                passed_oracle += 1
            if pll_pass:
                passed_pll += 1
            if specificity_pass:
                passed_specificity += 1
            if dms_epistatic:
                passed_epistatic_dms += 1

            status = "PASS" if (oracle_pass or pll_pass) else "FAIL"
            dms_str = f"{dms_z:+.2f}" if dms_z is not None else "N/A"

            print(
                f"  + {rescue_label:6s} ({ref:15s})  "
                f"Δoracle={oracle_delta:+.4f} {'✓' if oracle_pass else '✗'}  "
                f"ΔPLL@cancer={pll_delta:+.4f} {'✓' if pll_pass else '✗'}  "
                f"specificity={specificity:+.4f} {'✓' if specificity_pass else '✗'}  "
                f"DMS={dms_str} {'✓' if dms_epistatic else '✗'}  "
                f"[{status}]"
            )

            results_table.append({
                "cancer": cancer_mut,
                "rescue": rescue_label,
                "reference": ref,
                "oracle_delta": oracle_delta,
                "pll_delta": pll_delta,
                "rescue_site_pll": rescue_site_pll,
                "specificity": specificity,
                "dms_z": dms_z,
                "oracle_pass": oracle_pass,
                "pll_pass": pll_pass,
                "specificity_pass": specificity_pass,
                "dms_epistatic": dms_epistatic,
            })

        print()

    # Summary
    print("=" * 70)
    print("  VALIDATION SUMMARY")
    print("=" * 70)
    print(f"  Total rescue pairs tested:     {total_tests}")
    print(f"  Oracle Δscore > 0:             {passed_oracle}/{total_tests} ({100*passed_oracle/max(total_tests,1):.0f}%)")
    print(f"  PLL Δ@cancer_site > 0:         {passed_pll}/{total_tests} ({100*passed_pll/max(total_tests,1):.0f}%)")
    print(f"  Rescue specificity > 0:        {passed_specificity}/{total_tests} ({100*passed_specificity/max(total_tests,1):.0f}%)")
    print(f"  DMS epistatic (Z>0 on WT):     {passed_epistatic_dms}/{total_tests} ({100*passed_epistatic_dms/max(total_tests,1):.0f}%)")
    print()

    any_signal = sum(
        1 for r in results_table
        if r["oracle_pass"] or r["pll_pass"]
    )
    print(f"  Any signal (oracle OR PLL):    {any_signal}/{total_tests} ({100*any_signal/max(total_tests,1):.0f}%)")
    print()

    # Diagnostic notes
    print("─── Diagnostic Notes ───")
    oracle_deltas = [abs(r["oracle_delta"]) for r in results_table]
    max_delta = max(oracle_deltas) if oracle_deltas else 0
    if max_delta < 0.01:
        print(f"  ⚠ Oracle deltas are very small (max |Δ| = {max_delta:.6f}).")
        print("    The 8M-param ESM-2 + mean pooling produces nearly identical")
        print("    embeddings for single-mutation variants.  A larger ESM-2 model")
        print("    (e.g. esm2_t33_650M) or per-position scoring may help.")
    if passed_epistatic_dms == total_tests:
        print(f"  ✓ All known rescues are LoF on WT (DMS Z > 0), confirming they")
        print("    are compensatory mutations — beneficial only in cancer context.")
    print()

    # Exit code: fail if no signal at all
    if total_tests > 0 and any_signal / total_tests < 0.3:
        print("  [WARN] Less than 30% showing any signal — pipeline needs improvement")
        sys.exit(1)
    else:
        print("  Pipeline validation: ACCEPTABLE")
        sys.exit(0)


if __name__ == "__main__":
    main()
