#!/usr/bin/env python3
"""
Calculate rescue-drug synergy scores.

Usage:
    python src/docking/calculate_synergy.py

Requires:
    - Docking results from ensemble docking
    - Functional scores from Initiative 1
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List


class SynergyCalculator:
    """Calculate synergy scores for rescue-drug pairs."""

    def __init__(self):
        self.docking_dir = Path("Data/processed/docking/results")
        self.rescues_dir = Path("Data/processed/rescues")
        self.output_dir = Path("Data/processed/synergy/scores")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        print("🧬 Synergy Calculator")
        print(f"   Docking results: {self.docking_dir}")
        print(f"   Rescue data: {self.rescues_dir}")
        print(f"   Output: {self.output_dir}")

    def load_functional_scores(self) -> pd.DataFrame:
        """Load functional scores from Initiative 1."""

        print("\n📊 Loading functional scores from Initiative 1...")

        all_rescues = []
        targets = ["R175H", "R248Q", "R273H", "Y220C"]

        for target in targets:
            parquet_path = self.rescues_dir / target / "pareto.parquet"
            if parquet_path.exists():
                df = pd.read_parquet(parquet_path)
                df["target"] = target
                all_rescues.append(df)

        if not all_rescues:
            raise FileNotFoundError("No functional scores found!")

        combined = pd.concat(all_rescues, ignore_index=True)

        # Normalize rescue mutations format (replace commas with underscores)
        combined["rescue_normalized"] = combined["rescue_mutations"].str.replace(",", "_")

        print(f"   ✅ Loaded {len(combined)} rescues across {len(targets)} targets")

        return combined

    def load_docking_results(self) -> List[Dict]:
        """Load all docking results."""

        print("\n📊 Loading docking results...")

        docking_results = []
        for result_file in self.docking_dir.rglob("docking_results.json"):
            with open(result_file) as f:
                data = json.load(f)
                docking_results.append(data)

        print(f"   ✅ Loaded {len(docking_results)} docking results")

        return docking_results

    def normalize_score(self, value: float, scale: float = 1.0, offset: float = 0.0) -> float:
        """
        Normalize score to [0, 1] range using sigmoid.

        Args:
            value: Raw value
            scale: Scaling factor (higher = steeper sigmoid)
            offset: Offset (shift sigmoid center)

        Returns:
            Normalized value in [0, 1]
        """
        return 1.0 / (1.0 + np.exp(scale * (value - offset)))

    def calculate_synergy_score(
        self,
        functional_score: float,
        consensus_affinity: float,
        affinity_std: float,
        ddg_folding: float,
        ddg_binding: float,
        ddg_interface: float
    ) -> Dict:
        """
        Calculate composite synergy score.

        Components:
        1. Functional Enhancement Score (FES)
        2. Drug binding affinity
        3. Binding consistency (low std = good)
        4. DNA safety check

        Args:
            functional_score: From Initiative 1 (0-1 scale)
            consensus_affinity: Median ΔG_binding from docking (kcal/mol)
            affinity_std: Std dev of ΔG_binding across ensemble
            ddg_folding: Folding stability (kcal/mol)
            ddg_binding: DNA binding affinity (kcal/mol)
            ddg_interface: Interface stability (kcal/mol)

        Returns:
            Dict with synergy components and overall score
        """

        # Component 1: Functional Enhancement Score
        # Does drug enhance beyond rescue alone?
        # Good affinity: < -8 kcal/mol
        # Excellent affinity: < -10 kcal/mol
        drug_binding_norm = self.normalize_score(consensus_affinity, scale=0.5, offset=-8.0)

        # FES: rescue function + drug binding
        fes = 0.6 * functional_score + 0.4 * drug_binding_norm

        # Component 2: Binding strength
        # More negative = better binding
        binding_strength_norm = self.normalize_score(consensus_affinity, scale=0.3, offset=-9.0)

        # Component 3: Binding consistency
        # Low std = consistent binding across conformers (good)
        # Typical std: 0.5-2.0 kcal/mol
        # Good: < 1.0, Excellent: < 0.5
        consistency_norm = self.normalize_score(affinity_std, scale=2.0, offset=1.0)

        # Component 4: DNA safety
        # Drug should not disrupt DNA binding
        # ΔΔG_binding < 2.0 kcal/mol = safe
        dna_safe = 1.0 if ddg_binding <= 2.0 else 0.5

        # Component 5: Interface safety
        # Drug should not disrupt tetramerization
        # ΔΔG_interface < 1.0 kcal/mol = safe
        interface_safe = 1.0 if ddg_interface <= 1.0 else 0.5

        # Overall synergy score (weighted combination)
        synergy_score = (
            0.35 * fes +
            0.25 * binding_strength_norm +
            0.20 * consistency_norm +
            0.10 * dna_safe +
            0.10 * interface_safe
        )

        # Categorize
        if synergy_score > 0.8:
            category = "superstar"
        elif synergy_score > 0.6:
            category = "promising"
        elif synergy_score > 0.4:
            category = "additive"
        else:
            category = "poor"

        return {
            "synergy_score": float(synergy_score),
            "category": category,
            "components": {
                "fes": float(fes),
                "binding_strength": float(binding_strength_norm),
                "consistency": float(consistency_norm),
                "dna_safe": float(dna_safe),
                "interface_safe": float(interface_safe)
            },
            "raw_values": {
                "functional_score": float(functional_score),
                "consensus_affinity": float(consensus_affinity),
                "affinity_std": float(affinity_std),
                "ddg_folding": float(ddg_folding),
                "ddg_binding": float(ddg_binding),
                "ddg_interface": float(ddg_interface)
            }
        }

    def calculate_all_synergies(self) -> pd.DataFrame:
        """Calculate synergy scores for all rescue-drug pairs."""

        print("\n🧮 Calculating synergy scores...")

        # Load data
        functional_df = self.load_functional_scores()
        docking_results = self.load_docking_results()

        if not docking_results:
            print("   ⚠️  No docking results found!")
            print("      Run docking first: python src/docking/run_ensemble_docking.py")
            return pd.DataFrame()

        # Calculate synergy for each pair
        synergy_data = []

        for docking_result in docking_results:
            rescue_name = docking_result["rescue"]
            drug_name = docking_result["drug"]

            # Find functional score
            rescue_row = functional_df[functional_df["rescue_normalized"] == rescue_name]

            if rescue_row.empty:
                print(f"   ⚠️  Functional score not found for {rescue_name}")
                continue

            rescue_row = rescue_row.iloc[0]

            # Calculate synergy
            synergy = self.calculate_synergy_score(
                functional_score=rescue_row["functional_score"],
                consensus_affinity=docking_result["consensus_affinity"],
                affinity_std=docking_result["affinity_std"],
                ddg_folding=rescue_row["ddg_folding"],
                ddg_binding=rescue_row["ddg_binding"],
                ddg_interface=rescue_row["ddg_interface"]
            )

            # Combine data
            synergy_data.append({
                "rescue": rescue_name,
                "drug": drug_name,
                "target": rescue_row["target"],
                "synergy_score": synergy["synergy_score"],
                "category": synergy["category"],
                "functional_score": rescue_row["functional_score"],
                "consensus_affinity": docking_result["consensus_affinity"],
                "affinity_std": docking_result["affinity_std"],
                "ddg_folding": rescue_row["ddg_folding"],
                "ddg_binding": rescue_row["ddg_binding"],
                "ddg_interface": rescue_row["ddg_interface"],
                "overall_category": rescue_row["overall_category"],
                **synergy["components"]
            })

            print(f"   ✅ {rescue_name} + {drug_name}: {synergy['synergy_score']:.3f} ({synergy['category']})")

        # Convert to DataFrame
        synergy_df = pd.DataFrame(synergy_data)

        if synergy_df.empty:
            print("   ⚠️  No synergy scores calculated")
            return synergy_df

        # Sort by synergy score
        synergy_df = synergy_df.sort_values("synergy_score", ascending=False)

        # Save results
        output_file = self.output_dir / "all_pairs_ranked.csv"
        synergy_df.to_csv(output_file, index=False)

        print(f"\n   💾 Saved results to: {output_file}")

        # Save superstars separately
        superstars = synergy_df[synergy_df["category"] == "superstar"]
        if not superstars.empty:
            superstar_file = self.output_dir / "superstar_pairs.csv"
            superstars.to_csv(superstar_file, index=False)
            print(f"   ⭐ Saved {len(superstars)} superstar pairs to: {superstar_file}")

        return synergy_df

    def generate_summary(self, synergy_df: pd.DataFrame):
        """Generate summary statistics."""

        if synergy_df.empty:
            return

        print("\n" + "=" * 80)
        print("Synergy Analysis Summary")
        print("=" * 80)

        print(f"\nTotal rescue-drug pairs: {len(synergy_df)}")

        print(f"\nCategory distribution:")
        category_counts = synergy_df["category"].value_counts()
        for category, count in category_counts.items():
            pct = 100 * count / len(synergy_df)
            print(f"  - {category.capitalize()}: {count} ({pct:.1f}%)")

        print(f"\nTop 10 rescue-drug pairs:")
        top_10 = synergy_df.head(10)
        for i, row in top_10.iterrows():
            print(f"  {row['rescue']} + {row['drug']}: "
                  f"{row['synergy_score']:.3f} ({row['category']})")

        print(f"\nBest drug per target:")
        for target in synergy_df["target"].unique():
            target_df = synergy_df[synergy_df["target"] == target]
            best = target_df.iloc[0]
            print(f"  {target}: {best['drug']} (synergy: {best['synergy_score']:.3f})")

        print(f"\nMean synergy scores by drug:")
        drug_means = synergy_df.groupby("drug")["synergy_score"].mean().sort_values(ascending=False)
        for drug, mean_score in drug_means.items():
            print(f"  {drug}: {mean_score:.3f}")


def main():
    calculator = SynergyCalculator()

    try:
        # Calculate synergies
        synergy_df = calculator.calculate_all_synergies()

        if not synergy_df.empty:
            # Generate summary
            calculator.generate_summary(synergy_df)

            print("\n✅ Synergy analysis complete!")
            exit(0)
        else:
            print("\n⚠️  No synergy scores calculated")
            print("   Make sure you have run:")
            print("   1. python src/md/setup_md_simulation.py (for all rescues)")
            print("   2. python src/md/analyze_trajectory.py (for all rescues)")
            print("   3. python src/docking/run_ensemble_docking.py (for all rescue-drug pairs)")
            exit(1)

    except Exception as e:
        print(f"\n❌ Synergy calculation failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()
