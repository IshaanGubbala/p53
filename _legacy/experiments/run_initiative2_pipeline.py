#!/usr/bin/env python3
"""
Run full Initiative 2 pipeline: MD → Analysis → Docking → Synergy

Usage:
    # Full pipeline (all steps)
    python experiments/run_initiative2_pipeline.py --rescues top10 --drugs all

    # Single rescue + single drug (for testing)
    python experiments/run_initiative2_pipeline.py \
        --rescues "A189S_M133L_S95T" \
        --drugs "COTI-2"

    # Skip MD (use existing trajectories)
    python experiments/run_initiative2_pipeline.py \
        --rescues top10 \
        --drugs all \
        --skip-md
"""

import argparse
import json
import subprocess
from pathlib import Path
from typing import List
import time


# Top 10 rescues from Initiative 1
TOP_10_RESCUES = [
    "A189S_M133L_S95T",    # R175H, 0.895
    "M133L_R196Q_S95T",    # R248Q, 0.883
    "R196Q_S215A_S95T",    # R273H, 0.880
    "C229A_R196H_S95T",    # R273H, 0.879
    "R196H_S95T",          # R273H/Y220C, 0.874
    "M133L_S95T",          # R175H/Y220C, 0.858
    "M133L_S95A_T155A",    # R175H, 0.856
    "L145I_M133L_S95A",    # R175H, 0.855
    "A189S_M133L_S95A",    # R175H, 0.852
    "M133L_S95A",          # Multiple, 0.849
]

# Map rescues to targets (for MD simulation)
RESCUE_TO_TARGET = {
    "A189S_M133L_S95T": "R175H",
    "M133L_R196Q_S95T": "R248Q",
    "R196Q_S215A_S95T": "R273H",
    "C229A_R196H_S95T": "R273H",
    "R196H_S95T": "R273H",
    "M133L_S95T": "R175H",
    "M133L_S95A_T155A": "R175H",
    "L145I_M133L_S95A": "R175H",
    "A189S_M133L_S95A": "R175H",
    "M133L_S95A": "R175H",
}

# All drugs
ALL_DRUGS = [
    "COTI-2",
    "APR-246",
    "Nutlin-3a",
    "PHI-KAN-083",
    "PRIMA-1",
    "CP-31398"
]


class Initiative2Pipeline:
    """Run complete Initiative 2 pipeline."""

    def __init__(
        self,
        rescues: List[str],
        drugs: List[str],
        skip_md: bool = False,
        skip_analysis: bool = False,
        skip_docking: bool = False,
        md_duration: int = 50
    ):
        self.rescues = rescues
        self.drugs = drugs
        self.skip_md = skip_md
        self.skip_analysis = skip_analysis
        self.skip_docking = skip_docking
        self.md_duration = md_duration

        print("=" * 80)
        print("Initiative 2: Full Pipeline")
        print("=" * 80)
        print(f"Rescues: {len(rescues)}")
        print(f"Drugs: {len(drugs)}")
        print(f"Total pairs: {len(rescues) * len(drugs)}")
        print(f"MD duration: {md_duration} ns")
        print(f"\nPipeline steps:")
        print(f"  - MD simulations: {'SKIP' if skip_md else 'RUN'}")
        print(f"  - Trajectory analysis: {'SKIP' if skip_analysis else 'RUN'}")
        print(f"  - Ensemble docking: {'SKIP' if skip_docking else 'RUN'}")
        print(f"  - Synergy scoring: RUN (always)")
        print()

    def run_command(self, cmd: List[str], description: str) -> bool:
        """
        Run a command with progress tracking.

        Args:
            cmd: Command to run
            description: Description for logging

        Returns:
            True if successful, False otherwise
        """
        print(f"\n{'─' * 80}")
        print(f"Running: {description}")
        print(f"Command: {' '.join(cmd)}")
        print(f"{'─' * 80}\n")

        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                check=True,
                text=True
            )

            elapsed = time.time() - start_time
            print(f"\n✅ {description} completed in {elapsed / 60:.1f} minutes")
            return True

        except subprocess.CalledProcessError as e:
            elapsed = time.time() - start_time
            print(f"\n❌ {description} failed after {elapsed / 60:.1f} minutes")
            print(f"   Error: {e}")
            return False

    def phase1_md_simulations(self) -> int:
        """Phase 1: Run MD simulations for all rescues."""

        if self.skip_md:
            print("\n⏭️  Skipping MD simulations (using existing trajectories)")
            return len(self.rescues)

        print("\n" + "=" * 80)
        print("Phase 1: MD Simulations")
        print("=" * 80)

        successful = 0
        failed = []

        for i, rescue in enumerate(self.rescues, start=1):
            target = RESCUE_TO_TARGET.get(rescue, "R175H")

            print(f"\n[{i}/{len(self.rescues)}] Running MD for {rescue} (target: {target})")

            cmd = [
                "python", "src/md/setup_md_simulation.py",
                "--rescue", rescue.replace("_", ","),  # Convert back to comma format
                "--target", target,
                "--duration", str(self.md_duration)
            ]

            success = self.run_command(cmd, f"MD simulation: {rescue}")

            if success:
                successful += 1
            else:
                failed.append(rescue)
                print(f"   ⚠️  Continuing with remaining rescues...")

        print("\n" + "=" * 80)
        print(f"Phase 1 Complete: {successful}/{len(self.rescues)} successful")
        if failed:
            print(f"Failed: {', '.join(failed)}")
        print("=" * 80)

        return successful

    def phase2_trajectory_analysis(self) -> int:
        """Phase 2: Analyze trajectories and extract conformers."""

        if self.skip_analysis:
            print("\n⏭️  Skipping trajectory analysis (using existing conformers)")
            return len(self.rescues)

        print("\n" + "=" * 80)
        print("Phase 2: Trajectory Analysis")
        print("=" * 80)

        successful = 0
        failed = []

        for i, rescue in enumerate(self.rescues, start=1):
            print(f"\n[{i}/{len(self.rescues)}] Analyzing trajectory: {rescue}")

            cmd = [
                "python", "src/md/analyze_trajectory.py",
                "--rescue", rescue,
                "--n-clusters", "10"
            ]

            success = self.run_command(cmd, f"Trajectory analysis: {rescue}")

            if success:
                successful += 1
            else:
                failed.append(rescue)
                print(f"   ⚠️  Continuing with remaining rescues...")

        print("\n" + "=" * 80)
        print(f"Phase 2 Complete: {successful}/{len(self.rescues)} successful")
        if failed:
            print(f"Failed: {', '.join(failed)}")
        print("=" * 80)

        return successful

    def phase3_ensemble_docking(self) -> int:
        """Phase 3: Run ensemble docking for all rescue-drug pairs."""

        if self.skip_docking:
            print("\n⏭️  Skipping ensemble docking (using existing results)")
            return len(self.rescues) * len(self.drugs)

        print("\n" + "=" * 80)
        print("Phase 3: Ensemble Docking")
        print("=" * 80)

        total_pairs = len(self.rescues) * len(self.drugs)
        successful = 0
        failed = []

        pair_num = 0
        for rescue in self.rescues:
            for drug in self.drugs:
                pair_num += 1
                print(f"\n[{pair_num}/{total_pairs}] Docking {rescue} + {drug}")

                cmd = [
                    "python", "src/docking/run_ensemble_docking.py",
                    "--rescue", rescue,
                    "--drug", drug
                ]

                success = self.run_command(cmd, f"Docking: {rescue} + {drug}")

                if success:
                    successful += 1
                else:
                    failed.append(f"{rescue}+{drug}")
                    print(f"   ⚠️  Continuing with remaining pairs...")

        print("\n" + "=" * 80)
        print(f"Phase 3 Complete: {successful}/{total_pairs} successful")
        if failed:
            print(f"Failed: {', '.join(failed[:10])}" +
                  (f" and {len(failed)-10} more" if len(failed) > 10 else ""))
        print("=" * 80)

        return successful

    def phase4_synergy_scoring(self) -> bool:
        """Phase 4: Calculate synergy scores."""

        print("\n" + "=" * 80)
        print("Phase 4: Synergy Scoring")
        print("=" * 80)

        cmd = [
            "python", "src/docking/calculate_synergy.py"
        ]

        success = self.run_command(cmd, "Synergy scoring")

        print("\n" + "=" * 80)
        print(f"Phase 4 Complete: {'Success' if success else 'Failed'}")
        print("=" * 80)

        return success

    def run_full_pipeline(self):
        """Run complete pipeline."""

        start_time = time.time()

        print("\n🚀 Starting Initiative 2 pipeline...")
        print(f"   Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

        # Phase 1: MD Simulations
        md_success = self.phase1_md_simulations()

        # Phase 2: Trajectory Analysis
        analysis_success = self.phase2_trajectory_analysis()

        # Phase 3: Ensemble Docking
        docking_success = self.phase3_ensemble_docking()

        # Phase 4: Synergy Scoring
        synergy_success = self.phase4_synergy_scoring()

        # Final summary
        elapsed = time.time() - start_time

        print("\n" + "=" * 80)
        print("🎉 INITIATIVE 2 PIPELINE COMPLETE")
        print("=" * 80)
        print(f"Total time: {elapsed / 3600:.1f} hours")
        print(f"\nResults:")
        print(f"  - MD simulations: {md_success if not self.skip_md else 'SKIPPED'}")
        print(f"  - Trajectory analyses: {analysis_success if not self.skip_analysis else 'SKIPPED'}")
        print(f"  - Docking runs: {docking_success if not self.skip_docking else 'SKIPPED'}")
        print(f"  - Synergy scoring: {'Success' if synergy_success else 'Failed'}")
        print(f"\nOutput files:")
        print(f"  - Synergy scores: Data/processed/synergy/scores/all_pairs_ranked.csv")
        print(f"  - Superstar pairs: Data/processed/synergy/scores/superstar_pairs.csv")
        print(f"  - Docking results: Data/processed/docking/results/")
        print()

        if synergy_success:
            print("✅ Initiative 2 complete! Ready for 95%+ rigor validation.")
        else:
            print("⚠️  Pipeline completed with some errors. Check logs above.")


def main():
    parser = argparse.ArgumentParser(description="Run Initiative 2 full pipeline")

    parser.add_argument(
        "--rescues",
        default="top10",
        help="Rescues to run: 'top10', 'all', or comma-separated list (e.g., 'A189S_M133L_S95T,M133L_R196Q_S95T')"
    )

    parser.add_argument(
        "--drugs",
        default="all",
        help="Drugs to run: 'all' or comma-separated list (e.g., 'COTI-2,APR-246')"
    )

    parser.add_argument(
        "--skip-md",
        action="store_true",
        help="Skip MD simulations (use existing trajectories)"
    )

    parser.add_argument(
        "--skip-analysis",
        action="store_true",
        help="Skip trajectory analysis (use existing conformers)"
    )

    parser.add_argument(
        "--skip-docking",
        action="store_true",
        help="Skip docking (use existing results)"
    )

    parser.add_argument(
        "--md-duration",
        type=int,
        default=50,
        help="MD simulation duration in ns (default: 50)"
    )

    args = parser.parse_args()

    # Parse rescues
    if args.rescues == "top10":
        rescues = TOP_10_RESCUES
    elif args.rescues == "all":
        # TODO: Load all 28 unique rescues
        rescues = TOP_10_RESCUES
        print("⚠️  'all' not yet implemented, using top10")
    else:
        rescues = args.rescues.split(",")

    # Parse drugs
    if args.drugs == "all":
        drugs = ALL_DRUGS
    else:
        drugs = args.drugs.split(",")

    # Create pipeline
    pipeline = Initiative2Pipeline(
        rescues=rescues,
        drugs=drugs,
        skip_md=args.skip_md,
        skip_analysis=args.skip_analysis,
        skip_docking=args.skip_docking,
        md_duration=args.md_duration
    )

    # Run pipeline
    pipeline.run_full_pipeline()


if __name__ == "__main__":
    main()
