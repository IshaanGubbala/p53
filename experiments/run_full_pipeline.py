#!/usr/bin/env python3
"""
p53 StabiliMut: Complete Automated Pipeline
Runs all stages for all proteins of interest with all enhancements enabled.

Usage:
    python -m experiments.run_full_pipeline --targets R175H R248Q R273H
    python -m experiments.run_full_pipeline --all-hotspots
    python -m experiments.run_full_pipeline --skip-msa --skip-generative
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

# Add project root to path
sys.path.append(str(Path.cwd()))

from src.core.logging import get_logger, setup_logging

# Import experiment modules
from experiments import run_build_msa
from experiments import run_design_rescues
from experiments import run_generative_design
from experiments import run_score_variants
from experiments import run_make_report

logger = get_logger(__name__)

# Default hotspot mutations
DEFAULT_TARGETS = ["R175H", "R248Q", "R273H"]

# Extended hotspots (optional)
EXTENDED_TARGETS = DEFAULT_TARGETS + [
    "Y220C",    # Y220 pocket mutation (druggable)
    "G245S",    # Loop mutation
    "R249S",    # Contact mutation
    "R282W",    # DNA contact
]


def stage_1_build_msa(skip: bool = False) -> bool:
    """
    Stage 1: Build MSA and compute conservation scores.

    Returns:
        True if successful, False otherwise
    """
    if skip:
        logger.info("⏭️  Stage 1 (MSA) skipped per user request")
        return True

    logger.info("=" * 80)
    logger.info("STAGE 1: MSA Conservation Analysis")
    logger.info("=" * 80)

    try:
        start = time.time()

        # Check if MSA already exists
        msa_file = Path("Data/processed/msa/P04637_conservation.json")
        if msa_file.exists():
            logger.info(f"✅ MSA conservation file already exists: {msa_file}")
            logger.info("   Delete it to regenerate.")
            return True

        # Run MSA building
        logger.info("Fetching TP53 homologs and computing conservation...")
        result = run_build_msa.main()

        elapsed = time.time() - start
        logger.info(f"✅ Stage 1 complete in {elapsed:.1f}s")
        return result == 0

    except Exception as e:
        logger.error(f"❌ Stage 1 failed: {e}")
        return False


def stage_2_design_rescues(targets: list[str]) -> bool:
    """
    Stage 2: Beam search rescue design with all scoring enabled.

    Args:
        targets: List of target mutations (e.g., ["R175H", "R248Q"])

    Returns:
        True if successful, False otherwise
    """
    logger.info("=" * 80)
    logger.info(f"STAGE 2: Rescue Design (Beam Search)")
    logger.info(f"Targets: {', '.join(targets)}")
    logger.info("=" * 80)

    try:
        start = time.time()

        # Override sys.argv for run_design_rescues
        original_argv = sys.argv.copy()
        sys.argv = ["run_design_rescues", "--targets"] + targets

        result = run_design_rescues.main()

        # Restore argv
        sys.argv = original_argv

        elapsed = time.time() - start
        logger.info(f"✅ Stage 2 complete in {elapsed:.1f}s")
        return result == 0

    except Exception as e:
        logger.error(f"❌ Stage 2 failed: {e}")
        return False


def stage_3_generative_design(targets: list[str], n_samples: int = 100, skip: bool = False) -> bool:
    """
    Stage 3: Generative design with EvoDiff.

    Args:
        targets: List of target mutations
        n_samples: Number of sequences to generate per target
        skip: Skip this stage

    Returns:
        True if successful, False otherwise
    """
    if skip:
        logger.info("⏭️  Stage 3 (Generative Design) skipped per user request")
        return True

    logger.info("=" * 80)
    logger.info(f"STAGE 3: Generative Design (EvoDiff)")
    logger.info(f"Targets: {', '.join(targets)}, Samples: {n_samples}")
    logger.info("=" * 80)

    try:
        start = time.time()

        # Override sys.argv
        original_argv = sys.argv.copy()
        sys.argv = [
            "run_generative_design",
            "--targets", *targets,
            "--n_samples", str(n_samples)
        ]

        result = run_generative_design.main()

        # Restore argv
        sys.argv = original_argv

        elapsed = time.time() - start
        logger.info(f"✅ Stage 3 complete in {elapsed:.1f}s")
        return result == 0

    except Exception as e:
        logger.error(f"❌ Stage 3 failed: {e}")
        logger.warning("Continuing pipeline (generative design is optional)")
        return True  # Don't fail pipeline if generative fails


def stage_4_score_variants() -> bool:
    """
    Stage 4: Score ClinVar variants for validation.

    Returns:
        True if successful, False otherwise
    """
    logger.info("=" * 80)
    logger.info("STAGE 4: ClinVar Validation Scoring")
    logger.info("=" * 80)

    try:
        start = time.time()

        # Check if scoring already exists
        # (run_score_variants checks cache internally)

        result = run_score_variants.main()

        elapsed = time.time() - start
        logger.info(f"✅ Stage 4 complete in {elapsed:.1f}s")
        return result == 0

    except Exception as e:
        logger.error(f"❌ Stage 4 failed: {e}")
        return False


def stage_5_generate_report() -> bool:
    """
    Stage 5: Generate comprehensive PDF report.

    Returns:
        True if successful, False otherwise
    """
    logger.info("=" * 80)
    logger.info("STAGE 5: Report Generation")
    logger.info("=" * 80)

    try:
        start = time.time()

        result = run_make_report.main()

        elapsed = time.time() - start
        logger.info(f"✅ Stage 5 complete in {elapsed:.1f}s")
        return result == 0

    except Exception as e:
        logger.error(f"❌ Stage 5 failed: {e}")
        return False


def print_summary(targets: list[str], total_time: float):
    """Print pipeline completion summary."""
    logger.info("")
    logger.info("=" * 80)
    logger.info("🎉 PIPELINE COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total runtime: {total_time / 60:.1f} minutes")
    logger.info(f"Targets processed: {', '.join(targets)}")
    logger.info("")
    logger.info("📊 Generated Outputs:")
    logger.info("   • Data/processed/msa/P04637_conservation.json")

    for target in targets:
        logger.info(f"   • Data/processed/rescues/{target}/candidates.parquet")
        logger.info(f"   • Data/processed/rescues/{target}/pareto.parquet")
        logger.info(f"   • Data/processed/rescues/{target}/summary.json")

    logger.info("   • Data/processed/generative_candidates.parquet")
    logger.info("   • reports/guide.pdf")
    logger.info("   • reports/executive_summary.csv")
    logger.info("   • reports/executive_summary.tex")
    logger.info("   • reports/figures/pareto_3d_*.png")
    logger.info("")
    logger.info("📖 Next Steps:")
    logger.info("   1. Review Pareto fronts: reports/figures/pareto_3d_*.png")
    logger.info("   2. Read guide: reports/guide.pdf")
    logger.info("   3. Check top rescues: reports/executive_summary.csv")
    logger.info("   4. Visualize structures: PyMOL scripts in reports/pymol_scripts/")
    logger.info("")


def main():
    parser = argparse.ArgumentParser(
        description="Run complete p53 StabiliMut pipeline for all targets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run for default hotspots (R175H, R248Q, R273H)
  python -m experiments.run_full_pipeline

  # Run for specific targets
  python -m experiments.run_full_pipeline --targets R175H Y220C

  # Run for all extended hotspots
  python -m experiments.run_full_pipeline --all-hotspots

  # Skip optional stages
  python -m experiments.run_full_pipeline --skip-msa --skip-generative

  # Adjust generative sampling
  python -m experiments.run_full_pipeline --n-samples 200
        """
    )

    parser.add_argument(
        "--targets",
        nargs="+",
        help="Target mutations to process (default: R175H R248Q R273H)"
    )
    parser.add_argument(
        "--all-hotspots",
        action="store_true",
        help="Process all 7 extended hotspots (includes Y220C, G245S, R249S, R282W)"
    )
    parser.add_argument(
        "--skip-msa",
        action="store_true",
        help="Skip MSA building (use existing conservation scores)"
    )
    parser.add_argument(
        "--skip-generative",
        action="store_true",
        help="Skip generative design stage"
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of generative samples per target (default: 100)"
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip ClinVar validation scoring"
    )
    parser.add_argument(
        "--skip-report",
        action="store_true",
        help="Skip report generation"
    )

    args = parser.parse_args()

    # Setup logging
    setup_logging()

    # Determine targets
    if args.all_hotspots:
        targets = EXTENDED_TARGETS
    elif args.targets:
        targets = args.targets
    else:
        targets = DEFAULT_TARGETS

    logger.info("=" * 80)
    logger.info("p53 StabiliMut: Complete Automated Pipeline")
    logger.info("=" * 80)
    logger.info(f"Targets: {', '.join(targets)}")
    logger.info(f"Features enabled:")
    logger.info(f"  • MSA conservation: {'❌ SKIPPED' if args.skip_msa else '✅'}")
    logger.info(f"  • Multi-structure scoring: ✅")
    logger.info(f"  • Train/test validation: ✅")
    logger.info(f"  • Checkpointing: ✅")
    logger.info(f"  • RaSP deep learning: ✅")
    logger.info(f"  • Generative design: {'❌ SKIPPED' if args.skip_generative else '✅'}")
    logger.info(f"  • Allosteric targeting: ✅")
    logger.info(f"  • Docking: ✅ (mock mode)")
    logger.info("")

    pipeline_start = time.time()

    # Stage 1: MSA
    if not stage_1_build_msa(skip=args.skip_msa):
        logger.error("Pipeline failed at Stage 1 (MSA)")
        return 1

    # Stage 2: Rescue Design
    if not stage_2_design_rescues(targets):
        logger.error("Pipeline failed at Stage 2 (Rescue Design)")
        return 1

    # Stage 3: Generative Design
    if not stage_3_generative_design(targets, args.n_samples, skip=args.skip_generative):
        logger.warning("Stage 3 (Generative Design) had issues, continuing...")

    # Stage 4: Validation
    if not args.skip_validation:
        if not stage_4_score_variants():
            logger.error("Pipeline failed at Stage 4 (Validation)")
            return 1
    else:
        logger.info("⏭️  Stage 4 (Validation) skipped per user request")

    # Stage 5: Report
    if not args.skip_report:
        if not stage_5_generate_report():
            logger.error("Pipeline failed at Stage 5 (Report)")
            return 1
    else:
        logger.info("⏭️  Stage 5 (Report) skipped per user request")

    total_time = time.time() - pipeline_start
    print_summary(targets, total_time)

    return 0


if __name__ == "__main__":
    sys.exit(main())
