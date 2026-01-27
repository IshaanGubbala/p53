#!/usr/bin/env python3
"""
Build mutant PDB structure using EvoEF2 BuildMutant.

Usage:
    python src/md/build_mutant_structure.py \
        --mutations "A189S,M133L,S95T" \
        --output "Data/processed/md_simulations/structures/A189S_M133L_S95T.pdb"
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.scoring.evoef2_runner import build_mutant_model


def load_config():
    """Load EvoEF2 configuration."""
    import yaml

    config_path = Path("configs/scoring.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    return config.get("evoef2", {})


def main():
    parser = argparse.ArgumentParser(description="Build mutant structure using EvoEF2")
    parser.add_argument("--mutations", required=True,
                        help="Comma-separated mutations (e.g., 'A189S,M133L,S95T')")
    parser.add_argument("--output", required=True,
                        help="Output PDB file path")
    parser.add_argument("--base-pdb",
                        default="Data/processed/cache/evoef2_repair/5a7051e8f97848d43b3747268825b29cafb5fa99d2e3fff6f58efcdb74e74cc3/AF-P04637-F1-model_v6.pdb",
                        help="Base PDB structure (default: AlphaFold p53)")

    args = parser.parse_args()

    # Parse mutations
    mutation_list = args.mutations.split(",")
    print(f"Building mutant structure with mutations: {mutation_list}")
    print(f"Base structure: {args.base_pdb}")
    print(f"Output: {args.output}")

    # Load EvoEF2 config
    evoef2_cfg = load_config()

    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build mutant model
    try:
        mutant_pdb = build_mutant_model(
            mutations=mutation_list,
            pdb_path=Path(args.base_pdb),
            evoef2_cfg=evoef2_cfg,
            work_root=output_path.parent / "build_tmp",
            recompute=False
        )

        # Copy to output location
        import shutil
        shutil.copy2(mutant_pdb, output_path)

        print(f"\n✅ Mutant structure built successfully: {output_path}")

        # Verify output
        if output_path.exists():
            file_size = output_path.stat().st_size / 1024  # KB
            print(f"   File size: {file_size:.1f} KB")

            # Count atoms
            with open(output_path) as f:
                atom_count = sum(1 for line in f if line.startswith("ATOM"))
            print(f"   Atoms: {atom_count}")

        return 0

    except Exception as e:
        print(f"\n❌ Failed to build mutant structure: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
