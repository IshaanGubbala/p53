#!/usr/bin/env python3
"""Generate advanced visualizations - standalone script."""
import sys
from pathlib import Path

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

from visualization.create_advanced_visualizations import (
    create_mutation_network,
    create_radar_chart_comparison,
    create_pipeline_flow_diagram,
    create_domain_architecture_diagram,
    create_progressive_optimization_plot,
)

def main():
    processed_dir = repo_root / "Data" / "processed"
    figures_dir = repo_root / "reports" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Creating advanced visualizations...")

    # 1. Mutation network
    print("  1/5 Creating mutation network...")
    create_mutation_network(processed_dir, figures_dir / "mutation_network_advanced.png")

    # 2. Radar charts
    print("  2/5 Creating radar chart comparison...")
    create_radar_chart_comparison(processed_dir, figures_dir / "radar_comparison_advanced.png")

    # 3. Pipeline flow
    print("  3/5 Creating pipeline flow diagram...")
    create_pipeline_flow_diagram(figures_dir / "pipeline_flow_advanced.png")

    # 4. Domain architecture
    print("  4/5 Creating domain architecture...")
    create_domain_architecture_diagram(processed_dir, figures_dir / "domain_architecture_advanced.png")

    # 5. Progressive optimization
    print("  5/5 Creating progressive optimization plot...")
    create_progressive_optimization_plot(processed_dir, figures_dir / "progressive_optimization_advanced.png")

    print(f"\n✓ All advanced visualizations saved to {figures_dir}")
    print("\nGenerated files:")
    for f in sorted(figures_dir.glob("*advanced.png")):
        print(f"  - {f.name}")

if __name__ == "__main__":
    main()
