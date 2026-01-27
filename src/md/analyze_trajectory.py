#!/usr/bin/env python3
"""
Analyze MD trajectory and extract representative conformers.

Usage:
    python src/md/analyze_trajectory.py --rescue "A189S_M133L_S95T"

Requires:
    - MDAnalysis
    - matplotlib
    - scikit-learn (for clustering)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from sklearn.cluster import KMeans
except ImportError:
    print("❌ Missing dependencies. Install with:")
    print("   pip install MDAnalysis matplotlib scikit-learn")
    exit(1)


class TrajectoryAnalyzer:
    """Analyze MD trajectory and extract representative conformers."""

    def __init__(self, rescue_name: str, base_dir: Path):
        self.rescue_name = rescue_name
        self.work_dir = base_dir / rescue_name
        self.analysis_dir = self.work_dir / "analysis"
        self.analysis_dir.mkdir(exist_ok=True)

        print(f"🔬 Analyzing trajectory: {rescue_name}")
        print(f"   Work directory: {self.work_dir}")
        print(f"   Analysis directory: {self.analysis_dir}")

    def load_trajectory(self) -> mda.Universe:
        """Load trajectory using MDAnalysis."""

        # Try multiple possible file names
        topology_paths = [
            self.work_dir / "md.tpr",
            self.work_dir / "npt.gro",
            self.work_dir / "system.gro"
        ]

        trajectory_paths = [
            self.work_dir / "md.xtc",
            self.work_dir / "md.trr"
        ]

        topology = None
        for path in topology_paths:
            if path.exists():
                topology = path
                break

        trajectory = None
        for path in trajectory_paths:
            if path.exists():
                trajectory = path
                break

        if not topology or not trajectory:
            raise FileNotFoundError(
                f"Could not find topology or trajectory in {self.work_dir}\n"
                f"Looked for topology: {topology_paths}\n"
                f"Looked for trajectory: {trajectory_paths}"
            )

        print(f"\n📂 Loading trajectory...")
        print(f"   Topology: {topology}")
        print(f"   Trajectory: {trajectory}")

        u = mda.Universe(str(topology), str(trajectory))
        print(f"   ✅ Loaded {len(u.trajectory)} frames")
        print(f"   Duration: {u.trajectory.totaltime / 1000:.1f} ns")

        return u

    def calculate_rmsd(self, u: mda.Universe) -> np.ndarray:
        """Calculate backbone RMSD over trajectory."""

        print(f"\n📊 Calculating RMSD...")

        # Select backbone atoms
        backbone = u.select_atoms("backbone")

        # Align to first frame
        align.AlignTraj(u, u, select="backbone", in_memory=True).run()

        # Calculate RMSD
        rmsd_analysis = rms.RMSD(u, select="backbone")
        rmsd_analysis.run()

        rmsd_data = rmsd_analysis.results.rmsd
        times = rmsd_data[:, 1] / 1000  # Convert ps to ns
        rmsd_values = rmsd_data[:, 2]   # RMSD in Å

        # Save to file
        output_file = self.analysis_dir / "rmsd.txt"
        np.savetxt(
            output_file,
            np.column_stack([times, rmsd_values]),
            header="Time(ns) RMSD(Å)",
            fmt="%.3f"
        )

        print(f"   ✅ RMSD calculated")
        print(f"   Mean RMSD: {rmsd_values.mean():.2f} Å")
        print(f"   Std RMSD: {rmsd_values.std():.2f} Å")
        print(f"   Max RMSD: {rmsd_values.max():.2f} Å")
        print(f"   Saved to: {output_file}")

        # Plot RMSD
        self._plot_rmsd(times, rmsd_values)

        return rmsd_values

    def _plot_rmsd(self, times: np.ndarray, rmsd_values: np.ndarray):
        """Plot RMSD over time."""

        plt.figure(figsize=(10, 6))
        plt.plot(times, rmsd_values, linewidth=1.5)
        plt.xlabel("Time (ns)", fontsize=12)
        plt.ylabel("RMSD (Å)", fontsize=12)
        plt.title(f"Backbone RMSD: {self.rescue_name}", fontsize=14)
        plt.grid(True, alpha=0.3)

        # Add mean line
        mean_rmsd = rmsd_values.mean()
        plt.axhline(mean_rmsd, color='r', linestyle='--',
                   label=f'Mean: {mean_rmsd:.2f} Å')
        plt.legend()

        output_file = self.analysis_dir / "rmsd.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"   📊 RMSD plot saved to: {output_file}")

    def calculate_rmsf(self, u: mda.Universe) -> np.ndarray:
        """Calculate per-residue RMSF."""

        print(f"\n📊 Calculating RMSF...")

        # Select C-alpha atoms
        ca = u.select_atoms("name CA")

        # Calculate average positions
        avg_positions = np.zeros((len(ca), 3))
        for ts in u.trajectory:
            avg_positions += ca.positions
        avg_positions /= len(u.trajectory)

        # Calculate RMSF
        rmsf_values = np.zeros(len(ca))
        for ts in u.trajectory:
            diff = ca.positions - avg_positions
            rmsf_values += np.sum(diff**2, axis=1)

        rmsf_values = np.sqrt(rmsf_values / len(u.trajectory))

        # Get residue numbers
        residue_nums = [atom.resid for atom in ca]

        # Save to file
        output_file = self.analysis_dir / "rmsf.txt"
        np.savetxt(
            output_file,
            np.column_stack([residue_nums, rmsf_values]),
            header="Residue RMSF(Å)",
            fmt="%d %.3f"
        )

        print(f"   ✅ RMSF calculated")
        print(f"   Mean RMSF: {rmsf_values.mean():.2f} Å")
        print(f"   Max RMSF: {rmsf_values.max():.2f} Å (residue {residue_nums[np.argmax(rmsf_values)]})")
        print(f"   Saved to: {output_file}")

        # Plot RMSF
        self._plot_rmsf(residue_nums, rmsf_values)

        return rmsf_values

    def _plot_rmsf(self, residue_nums: List[int], rmsf_values: np.ndarray):
        """Plot RMSF per residue."""

        plt.figure(figsize=(12, 6))
        plt.plot(residue_nums, rmsf_values, linewidth=1.5)
        plt.xlabel("Residue Number", fontsize=12)
        plt.ylabel("RMSF (Å)", fontsize=12)
        plt.title(f"Per-Residue RMSF: {self.rescue_name}", fontsize=14)
        plt.grid(True, alpha=0.3)

        output_file = self.analysis_dir / "rmsf.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"   📊 RMSF plot saved to: {output_file}")

    def cluster_conformers(self, u: mda.Universe, n_clusters: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Cluster trajectory conformers using K-means.

        Args:
            u: MDAnalysis Universe
            n_clusters: Number of clusters (default: 10)

        Returns:
            cluster_labels: Cluster assignment for each frame
            cluster_centers: Indices of representative frames
        """

        print(f"\n🔍 Clustering conformers (K-means, n={n_clusters})...")

        # Select backbone atoms
        backbone = u.select_atoms("backbone")

        # Extract coordinates for all frames
        n_frames = len(u.trajectory)
        n_atoms = len(backbone)

        coords = np.zeros((n_frames, n_atoms * 3))
        for i, ts in enumerate(u.trajectory):
            coords[i] = backbone.positions.flatten()

        # K-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(coords)

        # Find representative frame for each cluster (closest to centroid)
        cluster_centers = []
        for i in range(n_clusters):
            cluster_frames = np.where(cluster_labels == i)[0]
            cluster_coords = coords[cluster_frames]
            centroid = kmeans.cluster_centers_[i]

            # Find closest frame to centroid
            distances = np.linalg.norm(cluster_coords - centroid, axis=1)
            closest_idx = cluster_frames[np.argmin(distances)]
            cluster_centers.append(closest_idx)

        cluster_centers = np.array(cluster_centers)

        print(f"   ✅ Clustering complete")
        print(f"   Cluster sizes: {np.bincount(cluster_labels)}")
        print(f"   Representative frames: {cluster_centers}")

        # Plot cluster distribution
        self._plot_clusters(cluster_labels, n_clusters)

        return cluster_labels, cluster_centers

    def _plot_clusters(self, cluster_labels: np.ndarray, n_clusters: int):
        """Plot cluster distribution over time."""

        plt.figure(figsize=(12, 6))

        # Plot cluster assignments over time
        times = np.arange(len(cluster_labels))
        plt.scatter(times, cluster_labels, c=cluster_labels,
                   cmap='tab10', s=10, alpha=0.6)

        plt.xlabel("Frame", fontsize=12)
        plt.ylabel("Cluster ID", fontsize=12)
        plt.title(f"Conformer Clustering: {self.rescue_name}", fontsize=14)
        plt.yticks(range(n_clusters))
        plt.grid(True, alpha=0.3)

        # Add cluster size histogram
        plt.figure(figsize=(8, 6))
        cluster_sizes = np.bincount(cluster_labels)
        plt.bar(range(n_clusters), cluster_sizes, color='steelblue', alpha=0.7)
        plt.xlabel("Cluster ID", fontsize=12)
        plt.ylabel("Number of Frames", fontsize=12)
        plt.title(f"Cluster Sizes: {self.rescue_name}", fontsize=14)
        plt.xticks(range(n_clusters))
        plt.grid(True, alpha=0.3)

        output_file = self.analysis_dir / "clusters.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close('all')

        print(f"   📊 Cluster plot saved to: {output_file}")

    def extract_representative_conformers(
        self,
        u: mda.Universe,
        cluster_centers: np.ndarray
    ) -> Path:
        """
        Extract representative conformers to multi-model PDB.

        Args:
            u: MDAnalysis Universe
            cluster_centers: Frame indices of representative conformers

        Returns:
            Path to multi-model PDB file
        """

        print(f"\n💾 Extracting {len(cluster_centers)} representative conformers...")

        # Select protein atoms
        protein = u.select_atoms("protein")

        # Output file
        output_file = self.analysis_dir / "representative_conformers.pdb"

        with mda.Writer(str(output_file), multiframe=True) as writer:
            for model_num, frame_idx in enumerate(cluster_centers, start=1):
                u.trajectory[frame_idx]
                writer.write(protein)
                print(f"   ✅ Extracted conformer {model_num} (frame {frame_idx})")

        print(f"   💾 Saved to: {output_file}")

        # Also save individual PDB files
        individual_dir = self.analysis_dir / "conformers"
        individual_dir.mkdir(exist_ok=True)

        for model_num, frame_idx in enumerate(cluster_centers, start=1):
            u.trajectory[frame_idx]
            individual_file = individual_dir / f"conformer_{model_num}.pdb"
            protein.write(str(individual_file))

        print(f"   💾 Individual PDBs saved to: {individual_dir}")

        return output_file

    def run_full_analysis(self, n_clusters: int = 10):
        """Run complete trajectory analysis pipeline."""

        try:
            # Load trajectory
            u = self.load_trajectory()

            # Calculate RMSD
            rmsd_values = self.calculate_rmsd(u)

            # Check convergence
            if rmsd_values.mean() > 5.0:
                print(f"\n   ⚠️  WARNING: High RMSD ({rmsd_values.mean():.2f} Å)")
                print(f"      Simulation may not be converged")
            elif rmsd_values.mean() > 3.0:
                print(f"\n   ⚠️  CAUTION: Moderate RMSD ({rmsd_values.mean():.2f} Å)")
                print(f"      Consider extending simulation")
            else:
                print(f"\n   ✅ RMSD converged ({rmsd_values.mean():.2f} Å)")

            # Calculate RMSF
            rmsf_values = self.calculate_rmsf(u)

            # Cluster conformers
            cluster_labels, cluster_centers = self.cluster_conformers(u, n_clusters)

            # Extract representative conformers
            output_pdb = self.extract_representative_conformers(u, cluster_centers)

            # Summary
            print("\n" + "=" * 80)
            print("✅ Trajectory Analysis Complete")
            print("=" * 80)
            print(f"Rescue: {self.rescue_name}")
            print(f"Frames analyzed: {len(u.trajectory)}")
            print(f"Duration: {u.trajectory.totaltime / 1000:.1f} ns")
            print(f"\nMetrics:")
            print(f"  - Mean RMSD: {rmsd_values.mean():.2f} Å")
            print(f"  - Mean RMSF: {rmsf_values.mean():.2f} Å")
            print(f"  - Clusters: {n_clusters}")
            print(f"\nOutput:")
            print(f"  - Analysis directory: {self.analysis_dir}")
            print(f"  - Representative conformers: {output_pdb}")
            print(f"  - Ready for docking!")

            return True

        except Exception as e:
            print(f"\n❌ Analysis failed: {e}")
            import traceback
            traceback.print_exc()
            return False


def main():
    parser = argparse.ArgumentParser(description="Analyze MD trajectory")
    parser.add_argument("--rescue", required=True, help="Rescue name (e.g., 'A189S_M133L_S95T')")
    parser.add_argument("--n-clusters", type=int, default=10, help="Number of clusters (default: 10)")
    parser.add_argument("--base-dir", type=Path, default=Path("Data/processed/md_simulations"),
                        help="Base directory for MD simulations")

    args = parser.parse_args()

    # Create analyzer
    analyzer = TrajectoryAnalyzer(args.rescue, args.base_dir)

    # Run analysis
    success = analyzer.run_full_analysis(n_clusters=args.n_clusters)

    if success:
        exit(0)
    else:
        exit(1)


if __name__ == "__main__":
    main()
