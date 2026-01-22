from __future__ import annotations

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Any, List

from src.core.logging import get_logger


def load_variant_scores(scores_file: Path) -> Dict[str, Any]:
    """Load variant scores from JSON file"""
    with scores_file.open('r') as f:
        return json.load(f)


def create_mutation_landscape_plot(scores: Dict[str, Any], output_dir: Path) -> Path:
    """Create a publication-ready mutation landscape plot"""
    logger = get_logger(__name__)
    
    # Prepare data for plotting
    mutations = []
    energies = []
    risks = []
    
    for mut_id, data in scores.items():
        if isinstance(data, dict):
            mutation = data.get('mutation', 'Unknown')
            energy = data.get('energy', 0.0)
            risk = data.get('risk_score', 0.0)
            
            mutations.append(mutation)
            energies.append(energy)
            risks.append(risk)
    
    if not mutations:
        logger.warning("No valid mutations found for plotting")
        return None
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('TP53 Mutation Landscape', fontsize=16, fontweight='bold')
    
    # 1. Energy distribution histogram
    ax1 = axes[0, 0]
    ax1.hist(energies, bins=30, alpha=0.7, color='steelblue', edgecolor='black')
    ax1.set_xlabel('Energy (ΔG)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Energy Distribution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Risk score scatter plot
    ax2 = axes[0, 1]
    scatter = ax2.scatter(range(len(mutations)), energies, c=risks, cmap='RdYlBu_r', s=60, alpha=0.7)
    ax2.set_xlabel('Mutation Index', fontsize=12)
    ax2.set_ylabel('Energy (ΔG)', fontsize=12)
    ax2.set_title('Mutation Energy Landscape', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add colorbar for risk
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label('Risk Score', fontsize=12)
    
    # Add mutation labels for top 10
    for i, (idx, mutation) in enumerate(list(scores.keys())[:10]):
        if idx < len(mutations) and i < len(energies):
            label_x = idx - 0.5
            label_y = energies[idx] + np.array(energies).std() * 0.1
            ax2.annotate(mutation, (label_x, label_y), fontsize=8, 
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    plt.tight_layout()
    
    # Save high-quality figure
    output_file = output_dir / "tp53_mutation_landscape.png"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved mutation landscape plot: {output_file}")
    
    return output_file


def create_mutation_summary_table(scores: Dict[str, Any], output_dir: Path) -> Path:
    """Create a publication-ready summary table"""
    logger = get_logger(__name__)
    
    # Sort mutations by risk score
    sorted_mutations = sorted(scores.items(), key=lambda x: x[1].get('risk_score', 0))
    
    # Create summary data
    summary_data = []
    for i, (mutation, data) in enumerate(sorted_mutations):
        energy = data.get('energy', 0.0)
        risk = data.get('risk_score', 0.0)
        engine = data.get('engine', 'fast')
        
        summary_data.append({
            'Rank': i + 1,
            'Mutation': mutation,
            'Energy (ΔG)': f"{energy:.3f}",
            'Risk Score': f"{risk:.3f}",
            'Confidence': data.get('energy_confidence', 1.0),
            'Engine': engine.capitalize()
        })
    
    # Create publication-ready DataFrame
    df_summary = pd.DataFrame(summary_data)
    
    # Save as both CSV and high-quality image
    csv_file = output_dir / "tp53_mutation_summary.csv"
    df_summary.to_csv(csv_file, index=False)
    
    # Create figure version of table
    fig_table, ax_table = plt.subplots(1, 1, figsize=(12, 8))
    ax_table.axis('tight')
    ax_table.axis('off')
    
    # Create table
    table_data = []
    headers = ['Rank', 'Mutation', 'Energy (ΔG)', 'Risk Score', 'Confidence', 'Engine']
    
    # Add header row
    table_data.append(headers)
    
    # Add data rows
    for i in range(min(10, len(sorted_mutations))):
        mutation = sorted_mutations[i][0]
        data = sorted_mutations[i][1]
        
        table_data.append([
            str(i + 1),
            mutation,
            f"{data.get('energy', 0.0):.3f}",
            f"{data.get('risk_score', 0.0):.3f}",
            str(data.get('energy_confidence', 1.0)),
            data.get('engine', 'fast').capitalize()
        ])
    
    table = ax_table.table(cellText=table_data, colLabels=headers, cellLoc='center', 
                        colWidths=[0.1, 0.25, 0.2, 0.15, 0.1, 0.15],
                        rowLocs=[len(table_data)-1] * [0.2, 0.1],
                        colColours=['#f0f0f0'] * len(headers),
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    # Style and save
    table.auto_set_fontsize(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)
    
    # Save table as high-quality image
    table_file = output_dir / "tp53_mutation_table.png"
    fig_table.savefig(table_file, dpi=300, bbox_inches='tight')
    
    logger.info(f"Saved mutation summary table: {table_file}")
    logger.info(f"Saved mutation summary CSV: {csv_file}")
    
    return csv_file, table_file


def generate_heatmap(scores: Dict[str, Any], output_dir: Path) -> Path:
    """Create a heatmap of mutation positions vs energy"""
    logger = get_logger(__name__)
    
    # Extract positions and energies
    position_energy = {}
    for mutation, data in scores.items():
        if isinstance(data, dict):
            # Extract position from mutation string
            pos_str = ''.join(filter(str.isdigit, mutation))
            if pos_str:
                position = int(pos_str)
                energy = data.get('energy', 0.0)
                position_energy[position] = min(position_energy.get(position, float('inf')), energy)
    
    if not position_energy:
        logger.warning("No valid mutation positions found for heatmap")
        return None
    
    # Create heatmap data
    positions = sorted(position_energy.keys())
    energy_matrix = []
    
    for pos in positions:
        row_energies = []
        for check_pos in range(1, 400):  # TP53 is ~393 aa
            if check_pos in position_energy:
                row_energies.append(position_energy[check_pos])
            else:
                row_energies.append(float('inf'))
        energy_matrix.append(row_energies)
    
    # Create heatmap
    fig_heatmap, ax_heatmap = plt.subplots(1, 1, figsize=(14, 8))
    
    im = ax_heatmap.imshow(energy_matrix, cmap='viridis_r', aspect='auto', 
                         interpolation='nearest', vmin=0, vmax=10)
    
    ax_heatmap.set_xlabel('Mutation Position', fontsize=12)
    ax_heatmap.set_ylabel('Check Position', fontsize=12)
    ax_heatmap.set_title('TP53 Mutation Energy Heatmap', fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt_heatmap.colorbar(im, ax=ax_heatmap)
    cbar.set_label('Energy (ΔG)', fontsize=12)
    
    heatmap_file = output_dir / "tp53_energy_heatmap.png"
    fig_heatmap.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    
    logger.info(f"Saved energy heatmap: {heatmap_file}")
    
    return heatmap_file


def create_visualizations(scores_file: Path, output_dir: Path) -> List[Path]:
    """Create all visualizations"""
    logger = get_logger(__name__)
    
    # Load scores
    scores = load_variant_scores(scores_file)
    
    if not scores:
        logger.warning("No scores found for visualization")
        return []
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create all visualizations
    files = []
    
    # 1. Mutation landscape
    landscape_file = create_mutation_landscape_plot(scores, output_dir)
    files.append(landscape_file)
    
    # 2. Summary table  
    csv_file, table_file = create_mutation_summary_table(scores, output_dir)
    files.extend([csv_file, table_file])
    
    # 3. Energy heatmap
    heatmap_file = generate_heatmap(scores, output_dir)
    if heatmap_file:
        files.append(heatmap_file)
    
    logger.info(f"Created {len(files)} visualizations in {output_dir}")
    return files