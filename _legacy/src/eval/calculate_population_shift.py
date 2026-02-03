
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

def calculate_folded_fraction(delta_g, temperature_k=310.15):
    """
    Calculate the fraction of folded protein using the Boltzmann distribution.
    Delta G here is the stability (kcal/mol), where POSITIVE means stable.
    """
    R = 0.001987  # Gas constant in kcal/(mol*K)
    
    exponent = -delta_g / (R * temperature_k)
    # Prevent overflow/underflow
    if exponent > 50:
        return 0.0
    if exponent < -50:
        return 1.0
        
    k_eq = np.exp(exponent)
    return 1.0 / (1.0 + k_eq)

def process_reports(summary_path):
    """
    Reads the executive summary and adds folded population metrics.
    """
    df = pd.read_csv(summary_path)
    
    # Constants for p53 (Approximate stability of core domain at 37C)
    # Literature suggests p53 is highly marginally stable.
    DG_WT = 3.0  # kcal/mol at 37C (Estimated)
    
    # In the reporter, ddg_seed is the DESTABILIZATION of the mutant (Energy_mut - Energy_wt)
    # ddg_gain is the STABILIZATION from the rescue (negative value)
    
    # % Folded WT
    folded_wt = calculate_folded_fraction(DG_WT)
    
    # % Folded Mutant
    # dG_mut = dG_wt - ddg_seed
    df['dG_mutant'] = DG_WT - df['ddg_seed']
    df['Folded_Mutant_%'] = df['dG_mutant'].apply(calculate_folded_fraction) * 100
    
    # % Folded Rescued
    # dG_rescue = dG_wt - (ddg_seed + ddg_gain)
    df['dG_rescued'] = DG_WT - (df['ddg_seed'] + df['ddg_gain'])
    df['Folded_Rescued_%'] = df['dG_rescued'].apply(calculate_folded_fraction) * 100
    
    # Improvement factor
    # Handle division by zero if mutant is 0%
    df['Rescue_Efficiency'] = df['Folded_Rescued_%'] / df['Folded_Mutant_%'].replace(0, 0.0001)
    
    # Format for readability
    df['Folded_Mutant_%'] = df['Folded_Mutant_%'].round(2)
    df['Folded_Rescued_%'] = df['Folded_Rescued_%'].round(2)
    df['Rescue_Efficiency'] = df['Rescue_Efficiency'].round(1)
    
    # Update the CSV
    df.to_csv(summary_path, index=False)
    print(f"Updated {summary_path}")
    print(f"WT Folded Population (Estimated): {folded_wt*100:.2f}%")
    print(df[['target', 'rescue_mutations', 'Folded_Mutant_%', 'Folded_Rescued_%', 'Rescue_Efficiency']])

if __name__ == "__main__":
    summary_file = Path("/Users/ishaangubbala/Documents/p53/reports/executive_summary.csv")
    if summary_file.exists():
        process_reports(summary_file)
    else:
        print(f"File not found: {summary_file}")
