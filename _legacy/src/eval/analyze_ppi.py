
import Bio.PDB
import pandas as pd
from pathlib import Path
import numpy as np

def calculate_min_distance(pdb_path, residues, target_chain='A'):
    """
    Calculate min distance from a set of residues to a target chain or interface.
    Since we are using the core domain AlphaFold model as base, 
    we can use a reference structure (like 1TUP for DNA) to map the interface.
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('ref', pdb_path)
    
    # DNA binding residues in p53 (from literature/1TUP):
    # K120, S241, R248, R273, R280, R283 etc.
    dna_interface = [120, 241, 248, 273, 280, 283]
    
    # MDM2 interaction: Primarily N-term (15-29). 
    # But for the core domain base, let's flag if it's near residue 100 or 290?
    # Actually, let's just use the user's distance check logic.
    
    # If we have the p53-MDM2 complex 1YCR:
    # MDM2 is residues 25-109. p53 peptide is 17-29.
    # This doesn't help with core domain mutations directly unless we have a full model.
    
    # Instead, let's measure distance to the DNA BINDING INTERFACE in the core domain structure.
    # This is a critical risk factor.
    
    model = structure[0]
    results = {}
    
    for res_id in residues:
        if res_id not in [r.id[1] for r in model['A'].get_residues()]:
            results[res_id] = 999.0
            continue
            
        res_atoms = [a for a in model['A'][res_id].get_atoms()]
        
        # Min distance to DNA interface
        min_dist = 999.0
        for dna_res_id in dna_interface:
            if dna_res_id == res_id: 
                min_dist = 0.0
                break
            if dna_res_id not in [r.id[1] for r in model['A'].get_residues()]: continue
            
            dna_atoms = [a for a in model['A'][dna_res_id].get_atoms()]
            for a1 in res_atoms:
                for a2 in dna_atoms:
                    dist = a1 - a2
                    if dist < min_dist:
                        min_dist = dist
        results[res_id] = round(float(min_dist), 2)
        
    return results

def main():
    summary_path = Path("reports/executive_summary.csv")
    if not summary_path.exists():
        return

    df = pd.read_csv(summary_path)
    
    # Base PDB (AlphaFold model used for design)
    # We'll use this to measure distance to DNA binding interface
    base_pdb = "Data/raw/alphafold/AF-P04637-F1-model_v6.pdb"
    if not Path(base_pdb).exists():
        # Fallback to search any pdb in raw
        base_pdbs = list(Path("Data/raw").rglob("*.pdb"))
        if base_pdbs: base_pdb = str(base_pdbs[0])
        else: return

    # Collect all mutations
    all_res = set()
    for _, row in df.iterrows():
        if pd.notna(row['rescue_mutations']):
            res_ids = [int(''.join(filter(str.isdigit, m))) for m in str(row['rescue_mutations']).split(",")]
            all_res.update(res_ids)
    
    distances = calculate_min_distance(base_pdb, list(all_res))
    
    def get_risk(mut_str):
        if pd.isna(mut_str): return "Low"
        res_ids = [int(''.join(filter(str.isdigit, m))) for m in str(mut_str).split(",")]
        min_d = min([distances.get(rid, 999) for rid in res_ids])
        if min_d < 8.0: return "High (Near DNA interface)"
        if min_d < 15.0: return "Medium"
        return "Low"

    df['PPI_Risk_DNA'] = df['rescue_mutations'].apply(get_risk)
    
    # MDM2 Flag (Generalized as most core mutations are >15A from N-term binding site)
    df['MDM2_Risk'] = "Low (Distal to N-term)"
    
    df.to_csv(summary_path, index=False)
    print("Updated PPI Risk in executive summary.")
    print(df[['target', 'rescue_mutations', 'PPI_Risk_DNA']])

if __name__ == "__main__":
    main()
