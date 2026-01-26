
import subprocess
import os
import pandas as pd
from pathlib import Path

def run_fpocket(pdb_path):
    """
    Run fpocket on a PDB and parse the druggability score of the top pocket.
    Generates {pdb}_out directory.
    """
    if not Path(pdb_path).exists():
        return None
    
    cmd = ["fpocket", "-f", str(pdb_path)]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # fpocket generates a folder with _out suffix containing a _info.txt
    pdb_name = Path(pdb_path).stem
    out_dir = Path(pdb_path).parent / f"{pdb_name}_out"
    info_file = out_dir / f"{pdb_name}_info.txt"
    
    if not info_file.exists():
        return None
        
    # Read the first pocket's druggability score
    # Usually: "Druggability Score : 0.456"
    try:
        with open(info_file, 'r') as f:
            for line in f:
                if "Druggability Score" in line:
                    score = float(line.split(":")[1].strip())
                    return score
    except:
        return None
    return 0.0

def main():
    summary_path = Path("reports/executive_summary.csv")
    if not summary_path.exists():
        return

    df = pd.read_csv(summary_path)
    
    # PDBs to screen
    # We screen the Cancer Mutant and the Top Rescue
    targets = ["R175H", "R248Q"]
    pocket_scores = {}
    
    for target in targets:
        # Cancer mutant base
        cancer_pdb = Path(f"Data/processed/modeled_structures/{target}/cancer_build/{target}_cancer.pdb")
        if not cancer_pdb.exists():
             cancer_pdb = Path(f"Data/processed/modeled_structures/{target}/wild_type.pdb")
             
        score_cancer = run_fpocket(cancer_pdb)
        
        # Rescued mutant
        rescue_pdb = Path(f"Data/processed/modeled_structures/{target}/rescue_build/{target}_rescued.pdb")
        score_rescue = run_fpocket(rescue_pdb)
        
        pocket_scores[target] = (score_cancer, score_rescue)
        print(f"Target {target}: Cancer Pocket Score: {score_cancer}, Rescue Pocket Score: {score_rescue}")

    # For now, let's just add a flag if the rescue creates a better pocket
    def get_pocket_flag(row):
        target = row['target']
        if target in pocket_scores:
            c, r = pocket_scores[target]
            if r and c and r > c + 0.1:
                return f"🚩 New Pocket (D={r:.2f})"
            if r and r > 0.5:
                return f"Druggable Pocket (D={r:.2f})"
        return "None"

    df['Druggable_Pockets'] = df.apply(get_pocket_flag, axis=1)
    df.to_csv(summary_path, index=False)
    print("Updated Druggable Pockets in executive summary.")

if __name__ == "__main__":
    main()
