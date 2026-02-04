from pathlib import Path
import pandas as pd
from p53cad.core.logging import get_logger

logger = get_logger(__name__)

class PyMolGenerator:
    """
    Generates PyMol scripts (.pml) to visualize p53 design candidates.
    Since we don't have PDBs for every mutant, we visualize the mutations 
    mapped onto the Wild-Type structure (P04637).
    """
    def __init__(self, pdb_path: str = "data/raw/p53_wt.pdb"):
        self.pdb_path = pdb_path
        from p53cad.data.dms import P53_WT
        self.wt_seq = P53_WT
        
    def generate_session(self, candidates_csv: Path, output_path: Path):
        """
        Creates a .pml script that:
        1. Loads the WT structure.
        2. Creates selections for each candidate's mutations.
        3. Colors them for checking.
        """
        if not candidates_csv.exists():
            logger.error(f"Candidates file not found: {candidates_csv}")
            return
            
        df = pd.read_csv(candidates_csv)
        
        # Start building the PML content
        pml_lines = []
        pml_lines.append(f"load {self.pdb_path}, wt_p53")
        pml_lines.append("hide everything")
        pml_lines.append("show cartoon, wt_p53")
        pml_lines.append("color white, wt_p53")
        pml_lines.append("set transparency, 0.3")
        
        # Highlight DNA binding domain (approx 100-300)
        pml_lines.append("color gray80, resi 94-292")
        
        # For each candidate, identify mutations and create a scene or selection
        # We need to know the mutation positions. 
        # Assuming we have a 'mutations' column (e.g. "R175H,A123V")
        # If not, we have to diff against WT again.
        
        # Simplification: Let's visualize the "Top 5" by score if available, or just the first few.
        if "predicted_function" in df.columns:
            df = df.sort_values("predicted_function", ascending=False)
            
        top_candidates = df.head(10)
        
        for i, row in top_candidates.iterrows():
            cand_id = f"cand_{i}"
            muts_str = row.get("mutations")
            if pd.isna(muts_str) or not muts_str:
                # Try to derive from sequence vs WT
                cand_seq = row.get("sequence")
                if isinstance(cand_seq, str) and len(cand_seq) == len(self.wt_seq):
                    muts = [f"{self.wt_seq[j]}{j+1}{cand_seq[j]}" for j in range(len(self.wt_seq)) if self.wt_seq[j] != cand_seq[j]]
                else:
                    continue
            else:
                muts = str(muts_str).split(",")
                muts = [m.strip() for m in muts if m and m != "nan"]
            
            if not muts:
                continue
                
            # Extract residue numbers
            # Format: A123V -> 123
            res_nums = []
            for m in muts:
                # Naive parse
                import re
                match = re.search(r"(\d+)", m)
                if match:
                    res_nums.append(match.group(1))
            
            if not res_nums:
                continue
                
            resi_str = "+".join(res_nums)
            
            # Create a selection
            pml_lines.append(f"select {cand_id}_muts, resi {resi_str}")
            pml_lines.append(f"show sticks, {cand_id}_muts")
            pml_lines.append(f"color red, {cand_id}_muts")
            
            # Label
            label_text = f"Rank {i+1}: {','.join(muts)}"
            # We can't easy label in 3D without an atom anchor, but we can print
            pml_lines.append(f"print 'Candidate {i+1}: {label_text}'")
            
        # Add some nice view settings
        pml_lines.append("bg_color white")
        pml_lines.append("set ray_shadows, on")
        pml_lines.append("zoom resi 100-300")
        
        with open(output_path, "w") as f:
            f.write("\n".join(pml_lines))
            
    def render_candidate_png(self, sequence: str, output_image: Path) -> bool:
        """
        Renders a specific protein sequence as a PNG by mapping mutations onto the WT structure.
        """
        import subprocess
        import re
        
        muts = [f"{self.wt_seq[j]}{j+1}{sequence[j]}" for j in range(len(self.wt_seq)) if self.wt_seq[j] != sequence[j]]
        res_nums = []
        for m in muts:
            match = re.search(r"(\d+)", m)
            if match:
                res_nums.append(match.group(1))
        
        pml_path = output_image.with_suffix(".pml")
        pml_lines = [
            f"load {self.pdb_path}, wt_p53",
            "hide everything",
            "show cartoon, wt_p53",
            "color gray70, wt_p53",
            "set transparency, 0.4",
            "bg_color white"
        ]
        
        if res_nums:
            resi_str = "+".join(res_nums)
            pml_lines.extend([
                f"select muts, resi {resi_str}",
                "show sticks, muts",
                "color red, muts",
                "set stick_radius, 0.3",
                "zoom muts, 20"
            ])
        else:
            pml_lines.append("zoom wt_p53")
            
        pml_lines.extend([
            "set ray_opaque_background, off",
            "set ray_shadows, on",
            f"png {output_image}, width=800, height=600, ray=1",
            "quit"
        ])
        
        with open(pml_path, "w") as f:
            f.write("\n".join(pml_lines))
            
        try:
            # Run PyMol in batch mode (no GUI)
            result = subprocess.run(["pymol", "-cq", str(pml_path)], capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"PyMol error: {result.stderr}")
                return False
            return output_image.exists()
        except Exception as e:
            logger.error(f"Failed to render PNG: {e}")
            return False
