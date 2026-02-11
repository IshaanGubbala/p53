from pathlib import Path
import re
from typing import List

import json
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

    # ------------------------------------------------------------------
    # New campaign-oriented helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_position(mutation_str: str) -> str | None:
        """Extract the numeric residue position from a mutation string like 'R175H'."""
        match = re.search(r"(\d+)", mutation_str)
        return match.group(1) if match else None

    def generate_rescue_pml(
        self,
        target_label: str,
        rescue_mutations: List[str],
        cancer_mutations: List[str],
        output_path: Path,
        delivery_method: str = "",
        evidence_score: float = 0.0,
    ) -> Path:
        """Generate a publication-quality PML script for a rescue candidate.

        Parameters
        ----------
        target_label : str
            Human-readable label for this candidate (used in the title).
        rescue_mutations : list[str]
            Rescue mutations such as ``["N239Y", "S240R"]``.
        cancer_mutations : list[str]
            Cancer (hotspot) mutations such as ``["R175H"]``.
        output_path : Path
            Where to write the ``.pml`` file.
        delivery_method : str, optional
            Delivery annotation to include in the header comment.
        evidence_score : float, optional
            Numeric evidence / confidence score.

        Returns
        -------
        Path
            The *output_path* that was written.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        pml: list[str] = []

        # --- Header comment ---------------------------------------------------
        pml.append(f"# PyMOL script for rescue candidate: {target_label}")
        pml.append(f"# Rescue mutations : {', '.join(rescue_mutations)}")
        pml.append(f"# Cancer mutations  : {', '.join(cancer_mutations)}")
        if delivery_method:
            pml.append(f"# Delivery method   : {delivery_method}")
        pml.append(f"# Evidence score    : {evidence_score:.3f}")
        pml.append("")

        # --- Load structure & base representation ----------------------------
        pml.append(f"load {self.pdb_path}, wt_p53")
        pml.append("hide everything")
        pml.append("show cartoon, wt_p53")
        pml.append("color gray80, wt_p53")

        # --- DNA-binding domain highlight ------------------------------------
        pml.append("color lightblue, wt_p53 and resi 94-292")

        # --- Cancer mutations (red sticks + labels) --------------------------
        cancer_positions = [
            self._extract_position(m) for m in cancer_mutations
        ]
        cancer_positions = [p for p in cancer_positions if p is not None]

        if cancer_positions:
            cancer_resi = "+".join(cancer_positions)
            pml.append(f"select cancer_muts, resi {cancer_resi}")
            pml.append("show sticks, cancer_muts")
            pml.append("color red, cancer_muts")
            pml.append("set stick_radius, 0.25, cancer_muts")
            for mut, pos in zip(cancer_mutations, cancer_positions):
                pml.append(
                    f'label resi {pos} and name CA, "{mut}"'
                )

        # --- Rescue mutations (green sticks + labels) ------------------------
        rescue_positions = [
            self._extract_position(m) for m in rescue_mutations
        ]
        rescue_positions = [p for p in rescue_positions if p is not None]

        if rescue_positions:
            rescue_resi = "+".join(rescue_positions)
            pml.append(f"select rescue_muts, resi {rescue_resi}")
            pml.append("show sticks, rescue_muts")
            pml.append("color green, rescue_muts")
            pml.append("set stick_radius, 0.25, rescue_muts")
            for mut, pos in zip(rescue_mutations, rescue_positions):
                pml.append(
                    f'label resi {pos} and name CA, "{mut}"'
                )

        # --- Title annotation -------------------------------------------------
        rescue_str = ", ".join(rescue_mutations)
        title_line = f"{target_label} | rescue: {rescue_str} | score: {evidence_score:.3f}"
        pml.append(f'set title, "{title_line}"')

        # --- Ray-tracing / display defaults -----------------------------------
        pml.append("bg_color white")
        pml.append("set ray_shadows, on")
        pml.append("set ray_opaque_background, on")
        pml.append("set label_size, 14")
        pml.append("set label_color, black")
        pml.append("set label_font_id, 7")
        pml.append("set antialias, 2")
        pml.append("set ray_trace_mode, 1")

        # --- Zoom to mutation region with 15 A buffer ------------------------
        all_positions = cancer_positions + rescue_positions
        if all_positions:
            zoom_resi = "+".join(all_positions)
            pml.append(f"zoom resi {zoom_resi}, 15")
        else:
            pml.append("zoom wt_p53")

        pml.append("deselect")
        pml.append("")

        with open(output_path, "w") as f:
            f.write("\n".join(pml))

        logger.info("Wrote rescue PML script to %s", output_path)
        return output_path

    def generate_campaign_pml(
        self,
        top_df: pd.DataFrame,
        output_dir: Path,
        top_n: int = 5,
    ) -> List[Path]:
        """Generate individual PML scripts for the top *top_n* candidates.

        Parameters
        ----------
        top_df : pd.DataFrame
            Campaign results dataframe.  Expected columns include
            ``mutations_json`` (JSON list of rescue mutations),
            ``targets_json`` (JSON list of cancer/target mutations),
            and optionally ``delivery_method`` and ``evidence_score``.
        output_dir : Path
            Directory in which to write the ``.pml`` files.
        top_n : int
            Number of candidates to visualize (default 5).

        Returns
        -------
        list[Path]
            Paths of the generated ``.pml`` files.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        generated: list[Path] = []
        subset = top_df.head(top_n)

        for rank, (idx, row) in enumerate(subset.iterrows(), start=1):
            # --- Parse mutation lists ----------------------------------------
            try:
                rescue_mutations = json.loads(row["mutations_json"]) if isinstance(row.get("mutations_json"), str) else []
            except (json.JSONDecodeError, TypeError):
                rescue_mutations = []

            try:
                cancer_mutations = json.loads(row["targets_json"]) if isinstance(row.get("targets_json"), str) else []
            except (json.JSONDecodeError, TypeError):
                cancer_mutations = []

            if not rescue_mutations and not cancer_mutations:
                logger.warning("Rank %d (index %s): no mutations found, skipping PML", rank, idx)
                continue

            target_label = f"Rank{rank}"
            delivery_method = str(row.get("delivery_method", "")) if pd.notna(row.get("delivery_method")) else ""
            evidence_score = float(row.get("evidence_score", 0.0)) if pd.notna(row.get("evidence_score")) else 0.0

            out_path = output_dir / f"candidate_rank{rank}.pml"
            self.generate_rescue_pml(
                target_label=target_label,
                rescue_mutations=rescue_mutations,
                cancer_mutations=cancer_mutations,
                output_path=out_path,
                delivery_method=delivery_method,
                evidence_score=evidence_score,
            )
            generated.append(out_path)

        logger.info("Generated %d PML scripts in %s", len(generated), output_dir)
        return generated
