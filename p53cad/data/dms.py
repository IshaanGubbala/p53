from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
from p53cad.core.logging import get_logger

logger = get_logger(__name__)

P53_WT = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

def apply_mutation(wt_seq: str, mutation: str) -> str | None:
    """
    Applies a mutation string (e.g. R175H or R175H, T211I) to WT sequence.
    Returns None if mutation is complex (fs, *, del).
    """
    if "fs" in mutation or "*" in mutation or "del" in mutation or "ins" in mutation:
        return None
        
    seq = list(wt_seq)
    parts = [p.strip() for p in mutation.split(",")]
    
    try:
        for p in parts:
            # Parse R175H -> (wt_aa, pos, mut_aa)
            import re
            match = re.match(r"([A-Z])(\d+)([A-Z])", p)
            if match:
                wt_aa, pos, mut_aa = match.groups()
                pos = int(pos) - 1 # 1-indexed to 0-indexed
                if pos < len(seq) and seq[pos] == wt_aa:
                    seq[pos] = mut_aa
                else:
                    return None # Mismatch or out of bounds
            else:
                return None
        return "".join(seq)
    except Exception:
        return None

def hydrate_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a 'sequence' column to a dataframe with a 'mutation' column.
    """
    hydrated = []
    for _, row in df.iterrows():
        seq = apply_mutation(P53_WT, str(row["mutation"]))
        if seq:
            new_row = row.to_dict()
            new_row["sequence"] = seq
            hydrated.append(new_row)
            
    return pd.DataFrame(hydrated)

def load_dms_data(file_path: Path | str) -> pd.DataFrame:
    """
    Load Giacomelli 2018 DMS data from CVS/Excel.
    Expected columns: 'mutation' (e.g. A123V), 'score' (functional score).
    """
    path = Path(file_path)
    if not path.exists():
        logger.warning(f"DMS data not found at {path}. Generating MOCK data for testing.")
        return generate_mock_dms_data()
    
    try:
        if path.suffix == '.xlsx':
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path)
            
        # Basic cleaning - adjust column names as needed based on actual file
        # Assuming standard format or user processed it
        if "mutation" not in df.columns or "score" not in df.columns:
            # Try to infer or return as is
            pass
            
        return df
    except Exception as e:
        logger.error(f"Failed to load DMS data: {e}. Returning mock.")
        return generate_mock_dms_data()

def generate_mock_dms_data(n_samples: int = 1000) -> pd.DataFrame:
    """
    Generates synthetic p53 DMS data for pipeline verification.
    Simulates correlation between structural constraint and function.
    """
    logger.info(f"Generating {n_samples} mock DMS records...")
    
    # Random mutations
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    data = []
    
    for _ in range(n_samples):
        pos = np.random.randint(1, 393)
        wt = np.random.choice(aa)
        mut = np.random.choice(aa)
        while mut == wt:
            mut = np.random.choice(aa)
            
        mutation = f"{wt}{pos}{mut}"
        
        # Mock Logic:
        # DNA binding domain (100-300) is sensitive -> low score
        # Hotspots (175, 248, 273) -> very low score
        # Surface -> neutral score
        
        score = np.random.normal(0, 0.5) # mostly neutral
        
        if 100 <= pos <= 300:
            score -= 1.0 # Generally deletrious
            
        if pos in [175, 248, 273, 220, 249, 282]:
            score -= 2.0 # Pathogenic
            
        # Add some noise
        score += np.random.normal(0, 0.2)
        
        # Normalize to roughly -1 to 1 range for 'fitness'
        # Giacomelli scores are often enriched/depleted logs.
        
        data.append({
            "mutation": mutation,
            "score": score,
            "pos": pos,
            "wt": wt,
            "alt": mut
        })
        
    return pd.DataFrame(data)
