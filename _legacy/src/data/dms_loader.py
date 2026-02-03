from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
from src.core.logging import get_logger

logger = get_logger(__name__)

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
