from __future__ import annotations

import json
import hashlib
import os
import shutil
import glob
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import numpy as np
import torch

from src.core.hashing import file_sha256
from src.core.logging import get_logger

# Import RaSP modules
RASP_DIR = Path(__file__).parent / "rasp_repo" / "src"
if str(RASP_DIR) not in sys.path:
    sys.path.append(str(RASP_DIR))

try:
    from cavity_model import CavityModel, DownstreamModel, ResidueEnvironmentsDataset
    from helpers import ds_pred, init_lin_weights
    from Bio.PDB.Polypeptide import index_to_one, one_to_index
    RaSP_OK = True
except ImportError as e:
    print(f"RaSP import failed: {e}")
    RaSP_OK = False

CACHE_PREFIX = "rasp"

class RaSPScorer:
    pass # Not used directly, function based interface preferred

def score_mutation_set_rasp(
    mutations: list[str],
    pdb_path: Path,
    cache_dir: Path,
    rasp_cfg: dict[str, Any],
    base_results: dict[str, float] | None = None,
) -> float:
    """
    Score mutations using RaSP. Requires RaSP repo set up in src/scoring/rasp_repo.
    """
    logger = get_logger(__name__)
    if not RaSP_OK:
        raise RuntimeError("RaSP modules missing.")

    mutations = sorted(list(mutations))
    pdb_hash = file_sha256(pdb_path)
    
    # Cache key for the entire PDB predictions
    full_cache_key = f"rasp_full_{pdb_hash}"
    full_cache_file = cache_dir / f"{full_cache_key}.json"
    
    predictions = {}
    
    if full_cache_file.exists():
        # logger.debug("Loading cached RaSP predictions")
        predictions = json.loads(full_cache_file.read_text())
    else:
        logger.info(f"Running RaSP inference on {pdb_path}")
        
        # Setup temporary directories for RaSP 
        # RaSP scripts are messy with files, so use a dedicated temp dir in cache
        rasp_work_dir = cache_dir / f"rasp_work_{pdb_hash}"
        rasp_work_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            predictions = _run_rasp_pipeline(pdb_path, rasp_work_dir, rasp_cfg)
            full_cache_file.write_text(json.dumps(predictions))
        except Exception as e:
            logger.error(f"RaSP pipeline failed: {e}")
            raise
        finally:
            # Cleanup work dir to save space (parquet/npz files)
            # shutil.rmtree(rasp_work_dir, ignore_errors=True)
            pass

    total_ddg = 0.0
    for mut in mutations:
        # Standardize mut: e.g. R175H
        if mut not in predictions:
            # Try to handle chain? RaSP usually predicts for chain A if we filtered.
            # predictions keys are likely "A15D" (WT+Pos+Mut) or "15D" if chain is dropped
            # My current prediction keys are "WT+Pos+Mut" (e.g. M1A)
            msg = f"Mutation {mut} not found in RaSP predictions."
            logger.warning(msg)
            continue
            
        total_ddg += predictions[mut]

    return total_ddg

def _run_rasp_pipeline(pdb_path: Path, work_dir: Path, config: dict) -> dict[str, float]:
    """
    Execute RaSP inference flow: clean -> extract -> predict.
    Returns dict { "A123C": 1.25 }
    """
    logger = get_logger(__name__)
    
    # 1. Paths
    raw_dir = work_dir / "raw"
    cleaned_dir = work_dir / "cleaned"
    parsed_dir = work_dir / "parsed"
    for d in [raw_dir, cleaned_dir, parsed_dir]:
        d.mkdir(exist_ok=True)
    
    # Copy PDB to raw
    local_pdb = raw_dir / "query.pdb"
    shutil.copy(pdb_path, local_pdb)
    
    # 2. Clean PDB
    # Use simplified clean_pdb_simple.py (doesn't require pdbfixer/openmm)
    clean_script = RASP_DIR / "pdb_parser_scripts" / "clean_pdb_simple.py"
    reduce_exe = RASP_DIR / "pdb_parser_scripts" / "reduce" / "reduce_src" / "reduce"

    if not reduce_exe.exists():
        raise RuntimeError(f"Reduce binary not found at {reduce_exe}")

    cmd_clean = [
        sys.executable, str(clean_script),
        "--pdb_file_in", str(local_pdb),
        "--out_dir", str(cleaned_dir),
        "--reduce_exe", str(reduce_exe)
    ]

    logger.info("Cleaning PDB with simplified cleaner...")
    result = subprocess.run(cmd_clean, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Clean PDB failed: {result.stderr}")
        raise RuntimeError(f"PDB cleaning failed: {result.stderr}")
    
    # Find cleaned file
    cleaned_pdb = next(cleaned_dir.glob("*_clean.pdb"))
    
    # 3. Extract Environments
    # python extract_environments.py --pdb_in ... --out_dir ...
    extract_script = RASP_DIR / "pdb_parser_scripts" / "extract_environments.py"
    cmd_extract = [
        sys.executable, str(extract_script),
        "--pdb_in", str(cleaned_pdb),
        "--out_dir", str(parsed_dir)
    ]
    
    logger.info("Extracting environments...")
    subprocess.run(cmd_extract, check=True, capture_output=True, text=True)
    
    # Find npz
    npz_file = next(parsed_dir.glob("*_coordinate_features.npz"))
    
    # 4. Load Data & Predict
    DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Running inference on {DEVICE}")
    
    # Dataset
    # ResidueEnvironmentsDataset expects a LIST of files
    dataset = ResidueEnvironmentsDataset([str(npz_file)], transformer=None)
    
    resenv_dataset = {}
    # Key format in dataset: {pdbid}{chain}{resno}{restype}
    # We need to reconstruct this carefully.
    
    # Replicate run_pipeline.py logic lines 360+
    # But simplified for single PDB
    
    for resenv in dataset:
        # We don't care about the key format as long as we map correctly later
        # RaSP pipeline logic:
        # key = f"--{PDB_ID}--{resenv.chain_id}--{...}"
        # We can just use resenv properties directly
        # But ds_pred expects a dataframe populated with resenvs.
        pass

    # Create DataFrame structure for saturation mutagenesis
    # Rows: 1 per residue * 20 AAs
    
    res_list = []
    for resenv in dataset:
        res_list.append({
            "chainid": resenv.chain_id,
            "pos": resenv.pdb_residue_number,
            "wt_AA": index_to_one(resenv.restype_index),
            "resenv": resenv
        })
        
    df_structure_no_mt = pd.DataFrame(res_list)
    df_structure_no_mt["pdbid"] = "query" 
    
    # Expand to 20 AAs
    aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    
    rows = []
    for _, row in df_structure_no_mt.iterrows():
        for mt in aa_list:
            new_row = row.copy()
            new_row["variant"] = f"{row['wt_AA']}{row['pos']}{mt}"
            rows.append(new_row)
            
    df_structure = pd.DataFrame(rows)
    
    # Load frequencies
    # Located at RASP_DIR/../data/pdb_frequencies.npz
    freq_path = RASP_DIR.parent / "data" / "pdb_frequencies.npz"
    pdb_nlfs = -np.log(np.load(str(freq_path))["frequencies"])
    
    df_structure["wt_idx"] = df_structure.apply(lambda row: one_to_index(row["variant"][0]), axis=1)
    df_structure["mt_idx"] = df_structure.apply(lambda row: one_to_index(row["variant"][-1]), axis=1)
    df_structure["wt_nlf"] = df_structure.apply(lambda row: pdb_nlfs[row["wt_idx"]], axis=1)
    df_structure["mt_nlf"] = df_structure.apply(lambda row: pdb_nlfs[row["mt_idx"]], axis=1)
    
    # Load Models
    repo_root = RASP_DIR.parent
    cavity_models_dir = repo_root / "output" / "cavity_models"
    ds_models_dir = repo_root / "output" / "ds_models"
    
    # Cavity Model
    # Determine best model path
    best_model_file = cavity_models_dir / "best_model_path.txt"
    best_model_rel = best_model_file.read_text().strip()
    # Strip paths if present
    best_model_rel = Path(best_model_rel).name
    best_model_path = cavity_models_dir / best_model_rel
    
    cavity_model = CavityModel(get_latent=True).to(DEVICE)
    cavity_model.load_state_dict(torch.load(str(best_model_path), map_location=torch.device(DEVICE)))
    cavity_model.eval()
    
    ds_model = DownstreamModel().to(DEVICE)
    ds_model.apply(init_lin_weights)
    ds_model.eval()
    
    # Predict
    # ds_pred takes NUM_ENSEMBLE.
    # It loads models internally using f"/content/..." hardcoded paths in helpers.py.
    # I PATCHED helpers.py to look in `repo_root`.
    # `helpers.py` calculates `repo_root` as `dirname(dirname(abspath(__file__)))`.
    # `__file__` is `helpers.py`. `dir` is `src`. `dir(dir)` is `rasp_repo`.
    # So `repo_root` in helpers.py = `.../rasp_repo`.
    # Which matches my structure.
    
    # However, ds_pred loads models by `model_idx` 0..9.
    
    logger.info("Running RaSP prediction loop...")
    df_ml = ds_pred(
        cavity_model,
        ds_model,
        df_structure,
        "predictions",
        10, # NUM_ENSEMBLE
        DEVICE
    )
    
    # df_ml has score_ml, pdbid, chainid, variant
    results = {}
    for _, row in df_ml.iterrows():
        # variant e.g. M1A
        var = row['variant']
        score = row['score_ml']
        results[var] = float(score)
        
    return results
