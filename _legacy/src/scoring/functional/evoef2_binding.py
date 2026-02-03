"""
EvoEF2-based binding affinity and interface stability calculations.

Uses EvoEF2 ComputeBinding to calculate:
1. DNA binding affinity (ΔΔG_binding)
2. Tetramer interface stability (ΔΔG_interface)

This replaces heuristic calculations with actual energy function computations.
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import yaml

from src.scoring.evoef2_runner import compute_binding, build_mutant_model
from src.core.logging import get_logger


@dataclass
class EvoEF2BindingResult:
    """Results from EvoEF2 binding energy calculation."""
    ddg_binding: float  # ΔΔG_binding (kcal/mol, positive = weaker binding)
    wt_binding: float  # WT binding energy
    mut_binding: float  # Mutant binding energy
    category: str  # "good", "acceptable", "bad"
    details: Dict


def load_config(config_path: str = "configs/functional_scoring.yaml") -> Dict:
    """Load functional scoring configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_evoef2_config(config_path: str = "configs/scoring.yaml") -> Dict:
    """Load EvoEF2 configuration."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    return cfg.get('evoef2', {})


def calculate_binding_evoef2(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],  # (position, wt_aa, mut_aa)
    split_chains: str,  # e.g., "ABC,EF" for protein vs DNA
    config: Optional[Dict] = None,
    evoef2_cfg: Optional[Dict] = None,
    cache_dir: str = "Data/processed/cache/functional_scoring"
) -> EvoEF2BindingResult:
    """
    Calculate ΔΔG_binding using EvoEF2 ComputeBinding.

    Args:
        pdb_path: Path to complex structure (protein+DNA or protein+protein)
        mutations: List of (position, wt_aa, mut_aa) tuples
        split_chains: Chain splitting (e.g., "ABC,EF" for protein vs DNA)
        config: Functional scoring config
        evoef2_cfg: EvoEF2 config
        cache_dir: Cache directory

    Returns:
        EvoEF2BindingResult with ΔΔG_binding
    """
    logger = get_logger(__name__)

    if config is None:
        config = load_config()
    if evoef2_cfg is None:
        evoef2_cfg = load_evoef2_config()

    # Create cache key
    mut_str = ",".join([f"{pos}{wt}{mut}" for pos, wt, mut in mutations])
    cache_key = hashlib.sha256(
        f"{pdb_path}:{split_chains}:{mut_str}:evoef2".encode()
    ).hexdigest()

    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    cache_file = Path(cache_dir) / f"evoef2_binding_{cache_key}.json"

    # Check cache
    if cache_file.exists():
        logger.debug("Loading cached binding result: %s", cache_key[:8])
        with open(cache_file) as f:
            cached = json.load(f)
        return EvoEF2BindingResult(**cached)

    logger.info("Calculating ΔΔG_binding for %s using EvoEF2", mut_str)

    # Create working directories
    work_base = Path(cache_dir) / "evoef2_work"
    wt_workdir = work_base / f"{cache_key}_wt"
    mut_workdir = work_base / f"{cache_key}_mut"
    wt_workdir.mkdir(parents=True, exist_ok=True)
    mut_workdir.mkdir(parents=True, exist_ok=True)

    pdb_path_obj = Path(pdb_path)

    # Step 1: Calculate WT binding energy
    logger.debug("Computing WT binding energy")
    wt_binding = compute_binding(
        pdb_path_obj,
        split_chains,
        evoef2_cfg,
        wt_workdir
    )

    # Step 2: Build mutant structure
    logger.debug("Building mutant structure")

    # Create mutation specification for EvoEF2
    # Format: wt_pos_mut (e.g., "M133L")
    mutation_specs = [
        f"{wt}{pos}{mut}"
        for pos, wt, mut in mutations
    ]

    # Determine chain for mutations (first protein chain from split)
    protein_chains = split_chains.split(',')[0]
    mutation_chain = protein_chains[0] if protein_chains else 'A'

    # Create EvoEF2 config with correct chain_id
    evoef2_cfg_mut = evoef2_cfg.copy()
    evoef2_cfg_mut['chain_id'] = mutation_chain
    # Don't use repaired_pdb for complex structures
    if 'repaired_pdb' in evoef2_cfg_mut:
        del evoef2_cfg_mut['repaired_pdb']

    # Build mutant using build_mutant_model
    mutant_pdb = build_mutant_model(
        mutations=mutation_specs,
        pdb_path=pdb_path_obj,
        evoef2_cfg=evoef2_cfg_mut,
        work_root=mut_workdir,
        recompute=False
    )

    # Step 3: Calculate mutant binding energy
    logger.debug("Computing mutant binding energy")
    mut_binding = compute_binding(
        mutant_pdb,
        split_chains,
        evoef2_cfg,
        mut_workdir
    )

    # Step 4: Calculate ΔΔG
    ddg_binding = mut_binding - wt_binding
    logger.info("ΔΔG_binding = %.4f kcal/mol", ddg_binding)

    # Categorize based on config thresholds
    thresholds = config['dna_binding']['thresholds']
    if ddg_binding <= thresholds['good']:
        category = "good"
    elif ddg_binding <= thresholds['acceptable']:
        category = "acceptable"
    else:
        category = "bad"

    # Build result
    result = EvoEF2BindingResult(
        ddg_binding=ddg_binding,
        wt_binding=wt_binding,
        mut_binding=mut_binding,
        category=category,
        details={
            'mutations': mutations,
            'split_chains': split_chains,
            'method': 'evoef2',
            'wt_pdb': str(pdb_path),
            'mut_pdb': str(mutant_pdb)
        }
    )

    # Cache result
    with open(cache_file, 'w') as f:
        result_dict = {
            'ddg_binding': result.ddg_binding,
            'wt_binding': result.wt_binding,
            'mut_binding': result.mut_binding,
            'category': result.category,
            'details': result.details
        }
        json.dump(result_dict, f, indent=2)

    return result


def calculate_interface_evoef2(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],
    split_chains: str = "A,B",  # Default: monomer A vs B
    config: Optional[Dict] = None,
    evoef2_cfg: Optional[Dict] = None,
    cache_dir: str = "Data/processed/cache/functional_scoring"
) -> EvoEF2BindingResult:
    """
    Calculate ΔΔG_interface using EvoEF2 ComputeBinding.

    This is identical to calculate_binding_evoef2 but with different default
    split_chains and uses interface thresholds for categorization.

    Args:
        pdb_path: Path to tetramer structure
        mutations: List of (position, wt_aa, mut_aa) tuples
        split_chains: Chain splitting (e.g., "A,B" for dimer interface)
        config: Functional scoring config
        evoef2_cfg: EvoEF2 config
        cache_dir: Cache directory

    Returns:
        EvoEF2BindingResult with ΔΔG_interface
    """
    logger = get_logger(__name__)

    if config is None:
        config = load_config()

    # Calculate binding energy (same method, different interpretation)
    result = calculate_binding_evoef2(
        pdb_path,
        mutations,
        split_chains,
        config,
        evoef2_cfg,
        cache_dir
    )

    # Re-categorize using interface thresholds instead of DNA binding thresholds
    thresholds = config['tetramer_interface']['thresholds']
    if result.ddg_binding <= thresholds['good']:
        category = "good"
    elif result.ddg_binding <= thresholds['acceptable']:
        category = "acceptable"
    else:
        category = "bad"

    # Update category
    result.category = category
    logger.info("ΔΔG_interface = %.4f kcal/mol (%s)", result.ddg_binding, category)

    return result


# Example usage and testing
if __name__ == "__main__":
    import sys
    from pathlib import Path

    # Test DNA binding affinity
    print("=== Testing EvoEF2 Binding Calculations ===\n")

    dna_pdb = "Data/raw/experimental_pdbs/1TSR.pdb"
    tetramer_pdb = "Data/raw/experimental_pdbs/3KMD.pdb"

    if Path(dna_pdb).exists():
        print("Test 1: DNA Binding - M133L (should be good)")
        try:
            result = calculate_binding_evoef2(
                dna_pdb,
                [(133, 'M', 'L')],
                split_chains="ABC,EF",  # Protein vs DNA
                cache_dir="/tmp/test_evoef2"
            )
            print(f"  WT binding: {result.wt_binding:.4f} kcal/mol")
            print(f"  Mutant binding: {result.mut_binding:.4f} kcal/mol")
            print(f"  ΔΔG_binding: {result.ddg_binding:.4f} kcal/mol")
            print(f"  Category: {result.category}")
            print()
        except Exception as e:
            print(f"  Error: {e}\n")

        print("Test 2: DNA Binding - R248Q (should be bad)")
        try:
            result = calculate_binding_evoef2(
                dna_pdb,
                [(248, 'R', 'Q')],
                split_chains="ABC,EF",
                cache_dir="/tmp/test_evoef2"
            )
            print(f"  WT binding: {result.wt_binding:.4f} kcal/mol")
            print(f"  Mutant binding: {result.mut_binding:.4f} kcal/mol")
            print(f"  ΔΔG_binding: {result.ddg_binding:.4f} kcal/mol")
            print(f"  Category: {result.category}")
            print()
        except Exception as e:
            print(f"  Error: {e}\n")
    else:
        print(f"Skipping DNA binding tests: {dna_pdb} not found\n")

    if Path(tetramer_pdb).exists():
        print("Test 3: Interface - M133L (should be good)")
        try:
            result = calculate_interface_evoef2(
                tetramer_pdb,
                [(133, 'M', 'L')],
                split_chains="A,B",
                cache_dir="/tmp/test_evoef2"
            )
            print(f"  ΔΔG_interface: {result.ddg_binding:.4f} kcal/mol")
            print(f"  Category: {result.category}")
            print()
        except Exception as e:
            print(f"  Error: {e}\n")
    else:
        print(f"Skipping interface tests: {tetramer_pdb} not found\n")

    print("✅ EvoEF2 binding tests complete")
