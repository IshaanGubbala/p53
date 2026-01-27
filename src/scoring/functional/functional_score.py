"""
Composite Functional Rescue Score

Combines multiple dimensions of rescue mutation quality:
1. Monomer stability (ΔΔG_folding from EvoEF2)
2. DNA binding affinity (ΔΔG_binding)
3. Tetramer interface stability (ΔΔG_interface)
4. Risk score (MSA conservation, burial, etc.)

Produces a multi-dimensional functional score for Pareto optimization.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import yaml
from dataclasses import dataclass, asdict

try:
    from .binding_affinity import calculate_binding_affinity_cached, BindingAffinityResult
    from .interface_stability import calculate_interface_stability_cached, InterfaceStabilityResult
except ImportError:
    # For standalone testing
    from binding_affinity import calculate_binding_affinity_cached, BindingAffinityResult
    from interface_stability import calculate_interface_stability_cached, InterfaceStabilityResult


@dataclass
class FunctionalScoreResult:
    """Complete functional scoring result."""
    # Component scores (raw)
    ddg_folding: float
    ddg_binding: float
    ddg_interface: float
    risk: float

    # Normalized scores (0-1, higher is better)
    folding_norm: float
    binding_norm: float
    interface_norm: float
    risk_norm: float

    # Composite score
    functional_score: float  # Weighted sum of normalized components

    # Categories
    folding_category: str
    binding_category: str
    interface_category: str
    overall_category: str  # "excellent", "good", "acceptable", "poor"

    # Sub-results
    binding_result: Optional[Dict] = None
    interface_result: Optional[Dict] = None


def load_config(config_path: str = "configs/functional_scoring.yaml") -> Dict:
    """Load functional scoring configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def normalize_score(value: float, min_val: float, max_val: float,
                    higher_is_better: bool = False) -> float:
    """
    Normalize a score to 0-1 range.

    Args:
        value: Raw score
        min_val: Minimum possible value (best)
        max_val: Maximum possible value (worst)
        higher_is_better: If True, higher values are better

    Returns:
        Normalized score (0 = worst, 1 = best)
    """
    if max_val == min_val:
        return 0.5

    # Clamp value to range
    value = np.clip(value, min_val, max_val)

    # Normalize to 0-1
    norm = (value - min_val) / (max_val - min_val)

    # Invert if lower is better (default for energies)
    if not higher_is_better:
        norm = 1.0 - norm

    return norm


def categorize_overall(folding_cat: str, binding_cat: str,
                       interface_cat: str) -> str:
    """
    Determine overall category from component categories.

    Rules:
    - Excellent: All components "good"
    - Good: ≥2 components "good", none "bad"
    - Acceptable: ≥1 component "good", or all "acceptable"
    - Poor: Any component "bad"
    """
    cats = [folding_cat, binding_cat, interface_cat]

    # Any bad component → poor overall
    if "bad" in cats:
        return "poor"

    good_count = sum(1 for c in cats if c == "good")
    acceptable_count = sum(1 for c in cats if c == "acceptable")

    if good_count == 3:
        return "excellent"
    elif good_count >= 2:
        return "good"
    elif good_count >= 1 or acceptable_count == 3:
        return "acceptable"
    else:
        return "poor"


def calculate_functional_score(
    mutations: List[Tuple[int, str, str]],  # (position, wt_aa, mut_aa)
    ddg_folding: float,  # From EvoEF2
    risk: float,  # From risk scoring
    config: Optional[Dict] = None,
    structure_paths: Optional[Dict[str, str]] = None,
    cache_dir: str = "Data/processed/cache/functional_scoring"
) -> FunctionalScoreResult:
    """
    Calculate comprehensive functional score for a rescue mutation.

    Args:
        mutations: List of (position, wt_aa, mut_aa) tuples
        ddg_folding: ΔΔG from monomer stability (EvoEF2)
        risk: Risk score (0-1, MSA conservation etc.)
        config: Configuration dict
        structure_paths: Override structure paths (for testing)
        cache_dir: Cache directory

    Returns:
        FunctionalScoreResult with all components
    """
    if config is None:
        config = load_config()

    func_config = config['functional_score']
    norm_config = func_config['normalization']
    weights = func_config['weights']

    # Get structure paths
    if structure_paths is None:
        structure_paths = {
            'dna_bound': config['structures']['dna_bound']['pdb'],
            'tetramer': config['structures']['tetramer']['pdb']
        }

    # Calculate DNA binding affinity
    binding_result = calculate_binding_affinity_cached(
        structure_paths['dna_bound'],
        mutations,
        config,
        cache_dir
    )
    ddg_binding = binding_result.ddg_binding

    # Calculate tetramer interface stability
    interface_result = calculate_interface_stability_cached(
        structure_paths['tetramer'],
        mutations,
        config,
        chain_a_id="A",
        chain_b_id="B",
        cache_dir=cache_dir
    )
    ddg_interface = interface_result.ddg_interface

    # Normalize all components to 0-1 (higher is better)
    folding_norm = normalize_score(
        ddg_folding,
        norm_config['folding']['min'],
        norm_config['folding']['max'],
        higher_is_better=False  # Lower ΔΔG is better
    )

    binding_norm = normalize_score(
        ddg_binding,
        norm_config['dna_binding']['min'],
        norm_config['dna_binding']['max'],
        higher_is_better=False
    )

    interface_norm = normalize_score(
        ddg_interface,
        norm_config['interface']['min'],
        norm_config['interface']['max'],
        higher_is_better=False
    )

    risk_norm = normalize_score(
        risk,
        norm_config['risk']['min'],
        norm_config['risk']['max'],
        higher_is_better=False  # Lower risk is better
    )

    # Calculate composite functional score
    functional_score = (
        weights['folding'] * folding_norm +
        weights['dna_binding'] * binding_norm +
        weights['interface'] * interface_norm +
        weights['risk'] * risk_norm
    )

    # Categorize folding (use typical thresholds)
    if ddg_folding <= -3.0:
        folding_cat = "good"
    elif ddg_folding <= 0.0:
        folding_cat = "acceptable"
    else:
        folding_cat = "bad"

    # Overall category
    overall_cat = categorize_overall(
        folding_cat,
        binding_result.category,
        interface_result.category
    )

    # Build result
    result = FunctionalScoreResult(
        ddg_folding=ddg_folding,
        ddg_binding=ddg_binding,
        ddg_interface=ddg_interface,
        risk=risk,
        folding_norm=folding_norm,
        binding_norm=binding_norm,
        interface_norm=interface_norm,
        risk_norm=risk_norm,
        functional_score=functional_score,
        folding_category=folding_cat,
        binding_category=binding_result.category,
        interface_category=interface_result.category,
        overall_category=overall_cat,
        binding_result=asdict(binding_result),
        interface_result=asdict(interface_result)
    )

    return result


def score_rescue_candidates(
    candidates_df: pd.DataFrame,
    config: Optional[Dict] = None,
    cache_dir: str = "Data/processed/cache/functional_scoring",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Score a dataframe of rescue candidates with functional scores.

    Args:
        candidates_df: DataFrame with columns: rescue_mutations, ddg_gain, risk
        config: Configuration
        cache_dir: Cache directory
        verbose: Print progress

    Returns:
        DataFrame with added functional score columns
    """
    if config is None:
        config = load_config()

    results = []

    for idx, row in candidates_df.iterrows():
        if verbose and idx % 10 == 0:
            print(f"Scoring {idx+1}/{len(candidates_df)}...")

        # Parse mutations from rescue_mutations string
        # Format: "M133L,T155A" → [(133, 'M', 'L'), (155, 'T', 'A')]
        mutations = []
        for mut in row['rescue_mutations'].split(','):
            if not mut:
                continue
            wt_aa = mut[0]
            mut_aa = mut[-1]
            pos = int(mut[1:-1])
            mutations.append((pos, wt_aa, mut_aa))

        # Calculate functional score
        try:
            func_result = calculate_functional_score(
                mutations,
                row['ddg_gain'],  # Use ddg_gain as folding score
                row['risk'],
                config,
                cache_dir=cache_dir
            )

            # Add to results
            results.append({
                'rescue_mutations': row['rescue_mutations'],
                'functional_score': func_result.functional_score,
                'ddg_folding': func_result.ddg_folding,
                'ddg_binding': func_result.ddg_binding,
                'ddg_interface': func_result.ddg_interface,
                'risk': func_result.risk,
                'folding_norm': func_result.folding_norm,
                'binding_norm': func_result.binding_norm,
                'interface_norm': func_result.interface_norm,
                'risk_norm': func_result.risk_norm,
                'overall_category': func_result.overall_category,
                'binding_category': func_result.binding_category,
                'interface_category': func_result.interface_category
            })
        except Exception as e:
            if verbose:
                print(f"  Error scoring {row['rescue_mutations']}: {e}")
            # Add placeholder with NaN
            results.append({
                'rescue_mutations': row['rescue_mutations'],
                'functional_score': np.nan,
                'ddg_folding': row['ddg_gain'],
                'ddg_binding': np.nan,
                'ddg_interface': np.nan,
                'risk': row['risk'],
                'overall_category': 'error'
            })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Merge with original
    output_df = candidates_df.merge(
        results_df,
        on='rescue_mutations',
        how='left',
        suffixes=('', '_func')
    )

    return output_df


# Example usage and testing
if __name__ == "__main__":
    import sys
    import os

    # Change to project root if needed
    if not os.path.exists('configs/functional_scoring.yaml'):
        os.chdir('../../../')  # Go up from src/scoring/functional to project root

    if not os.path.exists('configs/functional_scoring.yaml'):
        print("Error: Cannot find config file. Please run from project root.")
        sys.exit(1)

    print("=== Functional Score Calculator Test ===\n")

    # Test 1: Excellent rescue (stabilizing, safe for DNA and interface)
    print("Test 1: M133L (buried core stabilizer)")
    result = calculate_functional_score(
        mutations=[(133, 'M', 'L')],
        ddg_folding=-5.6,  # Strong stabilization
        risk=0.120
    )
    print(f"  Functional score: {result.functional_score:.3f}")
    print(f"  Overall category: {result.overall_category}")
    print(f"  Components:")
    print(f"    Folding: {result.ddg_folding:.2f} kcal/mol ({result.folding_category})")
    print(f"    DNA binding: {result.ddg_binding:.2f} kcal/mol ({result.binding_category})")
    print(f"    Interface: {result.ddg_interface:.2f} kcal/mol ({result.interface_category})")
    print(f"    Risk: {result.risk:.3f}")
    print()

    # Test 2: Bad rescue (disrupts DNA binding)
    print("Test 2: R248Q (DNA-binding hotspot)")
    result = calculate_functional_score(
        mutations=[(248, 'R', 'Q')],
        ddg_folding=-2.0,  # Modest stabilization
        risk=0.150
    )
    print(f"  Functional score: {result.functional_score:.3f}")
    print(f"  Overall category: {result.overall_category}")
    print(f"  Components:")
    print(f"    Folding: {result.ddg_folding:.2f} kcal/mol ({result.folding_category})")
    print(f"    DNA binding: {result.ddg_binding:.2f} kcal/mol ({result.binding_category})")
    print(f"    Interface: {result.ddg_interface:.2f} kcal/mol ({result.interface_category})")
    print(f"    Risk: {result.risk:.3f}")
    print()

    # Test 3: Interface-disrupting rescue
    print("Test 3: R181E (interface residue)")
    result = calculate_functional_score(
        mutations=[(181, 'R', 'E')],
        ddg_folding=-1.5,
        risk=0.100
    )
    print(f"  Functional score: {result.functional_score:.3f}")
    print(f"  Overall category: {result.overall_category}")
    print(f"  Components:")
    print(f"    Folding: {result.ddg_folding:.2f} kcal/mol ({result.folding_category})")
    print(f"    DNA binding: {result.ddg_binding:.2f} kcal/mol ({result.binding_category})")
    print(f"    Interface: {result.ddg_interface:.2f} kcal/mol ({result.interface_category})")
    print(f"    Risk: {result.risk:.3f}")
