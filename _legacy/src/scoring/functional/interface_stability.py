"""
Tetramer Interface Stability Calculator

Estimates ΔΔG_interface using structural analysis of protein-protein contacts.
This identifies mutations that may disrupt tetramerization (dominant-negative effect).

Method:
1. Identify interface residues in WT tetramer structure
2. Introduce mutation and check for contact changes
3. Score based on lost/gained interface contacts and buried surface area
4. Return estimated ΔΔG_interface (positive = destabilizing)

Note: This is an approximation. For final validation, use FoldX or experimental data.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from Bio.PDB import PDBParser, NeighborSearch, Selection
import yaml
import hashlib
import json
from dataclasses import dataclass


@dataclass
class InterfaceStabilityResult:
    """Results from tetramer interface stability calculation."""
    ddg_interface: float  # Estimated ΔΔG_interface (kcal/mol)
    wt_interface_contacts: int
    mut_interface_contacts: int
    lost_contacts: int
    gained_contacts: int
    wt_buried_area: float  # Approximate buried surface area (Å²)
    mut_buried_area: float
    clashes: int
    category: str  # "good", "acceptable", "bad"
    details: Dict


def load_config(config_path: str = "configs/functional_scoring.yaml") -> Dict:
    """Load functional scoring configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_interface_residues(structure, chain_a_id: str, chain_b_id: str,
                           cutoff: float = 8.0) -> Tuple[Set[int], Set[int]]:
    """
    Identify residues at the interface between two chains.

    Args:
        structure: BioPython structure object
        chain_a_id: First chain ID
        chain_b_id: Second chain ID
        cutoff: Interface distance cutoff in Angstroms

    Returns:
        (chain_a_interface_resnums, chain_b_interface_resnums)
    """
    chain_a = structure[0][chain_a_id]
    chain_b = structure[0][chain_b_id]

    # Get all atoms from each chain
    atoms_a = []
    for residue in chain_a:
        if residue.id[0] == ' ':  # Standard residue
            for atom in residue:
                atoms_a.append((atom, residue.id[1]))

    atoms_b = []
    for residue in chain_b:
        if residue.id[0] == ' ':
            for atom in residue:
                atoms_b.append((atom, residue.id[1]))

    # Use NeighborSearch to find interface
    ns = NeighborSearch([atom for atom, _ in atoms_a + atoms_b])

    interface_a = set()
    interface_b = set()

    # Check each atom in chain A for neighbors in chain B
    for atom_a, resnum_a in atoms_a:
        nearby = ns.search(atom_a.coord, cutoff)
        for neighbor in nearby:
            # Check if neighbor is from chain B
            if neighbor.get_parent().get_parent().id == chain_b_id:
                interface_a.add(resnum_a)
                break

    # Check each atom in chain B for neighbors in chain A
    for atom_b, resnum_b in atoms_b:
        nearby = ns.search(atom_b.coord, cutoff)
        for neighbor in nearby:
            if neighbor.get_parent().get_parent().id == chain_a_id:
                interface_b.add(resnum_b)
                break

    return interface_a, interface_b


def get_interface_contacts(structure, chain_a_id: str, chain_b_id: str,
                           contact_cutoff: float = 4.5) -> Dict[int, int]:
    """
    Count interface contacts for each residue.

    Args:
        structure: BioPython structure
        chain_a_id: First chain ID
        chain_b_id: Second chain ID
        contact_cutoff: Contact distance in Angstroms

    Returns:
        Dict mapping residue number to contact count
    """
    chain_a = structure[0][chain_a_id]
    chain_b = structure[0][chain_b_id]

    # Get atoms
    atoms_a = []
    for residue in chain_a:
        if residue.id[0] == ' ':
            for atom in residue:
                atoms_a.append((atom, residue.id[1]))

    atoms_b = []
    for residue in chain_b:
        if residue.id[0] == ' ':
            for atom in residue:
                atoms_b.append(atom)

    # Count contacts
    ns = NeighborSearch([atom for atom, _ in atoms_a] + atoms_b)

    contacts = {}
    for atom_a, resnum in atoms_a:
        nearby = ns.search(atom_a.coord, contact_cutoff)
        contact_count = sum(1 for atom in nearby if atom in atoms_b)

        if contact_count > 0:
            contacts[resnum] = contacts.get(resnum, 0) + contact_count

    return contacts


def estimate_buried_surface_area(structure, chain_a_id: str, chain_b_id: str,
                                 residue_num: int, probe_radius: float = 1.4) -> float:
    """
    Estimate buried surface area for a residue at the interface.

    Simplified calculation based on number of contacts and residue size.

    Args:
        structure: BioPython structure
        chain_a_id: Chain containing the residue
        chain_b_id: Partner chain
        residue_num: Residue number
        probe_radius: Probe radius (water = 1.4 Å)

    Returns:
        Approximate buried surface area in Ų
    """
    chain_a = structure[0][chain_a_id]

    # Get residue
    try:
        residue = chain_a[residue_num]
    except KeyError:
        return 0.0

    # Count heavy atoms as proxy for size
    num_atoms = len([atom for atom in residue if atom.element != 'H'])

    # Get interface contacts
    contacts = get_interface_contacts(structure, chain_a_id, chain_b_id, 4.5)
    num_contacts = contacts.get(residue_num, 0)

    # Empirical approximation: BSA ≈ 10-20 Ų per contact
    # Scale by residue size
    bsa = num_contacts * 15.0 * (num_atoms / 10.0)

    return bsa


def calculate_interface_stability(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],  # (position, wt_aa, mut_aa)
    config: Optional[Dict] = None,
    chain_a_id: str = "A",
    chain_b_id: str = "B"
) -> InterfaceStabilityResult:
    """
    Calculate estimated ΔΔG_interface for mutation(s).

    Args:
        pdb_path: Path to tetramer structure (e.g., 3KMD)
        mutations: List of (position, wt_aa, mut_aa) tuples
        config: Configuration dict
        chain_a_id: First chain in dimer
        chain_b_id: Second chain in dimer

    Returns:
        InterfaceStabilityResult
    """
    if config is None:
        config = load_config()

    params = config['tetramer_interface']

    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('tetramer', pdb_path)

    # Get WT interface residues and contacts
    interface_a, interface_b = get_interface_residues(
        structure, chain_a_id, chain_b_id, params['interface_cutoff']
    )

    wt_contacts = get_interface_contacts(
        structure, chain_a_id, chain_b_id, params['contact_cutoff']
    )

    # Calculate WT buried surface area (sum over all interface residues)
    wt_bsa = sum([
        estimate_buried_surface_area(structure, chain_a_id, chain_b_id, resnum)
        for resnum in interface_a
    ])

    # Analyze each mutation
    ddg_components = []
    total_clashes = 0
    total_bsa_loss = 0.0

    for position, wt_aa, mut_aa in mutations:
        # Check if this position is at the interface
        is_interface = position in interface_a or position in interface_b
        wt_num_contacts = wt_contacts.get(position, 0)

        if not is_interface:
            # Not at interface - minimal effect
            lost_contacts = 0
            gained_contacts = 0
            clashes = 0
            bsa_loss = 0.0
        else:
            # At interface - assess mutation impact (heuristic)

            # Contact changes based on residue properties
            if mut_aa in ['G', 'A']:  # Small residues - lose contacts
                lost_contacts = max(1, wt_num_contacts // 2)
                gained_contacts = 0
            elif mut_aa in ['W', 'F', 'Y', 'R', 'K']:  # Large residues - may clash
                lost_contacts = 0
                gained_contacts = 1  # May be unfavorable
                clashes = 1
            elif mut_aa == 'P':  # Proline - disruptive
                lost_contacts = wt_num_contacts
                gained_contacts = 0
                clashes = 0
            elif wt_aa in ['L', 'I', 'V'] and mut_aa in ['L', 'I', 'V']:  # Conservative
                lost_contacts = 0
                gained_contacts = 0
                clashes = 0
            else:  # Moderate change
                lost_contacts = max(0, wt_num_contacts // 3)
                gained_contacts = 0
                clashes = 0

            # Estimate buried surface area loss
            if lost_contacts > 0:
                res_bsa = estimate_buried_surface_area(
                    structure, chain_a_id, chain_b_id, position
                )
                bsa_loss = res_bsa * (lost_contacts / max(1, wt_num_contacts))
            else:
                bsa_loss = 0.0

            total_clashes += clashes
            total_bsa_loss += bsa_loss

        # Calculate ΔΔG component
        ddg = (lost_contacts * params['weights']['contact_loss'] +
               gained_contacts * params['weights']['contact_gain'] +
               clashes * params['weights']['clash'] +
               bsa_loss * params['weights']['buried_area_loss'])

        ddg_components.append(ddg)

    # Total ΔΔG_interface
    ddg_interface = sum(ddg_components)

    # Categorize
    thresholds = params['thresholds']
    if ddg_interface <= thresholds['good']:
        category = "good"
    elif ddg_interface <= thresholds['acceptable']:
        category = "acceptable"
    else:
        category = "bad"

    # Construct result
    result = InterfaceStabilityResult(
        ddg_interface=ddg_interface,
        wt_interface_contacts=sum(wt_contacts.values()),
        mut_interface_contacts=0,  # Placeholder - would need modeling
        lost_contacts=sum([1 for comp in ddg_components if comp > 0]),
        gained_contacts=0,
        wt_buried_area=wt_bsa,
        mut_buried_area=wt_bsa - total_bsa_loss,
        clashes=total_clashes,
        category=category,
        details={
            'mutations': mutations,
            'ddg_components': ddg_components,
            'interface_residues_chain_a': sorted(list(interface_a)),
            'interface_residues_chain_b': sorted(list(interface_b)),
            'method': 'structural_heuristic',
            'note': 'Heuristic estimate - validate with FoldX or experiment'
        }
    )

    return result


def calculate_interface_stability_cached(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],
    config: Optional[Dict] = None,
    chain_a_id: str = "A",
    chain_b_id: str = "B",
    cache_dir: str = "Data/processed/cache/functional_scoring"
) -> InterfaceStabilityResult:
    """Calculate interface stability with caching."""
    # Create cache key
    mut_str = ",".join([f"{pos}{wt}{mut}" for pos, wt, mut in mutations])
    cache_key = hashlib.sha256(
        f"{pdb_path}:{chain_a_id}:{chain_b_id}:{mut_str}".encode()
    ).hexdigest()

    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    cache_file = Path(cache_dir) / f"interface_{cache_key}.json"

    # Check cache
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)
        return InterfaceStabilityResult(**cached)

    # Calculate
    result = calculate_interface_stability(
        pdb_path, mutations, config, chain_a_id, chain_b_id
    )

    # Save to cache
    with open(cache_file, 'w') as f:
        # Convert sets to lists for JSON serialization
        result_dict = result.__dict__.copy()
        if 'interface_residues_chain_a' in result.details:
            result_dict['details']['interface_residues_chain_a'] = \
                list(result.details['interface_residues_chain_a'])
        if 'interface_residues_chain_b' in result.details:
            result_dict['details']['interface_residues_chain_b'] = \
                list(result.details['interface_residues_chain_b'])
        json.dump(result_dict, f, indent=2)

    return result


# Example usage and testing
if __name__ == "__main__":
    import sys

    pdb_path = "Data/raw/experimental_pdbs/3KMD.pdb"

    if Path(pdb_path).exists():
        # Test 1: S95A (near N-terminal interface - should have some effect)
        print("=== Test 1: S95A (near N-terminal interface) ===")
        result = calculate_interface_stability(
            pdb_path,
            [(95, 'S', 'A')],
            chain_a_id="A",
            chain_b_id="B"
        )
        print(f"ΔΔG_interface: {result.ddg_interface:.2f} kcal/mol")
        print(f"Category: {result.category}")
        print(f"WT interface contacts: {result.wt_interface_contacts}")
        print(f"WT buried surface area: {result.wt_buried_area:.1f} Ų")
        print(f"Interface residues (chain A): {len(result.details['interface_residues_chain_a'])}")
        print()

        # Test 2: M133L (buried core, away from interface - should be good)
        print("=== Test 2: M133L (buried core, not at interface) ===")
        result = calculate_interface_stability(
            pdb_path,
            [(133, 'M', 'L')],
            chain_a_id="A",
            chain_b_id="B"
        )
        print(f"ΔΔG_interface: {result.ddg_interface:.2f} kcal/mol")
        print(f"Category: {result.category}")
        print(f"Lost contacts: {result.lost_contacts}")
        print()

        # Test 3: Known interface residue (R181 - dimer interface)
        print("=== Test 3: R181E (interface residue - should be disruptive) ===")
        result = calculate_interface_stability(
            pdb_path,
            [(181, 'R', 'E')],
            chain_a_id="A",
            chain_b_id="B"
        )
        print(f"ΔΔG_interface: {result.ddg_interface:.2f} kcal/mol")
        print(f"Category: {result.category}")
        print(f"Lost contacts: {result.lost_contacts}")

    else:
        print(f"Error: {pdb_path} not found")
        print("Please ensure 3KMD.pdb is in Data/raw/experimental_pdbs/")
        sys.exit(1)
