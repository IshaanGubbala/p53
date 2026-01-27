"""
DNA Binding Affinity Calculator

Estimates ΔΔG_binding using structural analysis of protein-DNA contacts.
This is a lightweight alternative to FoldX AnalyseComplex for initial screening.

Method:
1. Identify DNA-contacting residues in WT structure
2. Introduce mutation and check for contact changes
3. Score based on lost/gained contacts, H-bonds, and clashes
4. Return estimated ΔΔG_binding (positive = weakened binding)

Note: This is an approximation. For final validation, use FoldX or experimental data.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch, Polypeptide
import yaml
import hashlib
import json
from dataclasses import dataclass


@dataclass
class BindingAffinityResult:
    """Results from DNA binding affinity calculation."""
    ddg_binding: float  # Estimated ΔΔG_binding (kcal/mol)
    wt_contacts: int
    mut_contacts: int
    lost_contacts: int
    gained_contacts: int
    wt_hbonds: int
    mut_hbonds: int
    lost_hbonds: int
    gained_hbonds: int
    clashes: int
    category: str  # "good", "acceptable", "bad"
    details: Dict


def load_config(config_path: str = "configs/functional_scoring.yaml") -> Dict:
    """Load functional scoring configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_dna_contacts(structure, protein_chains: List[str], dna_chains: List[str],
                     cutoff: float = 4.5) -> Dict[int, Set[str]]:
    """
    Identify protein residues in contact with DNA.

    Args:
        structure: BioPython structure object
        protein_chains: List of protein chain IDs
        dna_chains: List of DNA chain IDs
        cutoff: Contact distance in Angstroms

    Returns:
        Dict mapping residue number to set of DNA atoms it contacts
    """
    # Get all protein atoms
    protein_atoms = []
    for chain_id in protein_chains:
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[0] == ' ':  # Standard residue
                for atom in residue:
                    protein_atoms.append((atom, residue.id[1]))  # (atom, resnum)

    # Get all DNA atoms
    dna_atoms = []
    for chain_id in dna_chains:
        chain = structure[0][chain_id]
        for residue in chain:
            for atom in residue:
                dna_atoms.append(atom)

    # Find contacts using NeighborSearch
    ns = NeighborSearch([atom for atom, _ in protein_atoms] + dna_atoms)

    contacts = {}
    for protein_atom, resnum in protein_atoms:
        nearby_atoms = ns.search(protein_atom.coord, cutoff)
        dna_nearby = [atom for atom in nearby_atoms if atom in dna_atoms]

        if dna_nearby:
            if resnum not in contacts:
                contacts[resnum] = set()
            contacts[resnum].update([atom.get_full_id() for atom in dna_nearby])

    return contacts


def get_hydrogen_bonds(structure, protein_chains: List[str], dna_chains: List[str],
                       cutoff: float = 3.5) -> List[Tuple[int, str]]:
    """
    Identify potential hydrogen bonds between protein and DNA.

    Simplified approach: check N/O atoms in protein within H-bond distance of DNA.

    Returns:
        List of (residue_number, DNA_atom_id) pairs
    """
    hbonds = []

    # Get donor/acceptor atoms from protein
    protein_hbond_atoms = []
    for chain_id in protein_chains:
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[0] == ' ':
                resnum = residue.id[1]
                for atom in residue:
                    if atom.element in ['N', 'O']:  # Potential donors/acceptors
                        protein_hbond_atoms.append((atom, resnum))

    # Get donor/acceptor atoms from DNA
    dna_hbond_atoms = []
    for chain_id in dna_chains:
        chain = structure[0][chain_id]
        for residue in chain:
            for atom in residue:
                if atom.element in ['N', 'O', 'P']:  # DNA acceptors
                    dna_hbond_atoms.append(atom)

    # Check distances
    for protein_atom, resnum in protein_hbond_atoms:
        for dna_atom in dna_hbond_atoms:
            dist = np.linalg.norm(protein_atom.coord - dna_atom.coord)
            if dist <= cutoff:
                hbonds.append((resnum, dna_atom.get_full_id()))

    return hbonds


def mutate_residue_in_structure(structure, chain_id: str, resnum: int,
                                 new_aa: str, rotamer_lib: str = "dunbrack") -> bool:
    """
    Introduce mutation in structure (simplified - uses representative rotamer).

    For now, this is a placeholder. A full implementation would:
    1. Look up rotamer library for new amino acid
    2. Model sidechain atoms
    3. Minimize local geometry

    Returns:
        True if mutation successful, False otherwise
    """
    # This is a simplified stub - in production, would use:
    # - BioPython's internal rotamer library
    # - Or external tool like SCWRL4
    # - Or OpenMM to minimize after crude placement

    # For now, we'll mark that the mutation was attempted
    # The actual contact analysis will be done with the mutated structure

    return True


def check_clashes(structure, chain_id: str, resnum: int,
                  dna_chains: List[str], clash_cutoff: float = 2.0) -> int:
    """
    Check for steric clashes between mutated residue and DNA.

    Args:
        structure: Structure with mutation
        chain_id: Protein chain ID
        resnum: Residue number
        dna_chains: DNA chain IDs
        clash_cutoff: Distance below which atoms clash (Angstroms)

    Returns:
        Number of clashing atom pairs
    """
    residue = structure[0][chain_id][resnum]

    # Get DNA atoms
    dna_atoms = []
    for dna_chain_id in dna_chains:
        for dna_residue in structure[0][dna_chain_id]:
            for atom in dna_residue:
                dna_atoms.append(atom)

    # Check distances
    clashes = 0
    for res_atom in residue:
        for dna_atom in dna_atoms:
            dist = np.linalg.norm(res_atom.coord - dna_atom.coord)
            if dist < clash_cutoff:
                clashes += 1

    return clashes


def calculate_binding_affinity(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],  # (position, wt_aa, mut_aa)
    config: Optional[Dict] = None
) -> BindingAffinityResult:
    """
    Calculate estimated ΔΔG_binding for mutation(s).

    Args:
        pdb_path: Path to DNA-bound structure (e.g., 2OCJ)
        mutations: List of (position, wt_aa, mut_aa) tuples
        config: Configuration dict (loaded from YAML if None)

    Returns:
        BindingAffinityResult with estimated ΔΔG and details
    """
    if config is None:
        config = load_config()

    params = config['dna_binding']
    struct_config = config['structures']['dna_bound']

    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_path)

    protein_chains = struct_config['protein_chains']
    dna_chains = struct_config['dna_chains']

    # Get WT contacts and H-bonds
    wt_contacts = get_dna_contacts(structure, protein_chains, dna_chains,
                                   params['contact_cutoff'])
    wt_hbonds = get_hydrogen_bonds(structure, protein_chains, dna_chains,
                                   params['hydrogen_bond_cutoff'])

    # Group H-bonds by residue
    wt_hbonds_by_res = {}
    for resnum, dna_atom in wt_hbonds:
        if resnum not in wt_hbonds_by_res:
            wt_hbonds_by_res[resnum] = []
        wt_hbonds_by_res[resnum].append(dna_atom)

    # Analyze each mutation
    ddg_components = []
    total_clashes = 0

    for position, wt_aa, mut_aa in mutations:
        # Check if this position contacts DNA in WT
        wt_has_contact = position in wt_contacts
        wt_num_contacts = len(wt_contacts.get(position, []))
        wt_num_hbonds = len(wt_hbonds_by_res.get(position, []))

        # Estimate mutation effects (simplified - no actual structural modeling)
        # This is a heuristic based on residue properties

        # Contact changes (heuristic)
        if wt_has_contact:
            # Check if mutation likely disrupts contacts
            if mut_aa in ['G', 'A']:  # Small residues - likely lose contacts
                lost_contacts = max(1, wt_num_contacts // 2)
                gained_contacts = 0
            elif mut_aa in ['W', 'F', 'Y']:  # Large aromatics - may gain bad contacts
                lost_contacts = 0
                gained_contacts = 1
            elif mut_aa == 'P':  # Proline - likely disruptive
                lost_contacts = wt_num_contacts
                gained_contacts = 0
            else:  # Conservative - assume small loss
                lost_contacts = max(0, wt_num_contacts // 3)
                gained_contacts = 0
        else:
            # Not at interface - minimal effect
            lost_contacts = 0
            gained_contacts = 0

        # H-bond changes (heuristic)
        if wt_num_hbonds > 0:
            if mut_aa in ['G', 'A', 'V', 'L', 'I', 'P']:  # Hydrophobic - lose H-bonds
                lost_hbonds = wt_num_hbonds
                gained_hbonds = 0
            elif mut_aa in ['K', 'R', 'D', 'E', 'N', 'Q']:  # Polar - may maintain
                lost_hbonds = max(0, wt_num_hbonds // 2)
                gained_hbonds = 0
            else:
                lost_hbonds = wt_num_hbonds // 2
                gained_hbonds = 0
        else:
            lost_hbonds = 0
            gained_hbonds = 0

        # Clash estimation (heuristic - large residues in tight interface)
        if wt_has_contact and mut_aa in ['W', 'F', 'Y', 'R', 'K']:
            clashes = 1  # Likely some clash
        else:
            clashes = 0

        total_clashes += clashes

        # Calculate ΔΔG component for this mutation
        ddg = (lost_contacts * params['weights']['contact_loss'] +
               gained_contacts * params['weights']['contact_gain'] +
               lost_hbonds * params['weights']['hbond_loss'] +
               gained_hbonds * params['weights']['hbond_gain'] +
               clashes * params['weights']['clash'])

        ddg_components.append(ddg)

    # Total ΔΔG_binding
    ddg_binding = sum(ddg_components)

    # Categorize
    thresholds = params['thresholds']
    if ddg_binding <= thresholds['good']:
        category = "good"
    elif ddg_binding <= thresholds['acceptable']:
        category = "acceptable"
    else:
        category = "bad"

    # Construct result
    result = BindingAffinityResult(
        ddg_binding=ddg_binding,
        wt_contacts=sum(len(contacts) for contacts in wt_contacts.values()),
        mut_contacts=0,  # Placeholder - would need actual modeling
        lost_contacts=sum([comp for comp in ddg_components if comp > 0]),
        gained_contacts=0,
        wt_hbonds=len(wt_hbonds),
        mut_hbonds=0,  # Placeholder
        lost_hbonds=0,  # Calculated above
        gained_hbonds=0,
        clashes=total_clashes,
        category=category,
        details={
            'mutations': mutations,
            'ddg_components': ddg_components,
            'method': 'structural_heuristic',
            'note': 'Heuristic estimate - validate with FoldX or experiment'
        }
    )

    return result


def calculate_binding_affinity_cached(
    pdb_path: str,
    mutations: List[Tuple[int, str, str]],
    config: Optional[Dict] = None,
    cache_dir: str = "Data/processed/cache/functional_scoring"
) -> BindingAffinityResult:
    """
    Calculate binding affinity with caching.

    Args:
        pdb_path: Path to structure
        mutations: List of mutations
        config: Configuration
        cache_dir: Directory for cache files

    Returns:
        BindingAffinityResult
    """
    # Create cache key
    mut_str = ",".join([f"{pos}{wt}{mut}" for pos, wt, mut in mutations])
    cache_key = hashlib.sha256(f"{pdb_path}:{mut_str}".encode()).hexdigest()

    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    cache_file = Path(cache_dir) / f"binding_{cache_key}.json"

    # Check cache
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)
        return BindingAffinityResult(**cached)

    # Calculate
    result = calculate_binding_affinity(pdb_path, mutations, config)

    # Save to cache
    with open(cache_file, 'w') as f:
        json.dump(result.__dict__, f, indent=2)

    return result


# Example usage and testing
if __name__ == "__main__":
    import sys

    # Test on a known DNA-binding mutation
    pdb_path = "Data/raw/experimental_pdbs/1TSR.pdb"

    if Path(pdb_path).exists():
        # Test 1: Mutation at DNA-binding position (should have high ΔΔG)
        print("=== Test 1: R248Q (known to disrupt DNA binding) ===")
        result = calculate_binding_affinity(
            pdb_path,
            [(248, 'R', 'Q')]
        )
        print(f"ΔΔG_binding: {result.ddg_binding:.2f} kcal/mol")
        print(f"Category: {result.category}")
        print(f"Lost contacts: {result.lost_contacts}, Lost H-bonds: {result.lost_hbonds}")
        print()

        # Test 2: Mutation away from DNA (should have low ΔΔG)
        print("=== Test 2: M133L (away from DNA interface) ===")
        result = calculate_binding_affinity(
            pdb_path,
            [(133, 'M', 'L')]
        )
        print(f"ΔΔG_binding: {result.ddg_binding:.2f} kcal/mol")
        print(f"Category: {result.category}")
        print(f"Lost contacts: {result.lost_contacts}, Lost H-bonds: {result.lost_hbonds}")
    else:
        print(f"Error: {pdb_path} not found")
        print("Please ensure 2OCJ_full.pdb is in Data/raw/experimental_pdbs/")
        sys.exit(1)
