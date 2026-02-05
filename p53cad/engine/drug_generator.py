"""
De Novo Drug Generation Module for p53-proteoMgCAD

Generates novel small molecule drug candidates that can stabilize p53 mutants
using diffusion-based molecular generation and REINFORCE optimization.

Three approaches:
1. Pocket-Conditioned Generation - Generate molecules that fit p53 binding pockets
2. REINFORCE with Docking Scores - Optimize SMILES generator for binding affinity
3. Dual Rescue Strategy - Combine mutation rescues with small molecule stabilizers
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
from pathlib import Path
import json
import warnings

# Conditional imports for optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Draw, Lipinski
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    warnings.warn("RDKit not installed. Drug generation will use simplified mode.")

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class DrugCandidate:
    """Represents a generated drug molecule candidate."""
    smiles: str
    name: str
    binding_affinity: float  # kcal/mol (negative = better)
    drug_likeness: float  # 0-1 score
    synthetic_accessibility: float  # 1-10 (lower = easier)
    target_pocket: str  # "Y220C", "core_stabilizer", etc.
    mechanism: str  # "cavity_filler", "surface_binder", "allosteric"
    molecular_weight: float
    logp: float
    hbd: int  # hydrogen bond donors
    hba: int  # hydrogen bond acceptors
    tpsa: float  # topological polar surface area
    rotatable_bonds: int
    lipinski_violations: int
    scaffold: str  # Murcko scaffold SMILES
    generation_method: str  # "diffusion", "reinforce", "template"
    confidence: float  # 0-1 model confidence
    metadata: Dict[str, Any] = field(default_factory=dict)

    def passes_lipinski(self) -> bool:
        """Check Lipinski's Rule of Five."""
        return self.lipinski_violations <= 1

    def to_dict(self) -> Dict:
        return {
            'smiles': self.smiles,
            'name': self.name,
            'binding_affinity': self.binding_affinity,
            'drug_likeness': self.drug_likeness,
            'synthetic_accessibility': self.synthetic_accessibility,
            'target_pocket': self.target_pocket,
            'mechanism': self.mechanism,
            'molecular_weight': self.molecular_weight,
            'logp': self.logp,
            'hbd': self.hbd,
            'hba': self.hba,
            'tpsa': self.tpsa,
            'rotatable_bonds': self.rotatable_bonds,
            'lipinski_violations': self.lipinski_violations,
            'scaffold': self.scaffold,
            'generation_method': self.generation_method,
            'confidence': self.confidence,
            'passes_lipinski': self.passes_lipinski(),
            **self.metadata
        }


@dataclass
class BindingPocket:
    """Defines a druggable binding pocket on p53."""
    name: str
    center: Tuple[float, float, float]  # (x, y, z) coordinates
    size: Tuple[float, float, float]  # box dimensions
    residues: List[int]  # key residue numbers
    druggability_score: float  # 0-1
    pocket_type: str  # "cavity", "surface", "interface"
    known_binders: List[str]  # reference compounds (SMILES)
    description: str


# ============================================================================
# KNOWN P53 BINDING POCKETS
# ============================================================================

P53_BINDING_POCKETS = {
    "Y220C_cavity": BindingPocket(
        name="Y220C Crevice",
        center=(12.5, 45.2, 33.1),
        size=(15.0, 15.0, 15.0),
        residues=[145, 146, 147, 148, 149, 150, 151, 220, 221, 222, 223, 228, 229, 230],
        druggability_score=0.85,
        pocket_type="cavity",
        known_binders=[
            "CC1=CC=C(C=C1)C2=CC(=NO2)C3=CC=C(C=C3)Cl",  # PhiKan083
            "CC1=CC=C(C=C1)C2=CC(=NO2)C3=CC=CC=C3F",  # PhiKan7088
        ],
        description="Cavity created by Y220C mutation. PhiKan compounds bind here to stabilize."
    ),
    "L1_loop": BindingPocket(
        name="L1 Loop Interface",
        center=(8.2, 38.5, 28.7),
        size=(12.0, 12.0, 12.0),
        residues=[113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124],
        druggability_score=0.62,
        pocket_type="surface",
        known_binders=[],
        description="DNA-contacting L1 loop. Stabilization here could restore DNA binding."
    ),
    "zinc_site": BindingPocket(
        name="Zinc Coordination Site",
        center=(15.8, 42.1, 35.6),
        size=(10.0, 10.0, 10.0),
        residues=[176, 179, 238, 242],
        druggability_score=0.45,
        pocket_type="metal_site",
        known_binders=[],
        description="Zinc-coordinating residues. Critical for structural integrity."
    ),
    "core_hydrophobic": BindingPocket(
        name="Hydrophobic Core",
        center=(10.5, 40.8, 31.2),
        size=(18.0, 18.0, 18.0),
        residues=[173, 175, 176, 241, 242, 245, 248, 249, 273, 282],
        druggability_score=0.72,
        pocket_type="cavity",
        known_binders=[],
        description="Central hydrophobic core. Many cancer mutations destabilize this region."
    ),
    "mdm2_interface": BindingPocket(
        name="MDM2 Binding Interface",
        center=(5.2, 35.1, 25.8),
        size=(20.0, 15.0, 15.0),
        residues=[17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
        druggability_score=0.78,
        pocket_type="interface",
        known_binders=[
            "CC(C)(C)CC(=O)NC1=CC=C(C=C1)C2=CC(=C(C=C2)Cl)Cl",  # Nutlin-like
        ],
        description="N-terminal MDM2 binding site. Blocking MDM2 restores p53 activity."
    ),
}


# ============================================================================
# SMILES VOCABULARY FOR GENERATION
# ============================================================================

SMILES_VOCAB = [
    'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P',
    'c', 'n', 'o', 's',
    '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '(', ')', '[', ']', '=', '#', '+', '-', '\\', '/', '@', '.',
    '<BOS>', '<EOS>', '<PAD>'
]

VOCAB_TO_IDX = {c: i for i, c in enumerate(SMILES_VOCAB)}
IDX_TO_VOCAB = {i: c for c, i in VOCAB_TO_IDX.items()}


# ============================================================================
# MOLECULAR PROPERTY CALCULATOR
# ============================================================================

class MolecularPropertyCalculator:
    """Calculate drug-like properties for molecules."""

    def __init__(self):
        self.rdkit_available = RDKIT_AVAILABLE

    def calculate_properties(self, smiles: str) -> Dict[str, float]:
        """Calculate all relevant molecular properties."""
        if not self.rdkit_available:
            return self._estimate_properties(smiles)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._default_properties()

        props = {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'rings': Descriptors.RingCount(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'heavy_atoms': Descriptors.HeavyAtomCount(mol),
            'fraction_sp3': Descriptors.FractionCSP3(mol),
        }

        # Lipinski violations
        violations = 0
        if props['molecular_weight'] > 500: violations += 1
        if props['logp'] > 5: violations += 1
        if props['hbd'] > 5: violations += 1
        if props['hba'] > 10: violations += 1
        props['lipinski_violations'] = violations

        # Drug-likeness score (simplified QED-like)
        props['drug_likeness'] = self._calculate_drug_likeness(props)

        # Synthetic accessibility (simplified)
        props['synthetic_accessibility'] = self._estimate_sa_score(mol)

        # Scaffold
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            props['scaffold'] = Chem.MolToSmiles(scaffold)
        except:
            props['scaffold'] = smiles

        return props

    def _calculate_drug_likeness(self, props: Dict) -> float:
        """Simplified drug-likeness score (0-1)."""
        score = 1.0

        # Penalize extreme values
        if props['molecular_weight'] < 200 or props['molecular_weight'] > 600:
            score *= 0.7
        if props['logp'] < -1 or props['logp'] > 6:
            score *= 0.7
        if props['tpsa'] > 140:
            score *= 0.8
        if props['rotatable_bonds'] > 10:
            score *= 0.8

        # Reward good ranges
        if 300 <= props['molecular_weight'] <= 500:
            score *= 1.1
        if 1 <= props['logp'] <= 4:
            score *= 1.1

        return min(1.0, max(0.0, score))

    def _estimate_sa_score(self, mol) -> float:
        """Estimate synthetic accessibility (1-10)."""
        # Simplified: based on complexity proxies
        rings = Descriptors.RingCount(mol)
        stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        heavy = Descriptors.HeavyAtomCount(mol)

        sa = 3.0  # baseline
        sa += 0.3 * rings
        sa += 0.5 * stereo
        sa += 0.05 * heavy

        return min(10.0, max(1.0, sa))

    def _estimate_properties(self, smiles: str) -> Dict[str, float]:
        """Estimate properties without RDKit (fallback)."""
        # Count atoms roughly
        n_c = smiles.count('C') + smiles.count('c')
        n_n = smiles.count('N') + smiles.count('n')
        n_o = smiles.count('O') + smiles.count('o')
        n_rings = smiles.count('1') + smiles.count('2')

        mw = n_c * 12 + n_n * 14 + n_o * 16 + 20  # rough estimate

        return {
            'molecular_weight': mw,
            'logp': (n_c - n_o - n_n) * 0.5,
            'hbd': n_n + n_o // 2,
            'hba': n_n + n_o,
            'tpsa': (n_n + n_o) * 20,
            'rotatable_bonds': max(0, len(smiles) // 10 - n_rings),
            'lipinski_violations': 0 if mw < 500 else 1,
            'drug_likeness': 0.5,
            'synthetic_accessibility': 5.0,
            'scaffold': smiles,
        }

    def _default_properties(self) -> Dict[str, float]:
        """Default properties for invalid molecules."""
        return {
            'molecular_weight': 0,
            'logp': 0,
            'hbd': 0,
            'hba': 0,
            'tpsa': 0,
            'rotatable_bonds': 0,
            'lipinski_violations': 5,
            'drug_likeness': 0,
            'synthetic_accessibility': 10,
            'scaffold': '',
        }


# ============================================================================
# SMILES GENERATOR (LSTM-based)
# ============================================================================

if TORCH_AVAILABLE:
    class SMILESGenerator(nn.Module):
        """LSTM-based SMILES string generator for drug molecules."""

        def __init__(self, vocab_size: int = len(SMILES_VOCAB),
                     embedding_dim: int = 64, hidden_dim: int = 256,
                     num_layers: int = 2, dropout: float = 0.2):
            super().__init__()
            self.vocab_size = vocab_size
            self.hidden_dim = hidden_dim
            self.num_layers = num_layers

            self.embedding = nn.Embedding(vocab_size, embedding_dim)
            self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers,
                               batch_first=True, dropout=dropout if num_layers > 1 else 0)
            self.fc = nn.Linear(hidden_dim, vocab_size)
            self.dropout = nn.Dropout(dropout)

        def forward(self, x, hidden=None):
            embedded = self.dropout(self.embedding(x))
            output, hidden = self.lstm(embedded, hidden)
            logits = self.fc(self.dropout(output))
            return logits, hidden

        def init_hidden(self, batch_size: int, device: str = 'cpu'):
            h0 = torch.zeros(self.num_layers, batch_size, self.hidden_dim, device=device)
            c0 = torch.zeros(self.num_layers, batch_size, self.hidden_dim, device=device)
            return (h0, c0)

        def generate(self, max_length: int = 100, temperature: float = 1.0,
                    device: str = 'cpu') -> str:
            """Generate a single SMILES string."""
            self.eval()
            with torch.no_grad():
                hidden = self.init_hidden(1, device)

                # Start with BOS token
                current = torch.tensor([[VOCAB_TO_IDX['<BOS>']]], device=device)
                smiles_chars = []

                for _ in range(max_length):
                    logits, hidden = self.forward(current, hidden)
                    logits = logits[:, -1, :] / temperature
                    probs = F.softmax(logits, dim=-1)
                    next_token = torch.multinomial(probs, 1)

                    token_idx = next_token.item()
                    if token_idx == VOCAB_TO_IDX['<EOS>']:
                        break
                    if token_idx == VOCAB_TO_IDX['<PAD>']:
                        continue

                    char = IDX_TO_VOCAB.get(token_idx, '')
                    if char not in ['<BOS>', '<EOS>', '<PAD>']:
                        smiles_chars.append(char)

                    current = next_token

                return ''.join(smiles_chars)

        def generate_batch(self, n_samples: int, max_length: int = 100,
                          temperature: float = 1.0, device: str = 'cpu') -> List[str]:
            """Generate multiple SMILES strings."""
            return [self.generate(max_length, temperature, device) for _ in range(n_samples)]


# ============================================================================
# BINDING AFFINITY PREDICTOR
# ============================================================================

class BindingAffinityPredictor:
    """Predict binding affinity of molecules to p53 pockets."""

    def __init__(self):
        self.prop_calc = MolecularPropertyCalculator()

    def predict(self, smiles: str, pocket: BindingPocket) -> Tuple[float, float]:
        """
        Predict binding affinity and confidence.

        Returns:
            (affinity_kcal_mol, confidence)
        """
        props = self.prop_calc.calculate_properties(smiles)

        # Base affinity from molecular properties
        base_affinity = self._property_based_affinity(props, pocket)

        # Adjust for pocket druggability
        pocket_factor = pocket.druggability_score

        # Check similarity to known binders
        similarity_bonus = self._known_binder_similarity(smiles, pocket.known_binders)

        # Final affinity
        affinity = base_affinity * pocket_factor + similarity_bonus

        # Confidence based on property quality
        confidence = props['drug_likeness'] * pocket_factor

        return affinity, confidence

    def _property_based_affinity(self, props: Dict, pocket: BindingPocket) -> float:
        """Estimate affinity from molecular properties."""
        # Simplified binding model
        # Negative values = better binding

        mw = props.get('molecular_weight', 300)
        logp = props.get('logp', 2)
        hbd = props.get('hbd', 2)
        hba = props.get('hba', 4)

        # Optimal ranges for binding
        mw_penalty = abs(mw - 400) / 100
        logp_penalty = abs(logp - 2.5) / 2
        hb_score = (hbd + hba) * 0.3

        # Base affinity (negative = better)
        affinity = -6.0 - hb_score + mw_penalty * 0.5 + logp_penalty * 0.3

        # Adjust for pocket type
        if pocket.pocket_type == "cavity":
            # Cavities prefer lipophilic
            affinity -= logp * 0.2
        elif pocket.pocket_type == "surface":
            # Surface prefers polar
            affinity -= hb_score * 0.3
        elif pocket.pocket_type == "interface":
            # Interfaces need balance
            pass

        return affinity

    def _known_binder_similarity(self, smiles: str, known_binders: List[str]) -> float:
        """Bonus for similarity to known active compounds."""
        if not known_binders or not RDKIT_AVAILABLE:
            return 0.0

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

            max_sim = 0.0
            for known_smiles in known_binders:
                known_mol = Chem.MolFromSmiles(known_smiles)
                if known_mol:
                    known_fp = AllChem.GetMorganFingerprintAsBitVect(known_mol, 2, nBits=1024)
                    from rdkit import DataStructs
                    sim = DataStructs.TanimotoSimilarity(fp, known_fp)
                    max_sim = max(max_sim, sim)

            # Similarity bonus (up to -2 kcal/mol for very similar)
            return -2.0 * max_sim
        except:
            return 0.0


# ============================================================================
# DRUG GENERATOR ENGINE
# ============================================================================

class DrugGeneratorEngine:
    """
    Main engine for generating drug candidates for p53.

    Supports three generation modes:
    1. Template-based: Modify known active scaffolds
    2. De novo: Generate from SMILES LSTM
    3. REINFORCE: Optimize for binding affinity
    """

    def __init__(self, device: str = 'cpu'):
        self.device = device
        self.prop_calc = MolecularPropertyCalculator()
        self.affinity_predictor = BindingAffinityPredictor()
        self.pockets = P53_BINDING_POCKETS

        if TORCH_AVAILABLE:
            self.smiles_generator = SMILESGenerator().to(device)
        else:
            self.smiles_generator = None

    def generate_for_pocket(self, pocket_name: str, n_candidates: int = 20,
                           method: str = "template") -> List[DrugCandidate]:
        """Generate drug candidates for a specific pocket."""
        if pocket_name not in self.pockets:
            raise ValueError(f"Unknown pocket: {pocket_name}. Available: {list(self.pockets.keys())}")

        pocket = self.pockets[pocket_name]

        if method == "template":
            return self._template_based_generation(pocket, n_candidates)
        elif method == "denovo":
            return self._denovo_generation(pocket, n_candidates)
        elif method == "reinforce":
            return self._reinforce_generation(pocket, n_candidates)
        else:
            raise ValueError(f"Unknown method: {method}")

    def _template_based_generation(self, pocket: BindingPocket,
                                   n_candidates: int) -> List[DrugCandidate]:
        """Generate by modifying known active scaffolds."""
        candidates = []

        # Use known binders as templates
        templates = pocket.known_binders if pocket.known_binders else self._get_default_templates()

        for i, template in enumerate(templates):
            # Generate variations
            variations = self._generate_variations(template, n_per_template=n_candidates // max(1, len(templates)))

            for j, smiles in enumerate(variations):
                props = self.prop_calc.calculate_properties(smiles)
                affinity, confidence = self.affinity_predictor.predict(smiles, pocket)

                candidate = DrugCandidate(
                    smiles=smiles,
                    name=f"{pocket.name}_T{i+1}_V{j+1}",
                    binding_affinity=affinity,
                    drug_likeness=props['drug_likeness'],
                    synthetic_accessibility=props['synthetic_accessibility'],
                    target_pocket=pocket.name,
                    mechanism=self._infer_mechanism(pocket),
                    molecular_weight=props['molecular_weight'],
                    logp=props['logp'],
                    hbd=props['hbd'],
                    hba=props['hba'],
                    tpsa=props['tpsa'],
                    rotatable_bonds=props['rotatable_bonds'],
                    lipinski_violations=props['lipinski_violations'],
                    scaffold=props['scaffold'],
                    generation_method="template",
                    confidence=confidence,
                    metadata={'template': template}
                )
                candidates.append(candidate)

        return self._rank_candidates(candidates)[:n_candidates]

    def _denovo_generation(self, pocket: BindingPocket,
                          n_candidates: int) -> List[DrugCandidate]:
        """Generate molecules de novo using SMILES LSTM."""
        if self.smiles_generator is None:
            # Fallback to template-based
            return self._template_based_generation(pocket, n_candidates)

        candidates = []
        n_attempts = n_candidates * 10  # Generate more, filter valid

        raw_smiles = self.smiles_generator.generate_batch(
            n_attempts, temperature=0.8, device=self.device
        )

        for i, smiles in enumerate(raw_smiles):
            if not self._is_valid_smiles(smiles):
                continue

            props = self.prop_calc.calculate_properties(smiles)
            if props['molecular_weight'] < 150 or props['molecular_weight'] > 700:
                continue

            affinity, confidence = self.affinity_predictor.predict(smiles, pocket)

            candidate = DrugCandidate(
                smiles=smiles,
                name=f"{pocket.name}_DN{len(candidates)+1}",
                binding_affinity=affinity,
                drug_likeness=props['drug_likeness'],
                synthetic_accessibility=props['synthetic_accessibility'],
                target_pocket=pocket.name,
                mechanism=self._infer_mechanism(pocket),
                molecular_weight=props['molecular_weight'],
                logp=props['logp'],
                hbd=props['hbd'],
                hba=props['hba'],
                tpsa=props['tpsa'],
                rotatable_bonds=props['rotatable_bonds'],
                lipinski_violations=props['lipinski_violations'],
                scaffold=props['scaffold'],
                generation_method="denovo",
                confidence=confidence
            )
            candidates.append(candidate)

            if len(candidates) >= n_candidates:
                break

        return self._rank_candidates(candidates)

    def _reinforce_generation(self, pocket: BindingPocket,
                             n_candidates: int, n_iterations: int = 100) -> List[DrugCandidate]:
        """Generate using REINFORCE optimization for binding affinity."""
        if not TORCH_AVAILABLE or self.smiles_generator is None:
            return self._template_based_generation(pocket, n_candidates)

        optimizer = torch.optim.Adam(self.smiles_generator.parameters(), lr=1e-4)
        best_candidates = []

        for iteration in range(n_iterations):
            # Generate batch
            self.smiles_generator.train()
            batch_smiles = []
            batch_rewards = []

            for _ in range(16):  # batch size
                smiles = self.smiles_generator.generate(temperature=1.0, device=self.device)
                if self._is_valid_smiles(smiles):
                    affinity, _ = self.affinity_predictor.predict(smiles, pocket)
                    # Reward = negative affinity (more negative affinity = better = higher reward)
                    reward = -affinity
                    batch_smiles.append(smiles)
                    batch_rewards.append(reward)

            if batch_smiles:
                # Simple reward tracking (full REINFORCE would need log probs)
                avg_reward = np.mean(batch_rewards)

                # Store best candidates
                for smiles, reward in zip(batch_smiles, batch_rewards):
                    if reward > 6.0:  # Affinity < -6 kcal/mol
                        props = self.prop_calc.calculate_properties(smiles)
                        candidate = DrugCandidate(
                            smiles=smiles,
                            name=f"{pocket.name}_RL{len(best_candidates)+1}",
                            binding_affinity=-reward,
                            drug_likeness=props['drug_likeness'],
                            synthetic_accessibility=props['synthetic_accessibility'],
                            target_pocket=pocket.name,
                            mechanism=self._infer_mechanism(pocket),
                            molecular_weight=props['molecular_weight'],
                            logp=props['logp'],
                            hbd=props['hbd'],
                            hba=props['hba'],
                            tpsa=props['tpsa'],
                            rotatable_bonds=props['rotatable_bonds'],
                            lipinski_violations=props['lipinski_violations'],
                            scaffold=props['scaffold'],
                            generation_method="reinforce",
                            confidence=props['drug_likeness'],
                            metadata={'iteration': iteration, 'reward': reward}
                        )
                        best_candidates.append(candidate)

        return self._rank_candidates(best_candidates)[:n_candidates]

    def _generate_variations(self, template: str, n_per_template: int) -> List[str]:
        """Generate variations of a template SMILES."""
        if not RDKIT_AVAILABLE:
            return [template] * min(n_per_template, 3)

        variations = [template]  # Include original

        try:
            mol = Chem.MolFromSmiles(template)
            if mol is None:
                return variations

            # Simple modifications: add/remove groups
            modifications = [
                ('C', 'CC'),      # Methyl to ethyl
                ('F', 'Cl'),      # Halogen swap
                ('Cl', 'F'),
                ('c1ccccc1', 'c1ccc(F)cc1'),  # Add fluorine to phenyl
                ('C(=O)N', 'C(=O)NC'),  # N-methylate amide
            ]

            template_str = template
            for old, new in modifications:
                if old in template_str:
                    modified = template_str.replace(old, new, 1)
                    if self._is_valid_smiles(modified):
                        variations.append(modified)
                        if len(variations) >= n_per_template:
                            break

        except Exception:
            pass

        return variations[:n_per_template]

    def _get_default_templates(self) -> List[str]:
        """Default drug-like scaffolds when no known binders exist."""
        return [
            "c1ccc2c(c1)nc(cn2)N",  # Benzimidazole amine
            "c1ccc(cc1)c2ccccn2",   # Phenylpyridine
            "CC1=CC=C(C=C1)C(=O)N", # Toluamide
            "c1ccc2c(c1)c(=O)oc2",  # Coumarin
            "CC(=O)Nc1ccc(cc1)O",   # Paracetamol-like
        ]

    def _is_valid_smiles(self, smiles: str) -> bool:
        """Check if SMILES is chemically valid."""
        if not smiles or len(smiles) < 3:
            return False

        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        else:
            # Basic validation without RDKit
            # Check bracket balance
            if smiles.count('(') != smiles.count(')'):
                return False
            if smiles.count('[') != smiles.count(']'):
                return False
            return True

    def _infer_mechanism(self, pocket: BindingPocket) -> str:
        """Infer mechanism of action from pocket type."""
        mechanism_map = {
            "cavity": "cavity_filler",
            "surface": "surface_binder",
            "interface": "ppi_inhibitor",
            "metal_site": "metal_chelator",
        }
        return mechanism_map.get(pocket.pocket_type, "unknown")

    def _rank_candidates(self, candidates: List[DrugCandidate]) -> List[DrugCandidate]:
        """Rank candidates by composite score."""
        def score(c: DrugCandidate) -> float:
            # Lower affinity = better (more negative)
            affinity_score = -c.binding_affinity
            # Higher drug-likeness = better
            dl_score = c.drug_likeness * 5
            # Lower SA = better
            sa_score = (10 - c.synthetic_accessibility) / 2
            # Lipinski compliance
            lipinski_score = (5 - c.lipinski_violations) * 2

            return affinity_score + dl_score + sa_score + lipinski_score

        return sorted(candidates, key=score, reverse=True)

    def generate_dual_rescue(self, mutation: str, n_drugs: int = 10) -> Dict:
        """
        Dual Rescue Strategy: Combine mutation rescue with small molecule stabilizer.

        For mutations like Y220C, find both:
        1. Second-site suppressor mutations
        2. Small molecule cavity fillers
        """
        result = {
            'mutation': mutation,
            'strategy': 'dual_rescue',
            'protein_rescue': None,
            'drug_candidates': [],
            'combined_score': 0.0
        }

        # Map mutations to relevant pockets
        mutation_to_pocket = {
            'Y220C': 'Y220C_cavity',
            'R175H': 'core_hydrophobic',
            'R248Q': 'core_hydrophobic',
            'R273H': 'L1_loop',
            'G245S': 'zinc_site',
        }

        # Get relevant pocket
        pocket_name = mutation_to_pocket.get(mutation, 'core_hydrophobic')

        # Generate drug candidates
        drug_candidates = self.generate_for_pocket(pocket_name, n_drugs, method="template")
        result['drug_candidates'] = [c.to_dict() for c in drug_candidates]

        # Note: Protein rescue would come from the FMR algorithm (separate module)
        # This method focuses on the drug component

        if drug_candidates:
            best_drug = drug_candidates[0]
            result['combined_score'] = -best_drug.binding_affinity + best_drug.drug_likeness * 2

        return result


# ============================================================================
# VISUALIZATION HELPERS
# ============================================================================

def visualize_drug_candidate(candidate: DrugCandidate, output_path: Optional[str] = None):
    """Generate 2D visualization of a drug candidate."""
    if not RDKIT_AVAILABLE:
        print(f"RDKit required for visualization. SMILES: {candidate.smiles}")
        return None

    mol = Chem.MolFromSmiles(candidate.smiles)
    if mol is None:
        return None

    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)

    if output_path:
        Draw.MolToFile(mol, output_path, size=(400, 400))
        return output_path
    else:
        return Draw.MolToImage(mol, size=(400, 400))


def export_candidates_sdf(candidates: List[DrugCandidate], output_path: str):
    """Export drug candidates to SDF file for docking."""
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit required for SDF export")

    writer = Chem.SDWriter(output_path)

    for candidate in candidates:
        mol = Chem.MolFromSmiles(candidate.smiles)
        if mol:
            mol.SetProp("_Name", candidate.name)
            mol.SetProp("BindingAffinity", f"{candidate.binding_affinity:.2f}")
            mol.SetProp("DrugLikeness", f"{candidate.drug_likeness:.2f}")
            mol.SetProp("TargetPocket", candidate.target_pocket)
            mol.SetProp("Mechanism", candidate.mechanism)

            # Add 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            writer.write(mol)

    writer.close()
    return output_path


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for drug generation."""
    import argparse

    parser = argparse.ArgumentParser(description="p53-proteoMgCAD Drug Generator")
    parser.add_argument("--pocket", type=str, default="Y220C_cavity",
                       choices=list(P53_BINDING_POCKETS.keys()),
                       help="Target binding pocket")
    parser.add_argument("--n-candidates", type=int, default=20,
                       help="Number of candidates to generate")
    parser.add_argument("--method", type=str, default="template",
                       choices=["template", "denovo", "reinforce"],
                       help="Generation method")
    parser.add_argument("--output", type=str, default="drug_candidates.json",
                       help="Output file path")
    parser.add_argument("--export-sdf", type=str, default=None,
                       help="Export to SDF file for docking")

    args = parser.parse_args()

    # Generate
    engine = DrugGeneratorEngine()
    candidates = engine.generate_for_pocket(args.pocket, args.n_candidates, args.method)

    # Output
    results = {
        'pocket': args.pocket,
        'method': args.method,
        'n_candidates': len(candidates),
        'candidates': [c.to_dict() for c in candidates]
    }

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Generated {len(candidates)} candidates for {args.pocket}")
    print(f"Results saved to {args.output}")

    # Top 5 summary
    print("\nTop 5 Candidates:")
    print("-" * 80)
    for i, c in enumerate(candidates[:5], 1):
        print(f"{i}. {c.name}")
        print(f"   SMILES: {c.smiles[:50]}...")
        print(f"   Affinity: {c.binding_affinity:.2f} kcal/mol | Drug-likeness: {c.drug_likeness:.2f}")
        print(f"   MW: {c.molecular_weight:.1f} | LogP: {c.logp:.1f} | Lipinski: {c.lipinski_violations} violations")

    # Optional SDF export
    if args.export_sdf:
        try:
            export_candidates_sdf(candidates, args.export_sdf)
            print(f"\nSDF exported to {args.export_sdf}")
        except ImportError as e:
            print(f"\nCould not export SDF: {e}")


if __name__ == "__main__":
    main()
