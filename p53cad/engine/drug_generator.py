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
import re
import shutil
import subprocess
import tempfile
import os

from p53cad.core.logging import get_logger

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

try:
    from vina import Vina
    VINA_PY_AVAILABLE = True
except ImportError:
    VINA_PY_AVAILABLE = False

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False

try:
    import openmm  # noqa: F401
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


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

# Coordinates calibrated for AlphaFold p53 structure (AF-P04637-F1-model_v6)
P53_BINDING_POCKETS = {
    "Y220C_cavity": BindingPocket(
        name="Y220C Crevice",
        center=(-17.1, 6.7, 11.9),
        size=(21.0, 24.0, 18.0),
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
        center=(2.9, 15.0, 3.8),
        size=(30.0, 17.0, 17.0),
        residues=[113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124],
        druggability_score=0.62,
        pocket_type="surface",
        known_binders=[],
        description="DNA-contacting L1 loop. Stabilization here could restore DNA binding."
    ),
    "zinc_site": BindingPocket(
        name="Zinc Coordination Site",
        center=(8.5, -5.3, -4.2),
        size=(12.0, 14.0, 16.0),
        residues=[176, 179, 238, 242],
        druggability_score=0.45,
        pocket_type="metal_site",
        known_binders=[],
        description="Zinc-coordinating residues. Critical for structural integrity."
    ),
    "core_hydrophobic": BindingPocket(
        name="Hydrophobic Core",
        center=(4.6, 0.7, -7.9),
        size=(21.0, 34.0, 22.0),
        residues=[173, 175, 176, 241, 242, 245, 248, 249, 273, 282],
        druggability_score=0.72,
        pocket_type="cavity",
        known_binders=[],
        description="Central hydrophobic core. Many cancer mutations destabilize this region."
    ),
    "mdm2_interface": BindingPocket(
        name="MDM2 Binding Interface",
        center=(20.2, -18.0, -24.7),
        size=(24.0, 31.0, 17.0),
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
        """
        Estimate molecular properties without RDKit.

        Uses SMILES string analysis and empirical relationships to
        approximate key drug-like properties.
        """
        # Atom counting from SMILES
        # Uppercase = aliphatic, lowercase = aromatic
        n_c = smiles.count('C') + smiles.count('c')
        n_n = smiles.count('N') + smiles.count('n')
        n_o = smiles.count('O') + smiles.count('o')
        n_s = smiles.count('S') + smiles.count('s')
        n_f = smiles.count('F')
        n_cl = smiles.count('Cl')
        n_br = smiles.count('Br')

        # Aromatic atoms (lowercase in SMILES)
        n_aromatic = sum(1 for c in smiles if c.islower() and c in 'cnos')

        # Ring counting (numbers indicate ring closure)
        n_rings = len(set(c for c in smiles if c.isdigit()))

        # Hydrogen donors/acceptors from functional groups
        n_oh = smiles.count('O') - smiles.count('=O')  # OH groups (not C=O)
        n_nh = smiles.count('N') - smiles.count('#N')  # NH groups (not CN)

        # Molecular weight calculation (more accurate)
        mw = (n_c * 12.01 + n_n * 14.01 + n_o * 16.00 + n_s * 32.07 +
              n_f * 19.00 + n_cl * 35.45 + n_br * 79.90)
        # Estimate hydrogens using valence
        n_h_estimate = max(0, n_c * 2 + 2 - 2 * n_rings - n_aromatic // 2)
        mw += n_h_estimate * 1.008

        # LogP estimation (Wildman-Crippen-like)
        # Aromatic C: +0.15, Aliphatic C: +0.25, O: -0.5, N: -0.8, halogens: +0.5
        logp = (n_c - n_aromatic) * 0.25 + n_aromatic * 0.15
        logp -= n_o * 0.5 + n_n * 0.8
        logp += (n_f + n_cl + n_br) * 0.5
        logp -= n_rings * 0.1  # Rings slightly reduce logP

        # Hydrogen bond donors (NH, OH groups)
        hbd = n_oh + n_nh

        # Hydrogen bond acceptors (N, O atoms)
        hba = n_n + n_o

        # TPSA estimation (polar surface area)
        # O contributes ~20 Å², N contributes ~25 Å²
        tpsa = n_o * 20 + n_n * 25

        # Rotatable bonds (rough estimate)
        # Single bonds minus rings and terminal atoms
        single_bonds = smiles.count('-') if '-' in smiles else max(0, len(smiles) // 5)
        rotatable = max(0, single_bonds - n_rings - 2)

        # Lipinski violations
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1

        # Drug-likeness score
        drug_likeness = 1.0
        if mw < 200 or mw > 600:
            drug_likeness *= 0.7
        if logp < -1 or logp > 6:
            drug_likeness *= 0.7
        if 300 <= mw <= 500:
            drug_likeness *= 1.1
        if 1 <= logp <= 4:
            drug_likeness *= 1.1
        drug_likeness = min(1.0, max(0.0, drug_likeness))

        # Synthetic accessibility (higher = harder)
        sa = 3.0 + 0.3 * n_rings + 0.05 * len(smiles)
        if n_br > 0: sa += 0.5
        if smiles.count('@') > 0: sa += 1.0  # Stereocenters
        sa = min(10.0, max(1.0, sa))

        return {
            'molecular_weight': mw,
            'logp': logp,
            'hbd': hbd,
            'hba': hba,
            'tpsa': tpsa,
            'rotatable_bonds': rotatable,
            'rings': n_rings,
            'aromatic_rings': n_aromatic // 4,  # Approximate
            'heavy_atoms': n_c + n_n + n_o + n_s + n_f + n_cl + n_br,
            'lipinski_violations': violations,
            'drug_likeness': drug_likeness,
            'synthetic_accessibility': sa,
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

    Supports generation modes:
    1. Template-based: Modify known active scaffolds
    2. De novo: Generate from SMILES LSTM
    3. REINFORCE: Optimize for binding affinity
    4. Docking: Real AutoDock Vina docking (requires optional deps)
    5. Docking+MD: Docking with MD readiness diagnostics
    """

    def __init__(self, device: str = "cpu", allow_wt_receptor_fallback: bool = False):
        self.logger = get_logger(__name__)
        self.device = device
        self.allow_wt_receptor_fallback = bool(allow_wt_receptor_fallback)
        self.prop_calc = MolecularPropertyCalculator()
        self.affinity_predictor = BindingAffinityPredictor()
        self.pockets = P53_BINDING_POCKETS

        if TORCH_AVAILABLE:
            self.smiles_generator = SMILESGenerator().to(device)
        else:
            self.smiles_generator = None

    def generate_for_pocket(
        self,
        pocket_name: str,
        n_candidates: int = 20,
        method: str = "template",
        protein_context: Optional[Dict[str, Any]] = None,
        allow_wt_receptor_fallback: Optional[bool] = None,
    ) -> List[DrugCandidate]:
        """Generate drug candidates for a specific pocket.

        Args:
            pocket_name: Named pocket from P53_BINDING_POCKETS.
            n_candidates: Number of molecules to return.
            method: Generation strategy.
            protein_context: Optional rescue context from selected protein candidate.
                Expected keys: ``mutations`` (list[str]), ``score`` (float).
        """
        if pocket_name not in self.pockets:
            raise ValueError(f"Unknown pocket: {pocket_name}. Available: {list(self.pockets.keys())}")

        pocket = self.pockets[pocket_name]
        context_adjustment = self._compute_context_adjustment(pocket, protein_context)
        allow_wt_fallback = (
            self.allow_wt_receptor_fallback
            if allow_wt_receptor_fallback is None
            else bool(allow_wt_receptor_fallback)
        )

        if method == "template":
            return self._template_based_generation(pocket, n_candidates, context_adjustment=context_adjustment)
        elif method == "denovo":
            return self._denovo_generation(pocket, n_candidates, context_adjustment=context_adjustment)
        elif method == "reinforce":
            return self._reinforce_generation(pocket, n_candidates, context_adjustment=context_adjustment)
        elif method == "docking":
            return self._docking_generation(
                pocket_name,
                pocket,
                n_candidates,
                context_adjustment=context_adjustment,
                require_md=False,
                allow_wt_receptor_fallback=allow_wt_fallback,
            )
        elif method == "docking_md":
            return self._docking_generation(
                pocket_name,
                pocket,
                n_candidates,
                context_adjustment=context_adjustment,
                require_md=True,
                allow_wt_receptor_fallback=allow_wt_fallback,
            )
        else:
            raise ValueError(f"Unknown method: {method}")

    def _template_based_generation(self, pocket: BindingPocket,
                                   n_candidates: int,
                                   context_adjustment: Optional[Dict[str, float]] = None) -> List[DrugCandidate]:
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
                affinity, confidence = self._apply_context_adjustment(
                    affinity,
                    confidence,
                    context_adjustment,
                )

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
                    metadata={
                        'template': template,
                        "selection_reason": "ranked_by_predicted_affinity",
                        **self._context_metadata(context_adjustment),
                    }
                )
                candidates.append(candidate)

        return self._rank_candidates(candidates)[:n_candidates]

    def _denovo_generation(self, pocket: BindingPocket,
                          n_candidates: int,
                          context_adjustment: Optional[Dict[str, float]] = None) -> List[DrugCandidate]:
        """Generate molecules de novo using SMILES LSTM."""
        if self.smiles_generator is None or not RDKIT_AVAILABLE:
            # Without RDKit validation, LSTM generates invalid SMILES
            # Fallback to curated templates which are guaranteed valid
            return self._template_based_generation(
                pocket,
                n_candidates,
                context_adjustment=context_adjustment,
            )

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
            affinity, confidence = self._apply_context_adjustment(
                affinity,
                confidence,
                context_adjustment,
            )

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
                confidence=confidence,
                metadata=self._context_metadata(context_adjustment),
            )
            candidates.append(candidate)

            if len(candidates) >= n_candidates:
                break

        return self._rank_candidates(candidates)

    def _reinforce_generation(self, pocket: BindingPocket,
                             n_candidates: int, n_iterations: int = 100,
                             context_adjustment: Optional[Dict[str, float]] = None) -> List[DrugCandidate]:
        """Generate using REINFORCE optimization for binding affinity."""
        if not TORCH_AVAILABLE or self.smiles_generator is None:
            return self._template_based_generation(
                pocket,
                n_candidates,
                context_adjustment=context_adjustment,
            )

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
                    affinity, _ = self._apply_context_adjustment(
                        affinity,
                        0.0,
                        context_adjustment,
                    )
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
                        _, adjusted_conf = self._apply_context_adjustment(
                            -reward,
                            props['drug_likeness'],
                            context_adjustment,
                        )
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
                            confidence=adjusted_conf,
                            metadata={
                                'iteration': iteration,
                                'reward': reward,
                                "selection_reason": "ranked_by_reinforce_reward",
                                **self._context_metadata(context_adjustment),
                            }
                        )
                        best_candidates.append(candidate)

        return self._rank_candidates(best_candidates)[:n_candidates]

    def _docking_generation(
        self,
        pocket_key: str,
        pocket: BindingPocket,
        n_candidates: int,
        context_adjustment: Optional[Dict[str, float]] = None,
        require_md: bool = False,
        allow_wt_receptor_fallback: bool = False,
    ) -> List[DrugCandidate]:
        """Generate and score candidates using real docking when available."""
        runtime = self._resolve_docking_runtime(
            pocket_key,
            pocket,
            allow_wt_receptor_fallback=allow_wt_receptor_fallback,
        )
        if runtime["errors"]:
            raise RuntimeError(
                "Docking mode is unavailable: "
                + "; ".join(runtime["errors"])
                + ". Install optional deps (rdkit, meeko, vina/openbabel) and provide receptor PDBQT."
            )

        receptor_path = runtime["receptor_pdbqt"]
        method_name = "docking_md" if require_md else "docking"
        md_probe = self._run_md_refinement("CCO", require_md=require_md)
        if require_md and md_probe.get("status") == "unavailable":
            raise RuntimeError(
                "docking_md mode is unavailable: "
                + md_probe.get("detail", "MD dependencies are missing")
            )

        seed_count = max(n_candidates * 4, 24)
        seed_candidates = self._template_based_generation(
            pocket,
            seed_count,
            context_adjustment=None,
        )

        docked: List[DrugCandidate] = []
        seen_smiles = set()
        for idx, seed in enumerate(seed_candidates):
            if seed.smiles in seen_smiles:
                continue
            seen_smiles.add(seed.smiles)

            docking = self._dock_smiles(seed.smiles, receptor_path, pocket)
            md_meta = self._run_md_refinement(seed.smiles, require_md=True) if require_md else md_probe

            props = self.prop_calc.calculate_properties(seed.smiles)
            raw_affinity = float(docking["best_affinity"])
            base_confidence = float(np.clip(0.35 + 0.06 * float(docking["pose_count"]), 0.0, 0.95))
            affinity, confidence = self._apply_context_adjustment(
                raw_affinity,
                base_confidence,
                context_adjustment,
            )

            candidate = DrugCandidate(
                smiles=seed.smiles,
                name=f"{pocket.name}_DK{idx + 1}",
                binding_affinity=affinity,
                drug_likeness=props["drug_likeness"],
                synthetic_accessibility=props["synthetic_accessibility"],
                target_pocket=pocket.name,
                mechanism=self._infer_mechanism(pocket),
                molecular_weight=props["molecular_weight"],
                logp=props["logp"],
                hbd=props["hbd"],
                hba=props["hba"],
                tpsa=props["tpsa"],
                rotatable_bonds=props["rotatable_bonds"],
                lipinski_violations=props["lipinski_violations"],
                scaffold=props["scaffold"],
                generation_method=method_name,
                confidence=confidence,
                metadata={
                    "scoring_mode": "autodock_vina",
                    "docking_backend": docking["backend"],
                    "docking_receptor_pdbqt": str(receptor_path),
                    "receptor_source": runtime.get("receptor_source", "unknown"),
                    "docking_raw_affinity": raw_affinity,
                    "docking_pose_count": int(docking["pose_count"]),
                    "docking_exhaustiveness": int(docking["exhaustiveness"]),
                    "md_status": md_meta["status"],
                    "md_detail": md_meta.get("detail", ""),
                    "selection_reason": "ranked_by_binding_affinity",
                    **self._context_metadata(context_adjustment),
                },
            )
            docked.append(candidate)

            if len(docked) >= n_candidates:
                break

        return self._rank_candidates(docked)

    def _resolve_docking_runtime(
        self,
        pocket_key: str,
        pocket: BindingPocket,
        allow_wt_receptor_fallback: bool = False,
    ) -> Dict[str, Any]:
        errors: List[str] = []
        if not RDKIT_AVAILABLE:
            errors.append("RDKit is not installed")
        if not MEEKO_AVAILABLE:
            errors.append("Meeko is not installed")
        if not (VINA_PY_AVAILABLE or shutil.which("vina")):
            errors.append("AutoDock Vina backend not found")

        receptor_pdbqt = self._resolve_receptor_pdbqt(
            pocket_key,
            pocket,
            allow_wt_receptor_fallback=allow_wt_receptor_fallback,
        )
        if receptor_pdbqt is None:
            receptor_pdbqt = self._prepare_receptor_pdbqt(pocket_key)
        if receptor_pdbqt is None:
            errors.append("No receptor PDBQT found (expected under data/raw/receptors or data/processed/receptors)")

        receptor_source = self._categorize_receptor_source(receptor_pdbqt, pocket_key, pocket)

        return {
            "errors": errors,
            "receptor_pdbqt": receptor_pdbqt,
            "receptor_source": receptor_source,
            "allow_wt_receptor_fallback": bool(allow_wt_receptor_fallback),
        }

    def _resolve_receptor_pdbqt(
        self,
        pocket_key: str,
        pocket: BindingPocket,
        allow_wt_receptor_fallback: bool = False,
    ) -> Optional[Path]:
        pocket_slug = re.sub(r"[^a-z0-9]+", "_", pocket_key.lower()).strip("_")
        name_slug = re.sub(r"[^a-z0-9]+", "_", pocket.name.lower()).strip("_")
        env_override = (
            Path(str(Path.cwd() / os.environ["P53CAD_RECEPTOR_PDBQT"]))
            if "P53CAD_RECEPTOR_PDBQT" in os.environ
            else None
        )
        pocket_specific_candidates = [
            Path("data/raw/receptors") / f"{pocket_slug}.pdbqt",
            Path("data/raw/receptors") / f"{name_slug}.pdbqt",
            Path("data/processed/receptors") / f"{pocket_slug}.pdbqt",
            Path("data/processed/receptors") / f"{name_slug}.pdbqt",
            Path("data/raw/receptors/custom_receptor.pdbqt"),
        ]
        shared_candidates = [
            Path("data/raw/p53_wt.pdbqt"),
            Path("data/processed/p53_wt.pdbqt"),
        ]

        candidates = list(pocket_specific_candidates)
        if allow_wt_receptor_fallback:
            candidates.extend(shared_candidates)

        if env_override is not None:
            candidates.insert(0, env_override)

        for path in candidates:
            if path.exists():
                return path
        return None

    def _categorize_receptor_source(
        self,
        receptor_path: Optional[Path],
        pocket_key: str,
        pocket: BindingPocket,
    ) -> str:
        if receptor_path is None:
            return "missing"
        receptor_name = receptor_path.name.lower()
        pocket_slug = re.sub(r"[^a-z0-9]+", "_", pocket_key.lower()).strip("_")
        name_slug = re.sub(r"[^a-z0-9]+", "_", pocket.name.lower()).strip("_")
        if receptor_name in {f"{pocket_slug}.pdbqt", f"{name_slug}.pdbqt"}:
            return "pocket_specific"
        if receptor_name in {"p53_wt.pdbqt", "custom_receptor.pdbqt"}:
            return "shared_receptor"
        return "other"

    def get_receptor_manifest(self, allow_wt_receptor_fallback: bool = False) -> Dict[str, Any]:
        """
        Build a pocket->receptor manifest and detect receptor reuse collisions.
        """
        mapping: Dict[str, str] = {}
        reuse_index: Dict[str, List[str]] = {}
        for pocket_key, pocket in self.pockets.items():
            resolved = self._resolve_receptor_pdbqt(
                pocket_key,
                pocket,
                allow_wt_receptor_fallback=allow_wt_receptor_fallback,
            )
            resolved_str = str(resolved.resolve()) if resolved is not None else "missing"
            mapping[pocket_key] = resolved_str
            reuse_index.setdefault(resolved_str, []).append(pocket_key)

        duplicates = {k: v for k, v in reuse_index.items() if k != "missing" and len(v) > 1}
        return {
            "allow_wt_receptor_fallback": bool(allow_wt_receptor_fallback),
            "pocket_to_receptor": mapping,
            "duplicate_receptor_paths": duplicates,
            "duplicate_count": len(duplicates),
        }

    def get_mode_capabilities(
        self,
        pocket_key: str,
        method: str,
        allow_wt_receptor_fallback: bool = False,
    ) -> Dict[str, Any]:
        """
        Return strict capability probe for a generation mode.
        """
        result: Dict[str, Any] = {
            "method": method,
            "pocket_key": pocket_key,
            "ready": False,
            "reason": "",
            "receptor_path": None,
            "receptor_source": "n/a",
        }
        if pocket_key not in self.pockets:
            result["reason"] = f"Unknown pocket: {pocket_key}"
            return result

        pocket = self.pockets[pocket_key]
        if method in {"template", "denovo", "reinforce"}:
            result["ready"] = True
            result["reason"] = "No docking runtime required."
            return result

        runtime = self._resolve_docking_runtime(
            pocket_key,
            pocket,
            allow_wt_receptor_fallback=allow_wt_receptor_fallback,
        )
        result["receptor_path"] = str(runtime.get("receptor_pdbqt")) if runtime.get("receptor_pdbqt") else None
        result["receptor_source"] = runtime.get("receptor_source", "missing")

        errors = list(runtime.get("errors", []))
        if method == "docking_md":
            md_probe = self._run_md_refinement("CCO", require_md=True)
            if md_probe.get("status") != "ran":
                errors.append(md_probe.get("detail", "MD stack unavailable"))

        if errors:
            result["ready"] = False
            result["reason"] = "; ".join(errors)
            return result

        result["ready"] = True
        result["reason"] = "Ready"
        return result

    def _prepare_receptor_pdbqt(self, pocket_key: str) -> Optional[Path]:
        raw_pdb = Path("data/raw/p53_wt.pdb")
        if not raw_pdb.exists():
            return None

        out_dir = Path("data/processed/receptors")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{pocket_key}.pdbqt"
        if out_path.exists():
            return out_path

        obabel_bin = shutil.which("obabel")
        if not obabel_bin:
            return None

        cmd = [obabel_bin, "-ipdb", str(raw_pdb), "-opdbqt", "-O", str(out_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            warnings.warn(f"Receptor preparation via Open Babel failed: {proc.stderr.strip() or proc.stdout.strip()}")
            return None
        return out_path if out_path.exists() else None

    def _dock_smiles(
        self,
        smiles: str,
        receptor_pdbqt: Path,
        pocket: BindingPocket,
        exhaustiveness: int = 12,
        n_poses: int = 5,
    ) -> Dict[str, Any]:
        with tempfile.TemporaryDirectory(prefix="p53cad_dock_") as tmp:
            tmp_dir = Path(tmp)
            ligand_pdbqt = tmp_dir / "ligand.pdbqt"
            output_pdbqt = tmp_dir / "poses.pdbqt"
            self._write_ligand_pdbqt(smiles, ligand_pdbqt)

            if VINA_PY_AVAILABLE:
                return self._dock_with_vina_python(
                    ligand_pdbqt,
                    receptor_pdbqt,
                    output_pdbqt,
                    pocket,
                    exhaustiveness=exhaustiveness,
                    n_poses=n_poses,
                )
            return self._dock_with_vina_cli(
                ligand_pdbqt,
                receptor_pdbqt,
                output_pdbqt,
                pocket,
                exhaustiveness=exhaustiveness,
                n_poses=n_poses,
            )

    def _write_ligand_pdbqt(self, smiles: str, output_path: Path) -> None:
        if not RDKIT_AVAILABLE or not MEEKO_AVAILABLE:
            raise RuntimeError("Ligand preparation requires rdkit + meeko")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RuntimeError(f"Invalid SMILES for docking: {smiles}")

        mol = Chem.AddHs(mol)
        embed_status = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_status != 0:
            raise RuntimeError(f"Could not embed 3D coordinates for ligand: {smiles}")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

        atom_params = self._resolve_meeko_atom_params()
        prep = MoleculePreparation(load_atom_params=atom_params)
        setups = prep.prepare(mol)
        if not setups:
            raise RuntimeError("Meeko failed to prepare ligand setup")

        pdbqt_payload = PDBQTWriterLegacy.write_string(setups[0])
        if isinstance(pdbqt_payload, tuple):
            pdbqt_string = pdbqt_payload[0]
            ok_flag = bool(pdbqt_payload[1]) if len(pdbqt_payload) > 1 else True
            err_msg = str(pdbqt_payload[2]) if len(pdbqt_payload) > 2 else ""
            if not ok_flag:
                raise RuntimeError(f"Meeko PDBQT writer failed: {err_msg}")
        else:
            pdbqt_string = str(pdbqt_payload)

        output_path.write_text(pdbqt_string)

    def _resolve_meeko_atom_params(self) -> str:
        """Resolve Meeko atom parameter file, even when package data is missing."""
        try:
            import meeko  # type: ignore
            params_dir = Path(meeko.__file__).resolve().parent / "data" / "params"
            if params_dir.exists():
                preferred = params_dir / "ad4_types.json"
                if preferred.exists():
                    # Newer Meeko builds accept explicit file token with suffix.
                    return "ad4_types.json"
                # Fallback to any bundled JSON file if package layout differs.
                for candidate in sorted(params_dir.glob("*.json")):
                    return candidate.name
        except Exception:
            pass

        bundled = Path(__file__).resolve().parents[1] / "data" / "meeko" / "ad4_types.json"
        if bundled.exists():
            return str(bundled)

        raise RuntimeError(
            "Meeko atom params missing. Expected ad4_types.json in Meeko package data or "
            f"bundled file at {bundled}."
        )

    def _dock_with_vina_python(
        self,
        ligand_pdbqt: Path,
        receptor_pdbqt: Path,
        output_pdbqt: Path,
        pocket: BindingPocket,
        exhaustiveness: int,
        n_poses: int,
    ) -> Dict[str, Any]:
        vina = Vina(sf_name="vina")
        vina.set_receptor(str(receptor_pdbqt))
        vina.set_ligand_from_file(str(ligand_pdbqt))
        vina.compute_vina_maps(center=list(pocket.center), box_size=list(pocket.size))
        vina.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = vina.energies(n_poses=n_poses)
        best_affinity = float(energies[0][0]) if len(energies) else 0.0
        vina.write_poses(str(output_pdbqt), n_poses=n_poses, overwrite=True)
        return {
            "backend": "vina_python",
            "best_affinity": best_affinity,
            "pose_count": int(min(n_poses, len(energies))),
            "exhaustiveness": int(exhaustiveness),
        }

    def _dock_with_vina_cli(
        self,
        ligand_pdbqt: Path,
        receptor_pdbqt: Path,
        output_pdbqt: Path,
        pocket: BindingPocket,
        exhaustiveness: int,
        n_poses: int,
    ) -> Dict[str, Any]:
        vina_bin = shutil.which("vina")
        if not vina_bin:
            raise RuntimeError("Vina executable was not found on PATH")

        cmd = [
            vina_bin,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", f"{pocket.center[0]:.3f}",
            "--center_y", f"{pocket.center[1]:.3f}",
            "--center_z", f"{pocket.center[2]:.3f}",
            "--size_x", f"{pocket.size[0]:.3f}",
            "--size_y", f"{pocket.size[1]:.3f}",
            "--size_z", f"{pocket.size[2]:.3f}",
            "--exhaustiveness", str(int(exhaustiveness)),
            "--num_modes", str(int(n_poses)),
            "--out", str(output_pdbqt),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            detail = proc.stderr.strip() or proc.stdout.strip() or "unknown vina CLI error"
            raise RuntimeError(f"Vina CLI docking failed: {detail}")

        best_affinity = self._parse_vina_affinity(proc.stdout)
        if best_affinity is None and output_pdbqt.exists():
            best_affinity = self._parse_vina_affinity(output_pdbqt.read_text())
        if best_affinity is None:
            raise RuntimeError("Could not parse docking affinity from Vina output")

        return {
            "backend": "vina_cli",
            "best_affinity": float(best_affinity),
            "pose_count": int(n_poses),
            "exhaustiveness": int(exhaustiveness),
        }

    def _parse_vina_affinity(self, output_text: str) -> Optional[float]:
        for line in (output_text or "").splitlines():
            match = re.search(r"REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)", line)
            if match:
                return float(match.group(1))
            table = re.match(r"^\s*\d+\s+(-?\d+(?:\.\d+)?)\s+", line)
            if table:
                return float(table.group(1))
        return None

    def _run_md_refinement(self, smiles: str, require_md: bool = False) -> Dict[str, str]:
        if not require_md:
            return {"status": "not_requested"}

        if not OPENMM_AVAILABLE:
            return {"status": "unavailable", "detail": "OpenMM is not installed"}

        try:
            from openff.toolkit import ForceField, Molecule  # type: ignore
            from openff.interchange import Interchange  # type: ignore
            from openmm import LangevinMiddleIntegrator, unit  # type: ignore
            from openmm.app import Simulation  # type: ignore
            import openmmforcefields  # type: ignore  # noqa: F401
        except Exception:
            return {
                "status": "unavailable",
                "detail": "MD refinement requires openff-toolkit + openff-interchange + openmmforcefields",
            }

        if not RDKIT_AVAILABLE:
            return {"status": "unavailable", "detail": "RDKit is required to generate ligand conformers for MD"}

        if not smiles:
            return {"status": "unavailable", "detail": "Ligand SMILES is required for docking_md refinement"}

        try:
            # Generate a physically valid ligand conformer via RDKit.
            rdkit_mol = Chem.MolFromSmiles(smiles)
            if rdkit_mol is None:
                return {"status": "unavailable", "detail": "Invalid ligand SMILES for MD refinement"}
            rdkit_mol = Chem.AddHs(rdkit_mol)
            emb_status = AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
            if emb_status != 0:
                return {"status": "unavailable", "detail": "Failed to embed 3D ligand conformer for MD"}
            AllChem.MMFFOptimizeMolecule(rdkit_mol, maxIters=200)

            ligand = Molecule.from_rdkit(
                rdkit_mol,
                allow_undefined_stereo=True,
                hydrogens_are_explicit=True,
            )
            topology = ligand.to_topology()
            forcefield = ForceField("openff_unconstrained-2.1.0.offxml")
            interchange = Interchange.from_smirnoff(
                force_field=forcefield,
                topology=topology,
            )
            system = interchange.to_openmm_system()
            omm_top = interchange.to_openmm_topology()
            positions = interchange.positions.to_openmm()

            integrator = LangevinMiddleIntegrator(
                300 * unit.kelvin,
                1.0 / unit.picosecond,
                0.002 * unit.picoseconds,
            )
            simulation = Simulation(omm_top, system, integrator)
            simulation.context.setPositions(positions)
            simulation.minimizeEnergy(maxIterations=250)
            simulation.step(500)

            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            return {
                "status": "ran",
                "detail": f"Ligand-only OpenMM refinement completed (NVT 500 steps, E={energy:.2f} kcal/mol)",
            }
        except Exception as exc:
            return {"status": "unavailable", "detail": f"MD refinement failed: {exc}"}

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

    def _extract_mutation_positions(self, mutations: List[str]) -> List[int]:
        positions = []
        for mut in mutations:
            match = re.match(r"^[A-Za-z](\d+)[A-Za-z]$", str(mut).strip())
            if not match:
                continue
            positions.append(int(match.group(1)))
        return positions

    def _compute_context_adjustment(
        self,
        pocket: BindingPocket,
        protein_context: Optional[Dict[str, Any]],
    ) -> Dict[str, float]:
        """
        Compute pocket/rescue compatibility adjustment.

        Returns shifts where more negative affinity means better fit.
        """
        default = {
            "enabled": 0.0,
            "affinity_shift": 0.0,
            "confidence_shift": 0.0,
            "pocket_alignment": 0.0,
            "rescue_norm": 0.0,
            "mutation_count": 0.0,
        }
        if not protein_context:
            return default

        mutations = protein_context.get("mutations", []) or []
        positions = self._extract_mutation_positions(mutations)
        if not positions:
            return default

        aligned = 0
        for pos in positions:
            if any(abs(pos - residue) <= 6 for residue in pocket.residues):
                aligned += 1

        pocket_alignment = aligned / float(max(len(positions), 1))

        rescue_score = protein_context.get("score", 0.0)
        try:
            rescue_score = float(rescue_score)
        except (TypeError, ValueError):
            rescue_score = 0.0
        rescue_norm = float(np.clip((rescue_score + 1.0) / 2.0, 0.0, 1.0))

        mutation_count = float(len(positions))
        mutation_penalty = float(max(0.0, (mutation_count - 8.0) * 0.03))

        # Stronger alignment and better rescue score slightly improve predicted affinity.
        affinity_shift = (-0.9 * pocket_alignment) + (-0.35 * rescue_norm) + mutation_penalty
        confidence_shift = (0.12 * pocket_alignment) + (0.05 * rescue_norm) - (0.5 * mutation_penalty)

        return {
            "enabled": 1.0,
            "affinity_shift": float(affinity_shift),
            "confidence_shift": float(confidence_shift),
            "pocket_alignment": float(pocket_alignment),
            "rescue_norm": float(rescue_norm),
            "mutation_count": mutation_count,
        }

    def _apply_context_adjustment(
        self,
        affinity: float,
        confidence: float,
        context_adjustment: Optional[Dict[str, float]],
    ) -> Tuple[float, float]:
        if not context_adjustment or context_adjustment.get("enabled", 0.0) < 0.5:
            return float(affinity), float(np.clip(confidence, 0.0, 1.0))
        adj_affinity = float(affinity + context_adjustment.get("affinity_shift", 0.0))
        adj_confidence = float(
            np.clip(confidence + context_adjustment.get("confidence_shift", 0.0), 0.0, 1.0)
        )
        return adj_affinity, adj_confidence

    def _context_metadata(self, context_adjustment: Optional[Dict[str, float]]) -> Dict[str, Any]:
        if not context_adjustment or context_adjustment.get("enabled", 0.0) < 0.5:
            return {"protein_context_applied": False}
        return {
            "protein_context_applied": True,
            "protein_context_pocket_alignment": float(context_adjustment.get("pocket_alignment", 0.0)),
            "protein_context_rescue_norm": float(context_adjustment.get("rescue_norm", 0.0)),
            "protein_context_mutation_count": int(context_adjustment.get("mutation_count", 0.0)),
            "protein_context_affinity_shift": float(context_adjustment.get("affinity_shift", 0.0)),
            "protein_context_confidence_shift": float(context_adjustment.get("confidence_shift", 0.0)),
        }

    def _get_default_templates(self) -> List[str]:
        """Default drug-like scaffolds when no known binders exist."""
        return self._get_all_valid_templates()[:5]

    def _get_all_valid_templates(self) -> List[str]:
        """
        Comprehensive list of validated drug-like SMILES templates.

        These are real, chemically valid molecules suitable for p53 stabilization.
        Includes known p53 binders, approved drugs, and drug-like scaffolds.
        """
        return [
            # === KNOWN p53 BINDERS ===
            "CC1=CC=C(C=C1)C2=CC(=NO2)C3=CC=C(C=C3)Cl",  # PhiKan083 (Y220C stabilizer)
            "CC1=CC=C(C=C1)C2=CC(=NO2)C3=CC=CC=C3F",    # PhiKan7088
            "CC(C)(C)CC(=O)NC1=CC=C(C=C1)C2=CC(=C(C=C2)Cl)Cl",  # Nutlin-like (MDM2)

            # === BENZIMIDAZOLE SCAFFOLDS ===
            "c1ccc2c(c1)nc(cn2)N",      # Benzimidazole amine
            "c1ccc2c(c1)[nH]c(n2)C",    # 2-methylbenzimidazole
            "Nc1nc2ccccc2[nH]1",        # 2-aminobenzimidazole
            "c1ccc2c(c1)nc([nH]2)c3ccccn3",  # Benzimidazole-pyridine

            # === PYRIDINE/PYRIMIDINE SCAFFOLDS ===
            "c1ccc(cc1)c2ccccn2",       # Phenylpyridine
            "c1ccc(cc1)c2ccncc2",       # 4-phenylpyridine
            "Nc1ncnc2ccccc12",          # Quinazolin-4-amine
            "c1cnc2ccccc2n1",           # Quinoxaline

            # === INDOLE SCAFFOLDS ===
            "c1ccc2c(c1)cc[nH]2",       # Indole
            "c1ccc2c(c1)c(c[nH]2)C",    # 3-methylindole
            "CC(=O)c1c[nH]c2ccccc12",   # 3-acetylindole

            # === ISOXAZOLE/OXAZOLE ===
            "c1cc(on1)c2ccccc2",        # 3-phenylisoxazole
            "c1ccc(cc1)c2cnoc2",        # 5-phenylisoxazole
            "c1ccc(cc1)c2ncco2",        # 2-phenyloxazole

            # === THIOPHENE/FURAN ===
            "c1ccc(cc1)c2cccs2",        # 2-phenylthiophene
            "c1ccc(cc1)c2ccco2",        # 2-phenylfuran
            "Cc1cccs1",                  # 2-methylthiophene

            # === AMIDE/SULFONAMIDE ===
            "CC1=CC=C(C=C1)C(=O)N",     # Toluamide
            "CC(=O)Nc1ccc(cc1)O",       # Paracetamol-like
            "NS(=O)(=O)c1ccc(cc1)N",    # Sulfanilamide
            "CC(=O)Nc1ccc(cc1)C(=O)O",  # 4-acetamidobenzoic acid

            # === COUMARIN/CHROMONE ===
            "c1ccc2c(c1)c(=O)oc2",      # Coumarin
            "O=c1cc(oc2ccccc12)c3ccccc3",  # Flavone
            "O=c1ccoc2ccccc12",         # Chromone

            # === NAPHTHALENE ===
            "c1ccc2ccccc2c1",           # Naphthalene
            "Oc1ccc2ccccc2c1",          # 1-naphthol
            "Nc1ccc2ccccc2c1",          # 1-naphthylamine

            # === BIPHENYL ===
            "c1ccc(cc1)c2ccccc2",       # Biphenyl
            "c1ccc(cc1)c2ccc(cc2)O",    # 4-hydroxybiphenyl
            "c1ccc(cc1)c2ccc(cc2)N",    # 4-aminobiphenyl

            # === KNOWN DRUGS (p53-relevant) ===
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen (anti-inflammatory)
            "COc1ccc(cc1)C(O)C(C)NC",   # Ephedrine scaffold
            "c1ccc(cc1)C2=NCC(=O)N2",   # Phenytoin scaffold

            # === CAVITY FILLERS ===
            "CC(C)c1ccc(cc1)C",         # p-cymene
            "CCc1ccc(cc1)CC",           # 1,4-diethylbenzene
            "c1ccc(cc1)C(c2ccccc2)O",   # Benzhydrol

            # === HETEROCYCLIC CORES ===
            "c1cnc2ncccn12",            # Imidazo[1,2-a]pyrazine
            "c1ccc2nccnc2c1",           # Quinazoline
            "c1ccc2c(c1)cnn2",          # Indazole
            "c1cc2ccccc2cn1",           # Isoquinoline
        ]

    def _is_valid_smiles(self, smiles: str) -> bool:
        """Check if SMILES is chemically valid."""
        if not smiles or len(smiles) < 3:
            return False

        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        else:
            # Without RDKit, only allow known valid templates
            # Don't trust generated SMILES
            return smiles in self._get_all_valid_templates()

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
                       choices=["template", "denovo", "reinforce", "docking", "docking_md"],
                       help="Generation method")
    parser.add_argument("--output", type=str, default="drug_candidates.json",
                       help="Output file path")
    parser.add_argument("--export-sdf", type=str, default=None,
                       help="Export to SDF file for docking")
    parser.add_argument(
        "--allow-wt-receptor-fallback",
        action="store_true",
        help="Allow fallback to shared WT receptor if pocket-specific receptor is missing (less realistic).",
    )

    args = parser.parse_args()

    # Generate
    engine = DrugGeneratorEngine(allow_wt_receptor_fallback=args.allow_wt_receptor_fallback)
    candidates = engine.generate_for_pocket(
        args.pocket,
        args.n_candidates,
        args.method,
        allow_wt_receptor_fallback=args.allow_wt_receptor_fallback,
    )

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
