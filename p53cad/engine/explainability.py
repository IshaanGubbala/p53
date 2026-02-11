"""
Explainable AI Module for p53-proteoMgCAD

Provides mechanistic insights into WHY rescue mutations work through:
1. Attention Attribution Maps - Which residues ESM-2 focuses on
2. Counterfactual Analysis - What if we changed position X instead?
3. Energy Decomposition - Break down physics contributions
4. Rescue Mechanism Classification - Structural vs electrostatic vs allosteric
"""

import numpy as np
import os
import math
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
import warnings

try:
    import torch
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


class ExplainabilityDependencyError(RuntimeError):
    """Raised when explainability is requested without required real dependencies."""


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class ResidueAttribution:
    """Attribution score for a single residue."""
    position: int
    residue: str
    attention_score: float  # 0-1, how much model attends to this position
    gradient_score: float   # magnitude of gradient w.r.t. this position
    occlusion_score: float  # impact when masked
    integrated_gradient: float  # integrated gradients attribution
    composite_importance: float  # weighted combination

    # Interpretation
    functional_role: str  # "DNA_binding", "structural", "zinc", "unknown"
    is_hotspot: bool
    distance_to_mutation: float  # Angstroms to nearest cancer mutation


@dataclass
class CounterfactualResult:
    """Result of counterfactual analysis."""
    original_mutation: str
    alternative_mutation: str
    original_score: float
    alternative_score: float
    delta_score: float
    explanation: str
    confidence: float
    structural_impact: Dict[str, float]


@dataclass
class EnergyDecomposition:
    """Breakdown of energy contributions to rescue."""
    total_ddg: float  # Total ΔΔG

    # Component contributions
    electrostatic: float
    van_der_waals: float
    hydrogen_bonds: float
    solvation: float
    entropy: float

    # Structural contributions
    backbone_strain: float
    sidechain_packing: float
    cavity_formation: float

    # Functional contributions
    dna_binding_change: float
    zinc_coordination_change: float

    # Interpretation
    dominant_mechanism: str
    rescue_quality: str  # "excellent", "good", "marginal", "poor"

    def to_dict(self) -> Dict:
        return {
            'total_ddg': self.total_ddg,
            'electrostatic': self.electrostatic,
            'van_der_waals': self.van_der_waals,
            'hydrogen_bonds': self.hydrogen_bonds,
            'solvation': self.solvation,
            'entropy': self.entropy,
            'backbone_strain': self.backbone_strain,
            'sidechain_packing': self.sidechain_packing,
            'cavity_formation': self.cavity_formation,
            'dna_binding_change': self.dna_binding_change,
            'zinc_coordination_change': self.zinc_coordination_change,
            'dominant_mechanism': self.dominant_mechanism,
            'rescue_quality': self.rescue_quality
        }


@dataclass
class RescueMechanism:
    """Classification of rescue mechanism."""
    primary_mechanism: str  # "stability_restoration", "binding_enhancement", "allosteric", "compensatory"
    secondary_mechanism: Optional[str]
    confidence: float
    evidence: List[str]

    # Detailed breakdown
    stability_contribution: float  # -1 to 1
    binding_contribution: float
    structural_contribution: float
    electrostatic_contribution: float

    # Recommendations
    synergistic_positions: List[int]  # Other positions that might help
    risk_factors: List[str]


# ============================================================================
# FUNCTIONAL SITE DEFINITIONS
# ============================================================================

P53_FUNCTIONAL_SITES = {
    # DNA binding residues (direct contact)
    "DNA_binding_direct": [119, 120, 241, 243, 244, 245, 248, 273, 276, 277, 280, 281],

    # DNA binding residues (indirect/structural support)
    "DNA_binding_support": [174, 175, 176, 179, 239, 240, 242, 282, 283, 284],

    # Zinc coordination site
    "zinc_coordination": [176, 179, 238, 242],

    # L1 loop (DNA minor groove contact)
    "L1_loop": list(range(113, 125)),

    # L2 loop (structural/DNA contact)
    "L2_loop": list(range(163, 196)),

    # L3 loop (major groove contact)
    "L3_loop": list(range(236, 252)),

    # Hydrophobic core
    "hydrophobic_core": [132, 134, 138, 145, 146, 157, 158, 159, 173, 194, 195, 196, 214, 220, 229, 253, 254, 255, 256, 257],

    # Cancer hotspot positions
    "hotspots": [175, 220, 245, 248, 249, 273, 282],

    # Tetramerization domain contact
    "tetramer_interface": list(range(324, 356)),

    # MDM2 binding motif
    "mdm2_binding": list(range(17, 29)),
}


# ============================================================================
# ATTENTION ATTRIBUTION
# ============================================================================

class AttentionExtractor:
    """Extract and analyze attention patterns from ESM-2."""

    def __init__(self, model=None, tokenizer=None):
        self.model = model
        self.tokenizer = tokenizer
        self._attention_cache = {}

    def extract_attention(self, sequence: str, layer: int = -1) -> np.ndarray:
        """
        Extract attention weights from ESM-2 for a sequence.

        Args:
            sequence: Amino acid sequence
            layer: Which layer to extract (-1 = last layer)

        Returns:
            Attention matrix [n_heads, seq_len, seq_len]
        """
        if not TORCH_AVAILABLE:
            raise ExplainabilityDependencyError(
                "PyTorch is required for attention extraction."
            )
        if self.model is None or self.tokenizer is None:
            raise ExplainabilityDependencyError(
                "Attention extraction requires a loaded model and tokenizer."
            )

        cache_key = f"{sequence}_{layer}"
        if cache_key in self._attention_cache:
            return self._attention_cache[cache_key]

        # Tokenize
        inputs = self.tokenizer(sequence, return_tensors="pt")
        try:
            model_device = next(self.model.parameters()).device
        except StopIteration:
            model_device = torch.device("cpu")
        inputs = {k: v.to(model_device) for k, v in inputs.items()}

        # Forward with attention output
        try:
            with torch.no_grad():
                outputs = self.model(**inputs, output_attentions=True, return_dict=True)
        except RuntimeError as err:
            # MPS can intermittently fail on attention ops for this workload.
            # Retry once on CPU with the same real model and inputs.
            err_text = str(err)
            if model_device.type == "mps" and "Placeholder storage has not been allocated on MPS device" in err_text:
                warnings.warn(
                    "MPS attention execution failed; retrying explainability attention on CPU.",
                    RuntimeWarning,
                )
                self.model = self.model.to("cpu")
                inputs = {k: v.to("cpu") for k, v in inputs.items()}
                with torch.no_grad():
                    outputs = self.model(**inputs, output_attentions=True, return_dict=True)
            else:
                raise

        # Get attention from specified layer
        attentions = outputs.attentions
        if attentions is None or len(attentions) == 0:
            # Some HF/torch attention backends silently drop attention tensors.
            # Retry on CPU with explicit attention flags before failing closed.
            try:
                if hasattr(self.model, "config"):
                    self.model.config.output_attentions = True
                    self.model.config.return_dict = True
                esm_module = getattr(self.model, "esm", None)
                if esm_module is not None and hasattr(esm_module, "config"):
                    esm_module.config.output_attentions = True
                    esm_module.config.return_dict = True

                self.model = self.model.to("cpu")
                cpu_inputs = {k: v.to("cpu") for k, v in inputs.items()}
                with torch.no_grad():
                    retry_outputs = self.model(
                        **cpu_inputs,
                        output_attentions=True,
                        return_dict=True,
                    )
                attentions = retry_outputs.attentions
            except Exception as retry_err:
                raise ExplainabilityDependencyError(
                    "Model did not return attention tensors for explainability."
                ) from retry_err

        if attentions is None or len(attentions) == 0:
            raise ExplainabilityDependencyError(
                "Model did not return attention tensors for explainability."
            )
        layer_attention = attentions[layer].squeeze(0).detach().cpu().numpy()  # [n_heads, seq_len, seq_len]

        self._attention_cache[cache_key] = layer_attention
        return layer_attention

    def get_residue_attention_scores(self, sequence: str,
                                     focus_positions: Optional[List[int]] = None) -> Dict[int, float]:
        """
        Get attention scores for each residue position.

        Args:
            sequence: Amino acid sequence
            focus_positions: If provided, compute attention TO these positions

        Returns:
            Dict mapping position -> attention score
        """
        attention = self.extract_attention(sequence)

        # Average over heads
        avg_attention = attention.mean(axis=0)  # [seq_len, seq_len]

        scores = {}

        if focus_positions:
            # Attention flowing TO focus positions
            for pos in range(len(sequence)):
                score = sum(avg_attention[pos, fp] for fp in focus_positions if fp < len(sequence))
                scores[pos] = float(score / len(focus_positions))
        else:
            # Global attention received by each position
            for pos in range(len(sequence)):
                scores[pos] = float(avg_attention[:, pos].sum())

        # Normalize
        max_score = max(scores.values()) if scores else 1.0
        if max_score > 0:
            scores = {k: v / max_score for k, v in scores.items()}

        return scores

    def find_attention_paths(self, sequence: str, start_pos: int, end_pos: int,
                            threshold: float = 0.1) -> List[List[int]]:
        """
        Find attention paths between two positions (information flow).

        Returns list of paths (each path is list of positions).
        """
        attention = self.extract_attention(sequence)
        avg_attention = attention.mean(axis=0)

        # BFS to find paths above threshold
        paths = []
        queue = [[start_pos]]
        visited = {start_pos}

        while queue and len(paths) < 10:
            path = queue.pop(0)
            current = path[-1]

            if current == end_pos:
                paths.append(path)
                continue

            # Find high-attention neighbors
            for next_pos in range(len(sequence)):
                if next_pos not in visited and avg_attention[current, next_pos] > threshold:
                    new_path = path + [next_pos]
                    if len(new_path) <= 5:  # Max path length
                        queue.append(new_path)
                        visited.add(next_pos)

        return paths

    def _synthetic_attention(self, seq_len: int, sequence: str = None) -> np.ndarray:
        raise ExplainabilityDependencyError(
            "Synthetic attention fallback is disabled in strict mode."
        )


# ============================================================================
# GRADIENT-BASED ATTRIBUTION
# ============================================================================

class GradientAttributor:
    """Gradient-based attribution methods."""

    def __init__(self, model=None, oracle=None):
        self.model = model
        self.oracle = oracle

    def compute_gradients(self, embedding: np.ndarray,
                         target_fn: callable) -> np.ndarray:
        """
        Compute gradients of target function w.r.t. embedding.

        Args:
            embedding: [1, L, D] embedding tensor
            target_fn: Function that takes embedding, returns scalar

        Returns:
            Gradient array [L, D]
        """
        if not TORCH_AVAILABLE:
            raise ExplainabilityDependencyError(
                "Gradient attribution requires PyTorch."
            )

        emb_tensor = torch.tensor(embedding, requires_grad=True, dtype=torch.float32)

        # Forward
        score = target_fn(emb_tensor)

        # Backward
        score.backward()

        gradients = emb_tensor.grad.numpy().squeeze(0)  # [L, D]
        return gradients

    def _approximate_gradients(self, embedding: np.ndarray) -> np.ndarray:
        """
        Approximate gradients using embedding variance when autograd unavailable.

        Higher variance positions are likely more important for the model.
        """
        emb = embedding.squeeze(0)  # [L, D]
        L, D = emb.shape

        # Use embedding variance as importance proxy
        # Positions with distinct embeddings are likely more important
        mean_emb = emb.mean(axis=0, keepdims=True)
        deviation = emb - mean_emb

        # Magnitude of deviation indicates importance
        importance = np.linalg.norm(deviation, axis=1)  # [L]

        # Normalize and expand to embedding dimension
        importance = importance / (importance.max() + 1e-10)
        gradients = deviation * importance[:, np.newaxis]

        return gradients

    def compute_integrated_gradients(self, embedding: np.ndarray,
                                     baseline: np.ndarray,
                                     target_fn: callable,
                                     n_steps: int = 50) -> np.ndarray:
        """
        Compute integrated gradients attribution.

        IG = (x - baseline) * integral(gradient(baseline + alpha*(x-baseline)), alpha=0..1)
        """
        if not TORCH_AVAILABLE:
            raise ExplainabilityDependencyError(
                "Integrated gradients require PyTorch."
            )

        # Interpolate between baseline and input
        scaled_inputs = []
        for alpha in np.linspace(0, 1, n_steps):
            scaled = baseline + alpha * (embedding - baseline)
            scaled_inputs.append(scaled)
        # Compute gradients at each interpolation point
        gradients = []
        for scaled_input in scaled_inputs:
            grad = self.compute_gradients(scaled_input, target_fn)
            gradients.append(grad)

        # Average gradients
        avg_grad = np.mean(gradients, axis=0)

        # Scale by input - baseline
        ig = (embedding.squeeze(0) - baseline.squeeze(0)) * avg_grad

        return ig

    def get_position_importance(self, gradients: np.ndarray) -> Dict[int, float]:
        """
        Convert per-position gradient tensor to importance scores.

        Args:
            gradients: [L, D] gradient array

        Returns:
            Dict mapping position -> importance score
        """
        # L2 norm of gradients at each position
        importance = np.linalg.norm(gradients, axis=-1)

        # Normalize to 0-1
        max_imp = importance.max() if importance.max() > 0 else 1.0
        importance = importance / max_imp

        return {i: float(imp) for i, imp in enumerate(importance)}


# ============================================================================
# OCCLUSION-BASED ATTRIBUTION
# ============================================================================

class OcclusionAttributor:
    """Occlusion-based attribution (mask and measure impact)."""

    def __init__(self, oracle=None, embedder=None):
        self.oracle = oracle
        self.embedder = embedder

    def compute_occlusion_importance(self, sequence: str,
                                     base_score: float = None) -> Dict[int, float]:
        """
        Compute importance by masking each position and measuring score change.

        Returns:
            Dict mapping position -> importance (score drop when masked)
        """
        if base_score is None:
            base_score = self._get_score(sequence)

        importance = {}

        for pos in range(len(sequence)):
            # Mask position with 'X' (unknown)
            masked_seq = sequence[:pos] + 'X' + sequence[pos+1:]
            masked_score = self._get_score(masked_seq)

            # Importance = how much score drops
            importance[pos] = base_score - masked_score

        # Normalize
        max_imp = max(abs(v) for v in importance.values()) if importance else 1.0
        if max_imp > 0:
            importance = {k: v / max_imp for k, v in importance.items()}

        return importance

    def compute_ablation_importance(self, sequence: str,
                                   positions: List[int]) -> Dict[str, float]:
        """
        Measure importance of groups of positions by ablating them together.

        Returns:
            Dict mapping group_name -> importance
        """
        base_score = self._get_score(sequence)

        group_importance = {}

        for group_name, group_positions in P53_FUNCTIONAL_SITES.items():
            # Mask all positions in group
            masked_seq = list(sequence)
            for pos in group_positions:
                if pos < len(masked_seq):
                    masked_seq[pos] = 'X'
            masked_seq = ''.join(masked_seq)

            masked_score = self._get_score(masked_seq)
            group_importance[group_name] = base_score - masked_score

        return group_importance

    def _get_score(self, sequence: str) -> float:
        """Get functional score for sequence."""
        if self.oracle is None or self.embedder is None:
            raise ExplainabilityDependencyError(
                "Occlusion scoring requires both oracle and embedder."
            )

        # Use actual oracle
        embedding = self.embedder.encode(sequence)
        return float(self.oracle.predict(embedding.mean(axis=1)))

# ============================================================================
# COUNTERFACTUAL ANALYSIS
# ============================================================================

class CounterfactualAnalyzer:
    """Analyze what would happen with alternative mutations."""

    def __init__(self, oracle=None, embedder=None):
        self.oracle = oracle
        self.embedder = embedder
        self.amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

    def analyze_counterfactual(self, wt_sequence: str,
                              rescue_mutation: str,
                              n_alternatives: int = 5) -> List[CounterfactualResult]:
        """
        Analyze alternative mutations at the rescue position.

        Args:
            wt_sequence: Wild-type sequence
            rescue_mutation: Rescue mutation in format "A123B"
            n_alternatives: Number of alternatives to analyze

        Returns:
            List of CounterfactualResults
        """
        # Parse mutation
        original_aa = rescue_mutation[0]
        position = int(rescue_mutation[1:-1]) - 1  # 0-indexed
        rescue_aa = rescue_mutation[-1]

        # Get base score (with rescue)
        rescue_seq = wt_sequence[:position] + rescue_aa + wt_sequence[position+1:]
        rescue_score = self._get_score(rescue_seq)

        results = []

        # Test all other amino acids at this position
        for alt_aa in self.amino_acids:
            if alt_aa == rescue_aa:
                continue

            alt_seq = wt_sequence[:position] + alt_aa + wt_sequence[position+1:]
            alt_score = self._get_score(alt_seq)

            delta = alt_score - rescue_score

            # Generate explanation
            explanation = self._generate_explanation(
                rescue_aa, alt_aa, position, delta, wt_sequence
            )

            result = CounterfactualResult(
                original_mutation=rescue_mutation,
                alternative_mutation=f"{original_aa}{position+1}{alt_aa}",
                original_score=rescue_score,
                alternative_score=alt_score,
                delta_score=delta,
                explanation=explanation,
                confidence=abs(delta) / max(abs(rescue_score), 0.1),
                structural_impact=self._estimate_structural_impact(rescue_aa, alt_aa, position)
            )
            results.append(result)

        # Sort by score (best alternatives first)
        results.sort(key=lambda x: x.alternative_score, reverse=True)

        return results[:n_alternatives]

    def analyze_position_importance(self, wt_sequence: str,
                                   mutations: List[str]) -> Dict[int, float]:
        """
        Analyze which positions in a multi-mutation rescue are most important.

        Returns:
            Dict mapping position -> importance
        """
        # Apply all mutations
        full_rescue_seq = wt_sequence
        parsed_mutations = []

        for mut in mutations:
            original_aa = mut[0]
            position = int(mut[1:-1]) - 1
            new_aa = mut[-1]
            full_rescue_seq = full_rescue_seq[:position] + new_aa + full_rescue_seq[position+1:]
            parsed_mutations.append((position, original_aa, new_aa))

        full_score = self._get_score(full_rescue_seq)

        # Test removing each mutation
        position_importance = {}

        for pos, orig_aa, new_aa in parsed_mutations:
            # Remove this mutation (revert to original)
            partial_seq = full_rescue_seq[:pos] + orig_aa + full_rescue_seq[pos+1:]
            partial_score = self._get_score(partial_seq)

            # Importance = score drop when removed
            position_importance[pos] = full_score - partial_score

        return position_importance

    def _get_score(self, sequence: str) -> float:
        """Get functional score."""
        if self.oracle is None or self.embedder is None:
            raise ExplainabilityDependencyError(
                "Counterfactual scoring requires both oracle and embedder."
            )

        embedding = self.embedder.encode(sequence)
        return float(self.oracle.predict(embedding.mean(axis=1)))

    def _generate_explanation(self, original_aa: str, alternative_aa: str,
                             position: int, delta: float,
                             sequence: str) -> str:
        """Generate human-readable explanation for counterfactual."""

        # Classify amino acid types
        aa_properties = {
            'A': 'small hydrophobic', 'C': 'cysteine', 'D': 'acidic',
            'E': 'acidic', 'F': 'aromatic', 'G': 'glycine',
            'H': 'basic aromatic', 'I': 'hydrophobic', 'K': 'basic',
            'L': 'hydrophobic', 'M': 'hydrophobic', 'N': 'polar',
            'P': 'proline', 'Q': 'polar', 'R': 'basic',
            'S': 'small polar', 'T': 'polar', 'V': 'hydrophobic',
            'W': 'aromatic', 'Y': 'aromatic polar'
        }

        orig_prop = aa_properties.get(original_aa, 'unknown')
        alt_prop = aa_properties.get(alternative_aa, 'unknown')

        # Determine functional context
        position_1indexed = position + 1
        context = "unknown region"
        for site_name, positions in P53_FUNCTIONAL_SITES.items():
            if position in positions:
                context = site_name.replace('_', ' ')
                break

        # Build explanation
        if delta > 0.1:
            quality = "better"
        elif delta < -0.1:
            quality = "worse"
        else:
            quality = "similar"

        explanation = (
            f"Substituting {original_aa} ({orig_prop}) with {alternative_aa} ({alt_prop}) "
            f"at position {position_1indexed} ({context}) would result in {quality} function "
            f"(delta = {delta:+.3f})."
        )

        # Add mechanistic insight
        if 'hydrophobic' in orig_prop and 'polar' in alt_prop:
            explanation += " This would disrupt hydrophobic packing."
        elif 'basic' in orig_prop and 'acidic' in alt_prop:
            explanation += " This charge reversal may affect electrostatic interactions."
        elif original_aa == 'P' or alternative_aa == 'P':
            explanation += " Proline substitutions can affect backbone flexibility."
        elif original_aa == 'G' or alternative_aa == 'G':
            explanation += " Glycine substitutions can affect conformational freedom."

        return explanation

    def _estimate_structural_impact(self, original_aa: str,
                                   alternative_aa: str,
                                   position: int) -> Dict[str, float]:
        """Estimate structural impact of mutation."""

        # Amino acid properties
        size = {'G': 0, 'A': 1, 'S': 2, 'C': 2, 'P': 3, 'V': 3, 'T': 3,
                'I': 4, 'L': 4, 'N': 4, 'D': 4, 'M': 5, 'K': 5, 'E': 5,
                'Q': 5, 'H': 6, 'F': 7, 'R': 7, 'Y': 8, 'W': 9}

        hydrophobicity = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
                         'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
                         'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
                         'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}

        orig_size = size.get(original_aa, 5)
        alt_size = size.get(alternative_aa, 5)

        orig_hydro = hydrophobicity.get(original_aa, 0)
        alt_hydro = hydrophobicity.get(alternative_aa, 0)

        return {
            'size_change': alt_size - orig_size,
            'hydrophobicity_change': alt_hydro - orig_hydro,
            'steric_clash_risk': abs(alt_size - orig_size) / 9,
            'polarity_change': 1.0 if (orig_hydro > 0) != (alt_hydro > 0) else 0.0
        }


# ============================================================================
# ENERGY DECOMPOSITION
# ============================================================================

class EnergyDecomposer:
    """Decompose rescue energy into physical components."""

    def __init__(self):
        # Reference energy values (from literature, kcal/mol)
        self.reference_energies = {
            'hydrogen_bond': -2.0,
            'salt_bridge': -3.0,
            'hydrophobic_contact': -0.5,
            'aromatic_stacking': -1.0,
            'cavity_penalty': 2.0,
            'clash_penalty': 5.0,
        }

    def decompose_rescue(self, wt_sequence: str,
                        mutant_sequence: str,
                        rescue_sequence: str,
                        mutations: List[str]) -> EnergyDecomposition:
        """
        Decompose the energy contribution of a rescue.

        Args:
            wt_sequence: Wild-type p53
            mutant_sequence: Cancer mutant p53
            rescue_sequence: Rescued p53
            mutations: List of rescue mutations

        Returns:
            EnergyDecomposition object
        """
        # Parse mutations and analyze changes
        parsed = []
        for mut in mutations:
            orig_aa = mut[0]
            pos = int(mut[1:-1]) - 1
            new_aa = mut[-1]
            parsed.append((pos, orig_aa, new_aa))

        # Estimate energy components
        electrostatic = self._estimate_electrostatic(parsed)
        vdw = self._estimate_van_der_waals(parsed)
        hbonds = self._estimate_hydrogen_bonds(parsed)
        solvation = self._estimate_solvation(parsed)
        entropy = self._estimate_entropy(parsed)

        backbone_strain = self._estimate_backbone_strain(parsed)
        sidechain_packing = self._estimate_sidechain_packing(parsed)
        cavity = self._estimate_cavity_formation(parsed)

        dna_binding = self._estimate_dna_binding_change(parsed)
        zinc = self._estimate_zinc_coordination(parsed)

        # Apply rescue bonus: mutations positioned strategically should help
        rescue_bonus = 0.0
        for pos, orig, new in parsed:
            # Mutations in L2 loop or near zinc often compensate for structural mutants
            if pos in P53_FUNCTIONAL_SITES.get("L2_loop", []):
                rescue_bonus -= 0.8
            # DNA binding support mutations are known to help rescue function
            if pos in P53_FUNCTIONAL_SITES.get("DNA_binding_support", []):
                rescue_bonus -= 1.0
            # L3 loop modifications can restore DNA binding
            if pos in P53_FUNCTIONAL_SITES.get("L3_loop", []):
                rescue_bonus -= 0.7
            # Mutations that introduce favorable interactions
            if new in ['Y', 'F', 'W'] and orig not in ['Y', 'F', 'W']:
                rescue_bonus -= 0.6  # Aromatic packing improvement
            if new in ['R', 'K'] and pos not in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                rescue_bonus -= 0.5  # Surface positive charge often helps DNA binding
            # Known effective rescue positions (from literature)
            if pos in [239, 268, 284, 245]:  # Well-studied rescue positions
                rescue_bonus -= 0.8

        # Total with rescue bonus
        total_ddg = (electrostatic + vdw + hbonds + solvation + entropy +
                    backbone_strain + sidechain_packing + cavity +
                    dna_binding + zinc + rescue_bonus)

        # Determine dominant mechanism
        components = {
            'electrostatic': abs(electrostatic),
            'van_der_waals': abs(vdw),
            'hydrogen_bonds': abs(hbonds),
            'solvation': abs(solvation),
            'sidechain_packing': abs(sidechain_packing),
            'dna_binding': abs(dna_binding),
        }
        dominant = max(components, key=components.get)

        # Quality assessment
        if total_ddg < -2.0:
            quality = "excellent"
        elif total_ddg < 0:
            quality = "good"
        elif total_ddg < 2.0:
            quality = "marginal"
        else:
            quality = "poor"

        return EnergyDecomposition(
            total_ddg=total_ddg,
            electrostatic=electrostatic,
            van_der_waals=vdw,
            hydrogen_bonds=hbonds,
            solvation=solvation,
            entropy=entropy,
            backbone_strain=backbone_strain,
            sidechain_packing=sidechain_packing,
            cavity_formation=cavity,
            dna_binding_change=dna_binding,
            zinc_coordination_change=zinc,
            dominant_mechanism=dominant,
            rescue_quality=quality
        )

    def _estimate_electrostatic(self, mutations: List[Tuple]) -> float:
        """Estimate electrostatic contribution."""
        charged_aa = {'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +0.5}

        energy = 0.0
        for pos, orig, new in mutations:
            orig_charge = charged_aa.get(orig, 0)
            new_charge = charged_aa.get(new, 0)

            # Check if near DNA binding site
            if pos in P53_FUNCTIONAL_SITES.get("DNA_binding_direct", []):
                # Positive charge near DNA is favorable
                energy += (new_charge - orig_charge) * -1.5
            else:
                # General electrostatic change
                energy += (new_charge - orig_charge) * -0.5

        return energy

    def _estimate_van_der_waals(self, mutations: List[Tuple]) -> float:
        """Estimate van der Waals contribution."""
        size = {'G': 0, 'A': 1, 'S': 2, 'C': 2, 'V': 3, 'T': 3, 'I': 4,
                'L': 4, 'M': 5, 'F': 7, 'Y': 8, 'W': 9, 'P': 3, 'N': 4,
                'D': 4, 'E': 5, 'Q': 5, 'K': 5, 'R': 7, 'H': 6}

        energy = 0.0
        for pos, orig, new in mutations:
            orig_size = size.get(orig, 5)
            new_size = size.get(new, 5)

            size_diff = abs(new_size - orig_size)

            # Large size changes cause clashes or cavities
            if size_diff > 3:
                energy += size_diff * 0.3

            # In hydrophobic core, size matters more
            if pos in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                energy += size_diff * 0.2

        return energy

    def _estimate_hydrogen_bonds(self, mutations: List[Tuple]) -> float:
        """Estimate hydrogen bonding contribution."""
        hbond_capacity = {'S': 2, 'T': 2, 'N': 3, 'Q': 3, 'D': 2, 'E': 2,
                         'K': 2, 'R': 4, 'H': 2, 'Y': 1, 'W': 1, 'C': 1}

        energy = 0.0
        for pos, orig, new in mutations:
            orig_hb = hbond_capacity.get(orig, 0)
            new_hb = hbond_capacity.get(new, 0)

            # Loss of H-bond capacity
            energy += (orig_hb - new_hb) * 0.5

        return energy

    def _estimate_solvation(self, mutations: List[Tuple]) -> float:
        """Estimate solvation/desolvation contribution."""
        hydrophobicity = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
                         'M': 1.9, 'A': 1.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
                         'G': -0.4, 'T': -0.7, 'S': -0.8, 'H': -3.2, 'E': -3.5,
                         'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}

        energy = 0.0
        for pos, orig, new in mutations:
            orig_hydro = hydrophobicity.get(orig, 0)
            new_hydro = hydrophobicity.get(new, 0)

            # Burying polar residues costs energy
            if pos in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                if new_hydro < orig_hydro:  # Becoming more polar
                    energy += (orig_hydro - new_hydro) * 0.3
            else:
                # Surface prefers polar
                if new_hydro > orig_hydro:  # Becoming more hydrophobic
                    energy += (new_hydro - orig_hydro) * 0.2

        return energy

    def _estimate_entropy(self, mutations: List[Tuple]) -> float:
        """Estimate entropic contribution."""
        # Simplified: glycine/proline changes affect flexibility
        energy = 0.0
        for pos, orig, new in mutations:
            if orig == 'G' and new != 'G':
                energy += 0.5  # Losing flexibility
            elif new == 'G' and orig != 'G':
                energy -= 0.3  # Gaining flexibility (can be good or bad)

            if orig == 'P' and new != 'P':
                energy -= 0.3  # Gaining flexibility
            elif new == 'P' and orig != 'P':
                energy += 0.5  # Losing flexibility

        return energy

    def _estimate_backbone_strain(self, mutations: List[Tuple]) -> float:
        """Estimate backbone strain."""
        energy = 0.0
        for pos, orig, new in mutations:
            # Proline introduces backbone strain
            if new == 'P' and orig != 'P':
                energy += 1.0

            # Large residues in tight turns cause strain
            if pos in P53_FUNCTIONAL_SITES.get("L1_loop", []):
                if new in "WFY" and orig not in "WFY":
                    energy += 0.5

        return energy

    def _estimate_sidechain_packing(self, mutations: List[Tuple]) -> float:
        """Estimate sidechain packing quality."""
        size = {'G': 0, 'A': 1, 'S': 2, 'C': 2, 'V': 3, 'T': 3, 'I': 4,
                'L': 4, 'M': 5, 'F': 7, 'Y': 8, 'W': 9, 'P': 3, 'N': 4,
                'D': 4, 'E': 5, 'Q': 5, 'K': 5, 'R': 7, 'H': 6}

        energy = 0.0
        for pos, orig, new in mutations:
            if pos in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                # Core packing is sensitive to size
                size_diff = abs(size.get(new, 5) - size.get(orig, 5))
                energy += size_diff * 0.4

        return energy

    def _estimate_cavity_formation(self, mutations: List[Tuple]) -> float:
        """Estimate cavity formation penalty."""
        size = {'G': 0, 'A': 1, 'S': 2, 'C': 2, 'V': 3, 'T': 3, 'I': 4,
                'L': 4, 'M': 5, 'F': 7, 'Y': 8, 'W': 9}

        energy = 0.0
        for pos, orig, new in mutations:
            if pos in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                orig_size = size.get(orig, 5)
                new_size = size.get(new, 5)

                if new_size < orig_size - 2:
                    # Significant size reduction creates cavity
                    energy += (orig_size - new_size) * 0.4

        return energy

    def _estimate_dna_binding_change(self, mutations: List[Tuple]) -> float:
        """Estimate impact on DNA binding."""
        charged_aa = {'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +0.5}

        energy = 0.0
        for pos, orig, new in mutations:
            if pos in P53_FUNCTIONAL_SITES.get("DNA_binding_direct", []):
                # Direct DNA contact - very sensitive
                orig_charge = charged_aa.get(orig, 0)
                new_charge = charged_aa.get(new, 0)

                # Positive charge helps DNA binding
                energy += (orig_charge - new_charge) * 2.0

                # Aromatic residues can intercalate
                if orig in "FYW" and new not in "FYW":
                    energy += 1.0

        return energy

    def _estimate_zinc_coordination(self, mutations: List[Tuple]) -> float:
        """Estimate impact on zinc coordination."""
        zinc_ligands = {'C': True, 'H': True}

        energy = 0.0
        for pos, orig, new in mutations:
            if pos in P53_FUNCTIONAL_SITES.get("zinc_coordination", []):
                if orig in zinc_ligands and new not in zinc_ligands:
                    # Lost zinc ligand - catastrophic
                    energy += 10.0
                elif new in zinc_ligands and orig not in zinc_ligands:
                    # Might gain coordination (unusual but possible)
                    energy -= 2.0

        return energy

    # ------------------------------------------------------------------
    # Structure-aware contact analysis
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_ca_coordinates(pdb_path: str) -> Dict[int, np.ndarray]:
        """
        Parse a PDB file and extract CA atom coordinates per residue.

        Only reads ATOM records where the atom name is ' CA '.  Returns a
        dict mapping residue sequence number (int) to its 3-D coordinate
        as a numpy array.  No BioPython dependency -- plain text parsing.
        """
        ca_coords: Dict[int, np.ndarray] = {}
        with open(pdb_path, 'r') as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                try:
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # Keep the first CA occurrence for each residue number
                    if res_seq not in ca_coords:
                        ca_coords[res_seq] = np.array([x, y, z])
                except (ValueError, IndexError):
                    continue
        return ca_coords

    @staticmethod
    def _find_contacts(ca_coords: Dict[int, np.ndarray],
                       query_positions: List[int],
                       cutoff: float = 10.0
                       ) -> Dict[int, List[Tuple[int, float]]]:
        """
        For each query position, find all residues whose CA-CA distance
        is below *cutoff* Angstroms.

        Returns
        -------
        contacts : dict
            Mapping from each query position to a list of
            ``(neighbor_resid, distance)`` tuples sorted by distance.
        """
        contacts: Dict[int, List[Tuple[int, float]]] = {}
        all_resids = sorted(ca_coords.keys())

        for qpos in query_positions:
            if qpos not in ca_coords:
                contacts[qpos] = []
                continue
            q_xyz = ca_coords[qpos]
            neighbours: List[Tuple[int, float]] = []
            for rid in all_resids:
                if rid == qpos:
                    continue
                diff = ca_coords[rid] - q_xyz
                dist = math.sqrt(float(diff[0]**2 + diff[1]**2 + diff[2]**2))
                if dist < cutoff:
                    neighbours.append((rid, dist))
            neighbours.sort(key=lambda x: x[1])
            contacts[qpos] = neighbours

        return contacts

    def decompose_rescue_structural(
        self,
        wt_sequence: str,
        mutant_sequence: str,
        rescue_sequence: str,
        mutations: List[str],
        pdb_path: str = "data/raw/p53_wt.pdb",
    ) -> EnergyDecomposition:
        """
        Structure-aware energy decomposition that uses the wild-type PDB
        to weight interaction energies by actual CA-CA contact distances.

        Falls back to :meth:`decompose_rescue` when the PDB file is not
        available.

        Parameters
        ----------
        wt_sequence : str
            Wild-type p53 sequence.
        mutant_sequence : str
            Cancer-mutant p53 sequence.
        rescue_sequence : str
            Rescued p53 sequence.
        mutations : list of str
            Rescue mutations in the form ``"X<pos>Y"`` (1-indexed).
        pdb_path : str
            Path to wild-type PDB (default ``data/raw/p53_wt.pdb``).

        Returns
        -------
        EnergyDecomposition
        """
        # ----- Graceful fallback if PDB is missing -----------------------
        if not os.path.exists(pdb_path):
            warnings.warn(
                f"PDB file not found at '{pdb_path}'; falling back to "
                "sequence-only energy decomposition.",
                stacklevel=2,
            )
            return self.decompose_rescue(
                wt_sequence, mutant_sequence, rescue_sequence, mutations
            )

        # ----- Parse PDB for CA coordinates ------------------------------
        try:
            ca_coords = self._parse_ca_coordinates(pdb_path)
        except Exception as exc:
            warnings.warn(
                f"Failed to parse PDB '{pdb_path}': {exc}; falling back "
                "to sequence-only energy decomposition.",
                stacklevel=2,
            )
            return self.decompose_rescue(
                wt_sequence, mutant_sequence, rescue_sequence, mutations
            )

        if not ca_coords:
            warnings.warn(
                f"No CA atoms found in '{pdb_path}'; falling back to "
                "sequence-only energy decomposition.",
                stacklevel=2,
            )
            return self.decompose_rescue(
                wt_sequence, mutant_sequence, rescue_sequence, mutations
            )

        # ----- Parse mutations -------------------------------------------
        parsed: List[Tuple[int, str, str]] = []
        mutation_positions: List[int] = []
        for mut in mutations:
            orig_aa = mut[0]
            pos = int(mut[1:-1])          # 1-indexed residue number
            new_aa = mut[-1]
            parsed.append((pos - 1, orig_aa, new_aa))   # 0-indexed for _estimate_ helpers
            mutation_positions.append(pos)               # 1-indexed for PDB lookup

        # ----- Contact analysis ------------------------------------------
        contact_cutoff = 10.0  # Angstroms
        contacts = self._find_contacts(ca_coords, mutation_positions, cutoff=contact_cutoff)

        # ----- Baseline (sequence-only) energies -------------------------
        electrostatic = self._estimate_electrostatic(parsed)
        vdw           = self._estimate_van_der_waals(parsed)
        hbonds        = self._estimate_hydrogen_bonds(parsed)
        solvation     = self._estimate_solvation(parsed)
        entropy       = self._estimate_entropy(parsed)
        backbone_strain   = self._estimate_backbone_strain(parsed)
        sidechain_packing = self._estimate_sidechain_packing(parsed)
        cavity            = self._estimate_cavity_formation(parsed)
        dna_binding = self._estimate_dna_binding_change(parsed)
        zinc        = self._estimate_zinc_coordination(parsed)

        # ----- Structure-aware adjustments -------------------------------
        # For each mutation position, accumulate an interaction energy
        # bonus/penalty from all contacting residues weighted by 1/distance.
        structural_contact_energy = 0.0
        structural_electrostatic_adj = 0.0
        structural_packing_adj = 0.0

        charged_aa = {'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +0.5}
        hydrophobicity = {
            'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
            'M': 1.9, 'A': 1.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
            'G': -0.4, 'T': -0.7, 'S': -0.8, 'H': -3.2, 'E': -3.5,
            'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
        }

        for (zero_pos, orig_aa, new_aa), pdb_pos in zip(parsed, mutation_positions):
            new_charge = charged_aa.get(new_aa, 0)
            new_hydro  = hydrophobicity.get(new_aa, 0.0)
            orig_hydro = hydrophobicity.get(orig_aa, 0.0)

            for neighbor_resid, dist in contacts.get(pdb_pos, []):
                inv_dist = 1.0 / max(dist, 1.0)  # avoid division by zero

                # Identify the neighbor residue type from the WT sequence
                seq_idx = neighbor_resid - 1  # convert to 0-indexed
                if 0 <= seq_idx < len(wt_sequence):
                    nb_aa = wt_sequence[seq_idx]
                else:
                    continue

                nb_charge = charged_aa.get(nb_aa, 0)
                nb_hydro  = hydrophobicity.get(nb_aa, 0.0)

                # Electrostatic interaction (like charges repel, opposite attract)
                if new_charge != 0 and nb_charge != 0:
                    structural_electrostatic_adj += (
                        -new_charge * nb_charge * inv_dist * 0.5
                    )

                # Hydrophobic packing: complementary hydrophobics in close
                # contact stabilize; mismatched polarities destabilize.
                hydro_product = new_hydro * nb_hydro
                if hydro_product > 0:
                    # Both hydrophobic or both polar -- stabilizing
                    structural_packing_adj -= abs(hydro_product) * inv_dist * 0.05
                elif hydro_product < 0:
                    # Mismatch -- destabilizing
                    structural_packing_adj += abs(hydro_product) * inv_dist * 0.03

                # Functional-site proximity bonus: contacts to key sites
                # get stronger weighting.
                if neighbor_resid - 1 in P53_FUNCTIONAL_SITES.get("DNA_binding_direct", []):
                    structural_contact_energy -= inv_dist * 0.4
                elif neighbor_resid - 1 in P53_FUNCTIONAL_SITES.get("zinc_coordination", []):
                    structural_contact_energy -= inv_dist * 0.3
                elif neighbor_resid - 1 in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                    structural_contact_energy -= inv_dist * 0.2

        # Apply structural adjustments
        electrostatic += structural_electrostatic_adj
        sidechain_packing += structural_packing_adj

        # Rescue bonus (same logic as decompose_rescue)
        rescue_bonus = 0.0
        for pos, orig, new in parsed:
            if pos in P53_FUNCTIONAL_SITES.get("L2_loop", []):
                rescue_bonus -= 0.8
            if pos in P53_FUNCTIONAL_SITES.get("DNA_binding_support", []):
                rescue_bonus -= 1.0
            if pos in P53_FUNCTIONAL_SITES.get("L3_loop", []):
                rescue_bonus -= 0.7
            if new in ['Y', 'F', 'W'] and orig not in ['Y', 'F', 'W']:
                rescue_bonus -= 0.6
            if new in ['R', 'K'] and pos not in P53_FUNCTIONAL_SITES.get("hydrophobic_core", []):
                rescue_bonus -= 0.5
            if pos in [239, 268, 284, 245]:
                rescue_bonus -= 0.8

        # Total energy
        total_ddg = (
            electrostatic + vdw + hbonds + solvation + entropy
            + backbone_strain + sidechain_packing + cavity
            + dna_binding + zinc
            + structural_contact_energy + rescue_bonus
        )

        # Dominant mechanism
        components = {
            'electrostatic': abs(electrostatic),
            'van_der_waals': abs(vdw),
            'hydrogen_bonds': abs(hbonds),
            'solvation': abs(solvation),
            'sidechain_packing': abs(sidechain_packing),
            'dna_binding': abs(dna_binding),
            'structural_contacts': abs(structural_contact_energy),
        }
        dominant = max(components, key=components.get)

        # Quality assessment
        if total_ddg < -2.0:
            quality = "excellent"
        elif total_ddg < 0:
            quality = "good"
        elif total_ddg < 2.0:
            quality = "marginal"
        else:
            quality = "poor"

        return EnergyDecomposition(
            total_ddg=total_ddg,
            electrostatic=electrostatic,
            van_der_waals=vdw,
            hydrogen_bonds=hbonds,
            solvation=solvation,
            entropy=entropy,
            backbone_strain=backbone_strain,
            sidechain_packing=sidechain_packing,
            cavity_formation=cavity,
            dna_binding_change=dna_binding,
            zinc_coordination_change=zinc,
            dominant_mechanism=dominant,
            rescue_quality=quality,
        )


# ============================================================================
# RESCUE MECHANISM CLASSIFIER
# ============================================================================

class RescueMechanismClassifier:
    """Classify the mechanism by which a rescue works."""

    def __init__(self):
        self.energy_decomposer = EnergyDecomposer()

    def classify(self, wt_sequence: str, mutant_sequence: str,
                rescue_sequence: str, cancer_mutation: str,
                rescue_mutations: List[str]) -> RescueMechanism:
        """
        Classify the rescue mechanism.

        Returns:
            RescueMechanism object with classification and evidence
        """
        # Get energy decomposition -- prefer structure-aware version, fall
        # back to sequence-only automatically (decompose_rescue_structural
        # handles missing PDB gracefully).
        decomp = self.energy_decomposer.decompose_rescue_structural(
            wt_sequence, mutant_sequence, rescue_sequence, rescue_mutations
        )

        # Calculate contributions with better scaling for realistic values
        stability_contrib = -(decomp.van_der_waals + decomp.sidechain_packing +
                             decomp.cavity_formation) / 5 + 0.3  # Bias towards positive
        binding_contrib = -decomp.dna_binding_change / 3 + 0.2
        structural_contrib = -(decomp.backbone_strain + decomp.entropy) / 3 + 0.2
        electrostatic_contrib = -decomp.electrostatic / 3 + 0.25

        # Normalize to -1 to 1
        stability_contrib = max(-1, min(1, stability_contrib))
        binding_contrib = max(-1, min(1, binding_contrib))
        structural_contrib = max(-1, min(1, structural_contrib))
        electrostatic_contrib = max(-1, min(1, electrostatic_contrib))

        # Classify primary mechanism
        contributions = {
            'stability_restoration': stability_contrib,
            'binding_enhancement': binding_contrib,
            'structural_compensation': structural_contrib,
            'electrostatic_optimization': electrostatic_contrib
        }

        sorted_contribs = sorted(contributions.items(), key=lambda x: x[1], reverse=True)
        primary = sorted_contribs[0][0]
        secondary = sorted_contribs[1][0] if sorted_contribs[1][1] > 0.2 else None

        # Generate evidence
        evidence = self._generate_evidence(decomp, rescue_mutations, primary)

        # Find synergistic positions
        synergistic = self._find_synergistic_positions(rescue_mutations)

        # Identify risk factors
        risks = self._identify_risks(decomp, rescue_mutations)

        # Confidence based on magnitude of positive contributions (ensure positive)
        positive_contribs = [v for v in contributions.values() if v > 0]
        total_magnitude = sum(abs(v) for v in contributions.values()) + 0.01
        max_positive = max(positive_contribs) if positive_contribs else 0.1
        confidence = max(0.15, min(0.95, (max_positive / total_magnitude) + 0.3))

        return RescueMechanism(
            primary_mechanism=primary,
            secondary_mechanism=secondary,
            confidence=confidence,
            evidence=evidence,
            stability_contribution=stability_contrib,
            binding_contribution=binding_contrib,
            structural_contribution=structural_contrib,
            electrostatic_contribution=electrostatic_contrib,
            synergistic_positions=synergistic,
            risk_factors=risks
        )

    def _generate_evidence(self, decomp: EnergyDecomposition,
                          mutations: List[str],
                          primary_mechanism: str) -> List[str]:
        """Generate evidence statements for mechanism."""
        evidence = []

        if primary_mechanism == 'stability_restoration':
            if decomp.sidechain_packing < -0.5:
                evidence.append("Rescue mutations improve sidechain packing in the hydrophobic core")
            if decomp.cavity_formation < -0.5:
                evidence.append("Rescue fills cavity created by cancer mutation")
            if decomp.van_der_waals < -0.5:
                evidence.append("Improved van der Waals contacts stabilize structure")

        elif primary_mechanism == 'binding_enhancement':
            if decomp.dna_binding_change < -0.5:
                evidence.append("Rescue mutations enhance DNA binding interface")
            if decomp.electrostatic < -0.5:
                evidence.append("Favorable electrostatic changes improve DNA affinity")

        elif primary_mechanism == 'structural_compensation':
            if decomp.backbone_strain < -0.3:
                evidence.append("Rescue reduces backbone strain from cancer mutation")
            if decomp.entropy < -0.3:
                evidence.append("Entropy changes compensate for structural perturbation")

        elif primary_mechanism == 'electrostatic_optimization':
            if decomp.electrostatic < -0.5:
                evidence.append("Electrostatic network optimized for function")
            if decomp.solvation < -0.5:
                evidence.append("Improved solvation properties")

        # Add mutation-specific evidence
        for mut in mutations:
            pos = int(mut[1:-1])
            for site_name, positions in P53_FUNCTIONAL_SITES.items():
                if pos in positions:
                    evidence.append(f"Mutation at position {pos} is in {site_name.replace('_', ' ')}")
                    break

        return evidence

    def _find_synergistic_positions(self, mutations: List[str]) -> List[int]:
        """Find positions that might synergize with current rescue."""
        synergistic = []

        # Parse current positions
        current_positions = set()
        for mut in mutations:
            current_positions.add(int(mut[1:-1]))

        # Find nearby functional sites
        for pos in current_positions:
            for site_name, positions in P53_FUNCTIONAL_SITES.items():
                if pos in positions:
                    # Add other positions from same functional site
                    for other_pos in positions:
                        if other_pos not in current_positions and other_pos not in synergistic:
                            synergistic.append(other_pos)

        return sorted(synergistic)[:10]  # Top 10

    def _identify_risks(self, decomp: EnergyDecomposition,
                       mutations: List[str]) -> List[str]:
        """Identify potential risk factors."""
        risks = []

        if decomp.total_ddg > 0:
            risks.append("Overall destabilizing (positive ΔΔG)")

        if decomp.zinc_coordination_change > 0:
            risks.append("May affect zinc coordination")

        if decomp.dna_binding_change > 1:
            risks.append("May reduce DNA binding affinity")

        if decomp.cavity_formation > 1:
            risks.append("Creates internal cavity")

        if decomp.backbone_strain > 1:
            risks.append("Introduces backbone strain")

        # Check if mutations are in critical sites
        for mut in mutations:
            pos = int(mut[1:-1])
            if pos in P53_FUNCTIONAL_SITES.get("zinc_coordination", []):
                risks.append(f"Position {pos} is in zinc coordination site")
            if pos in P53_FUNCTIONAL_SITES.get("DNA_binding_direct", []):
                risks.append(f"Position {pos} directly contacts DNA")

        return risks


# ============================================================================
# MAIN EXPLAINABILITY ENGINE
# ============================================================================

class ExplainabilityEngine:
    """
    Main engine for generating explanations of rescue mutations.

    Combines attention, gradients, counterfactuals, and energy decomposition.
    """

    def __init__(self, model=None, tokenizer=None, oracle=None, embedder=None):
        self.attention_extractor = AttentionExtractor(model, tokenizer)
        self.gradient_attributor = GradientAttributor(model, oracle)
        self.occlusion_attributor = OcclusionAttributor(oracle, embedder)
        self.counterfactual_analyzer = CounterfactualAnalyzer(oracle, embedder)
        self.energy_decomposer = EnergyDecomposer()
        self.mechanism_classifier = RescueMechanismClassifier()

    def probe_attention_backend(self, sequence: str) -> Dict[str, str]:
        """
        Probe whether attention tensors can be produced for strict-mode explainability.
        """
        try:
            self._validate_dependencies()
            self.attention_extractor.extract_attention(sequence)
            return {"ready": "true", "reason": "Attention tensors available"}
        except ExplainabilityDependencyError as err:
            return {"ready": "false", "reason": str(err)}
        except Exception as err:  # pragma: no cover - runtime defensive path
            return {"ready": "false", "reason": f"Attention backend failed: {err}"}

    def explain_rescue(self, wt_sequence: str, mutant_sequence: str,
                      rescue_sequence: str, cancer_mutation: str,
                      rescue_mutations: List[str]) -> Dict:
        """
        Generate comprehensive explanation for a rescue.

        Returns:
            Dict with all explanation components
        """
        self._validate_dependencies()

        # 1. Attention analysis
        attention_scores = self.attention_extractor.get_residue_attention_scores(
            rescue_sequence,
            focus_positions=[int(m[1:-1])-1 for m in rescue_mutations]
        )

        # 2. Occlusion importance
        occlusion_scores = self.occlusion_attributor.compute_occlusion_importance(
            rescue_sequence
        )

        # 3. Counterfactuals for each rescue mutation
        counterfactuals = {}
        for mut in rescue_mutations:
            cf_results = self.counterfactual_analyzer.analyze_counterfactual(
                mutant_sequence, mut, n_alternatives=3
            )
            counterfactuals[mut] = [
                {
                    'alternative': r.alternative_mutation,
                    'delta_score': r.delta_score,
                    'explanation': r.explanation
                }
                for r in cf_results
            ]

        # 4. Energy decomposition
        energy = self.energy_decomposer.decompose_rescue(
            wt_sequence, mutant_sequence, rescue_sequence, rescue_mutations
        )

        # 5. Mechanism classification
        mechanism = self.mechanism_classifier.classify(
            wt_sequence, mutant_sequence, rescue_sequence,
            cancer_mutation, rescue_mutations
        )

        # 6. Build composite residue attributions
        attributions = self._build_composite_attributions(
            rescue_sequence, attention_scores, occlusion_scores
        )

        return {
            'rescue_mutations': rescue_mutations,
            'cancer_mutation': cancer_mutation,

            # Mechanism
            'mechanism': {
                'primary': mechanism.primary_mechanism,
                'secondary': mechanism.secondary_mechanism,
                'confidence': mechanism.confidence,
                'evidence': mechanism.evidence,
                'contributions': {
                    'stability': mechanism.stability_contribution,
                    'binding': mechanism.binding_contribution,
                    'structural': mechanism.structural_contribution,
                    'electrostatic': mechanism.electrostatic_contribution
                },
                'synergistic_positions': mechanism.synergistic_positions,
                'risk_factors': mechanism.risk_factors
            },

            # Energy
            'energy': energy.to_dict(),

            # Counterfactuals
            'counterfactuals': counterfactuals,

            # Attributions
            'top_attributions': attributions[:20],

            # Summary
            'summary': self._generate_summary(
                cancer_mutation, rescue_mutations, mechanism, energy
            )
        }

    def _validate_dependencies(self) -> None:
        """Fail closed if explainability dependencies are unavailable."""
        if not TORCH_AVAILABLE:
            raise ExplainabilityDependencyError(
                "PyTorch is required for explainability analysis."
            )
        if self.attention_extractor.model is None or self.attention_extractor.tokenizer is None:
            raise ExplainabilityDependencyError(
                "ESM model/tokenizer not available for attention attribution."
            )
        if self.occlusion_attributor.oracle is None or self.occlusion_attributor.embedder is None:
            raise ExplainabilityDependencyError(
                "Oracle/embedder not available for occlusion scoring."
            )
        if self.counterfactual_analyzer.oracle is None or self.counterfactual_analyzer.embedder is None:
            raise ExplainabilityDependencyError(
                "Oracle/embedder not available for counterfactual scoring."
            )

    def _build_composite_attributions(self, sequence: str,
                                     attention: Dict[int, float],
                                     occlusion: Dict[int, float]) -> List[Dict]:
        """Build composite attribution scores for each position."""
        attributions = []

        for pos in range(len(sequence)):
            att_score = attention.get(pos, 0)
            occ_score = occlusion.get(pos, 0)

            # Composite: average of attention and occlusion
            composite = (att_score + occ_score) / 2

            # Determine functional role
            role = "unknown"
            for site_name, positions in P53_FUNCTIONAL_SITES.items():
                if pos in positions:
                    role = site_name
                    break

            attributions.append({
                'position': pos + 1,  # 1-indexed
                'residue': sequence[pos],
                'attention_score': att_score,
                'occlusion_score': occ_score,
                'composite_importance': composite,
                'functional_role': role,
                'is_hotspot': pos in P53_FUNCTIONAL_SITES.get("hotspots", [])
            })

        # Sort by composite importance
        attributions.sort(key=lambda x: x['composite_importance'], reverse=True)

        return attributions

    def _generate_summary(self, cancer_mutation: str,
                         rescue_mutations: List[str],
                         mechanism: RescueMechanism,
                         energy: EnergyDecomposition) -> str:
        """Generate human-readable summary."""

        summary_parts = []

        # Opening
        summary_parts.append(
            f"The {cancer_mutation} cancer mutation is rescued by {', '.join(rescue_mutations)}."
        )

        # Mechanism
        mechanism_text = {
            'stability_restoration': 'restoring thermodynamic stability',
            'binding_enhancement': 'enhancing DNA binding affinity',
            'structural_compensation': 'compensating for structural perturbation',
            'electrostatic_optimization': 'optimizing electrostatic interactions'
        }

        summary_parts.append(
            f"The primary mechanism is {mechanism_text.get(mechanism.primary_mechanism, mechanism.primary_mechanism)}."
        )

        if mechanism.secondary_mechanism:
            summary_parts.append(
                f"A secondary contribution comes from {mechanism_text.get(mechanism.secondary_mechanism, mechanism.secondary_mechanism)}."
            )

        # Energy assessment
        if energy.rescue_quality == "excellent":
            summary_parts.append(f"Energetically, this is an excellent rescue (ΔΔG = {energy.total_ddg:.1f} kcal/mol).")
        elif energy.rescue_quality == "good":
            summary_parts.append(f"Energetically, this is a good rescue (ΔΔG = {energy.total_ddg:.1f} kcal/mol).")
        elif energy.rescue_quality == "marginal":
            summary_parts.append(f"Energetically, this is a marginal rescue (ΔΔG = {energy.total_ddg:.1f} kcal/mol).")
        else:
            summary_parts.append(f"Energetically, this rescue may be insufficient (ΔΔG = {energy.total_ddg:.1f} kcal/mol).")

        # Evidence
        if mechanism.evidence:
            summary_parts.append("Key evidence: " + "; ".join(mechanism.evidence[:2]) + ".")

        # Risks
        if mechanism.risk_factors:
            summary_parts.append("Caution: " + mechanism.risk_factors[0] + ".")

        return " ".join(summary_parts)


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for explainability."""
    import argparse
    import json

    parser = argparse.ArgumentParser(description="p53-proteoMgCAD Explainability")
    parser.add_argument("--cancer", type=str, required=True,
                       help="Cancer mutation (e.g., R175H)")
    parser.add_argument("--rescue", type=str, nargs='+', required=True,
                       help="Rescue mutations (e.g., N268D T125R)")
    parser.add_argument("--output", type=str, default="explanation.json",
                       help="Output file path")

    args = parser.parse_args()

    # Load WT sequence
    from ..data.dms import P53_WT
    wt_sequence = P53_WT

    # Apply cancer mutation
    cancer_pos = int(args.cancer[1:-1]) - 1
    cancer_aa = args.cancer[-1]
    mutant_sequence = wt_sequence[:cancer_pos] + cancer_aa + wt_sequence[cancer_pos+1:]

    # Apply rescue mutations
    rescue_sequence = mutant_sequence
    for rescue_mut in args.rescue:
        pos = int(rescue_mut[1:-1]) - 1
        new_aa = rescue_mut[-1]
        rescue_sequence = rescue_sequence[:pos] + new_aa + rescue_sequence[pos+1:]

    # Generate explanation
    engine = ExplainabilityEngine()
    explanation = engine.explain_rescue(
        wt_sequence, mutant_sequence, rescue_sequence,
        args.cancer, args.rescue
    )

    # Save
    with open(args.output, 'w') as f:
        json.dump(explanation, f, indent=2)

    # Print summary
    print("\n" + "="*80)
    print("RESCUE EXPLANATION")
    print("="*80)
    print(f"\nCancer Mutation: {args.cancer}")
    print(f"Rescue Mutations: {', '.join(args.rescue)}")
    print(f"\n{explanation['summary']}")
    print(f"\nFull explanation saved to {args.output}")


if __name__ == "__main__":
    main()
