"""
Multi-Target Generalization Module for p53-proteoMgCAD

Expands the platform beyond p53 to other tumor suppressors:
1. Other Tumor Suppressors - BRCA1, PTEN, RB1, APC, VHL
2. Universal Rescue Scaffold - Mutations that rescue multiple variants
3. Transfer Learning Analysis - Does p53 training help other proteins?
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any, Set
from dataclasses import dataclass, field
from collections import defaultdict
import json


# ============================================================================
# TUMOR SUPPRESSOR DATABASE
# ============================================================================

TUMOR_SUPPRESSORS = {
    "TP53": {
        "name": "Tumor Protein p53",
        "uniprot": "P04637",
        "length": 393,
        "domains": {
            "TAD1": (1, 40),      # Transactivation domain 1
            "TAD2": (40, 60),     # Transactivation domain 2
            "PRD": (60, 92),      # Proline-rich domain
            "DBD": (102, 292),    # DNA-binding domain
            "TD": (323, 356),     # Tetramerization domain
            "CTD": (363, 393),    # C-terminal domain
        },
        "hotspots": [175, 220, 245, 248, 249, 273, 282],
        "mutation_frequency_cancer": 50,  # % of cancers
        "esm_model": "facebook/esm2_t6_8M_UR50D",
        "known_rescues": {
            "R175H": ["N268D", "H178Y", "N239Y"],
            "Y220C": ["PhiKan083 (drug)"],
            "R248Q": ["T123A", "H178Y"],
            "R273H": ["T284R"],
        },
        "structure_pdb": "1TUP",
    },

    "BRCA1": {
        "name": "BRCA1 DNA Repair Associated",
        "uniprot": "P38398",
        "length": 1863,
        "domains": {
            "RING": (1, 109),         # RING finger domain
            "CC": (1364, 1437),       # Coiled-coil
            "BRCT1": (1642, 1736),    # BRCT repeat 1
            "BRCT2": (1756, 1855),    # BRCT repeat 2
        },
        "hotspots": [1699, 1708, 1775, 1853],
        "mutation_frequency_cancer": 15,  # breast/ovarian
        "esm_model": "facebook/esm2_t33_650M_UR50D",  # Needs larger model
        "known_rescues": {},
        "structure_pdb": "1JM7",
    },

    "PTEN": {
        "name": "Phosphatase and Tensin Homolog",
        "uniprot": "P60484",
        "length": 403,
        "domains": {
            "PBD": (1, 13),           # PIP2-binding domain
            "Phosphatase": (14, 185), # Phosphatase domain
            "C2": (186, 351),         # C2 domain
            "CTail": (352, 403),      # C-terminal tail
        },
        "hotspots": [124, 130, 136, 173, 233, 267, 289],
        "mutation_frequency_cancer": 30,  # prostate, glioblastoma
        "esm_model": "facebook/esm2_t6_8M_UR50D",
        "known_rescues": {},
        "structure_pdb": "1D5R",
    },

    "RB1": {
        "name": "RB Transcriptional Corepressor 1",
        "uniprot": "P06400",
        "length": 928,
        "domains": {
            "RbN": (52, 355),         # N-terminal domain
            "PocketA": (380, 579),    # Pocket domain A
            "PocketB": (640, 771),    # Pocket domain B
            "RbC": (772, 928),        # C-terminal domain
        },
        "hotspots": [455, 567, 661, 698],
        "mutation_frequency_cancer": 25,  # retinoblastoma, osteosarcoma
        "esm_model": "facebook/esm2_t33_650M_UR50D",
        "known_rescues": {},
        "structure_pdb": "1GUX",
    },

    "APC": {
        "name": "APC Regulator of WNT Signaling Pathway",
        "uniprot": "P25054",
        "length": 2843,
        "domains": {
            "Oligomerization": (1, 55),
            "Armadillo": (453, 767),
            "15R": (1020, 1170),       # 15-aa repeats
            "20R": (1262, 2033),       # 20-aa repeats
            "SAMP": (1562, 2053),      # Axin binding
            "Basic": (2200, 2400),     # Microtubule binding
        },
        "hotspots": [1309, 1450],      # Mutation cluster region
        "mutation_frequency_cancer": 80,  # colorectal
        "esm_model": "facebook/esm2_t33_650M_UR50D",
        "known_rescues": {},
        "structure_pdb": "1DEB",
    },

    "VHL": {
        "name": "Von Hippel-Lindau Tumor Suppressor",
        "uniprot": "P40337",
        "length": 213,
        "domains": {
            "Beta": (63, 154),         # Beta domain (HIF binding)
            "Alpha": (155, 204),       # Alpha domain (Elongin binding)
        },
        "hotspots": [98, 117, 161, 167],
        "mutation_frequency_cancer": 70,  # clear cell renal
        "esm_model": "facebook/esm2_t6_8M_UR50D",
        "known_rescues": {},
        "structure_pdb": "1LM8",
    },
}


# Amino acid sequences for tumor suppressors (truncated for readability)
SEQUENCES = {
    "TP53": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",

    "PTEN": "MTAIIKEIVSRNKRRYQEDGFDLDLTYIYPNIIAMGFPAERLEGVYRNNIDDVVRFLDSKHKNHYKIYNLCAERHYDTAKFNCRVAQYPFEDHNPPQLELIKPFCEDLDQWLSEDDNHVAAIHCKAGKGRTGVMICAYLLHRGKFLKAQEALDFYGEVRTRDKKGVTIPSQRRYVYYYSYLLKNHLDYRPVALLFHKMMFETIPMFSGGTCNPQFVVCQLKVKIYSSNSGPTRREDKFMYFEFPQPLPVCGDIKVEFFHKQNKMLKKDKMFHFWVNTFFIPGPEETSEKVENGSLCDQEIDSICSIERADNDKEYLVLTLTKNDLDKANKDKANRYFSPNFKVKLYFTKTVEEPSNPEASSSTSVTPDVSDNEPDHYRYSDTTDSDPENEPFDEDQHTQITKV",

    "VHL": "MPRRAENWDEAEVGAEEAGVEEYGPEEDGGEESGAEESGPEESGPEELGAEEEMEAGRPRPVLRSVNSREPSQVIFCNRSPRPVVLPVWLNFDGEPQPYPTLPPGTGRRIHSYRGHLWLFRDAGTHDGLLVNQTELFVPSLNVDGQPIFANITLPVYTLKERCLQVVRSLVKPENYRRLDIVRSLYEDLEDHPNVQKDLERLTQERIAHQRMGD",
}


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class TumorSuppressorProfile:
    """Profile of a tumor suppressor gene."""
    gene: str
    name: str
    length: int
    domains: Dict[str, Tuple[int, int]]
    hotspots: List[int]
    cancer_frequency: int
    sequence: str
    known_rescues: Dict[str, List[str]]


@dataclass
class UniversalRescue:
    """A rescue mutation that works across multiple variants."""
    position: int
    original_aa: str
    rescue_aa: str
    rescued_variants: List[str]
    rescue_efficiency: Dict[str, float]  # variant -> efficiency score
    mechanism: str
    universality_score: float  # 0-1, how broadly applicable


@dataclass
class TransferLearningResult:
    """Results from transfer learning analysis."""
    source_protein: str
    target_protein: str
    embedding_correlation: float
    attention_pattern_similarity: float
    rescue_prediction_accuracy: float
    transferable_features: List[str]
    recommendations: List[str]


@dataclass
class MultiTargetDesign:
    """Design that targets multiple tumor suppressors."""
    name: str
    targets: List[str]  # Gene names
    mutations_per_target: Dict[str, List[str]]
    combined_patient_impact: int
    development_complexity: str
    synergy_score: float


# ============================================================================
# UNIVERSAL RESCUE FINDER
# ============================================================================

class UniversalRescueFinder:
    """Find rescue mutations that work across multiple p53 variants."""

    def __init__(self):
        self.p53_sequence = SEQUENCES.get("TP53", "")
        self.known_rescues = TUMOR_SUPPRESSORS["TP53"]["known_rescues"]

    def find_universal_rescues(self, cancer_mutations: List[str],
                               min_coverage: float = 0.5) -> List[UniversalRescue]:
        """
        Find positions where a single mutation rescues multiple cancer variants.

        Args:
            cancer_mutations: List of cancer mutations to analyze
            min_coverage: Minimum fraction of mutations that must be rescued

        Returns:
            List of UniversalRescue objects
        """
        # Analyze each potential rescue position
        universal_rescues = []

        # Known rescue positions from literature
        rescue_positions = self._get_known_rescue_positions()

        # Add positions near hotspots (compensatory)
        hotspots = TUMOR_SUPPRESSORS["TP53"]["hotspots"]
        for hotspot in hotspots:
            # Positions within 10-50 residues often compensate
            for offset in [-30, -20, -10, 10, 20, 30]:
                pos = hotspot + offset
                if 1 <= pos <= 393:
                    rescue_positions.add(pos)

        # Test each position
        for pos in sorted(rescue_positions):
            result = self._evaluate_position(pos, cancer_mutations)
            if result and result.universality_score >= min_coverage:
                universal_rescues.append(result)

        # Sort by universality score
        universal_rescues.sort(key=lambda x: x.universality_score, reverse=True)

        return universal_rescues

    def _get_known_rescue_positions(self) -> Set[int]:
        """Extract positions of known rescue mutations."""
        positions = set()
        for variant, rescues in self.known_rescues.items():
            for rescue in rescues:
                # Parse mutation string like "N268D"
                if rescue and rescue[0].isalpha() and rescue[-1].isalpha():
                    try:
                        pos = int(rescue[1:-1])
                        positions.add(pos)
                    except ValueError:
                        pass
        return positions

    def _evaluate_position(self, position: int,
                          cancer_mutations: List[str]) -> Optional[UniversalRescue]:
        """Evaluate if a position can serve as universal rescue."""

        if position < 1 or position > len(self.p53_sequence):
            return None

        original_aa = self.p53_sequence[position - 1]

        # Find best rescue amino acid for this position
        best_aa = None
        best_efficiency = {}
        best_count = 0

        for test_aa in "ACDEFGHIKLMNPQRSTVWY":
            if test_aa == original_aa:
                continue

            efficiency = {}
            rescued_count = 0

            for cancer_mut in cancer_mutations:
                # Simulate rescue efficiency (placeholder)
                eff = self._simulate_rescue_efficiency(position, test_aa, cancer_mut)
                if eff > 0.3:  # Threshold for "rescued"
                    rescued_count += 1
                efficiency[cancer_mut] = eff

            if rescued_count > best_count:
                best_count = rescued_count
                best_aa = test_aa
                best_efficiency = efficiency

        if best_aa and best_count > 0:
            universality = best_count / len(cancer_mutations)

            return UniversalRescue(
                position=position,
                original_aa=original_aa,
                rescue_aa=best_aa,
                rescued_variants=[m for m, e in best_efficiency.items() if e > 0.3],
                rescue_efficiency=best_efficiency,
                mechanism=self._infer_mechanism(position),
                universality_score=universality
            )

        return None

    def _simulate_rescue_efficiency(self, rescue_pos: int, rescue_aa: str,
                                   cancer_mutation: str) -> float:
        """
        Simulate rescue efficiency (placeholder for real scoring).

        In production, would use the FMR oracle or physics-based scoring.
        """
        cancer_pos = int(cancer_mutation[1:-1])

        # Distance-based heuristic
        distance = abs(rescue_pos - cancer_pos)

        # Positions 10-50 residues away often compensate
        if 10 <= distance <= 50:
            base_score = 0.5
        elif distance < 10:
            base_score = 0.3  # Too close might cause issues
        else:
            base_score = 0.2  # Too far for direct compensation

        # Domain-based modulation
        dbd_start, dbd_end = TUMOR_SUPPRESSORS["TP53"]["domains"]["DBD"]
        if dbd_start <= rescue_pos <= dbd_end and dbd_start <= cancer_pos <= dbd_end:
            base_score += 0.2  # Both in DBD

        # Random variation for simulation
        noise = np.random.normal(0, 0.1)

        return max(0, min(1, base_score + noise))

    def _infer_mechanism(self, position: int) -> str:
        """Infer rescue mechanism from position."""
        domains = TUMOR_SUPPRESSORS["TP53"]["domains"]

        for domain_name, (start, end) in domains.items():
            if start <= position <= end:
                if domain_name == "DBD":
                    return "DNA_binding_compensation"
                elif domain_name == "TD":
                    return "tetramerization_stabilization"
                else:
                    return f"{domain_name}_stabilization"

        return "allosteric_compensation"

    def design_universal_scaffold(self, target_mutations: List[str],
                                  max_rescue_mutations: int = 5) -> Dict:
        """
        Design a universal p53 scaffold that rescues multiple variants.

        Returns a scaffold with positions to mutate for each target.
        """
        # Find universal rescues
        universal = self.find_universal_rescues(target_mutations)

        # Greedy selection to cover all targets
        selected_rescues = []
        covered = set()

        for rescue in universal:
            if len(selected_rescues) >= max_rescue_mutations:
                break

            new_coverage = set(rescue.rescued_variants) - covered
            if new_coverage:
                selected_rescues.append(rescue)
                covered.update(rescue.rescued_variants)

        return {
            "scaffold_name": f"Universal_p53_rescue_v{len(selected_rescues)}",
            "target_mutations": target_mutations,
            "rescue_mutations": [
                f"{r.original_aa}{r.position}{r.rescue_aa}"
                for r in selected_rescues
            ],
            "coverage": len(covered) / len(target_mutations) if target_mutations else 0,
            "covered_variants": list(covered),
            "uncovered_variants": list(set(target_mutations) - covered),
            "mechanism_summary": [r.mechanism for r in selected_rescues],
        }


# ============================================================================
# MULTI-TARGET DESIGN ENGINE
# ============================================================================

class MultiTargetDesignEngine:
    """Design rescue strategies across multiple tumor suppressors."""

    def __init__(self):
        self.proteins = TUMOR_SUPPRESSORS
        self.sequences = SEQUENCES

    def get_protein_profile(self, gene: str) -> Optional[TumorSuppressorProfile]:
        """Get profile for a tumor suppressor."""
        if gene not in self.proteins:
            return None

        info = self.proteins[gene]

        return TumorSuppressorProfile(
            gene=gene,
            name=info["name"],
            length=info["length"],
            domains=info["domains"],
            hotspots=info["hotspots"],
            cancer_frequency=info["mutation_frequency_cancer"],
            sequence=self.sequences.get(gene, ""),
            known_rescues=info.get("known_rescues", {})
        )

    def compare_proteins(self, gene1: str, gene2: str) -> Dict:
        """Compare two tumor suppressors for transferability."""

        prof1 = self.get_protein_profile(gene1)
        prof2 = self.get_protein_profile(gene2)

        if not prof1 or not prof2:
            return {"error": "Unknown protein"}

        # Compare domain architecture
        domain_similarity = self._compare_domains(prof1.domains, prof2.domains)

        # Compare hotspot patterns
        hotspot_similarity = self._compare_hotspots(
            prof1.hotspots, prof1.length,
            prof2.hotspots, prof2.length
        )

        # Sequence similarity (if sequences available)
        if prof1.sequence and prof2.sequence:
            seq_identity = self._calculate_sequence_identity(
                prof1.sequence, prof2.sequence
            )
        else:
            seq_identity = 0

        return {
            "gene1": gene1,
            "gene2": gene2,
            "domain_architecture_similarity": domain_similarity,
            "hotspot_pattern_similarity": hotspot_similarity,
            "sequence_identity": seq_identity,
            "transfer_learning_potential": (
                domain_similarity * 0.3 +
                hotspot_similarity * 0.3 +
                seq_identity * 0.4
            ),
            "recommendations": self._generate_transfer_recommendations(
                gene1, gene2, domain_similarity, hotspot_similarity
            )
        }

    def _compare_domains(self, domains1: Dict, domains2: Dict) -> float:
        """Compare domain architectures."""
        # Count similar domain types
        domain_types1 = set()
        domain_types2 = set()

        type_keywords = ["DBD", "DNA", "binding", "phosphatase", "kinase",
                        "RING", "BRCT", "C2", "coil", "tetramer"]

        for name in domains1:
            for kw in type_keywords:
                if kw.lower() in name.lower():
                    domain_types1.add(kw.lower())

        for name in domains2:
            for kw in type_keywords:
                if kw.lower() in name.lower():
                    domain_types2.add(kw.lower())

        if not domain_types1 or not domain_types2:
            return 0.1

        overlap = len(domain_types1 & domain_types2)
        union = len(domain_types1 | domain_types2)

        return overlap / union if union > 0 else 0

    def _compare_hotspots(self, hotspots1: List[int], length1: int,
                         hotspots2: List[int], length2: int) -> float:
        """Compare hotspot distribution patterns."""
        # Normalize positions to 0-1 scale
        norm1 = [h / length1 for h in hotspots1]
        norm2 = [h / length2 for h in hotspots2]

        # Compare distributions
        if not norm1 or not norm2:
            return 0.1

        # Simple comparison: how similar are the normalized positions?
        min_dist = 0
        for n1 in norm1:
            dists = [abs(n1 - n2) for n2 in norm2]
            min_dist += min(dists)

        avg_dist = min_dist / len(norm1)

        # Convert distance to similarity (closer = more similar)
        similarity = max(0, 1 - avg_dist * 2)

        return similarity

    def _calculate_sequence_identity(self, seq1: str, seq2: str) -> float:
        """Calculate global sequence identity (simplified)."""
        # For very different length sequences, use shorter alignment
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))

        if max_len == 0:
            return 0

        # Simple identity without alignment
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])

        return matches / max_len

    def _generate_transfer_recommendations(self, gene1: str, gene2: str,
                                          domain_sim: float,
                                          hotspot_sim: float) -> List[str]:
        """Generate recommendations for transfer learning."""
        recommendations = []

        if domain_sim > 0.5:
            recommendations.append(
                f"Similar domain architecture suggests {gene1} rescue strategies "
                f"may inform {gene2} design"
            )

        if hotspot_sim > 0.5:
            recommendations.append(
                f"Similar hotspot patterns - consider analogous positions for rescue"
            )

        if gene1 == "TP53" and gene2 in ["PTEN", "VHL"]:
            recommendations.append(
                f"Both are relatively short proteins - ESM-2 t6 model sufficient"
            )
        elif gene2 in ["BRCA1", "RB1", "APC"]:
            recommendations.append(
                f"{gene2} requires larger ESM-2 model (t33) due to length"
            )

        if not recommendations:
            recommendations.append(
                f"Limited structural similarity - independent optimization recommended"
            )

        return recommendations

    def design_multi_target_therapy(self, targets: List[str]) -> MultiTargetDesign:
        """
        Design therapy targeting multiple tumor suppressors.

        This is conceptual - in practice would require separate delivery
        or a very sophisticated multi-gene construct.
        """
        mutations_per_target = {}
        total_patients = 0

        for gene in targets:
            if gene not in self.proteins:
                continue

            info = self.proteins[gene]

            # Get top hotspots
            hotspots = info["hotspots"][:3]  # Top 3 mutations
            mutations_per_target[gene] = [f"Position {h}" for h in hotspots]

            # Estimate patients
            # Rough: cancer_frequency * relevant cancer incidence
            total_patients += info["mutation_frequency_cancer"] * 1000

        # Complexity assessment
        if len(targets) == 1:
            complexity = "low"
        elif len(targets) == 2:
            complexity = "medium"
        else:
            complexity = "high"

        return MultiTargetDesign(
            name=f"Multi-TSG_{'-'.join(targets)}",
            targets=targets,
            mutations_per_target=mutations_per_target,
            combined_patient_impact=total_patients,
            development_complexity=complexity,
            synergy_score=0.7 if len(targets) > 1 else 1.0  # Discount for complexity
        )


# ============================================================================
# TRANSFER LEARNING ANALYZER
# ============================================================================

class TransferLearningAnalyzer:
    """Analyze transfer learning potential between proteins."""

    def __init__(self):
        self.engine = MultiTargetDesignEngine()

    def analyze_transfer(self, source: str, target: str) -> TransferLearningResult:
        """
        Analyze how well rescue strategies transfer between proteins.
        """
        comparison = self.engine.compare_proteins(source, target)

        # Embedding correlation (placeholder - would need actual embeddings)
        embed_corr = self._estimate_embedding_correlation(source, target)

        # Attention pattern similarity
        attention_sim = comparison.get("domain_architecture_similarity", 0)

        # Rescue prediction accuracy (would need actual validation)
        rescue_accuracy = comparison.get("transfer_learning_potential", 0)

        # Transferable features
        transferable = self._identify_transferable_features(source, target)

        # Recommendations
        recommendations = comparison.get("recommendations", [])

        return TransferLearningResult(
            source_protein=source,
            target_protein=target,
            embedding_correlation=embed_corr,
            attention_pattern_similarity=attention_sim,
            rescue_prediction_accuracy=rescue_accuracy,
            transferable_features=transferable,
            recommendations=recommendations
        )

    def _estimate_embedding_correlation(self, source: str, target: str) -> float:
        """Estimate ESM-2 embedding correlation."""
        # Placeholder - in production would compute actual embeddings
        # Correlation tends to be higher for similar length/function proteins

        source_len = TUMOR_SUPPRESSORS.get(source, {}).get("length", 0)
        target_len = TUMOR_SUPPRESSORS.get(target, {}).get("length", 0)

        if source_len == 0 or target_len == 0:
            return 0

        # Length similarity impacts embedding similarity
        len_ratio = min(source_len, target_len) / max(source_len, target_len)

        # Base correlation from length similarity
        base_corr = len_ratio * 0.6

        # Bonus for functional similarity (both are tumor suppressors)
        functional_bonus = 0.2

        return min(1.0, base_corr + functional_bonus)

    def _identify_transferable_features(self, source: str, target: str) -> List[str]:
        """Identify features that transfer between proteins."""
        features = []

        source_info = TUMOR_SUPPRESSORS.get(source, {})
        target_info = TUMOR_SUPPRESSORS.get(target, {})

        # Check for similar domain types
        source_domains = set(source_info.get("domains", {}).keys())
        target_domains = set(target_info.get("domains", {}).keys())

        if "DBD" in source_domains and any("DNA" in d or "binding" in d.lower()
                                           for d in target_domains):
            features.append("DNA binding domain architecture")

        if len(source_domains) > 2 and len(target_domains) > 2:
            features.append("Multi-domain organization")

        # Check for similar size
        source_len = source_info.get("length", 0)
        target_len = target_info.get("length", 0)
        if 0.5 <= source_len / max(target_len, 1) <= 2.0:
            features.append("Similar protein size")

        # Universal features
        features.append("ESM-2 latent space representation")
        features.append("Secondary structure patterns")
        features.append("Hydrophobicity profiles")

        return features

    def recommend_model_for_protein(self, gene: str) -> Dict:
        """Recommend ESM-2 model size for a protein."""
        info = TUMOR_SUPPRESSORS.get(gene, {})
        length = info.get("length", 500)

        if length <= 500:
            model = "facebook/esm2_t6_8M_UR50D"
            reason = "Short protein, small model sufficient"
        elif length <= 1000:
            model = "facebook/esm2_t12_35M_UR50D"
            reason = "Medium protein, medium model recommended"
        else:
            model = "facebook/esm2_t33_650M_UR50D"
            reason = "Long protein, large model required for accurate embeddings"

        return {
            "gene": gene,
            "length": length,
            "recommended_model": model,
            "reason": reason,
            "alternative_models": [
                "facebook/esm2_t6_8M_UR50D",
                "facebook/esm2_t12_35M_UR50D",
                "facebook/esm2_t33_650M_UR50D",
            ]
        }


# ============================================================================
# MAIN INTERFACE
# ============================================================================

class MultiTargetPlatform:
    """
    Main platform for multi-target generalization.
    """

    def __init__(self):
        self.universal_finder = UniversalRescueFinder()
        self.design_engine = MultiTargetDesignEngine()
        self.transfer_analyzer = TransferLearningAnalyzer()

    def get_available_targets(self) -> List[Dict]:
        """Get list of available tumor suppressor targets."""
        targets = []
        for gene, info in TUMOR_SUPPRESSORS.items():
            targets.append({
                "gene": gene,
                "name": info["name"],
                "length": info["length"],
                "cancer_frequency": info["mutation_frequency_cancer"],
                "n_hotspots": len(info["hotspots"]),
                "has_known_rescues": bool(info.get("known_rescues")),
            })
        return targets

    def analyze_protein(self, gene: str) -> Dict:
        """Get comprehensive analysis for a tumor suppressor."""
        profile = self.design_engine.get_protein_profile(gene)

        if not profile:
            return {"error": f"Unknown protein: {gene}"}

        # Model recommendation
        model_rec = self.transfer_analyzer.recommend_model_for_protein(gene)

        # Transfer from p53
        if gene != "TP53":
            transfer = self.transfer_analyzer.analyze_transfer("TP53", gene)
            transfer_info = {
                "from_p53": {
                    "embedding_correlation": transfer.embedding_correlation,
                    "transfer_potential": transfer.rescue_prediction_accuracy,
                    "transferable_features": transfer.transferable_features,
                    "recommendations": transfer.recommendations,
                }
            }
        else:
            transfer_info = {"from_p53": "N/A (source protein)"}

        return {
            "gene": gene,
            "name": profile.name,
            "length": profile.length,
            "domains": profile.domains,
            "hotspots": profile.hotspots,
            "cancer_frequency": f"{profile.cancer_frequency}% of relevant cancers",
            "known_rescues": profile.known_rescues,
            "recommended_model": model_rec["recommended_model"],
            "model_reason": model_rec["reason"],
            "transfer_learning": transfer_info,
        }

    def find_universal_p53_rescues(self, mutations: Optional[List[str]] = None) -> Dict:
        """Find mutations that rescue multiple p53 variants."""
        if mutations is None:
            # Use common hotspot mutations
            mutations = ["R175H", "Y220C", "G245S", "R248Q", "R248W",
                        "R249S", "R273H", "R273C", "R282W"]

        universal = self.universal_finder.find_universal_rescues(mutations)

        return {
            "target_mutations": mutations,
            "universal_rescues": [
                {
                    "mutation": f"{r.original_aa}{r.position}{r.rescue_aa}",
                    "rescued_variants": r.rescued_variants,
                    "universality_score": r.universality_score,
                    "mechanism": r.mechanism,
                }
                for r in universal[:10]  # Top 10
            ],
            "scaffold_design": self.universal_finder.design_universal_scaffold(mutations)
        }

    def export_multi_target_report(self, output_path: str) -> str:
        """Export comprehensive multi-target report."""

        lines = []
        lines.append("# Multi-Target Tumor Suppressor Platform Report")
        lines.append(f"\n## Available Targets")

        for target in self.get_available_targets():
            lines.append(f"\n### {target['gene']} - {target['name']}")
            lines.append(f"- Length: {target['length']} aa")
            lines.append(f"- Cancer frequency: {target['cancer_frequency']}%")
            lines.append(f"- Known hotspots: {target['n_hotspots']}")
            lines.append(f"- Has known rescues: {'Yes' if target['has_known_rescues'] else 'No'}")

        lines.append(f"\n## Transfer Learning Analysis")
        lines.append("\nHow well do p53 rescue strategies transfer to other proteins?")

        for gene in ["PTEN", "VHL", "BRCA1"]:
            transfer = self.transfer_analyzer.analyze_transfer("TP53", gene)
            lines.append(f"\n### TP53 → {gene}")
            lines.append(f"- Embedding correlation: {transfer.embedding_correlation:.2f}")
            lines.append(f"- Transfer potential: {transfer.rescue_prediction_accuracy:.2f}")
            lines.append(f"- Transferable features: {', '.join(transfer.transferable_features[:3])}")

        lines.append(f"\n## Universal p53 Rescue Scaffold")
        universal = self.find_universal_p53_rescues()
        lines.append(f"\nTarget mutations: {', '.join(universal['target_mutations'])}")
        lines.append(f"\nTop universal rescues:")
        for rescue in universal['universal_rescues'][:5]:
            lines.append(f"- {rescue['mutation']}: rescues {len(rescue['rescued_variants'])} variants "
                        f"(score: {rescue['universality_score']:.2f})")

        content = "\n".join(lines)

        with open(output_path, 'w') as f:
            f.write(content)

        return output_path


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for multi-target platform."""
    import argparse

    parser = argparse.ArgumentParser(description="p53-proteoMgCAD Multi-Target Platform")
    parser.add_argument("--analyze", type=str, help="Analyze a specific tumor suppressor")
    parser.add_argument("--compare", type=str, nargs=2, help="Compare two proteins")
    parser.add_argument("--universal", action="store_true", help="Find universal p53 rescues")
    parser.add_argument("--report", type=str, help="Generate full report to file")

    args = parser.parse_args()

    platform = MultiTargetPlatform()

    if args.analyze:
        result = platform.analyze_protein(args.analyze)
        print(json.dumps(result, indent=2))

    elif args.compare:
        result = platform.design_engine.compare_proteins(args.compare[0], args.compare[1])
        print(json.dumps(result, indent=2))

    elif args.universal:
        result = platform.find_universal_p53_rescues()
        print(json.dumps(result, indent=2))

    elif args.report:
        path = platform.export_multi_target_report(args.report)
        print(f"Report saved to: {path}")

    else:
        # Default: list available targets
        print("\nAvailable Tumor Suppressor Targets:")
        print("-" * 60)
        for target in platform.get_available_targets():
            print(f"  {target['gene']:8s} | {target['name'][:30]:30s} | "
                  f"{target['cancer_frequency']:2d}% cancers | "
                  f"{target['length']:4d} aa")


if __name__ == "__main__":
    main()
