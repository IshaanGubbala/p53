"""
Clinical Impact Quantification Module for p53-proteoMgCAD

Quantifies the clinical relevance of rescue mutations:
1. Patient Stratification - How many patients would benefit (TCGA data)
2. Therapeutic Index - Safety margins based on sequence identity
3. Delivery Feasibility - Scoring for different delivery methods
4. Immunogenicity Risk - Assessment of immune response potential
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
import json


# ============================================================================
# TCGA MUTATION FREQUENCY DATA
# ============================================================================

# From TCGA pan-cancer analysis (approximate frequencies)
# Based on cBioPortal and COSMIC database
TCGA_P53_MUTATIONS = {
    # Hotspot mutations (% of all p53 mutations)
    "R175H": {"frequency": 5.6, "cancer_types": ["breast", "colon", "ovarian", "lung", "pancreatic"]},
    "R248Q": {"frequency": 4.4, "cancer_types": ["colon", "breast", "ovarian", "lung", "gastric"]},
    "R273H": {"frequency": 4.2, "cancer_types": ["colon", "breast", "lung", "head_neck", "bladder"]},
    "R248W": {"frequency": 3.8, "cancer_types": ["colon", "breast", "ovarian", "lung"]},
    "R273C": {"frequency": 2.1, "cancer_types": ["colon", "breast", "lung"]},
    "Y220C": {"frequency": 1.8, "cancer_types": ["breast", "lung", "colon", "ovarian"]},
    "G245S": {"frequency": 2.4, "cancer_types": ["breast", "colon", "lung", "ovarian"]},
    "R249S": {"frequency": 2.9, "cancer_types": ["liver", "lung", "esophageal"]},  # Aflatoxin-associated
    "R282W": {"frequency": 2.0, "cancer_types": ["breast", "colon", "ovarian", "lung"]},
    "G245D": {"frequency": 1.5, "cancer_types": ["colon", "breast", "lung"]},
    "R175G": {"frequency": 0.8, "cancer_types": ["breast", "colon"]},
    "C176F": {"frequency": 0.7, "cancer_types": ["breast", "ovarian"]},
    "H179R": {"frequency": 0.9, "cancer_types": ["breast", "colon", "lung"]},
    "H179Y": {"frequency": 0.6, "cancer_types": ["breast", "colon"]},
    "C242F": {"frequency": 0.8, "cancer_types": ["lung", "colon"]},
    "C242S": {"frequency": 0.5, "cancer_types": ["lung", "colon"]},
    "G245C": {"frequency": 0.7, "cancer_types": ["breast", "colon"]},
    "R280K": {"frequency": 0.6, "cancer_types": ["breast", "colon", "lung"]},
    "D281G": {"frequency": 0.4, "cancer_types": ["breast", "colon"]},
    "E285K": {"frequency": 0.5, "cancer_types": ["colon", "breast"]},

    # Contact mutations
    "R273L": {"frequency": 0.9, "cancer_types": ["colon", "breast", "lung"]},
    "S241F": {"frequency": 0.6, "cancer_types": ["breast", "colon"]},
    "R248L": {"frequency": 0.5, "cancer_types": ["colon", "breast"]},

    # Structural mutations
    "V143A": {"frequency": 0.4, "cancer_types": ["breast", "colon"]},
    "V157F": {"frequency": 0.5, "cancer_types": ["lung", "colon"]},
    "Y163C": {"frequency": 0.4, "cancer_types": ["breast", "lung"]},
    "V173L": {"frequency": 0.3, "cancer_types": ["colon", "breast"]},
    "I195T": {"frequency": 0.3, "cancer_types": ["breast", "colon"]},
    "I255F": {"frequency": 0.3, "cancer_types": ["lung", "colon"]},
}

# Cancer incidence data (new cases per year, US)
CANCER_INCIDENCE = {
    "breast": 290000,
    "lung": 235000,
    "colon": 150000,
    "ovarian": 20000,
    "pancreatic": 62000,
    "gastric": 27000,
    "liver": 42000,
    "bladder": 83000,
    "head_neck": 66000,
    "esophageal": 21000,
}

# p53 mutation rate by cancer type (% of cases with p53 mutation)
P53_MUTATION_RATE = {
    "ovarian": 96,  # High-grade serous
    "esophageal": 83,
    "lung": 75,  # NSCLC
    "colon": 60,
    "pancreatic": 70,
    "head_neck": 85,
    "gastric": 50,
    "bladder": 50,
    "breast": 30,  # Higher in triple-negative
    "liver": 30,
}


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class PatientPopulation:
    """Estimated patient population for a specific mutation."""
    mutation: str
    mutation_frequency: float  # % of p53 mutations
    total_patients_per_year: int  # US estimate
    by_cancer_type: Dict[str, int]  # patients per cancer type
    global_estimate: int  # worldwide estimate (3x US)
    five_year_cumulative: int


@dataclass
class TherapeuticIndex:
    """Therapeutic index assessment."""
    sequence_identity: float  # % identity to WT
    mutation_count: int
    immunogenicity_score: float  # 0-1 (lower is better)
    off_target_risk: float  # 0-1
    therapeutic_window: str  # "narrow", "moderate", "wide"
    fda_pathway: str  # "gene_therapy", "cell_therapy", "protein"
    regulatory_complexity: str  # "low", "medium", "high"


@dataclass
class DeliveryFeasibility:
    """Delivery method feasibility assessment."""
    method: str
    feasibility_score: float  # 0-1
    advantages: List[str]
    challenges: List[str]
    estimated_cost_per_dose: str
    regulatory_status: str
    clinical_examples: List[str]
    identity_requirement: float  # minimum % identity for this method


@dataclass
class ImmunogenicityAssessment:
    """Immunogenicity risk assessment."""
    overall_risk: str  # "low", "moderate", "high"
    risk_score: float  # 0-1
    neoantigens_predicted: int
    hla_binding_peptides: int
    risk_factors: List[str]
    mitigation_strategies: List[str]


@dataclass
class ClinicalImpactReport:
    """Complete clinical impact report."""
    rescue_name: str
    cancer_mutation: str
    rescue_mutations: List[str]

    patient_population: PatientPopulation
    therapeutic_index: TherapeuticIndex
    delivery_options: List[DeliveryFeasibility]
    immunogenicity: ImmunogenicityAssessment

    overall_clinical_score: float  # 0-100
    clinical_viability: str  # "high", "moderate", "low", "not_viable"
    key_advantages: List[str]
    key_challenges: List[str]
    recommended_development_path: str


# ============================================================================
# PATIENT STRATIFICATION
# ============================================================================

class PatientStratifier:
    """Estimate patient populations for p53 mutations."""

    def __init__(self):
        self.mutation_data = TCGA_P53_MUTATIONS
        self.incidence = CANCER_INCIDENCE
        self.p53_rate = P53_MUTATION_RATE

    def estimate_population(self, mutation: str) -> PatientPopulation:
        """Estimate patient population for a specific mutation."""

        if mutation not in self.mutation_data:
            # Unknown mutation - estimate conservatively
            return self._estimate_unknown_mutation(mutation)

        mut_info = self.mutation_data[mutation]
        frequency = mut_info["frequency"]
        cancer_types = mut_info["cancer_types"]

        by_cancer = {}
        total = 0

        for cancer in cancer_types:
            if cancer in self.incidence and cancer in self.p53_rate:
                # Patients with this cancer AND p53 mutation AND this specific variant
                cases = self.incidence[cancer]
                p53_cases = cases * (self.p53_rate[cancer] / 100)
                mutation_cases = p53_cases * (frequency / 100)
                by_cancer[cancer] = int(mutation_cases)
                total += int(mutation_cases)

        return PatientPopulation(
            mutation=mutation,
            mutation_frequency=frequency,
            total_patients_per_year=total,
            by_cancer_type=by_cancer,
            global_estimate=total * 3,  # Rough global multiplier
            five_year_cumulative=total * 5
        )

    def _estimate_unknown_mutation(self, mutation: str) -> PatientPopulation:
        """Estimate for unknown mutations."""
        # Assume rare mutation (~0.1% of p53 mutations)
        frequency = 0.1

        # Estimate across common p53-mutated cancers
        total = 0
        by_cancer = {}

        for cancer in ["breast", "colon", "lung", "ovarian"]:
            cases = self.incidence.get(cancer, 0)
            p53_cases = cases * (self.p53_rate.get(cancer, 50) / 100)
            mutation_cases = p53_cases * (frequency / 100)
            by_cancer[cancer] = int(mutation_cases)
            total += int(mutation_cases)

        return PatientPopulation(
            mutation=mutation,
            mutation_frequency=frequency,
            total_patients_per_year=total,
            by_cancer_type=by_cancer,
            global_estimate=total * 3,
            five_year_cumulative=total * 5
        )

    def rank_mutations_by_impact(self, top_n: int = 20) -> List[Tuple[str, int]]:
        """Rank mutations by patient impact."""
        populations = []
        for mutation in self.mutation_data:
            pop = self.estimate_population(mutation)
            populations.append((mutation, pop.total_patients_per_year))

        populations.sort(key=lambda x: x[1], reverse=True)
        return populations[:top_n]

    def get_addressable_market(self, mutations: List[str]) -> Dict:
        """Calculate total addressable market for a set of mutations."""
        total_patients = 0
        by_cancer = defaultdict(int)

        for mutation in mutations:
            pop = self.estimate_population(mutation)
            total_patients += pop.total_patients_per_year
            for cancer, count in pop.by_cancer_type.items():
                by_cancer[cancer] += count

        return {
            "mutations_targeted": mutations,
            "total_patients_per_year": total_patients,
            "global_estimate": total_patients * 3,
            "by_cancer_type": dict(by_cancer),
            "five_year_cumulative": total_patients * 5,
            "market_potential": self._estimate_market_value(total_patients)
        }

    def _estimate_market_value(self, patients: int) -> str:
        """Rough market value estimate."""
        # Assume $100-500k per patient for gene therapy
        low = patients * 100000
        high = patients * 500000

        if high > 1e9:
            return f"${low/1e9:.1f}B - ${high/1e9:.1f}B"
        else:
            return f"${low/1e6:.0f}M - ${high/1e6:.0f}M"


# ============================================================================
# THERAPEUTIC INDEX CALCULATOR
# ============================================================================

class TherapeuticIndexCalculator:
    """Calculate therapeutic index metrics."""

    def __init__(self):
        self.fda_identity_thresholds = {
            "gene_therapy": 90,  # Lower threshold acceptable
            "mrna_therapy": 92,
            "protein_therapy": 95,  # Highest requirement
            "cell_therapy": 90,
        }

    def calculate(self, wt_sequence: str, rescue_sequence: str,
                 mutations: List[str]) -> TherapeuticIndex:
        """Calculate therapeutic index for a rescue design."""

        # Sequence identity
        identity = self._calculate_identity(wt_sequence, rescue_sequence)

        # Mutation count
        mutation_count = len(mutations)

        # Immunogenicity score
        immuno_score = self._estimate_immunogenicity(mutations, identity)

        # Off-target risk
        off_target = self._estimate_off_target_risk(mutations)

        # Determine therapeutic window
        if identity >= 98 and immuno_score < 0.2:
            window = "wide"
        elif identity >= 95 and immuno_score < 0.4:
            window = "moderate"
        else:
            window = "narrow"

        # Determine best FDA pathway
        fda_pathway = self._recommend_fda_pathway(identity, mutation_count)

        # Regulatory complexity
        if mutation_count <= 3 and identity >= 97:
            complexity = "low"
        elif mutation_count <= 6 and identity >= 94:
            complexity = "medium"
        else:
            complexity = "high"

        return TherapeuticIndex(
            sequence_identity=identity,
            mutation_count=mutation_count,
            immunogenicity_score=immuno_score,
            off_target_risk=off_target,
            therapeutic_window=window,
            fda_pathway=fda_pathway,
            regulatory_complexity=complexity
        )

    def _calculate_identity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence identity percentage."""
        if len(seq1) != len(seq2):
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return (matches / len(seq1)) * 100

    def _estimate_immunogenicity(self, mutations: List[str],
                                identity: float) -> float:
        """Estimate immunogenicity risk (0-1, lower is better)."""
        base_risk = 0.1

        # Add risk per mutation
        mutation_risk = len(mutations) * 0.05

        # Identity factor
        identity_risk = max(0, (100 - identity) * 0.02)

        # Check for known immunogenic positions
        immunogenic_positions = [17, 18, 19, 20, 21, 22, 23, 24, 25]  # MDM2 binding / exposed
        for mut in mutations:
            pos = int(mut[1:-1])
            if pos in immunogenic_positions:
                mutation_risk += 0.1

        total_risk = base_risk + mutation_risk + identity_risk
        return min(1.0, total_risk)

    def _estimate_off_target_risk(self, mutations: List[str]) -> float:
        """Estimate off-target effects risk."""
        # Based on how close mutations are to critical functional sites
        critical_positions = {
            # Zinc coordination
            176: 0.3, 179: 0.3, 238: 0.3, 242: 0.3,
            # DNA binding
            248: 0.2, 273: 0.2, 280: 0.2,
            # Tetramerization
            341: 0.2, 344: 0.2, 348: 0.2,
        }

        risk = 0.05  # Base risk
        for mut in mutations:
            pos = int(mut[1:-1])
            if pos in critical_positions:
                risk += critical_positions[pos]

        return min(1.0, risk)

    def _recommend_fda_pathway(self, identity: float,
                              mutation_count: int) -> str:
        """Recommend FDA regulatory pathway."""
        if identity >= 97 and mutation_count <= 3:
            return "gene_therapy"  # AAV delivery feasible
        elif identity >= 94:
            return "mrna_therapy"  # mRNA with lipid nanoparticles
        elif identity >= 90:
            return "cell_therapy"  # Ex vivo modification
        else:
            return "protein_therapy"  # Purified protein (hardest)


# ============================================================================
# DELIVERY FEASIBILITY ASSESSOR
# ============================================================================

class DeliveryFeasibilityAssessor:
    """Assess feasibility of different delivery methods."""

    def __init__(self):
        self.delivery_methods = self._define_delivery_methods()

    def _define_delivery_methods(self) -> List[Dict]:
        """Define available delivery methods."""
        return [
            {
                "method": "AAV Gene Therapy",
                "description": "Adeno-associated virus delivery of p53 gene",
                "identity_requirement": 90,
                "advantages": [
                    "Long-lasting expression",
                    "Established clinical precedent (Luxturna, Zolgensma)",
                    "Tissue-specific targeting possible",
                    "Single dose potential"
                ],
                "challenges": [
                    "Pre-existing immunity (~30% population)",
                    "Limited cargo capacity (4.7kb)",
                    "Manufacturing complexity",
                    "High cost per dose"
                ],
                "cost": "$100k-$500k per dose",
                "regulatory_status": "Multiple approved products",
                "examples": ["Luxturna", "Zolgensma", "Hemgenix"]
            },
            {
                "method": "Lentiviral Gene Therapy",
                "description": "Lentiviral vector for stable integration",
                "identity_requirement": 90,
                "advantages": [
                    "Larger cargo capacity",
                    "Stable integration",
                    "Works in dividing cells",
                    "Ex vivo modification possible"
                ],
                "challenges": [
                    "Insertional mutagenesis risk",
                    "Ex vivo manufacturing required",
                    "Cannot target all tissues"
                ],
                "cost": "$300k-$1M per treatment",
                "regulatory_status": "Approved for blood disorders",
                "examples": ["Kymriah", "Yescarta", "Breyanzi"]
            },
            {
                "method": "mRNA-LNP Therapy",
                "description": "mRNA encapsulated in lipid nanoparticles",
                "identity_requirement": 92,
                "advantages": [
                    "Transient expression (safety)",
                    "Rapid manufacturing",
                    "Lower cost than viral vectors",
                    "Repeat dosing possible"
                ],
                "challenges": [
                    "Short duration of expression",
                    "Requires repeated dosing",
                    "Liver tropism of LNPs",
                    "Cold chain requirements"
                ],
                "cost": "$5k-$50k per dose course",
                "regulatory_status": "COVID vaccines approved",
                "examples": ["Comirnaty", "Spikevax"]
            },
            {
                "method": "Purified Protein Therapy",
                "description": "Direct delivery of purified p53 protein",
                "identity_requirement": 95,
                "advantages": [
                    "Precise dosing",
                    "Immediate effect",
                    "Well-understood PK/PD"
                ],
                "challenges": [
                    "Cell penetration barrier",
                    "Short half-life",
                    "Immunogenicity concerns",
                    "Manufacturing complexity"
                ],
                "cost": "$10k-$100k per course",
                "regulatory_status": "Challenging for transcription factors",
                "examples": ["Insulin (proof of concept)"]
            },
            {
                "method": "CRISPR Prime Editing",
                "description": "Direct genome editing to introduce rescue mutations",
                "identity_requirement": 100,  # Creates WT-like sequence
                "advantages": [
                    "Permanent correction",
                    "Creates near-WT sequence",
                    "No foreign protein expression",
                    "Precise edits possible"
                ],
                "challenges": [
                    "Delivery to solid tumors difficult",
                    "Off-target editing risk",
                    "Low efficiency in vivo",
                    "Regulatory uncertainty"
                ],
                "cost": "Uncertain ($100k+)",
                "regulatory_status": "Early clinical trials",
                "examples": ["EDIT-101 (Editas)", "CTX001 (CRISPR Therapeutics)"]
            },
            {
                "method": "Oncolytic Virus + p53",
                "description": "Armed oncolytic virus expressing p53",
                "identity_requirement": 90,
                "advantages": [
                    "Tumor selectivity",
                    "Immune activation",
                    "Bystander effect",
                    "Clinical precedent (IMLYGIC)"
                ],
                "challenges": [
                    "Systemic delivery limited",
                    "Anti-vector immunity",
                    "Variable tumor penetration"
                ],
                "cost": "$50k-$200k per course",
                "regulatory_status": "IMLYGIC approved (melanoma)",
                "examples": ["IMLYGIC", "ONYX-015 (historical)"]
            }
        ]

    def assess_all_methods(self, sequence_identity: float,
                          mutation_count: int) -> List[DeliveryFeasibility]:
        """Assess feasibility of all delivery methods."""
        assessments = []

        for method_info in self.delivery_methods:
            score = self._calculate_feasibility(
                method_info, sequence_identity, mutation_count
            )

            assessment = DeliveryFeasibility(
                method=method_info["method"],
                feasibility_score=score,
                advantages=method_info["advantages"],
                challenges=method_info["challenges"],
                estimated_cost_per_dose=method_info["cost"],
                regulatory_status=method_info["regulatory_status"],
                clinical_examples=method_info["examples"],
                identity_requirement=method_info["identity_requirement"]
            )
            assessments.append(assessment)

        # Sort by feasibility score
        assessments.sort(key=lambda x: x.feasibility_score, reverse=True)
        return assessments

    def _calculate_feasibility(self, method_info: Dict,
                              identity: float,
                              mutation_count: int) -> float:
        """Calculate feasibility score for a delivery method."""
        score = 0.5  # Base score

        # Identity requirement check
        req_identity = method_info["identity_requirement"]
        if identity >= req_identity:
            score += 0.3
        elif identity >= req_identity - 3:
            score += 0.1
        else:
            score -= 0.3

        # Mutation count factor (fewer = better)
        if mutation_count <= 3:
            score += 0.1
        elif mutation_count <= 5:
            score += 0.05
        elif mutation_count > 10:
            score -= 0.1

        # Regulatory status bonus
        if "approved" in method_info["regulatory_status"].lower():
            score += 0.1

        return max(0, min(1, score))


# ============================================================================
# IMMUNOGENICITY ASSESSOR
# ============================================================================

class ImmunogenicityAssessor:
    """Assess immunogenicity risk of rescue designs."""

    def __init__(self):
        # Common HLA alleles for neoantigen prediction (simplified)
        self.common_hlas = ["A*02:01", "A*01:01", "B*07:02", "B*08:01"]

    def assess(self, wt_sequence: str, rescue_sequence: str,
              mutations: List[str]) -> ImmunogenicityAssessment:
        """Assess immunogenicity risk."""

        # Count predicted neoantigens
        neoantigens = self._predict_neoantigens(mutations, rescue_sequence)

        # Estimate HLA binding peptides
        hla_peptides = self._estimate_hla_binders(mutations, rescue_sequence)

        # Identify risk factors
        risk_factors = self._identify_risk_factors(mutations, rescue_sequence)

        # Calculate risk score
        risk_score = self._calculate_risk_score(
            len(mutations), neoantigens, hla_peptides, risk_factors
        )

        # Determine overall risk
        if risk_score < 0.3:
            overall_risk = "low"
        elif risk_score < 0.6:
            overall_risk = "moderate"
        else:
            overall_risk = "high"

        # Mitigation strategies
        mitigations = self._suggest_mitigations(overall_risk, risk_factors)

        return ImmunogenicityAssessment(
            overall_risk=overall_risk,
            risk_score=risk_score,
            neoantigens_predicted=neoantigens,
            hla_binding_peptides=hla_peptides,
            risk_factors=risk_factors,
            mitigation_strategies=mitigations
        )

    def _predict_neoantigens(self, mutations: List[str],
                            sequence: str) -> int:
        """Predict number of potential neoantigens (simplified)."""
        # Each mutation creates potential neoantigen
        # Surface-exposed mutations have higher risk
        exposed_positions = set(range(1, 50)) | set(range(290, 393))  # N/C termini

        count = 0
        for mut in mutations:
            pos = int(mut[1:-1])
            if pos in exposed_positions:
                count += 2  # Higher risk
            else:
                count += 1

        return count

    def _estimate_hla_binders(self, mutations: List[str],
                             sequence: str) -> int:
        """Estimate HLA-binding peptides from mutations."""
        # Simplified: each mutation can generate ~4 potential 9-mer peptides
        # Real implementation would use NetMHC or similar
        return len(mutations) * 4

    def _identify_risk_factors(self, mutations: List[str],
                              sequence: str) -> List[str]:
        """Identify specific immunogenicity risk factors."""
        risks = []

        # High mutation count
        if len(mutations) > 5:
            risks.append(f"High mutation count ({len(mutations)} mutations)")

        # Check for immunogenic amino acid changes
        immunogenic_changes = {'A': 'K', 'A': 'R', 'A': 'E', 'A': 'D'}  # Non-self-like
        for mut in mutations:
            new_aa = mut[-1]
            if new_aa in 'KRDE':  # Charged residues often immunogenic
                risks.append(f"Charged residue introduction at {mut}")

        # Surface-exposed mutations
        exposed = set(range(1, 50)) | set(range(290, 393))
        for mut in mutations:
            pos = int(mut[1:-1])
            if pos in exposed:
                risks.append(f"Surface-exposed mutation {mut}")

        return risks[:5]  # Top 5 risks

    def _calculate_risk_score(self, mutation_count: int,
                             neoantigens: int, hla_peptides: int,
                             risk_factors: List[str]) -> float:
        """Calculate overall immunogenicity risk score."""
        score = 0.1  # Base risk

        # Mutation count contribution
        score += mutation_count * 0.05

        # Neoantigen contribution
        score += neoantigens * 0.03

        # Risk factor contribution
        score += len(risk_factors) * 0.05

        return min(1.0, score)

    def _suggest_mitigations(self, risk_level: str,
                            risk_factors: List[str]) -> List[str]:
        """Suggest immunogenicity mitigation strategies."""
        mitigations = []

        mitigations.append("Use immunosuppression during initial dosing")
        mitigations.append("Consider tolerization protocol before treatment")

        if risk_level in ["moderate", "high"]:
            mitigations.append("Prioritize mRNA delivery for transient expression")
            mitigations.append("Monitor for anti-drug antibodies")
            mitigations.append("Consider dose titration strategy")

        if "Surface-exposed" in str(risk_factors):
            mitigations.append("Evaluate alternative mutations at buried positions")

        return mitigations


# ============================================================================
# CLINICAL IMPACT ENGINE
# ============================================================================

class ClinicalImpactEngine:
    """
    Main engine for clinical impact assessment.
    """

    def __init__(self):
        self.patient_stratifier = PatientStratifier()
        self.ti_calculator = TherapeuticIndexCalculator()
        self.delivery_assessor = DeliveryFeasibilityAssessor()
        self.immuno_assessor = ImmunogenicityAssessor()

    def generate_report(self, name: str, wt_sequence: str,
                       rescue_sequence: str, cancer_mutation: str,
                       rescue_mutations: List[str]) -> ClinicalImpactReport:
        """Generate comprehensive clinical impact report."""

        # Patient population
        patient_pop = self.patient_stratifier.estimate_population(cancer_mutation)

        # Therapeutic index
        ti = self.ti_calculator.calculate(wt_sequence, rescue_sequence, rescue_mutations)

        # Delivery options
        delivery_options = self.delivery_assessor.assess_all_methods(
            ti.sequence_identity, ti.mutation_count
        )

        # Immunogenicity
        immuno = self.immuno_assessor.assess(wt_sequence, rescue_sequence, rescue_mutations)

        # Calculate overall clinical score
        clinical_score = self._calculate_clinical_score(
            patient_pop, ti, delivery_options, immuno
        )

        # Determine viability
        if clinical_score >= 70:
            viability = "high"
        elif clinical_score >= 50:
            viability = "moderate"
        elif clinical_score >= 30:
            viability = "low"
        else:
            viability = "not_viable"

        # Key advantages and challenges
        advantages = self._identify_advantages(patient_pop, ti, delivery_options)
        challenges = self._identify_challenges(ti, immuno, delivery_options)

        # Recommended path
        recommended_path = self._recommend_development_path(
            ti, delivery_options, patient_pop
        )

        return ClinicalImpactReport(
            rescue_name=name,
            cancer_mutation=cancer_mutation,
            rescue_mutations=rescue_mutations,
            patient_population=patient_pop,
            therapeutic_index=ti,
            delivery_options=delivery_options,
            immunogenicity=immuno,
            overall_clinical_score=clinical_score,
            clinical_viability=viability,
            key_advantages=advantages,
            key_challenges=challenges,
            recommended_development_path=recommended_path
        )

    def _calculate_clinical_score(self, patient_pop: PatientPopulation,
                                 ti: TherapeuticIndex,
                                 delivery_options: List[DeliveryFeasibility],
                                 immuno: ImmunogenicityAssessment) -> float:
        """Calculate overall clinical score (0-100)."""
        score = 0

        # Patient population (25 points max)
        if patient_pop.total_patients_per_year >= 10000:
            score += 25
        elif patient_pop.total_patients_per_year >= 5000:
            score += 20
        elif patient_pop.total_patients_per_year >= 1000:
            score += 15
        elif patient_pop.total_patients_per_year >= 100:
            score += 10
        else:
            score += 5

        # Therapeutic index (25 points max)
        if ti.sequence_identity >= 98:
            score += 15
        elif ti.sequence_identity >= 95:
            score += 10
        elif ti.sequence_identity >= 92:
            score += 5

        if ti.therapeutic_window == "wide":
            score += 10
        elif ti.therapeutic_window == "moderate":
            score += 5

        # Delivery feasibility (25 points max)
        best_delivery = delivery_options[0] if delivery_options else None
        if best_delivery:
            score += best_delivery.feasibility_score * 25

        # Immunogenicity (25 points max, inverted)
        immuno_penalty = immuno.risk_score * 25
        score += (25 - immuno_penalty)

        return min(100, max(0, score))

    def _identify_advantages(self, patient_pop: PatientPopulation,
                            ti: TherapeuticIndex,
                            delivery_options: List[DeliveryFeasibility]) -> List[str]:
        """Identify key advantages."""
        advantages = []

        if patient_pop.total_patients_per_year >= 5000:
            advantages.append(f"Large addressable patient population ({patient_pop.total_patients_per_year:,} patients/year in US)")

        if ti.sequence_identity >= 97:
            advantages.append(f"High sequence identity ({ti.sequence_identity:.1f}%) minimizes immunogenicity")

        if ti.mutation_count <= 3:
            advantages.append(f"Minimal mutations ({ti.mutation_count}) for regulatory simplicity")

        if delivery_options and delivery_options[0].feasibility_score >= 0.7:
            advantages.append(f"Favorable delivery via {delivery_options[0].method}")

        if ti.therapeutic_window == "wide":
            advantages.append("Wide therapeutic window provides safety margin")

        return advantages[:5]

    def _identify_challenges(self, ti: TherapeuticIndex,
                            immuno: ImmunogenicityAssessment,
                            delivery_options: List[DeliveryFeasibility]) -> List[str]:
        """Identify key challenges."""
        challenges = []

        if ti.sequence_identity < 95:
            challenges.append(f"Lower sequence identity ({ti.sequence_identity:.1f}%) may increase immunogenicity")

        if ti.mutation_count > 5:
            challenges.append(f"Multiple mutations ({ti.mutation_count}) increase regulatory complexity")

        if immuno.overall_risk in ["moderate", "high"]:
            challenges.append(f"{immuno.overall_risk.capitalize()} immunogenicity risk requires monitoring")

        if ti.regulatory_complexity == "high":
            challenges.append("Complex regulatory pathway expected")

        if delivery_options:
            challenges.extend(delivery_options[0].challenges[:2])

        return challenges[:5]

    def _recommend_development_path(self, ti: TherapeuticIndex,
                                   delivery_options: List[DeliveryFeasibility],
                                   patient_pop: PatientPopulation) -> str:
        """Recommend development path."""
        if ti.sequence_identity >= 97 and ti.mutation_count <= 3:
            best_delivery = delivery_options[0].method if delivery_options else "AAV"
            return f"Proceed with {best_delivery} delivery. Fast-track potential for orphan indication. Target initial trial in {patient_pop.by_cancer_type.get(list(patient_pop.by_cancer_type.keys())[0] if patient_pop.by_cancer_type else 'solid tumors', 'solid tumors')} patients."

        elif ti.sequence_identity >= 94:
            return "Consider mRNA-LNP delivery for transient expression. Requires repeat dosing but lower immunogenicity risk. Phase 1 dose-escalation recommended."

        elif ti.sequence_identity >= 90:
            return "Ex vivo cell therapy approach recommended. Modify patient cells, validate correction, then reinfuse. Higher complexity but better safety profile."

        else:
            return "Consider alternative rescue designs with fewer mutations. Current design may face significant immunogenicity and regulatory challenges."


# ============================================================================
# EXPORT FUNCTIONS
# ============================================================================

def export_clinical_report(report: ClinicalImpactReport, output_path: str) -> str:
    """Export clinical impact report to markdown."""

    lines = []
    lines.append(f"# Clinical Impact Report: {report.rescue_name}")
    lines.append(f"\nGenerated for ISEF 2026")

    lines.append(f"\n## Summary")
    lines.append(f"- **Cancer Mutation:** {report.cancer_mutation}")
    lines.append(f"- **Rescue Mutations:** {', '.join(report.rescue_mutations)}")
    lines.append(f"- **Clinical Viability:** {report.clinical_viability.upper()}")
    lines.append(f"- **Overall Score:** {report.overall_clinical_score:.0f}/100")

    lines.append(f"\n## Patient Population")
    lines.append(f"- **Annual US Patients:** {report.patient_population.total_patients_per_year:,}")
    lines.append(f"- **Global Estimate:** {report.patient_population.global_estimate:,}")
    lines.append(f"- **5-Year Cumulative:** {report.patient_population.five_year_cumulative:,}")
    lines.append(f"\n### By Cancer Type")
    for cancer, count in report.patient_population.by_cancer_type.items():
        lines.append(f"  - {cancer.title()}: {count:,}")

    lines.append(f"\n## Therapeutic Index")
    lines.append(f"- **Sequence Identity:** {report.therapeutic_index.sequence_identity:.1f}%")
    lines.append(f"- **Mutation Count:** {report.therapeutic_index.mutation_count}")
    lines.append(f"- **Therapeutic Window:** {report.therapeutic_index.therapeutic_window}")
    lines.append(f"- **Recommended FDA Pathway:** {report.therapeutic_index.fda_pathway}")
    lines.append(f"- **Regulatory Complexity:** {report.therapeutic_index.regulatory_complexity}")

    lines.append(f"\n## Delivery Options")
    for i, delivery in enumerate(report.delivery_options[:3], 1):
        lines.append(f"\n### {i}. {delivery.method}")
        lines.append(f"- **Feasibility Score:** {delivery.feasibility_score:.2f}")
        lines.append(f"- **Cost:** {delivery.estimated_cost_per_dose}")
        lines.append(f"- **Regulatory Status:** {delivery.regulatory_status}")
        lines.append(f"- **Advantages:** {', '.join(delivery.advantages[:2])}")
        lines.append(f"- **Challenges:** {', '.join(delivery.challenges[:2])}")

    lines.append(f"\n## Immunogenicity Assessment")
    lines.append(f"- **Overall Risk:** {report.immunogenicity.overall_risk.upper()}")
    lines.append(f"- **Risk Score:** {report.immunogenicity.risk_score:.2f}")
    lines.append(f"- **Predicted Neoantigens:** {report.immunogenicity.neoantigens_predicted}")
    if report.immunogenicity.risk_factors:
        lines.append(f"- **Risk Factors:**")
        for rf in report.immunogenicity.risk_factors[:3]:
            lines.append(f"  - {rf}")
    if report.immunogenicity.mitigation_strategies:
        lines.append(f"- **Mitigation Strategies:**")
        for ms in report.immunogenicity.mitigation_strategies[:3]:
            lines.append(f"  - {ms}")

    lines.append(f"\n## Key Advantages")
    for adv in report.key_advantages:
        lines.append(f"- {adv}")

    lines.append(f"\n## Key Challenges")
    for chal in report.key_challenges:
        lines.append(f"- {chal}")

    lines.append(f"\n## Recommended Development Path")
    lines.append(f"\n{report.recommended_development_path}")

    content = "\n".join(lines)

    with open(output_path, 'w') as f:
        f.write(content)

    return output_path


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for clinical impact assessment."""
    import argparse

    parser = argparse.ArgumentParser(description="p53-proteoMgCAD Clinical Impact")
    parser.add_argument("--cancer", type=str, required=True,
                       help="Cancer mutation (e.g., R175H)")
    parser.add_argument("--rescue", type=str, nargs='+', required=True,
                       help="Rescue mutations (e.g., N268D T125R)")
    parser.add_argument("--output", type=str, default="clinical_impact_report.md",
                       help="Output file path")

    args = parser.parse_args()

    # Load WT sequence
    from ..data.dms import P53_WT
    wt_sequence = P53_WT

    # Apply rescue mutations
    rescue_sequence = wt_sequence
    for mut in args.rescue:
        pos = int(mut[1:-1]) - 1
        new_aa = mut[-1]
        rescue_sequence = rescue_sequence[:pos] + new_aa + rescue_sequence[pos+1:]

    # Generate report
    engine = ClinicalImpactEngine()
    report = engine.generate_report(
        name=f"p53_{args.cancer}_rescue",
        wt_sequence=wt_sequence,
        rescue_sequence=rescue_sequence,
        cancer_mutation=args.cancer,
        rescue_mutations=args.rescue
    )

    # Export
    export_clinical_report(report, args.output)

    # Print summary
    print("\n" + "="*60)
    print("CLINICAL IMPACT ASSESSMENT")
    print("="*60)
    print(f"\nCancer Mutation: {args.cancer}")
    print(f"Rescue Mutations: {', '.join(args.rescue)}")
    print(f"\nClinical Viability: {report.clinical_viability.upper()}")
    print(f"Overall Score: {report.overall_clinical_score:.0f}/100")
    print(f"Patient Population: {report.patient_population.total_patients_per_year:,} patients/year (US)")
    print(f"\nFull report saved to: {args.output}")


if __name__ == "__main__":
    main()
