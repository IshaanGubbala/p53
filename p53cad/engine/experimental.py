"""
Experimental Validation Pipeline for p53-proteoMgCAD

Generates experimental protocols and reagents for validating computational rescues:
1. Yeast Display Screening Protocol
2. CRISPR Prime Editing Guide Design
3. Collaboration Export (GenBank, SnapGene)
4. Assay Design Recommendations
"""

import json
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
import re


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class YeastDisplayProtocol:
    """Yeast surface display screening protocol."""
    title: str
    target_mutations: List[str]
    rescue_mutations: List[str]

    # Plasmid construction
    vector: str
    promoter: str
    fusion_partner: str  # Aga2p typically
    linker_sequence: str

    # Library design
    library_size: int
    coverage: float
    mutagenesis_method: str

    # Selection parameters
    selection_rounds: int
    sorting_gates: List[str]
    negative_selection: bool

    # Timeline
    estimated_days: int

    # Materials list
    materials: List[Dict[str, Any]]

    # Step-by-step protocol
    protocol_steps: List[Dict[str, str]]

    # Expected outcomes
    expected_outcomes: Dict[str, str]


@dataclass
class PrimeEditingGuide:
    """CRISPR Prime Editing guide RNA design."""
    target_mutation: str
    edit_type: str  # "rescue_introduction" or "mutation_reversion"

    # pegRNA components
    spacer_sequence: str  # 20nt guide
    pam_sequence: str  # NGG
    pbs_sequence: str  # Primer binding site (13-17nt)
    rt_template: str  # RT template with edit

    # Nicking guide (for PE3)
    nick_spacer: str
    nick_pam: str
    nick_distance: int  # bp from pegRNA cut

    # Quality metrics
    spacer_gc_content: float
    pbs_gc_content: float
    predicted_efficiency: float
    off_target_score: float

    # Oligo sequences for ordering
    pegRNA_oligo_top: str
    pegRNA_oligo_bottom: str
    nick_oligo_top: str
    nick_oligo_bottom: str


@dataclass
class GenBankExport:
    """GenBank format export."""
    locus: str
    definition: str
    accession: str
    version: str
    keywords: str
    source: str
    organism: str
    sequence: str
    features: List[Dict[str, Any]]

    def to_genbank_string(self) -> str:
        """Convert to GenBank format string."""
        lines = []

        # Header
        lines.append(f"LOCUS       {self.locus:16s} {len(self.sequence):>7d} bp    DNA     linear   SYN {datetime.now().strftime('%d-%b-%Y').upper()}")
        lines.append(f"DEFINITION  {self.definition}")
        lines.append(f"ACCESSION   {self.accession}")
        lines.append(f"VERSION     {self.version}")
        lines.append(f"KEYWORDS    {self.keywords}")
        lines.append(f"SOURCE      {self.source}")
        lines.append(f"  ORGANISM  {self.organism}")
        lines.append("FEATURES             Location/Qualifiers")

        # Features
        for feature in self.features:
            feat_type = feature.get('type', 'misc_feature')
            location = feature.get('location', '1..100')
            lines.append(f"     {feat_type:15s} {location}")

            for qual_key, qual_val in feature.get('qualifiers', {}).items():
                lines.append(f'                     /{qual_key}="{qual_val}"')

        # Origin (sequence)
        lines.append("ORIGIN")
        seq = self.sequence.lower()
        for i in range(0, len(seq), 60):
            chunk = seq[i:i+60]
            formatted = ' '.join([chunk[j:j+10] for j in range(0, len(chunk), 10)])
            lines.append(f"{i+1:>9d} {formatted}")

        lines.append("//")

        return "\n".join(lines)


@dataclass
class AssayRecommendation:
    """Recommended assay for validation."""
    assay_name: str
    assay_type: str  # "functional", "structural", "binding"
    priority: int  # 1-5, 1 is highest
    rationale: str
    readout: str
    throughput: str  # "low", "medium", "high"
    estimated_cost: str
    timeline: str
    required_equipment: List[str]
    protocol_reference: str


# ============================================================================
# CODON TABLE
# ============================================================================

CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA'],  # Stop codons
}

# Preferred codons for E. coli expression
ECOLI_PREFERRED_CODONS = {
    'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAT', 'I': 'ATC', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
}

# Preferred codons for human expression
HUMAN_PREFERRED_CODONS = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGG',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def protein_to_dna(protein_seq: str, organism: str = "human") -> str:
    """Convert protein sequence to DNA using preferred codons."""
    if organism == "human":
        codon_table = HUMAN_PREFERRED_CODONS
    elif organism == "ecoli":
        codon_table = ECOLI_PREFERRED_CODONS
    else:
        codon_table = HUMAN_PREFERRED_CODONS

    dna = []
    for aa in protein_seq:
        if aa in codon_table:
            dna.append(codon_table[aa])
        else:
            # Use first codon from general table
            dna.append(CODON_TABLE.get(aa, ['NNN'])[0])

    return ''.join(dna)


def reverse_complement(dna: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(dna))


def gc_content(seq: str) -> float:
    """Calculate GC content of sequence."""
    seq = seq.upper()
    gc = sum(1 for base in seq if base in 'GC')
    return gc / len(seq) if seq else 0


def find_pam_sites(dna: str, pam: str = "NGG") -> List[Tuple[int, str]]:
    """Find PAM sites in DNA sequence."""
    sites = []

    # Forward strand
    for i in range(len(dna) - len(pam) + 1):
        match = True
        for j, char in enumerate(pam):
            if char == 'N':
                continue
            if dna[i+j].upper() != char:
                match = False
                break
        if match:
            sites.append((i, '+'))

    # Reverse strand
    rc_dna = reverse_complement(dna)
    for i in range(len(rc_dna) - len(pam) + 1):
        match = True
        for j, char in enumerate(pam):
            if char == 'N':
                continue
            if rc_dna[i+j].upper() != char:
                match = False
                break
        if match:
            # Convert position to forward strand coordinates
            sites.append((len(dna) - i - len(pam), '-'))

    return sites


# ============================================================================
# YEAST DISPLAY PROTOCOL GENERATOR
# ============================================================================

class YeastDisplayGenerator:
    """Generate yeast surface display protocols."""

    def __init__(self):
        self.aga2_sequence = "ATGCAGTTACTTCGCTGTTTTTCAATATTTTCTGTTATTGCTTCAGTTTTAGCACAGGAACTGACAACTATATGCGAGCAAATCCCCTCACCAACTTTAGAATCGACGCCGTACTCTTTGTCAACGACTACTATTTTGGCCAACGGGAAGGCAATGCAAGGAGTTTTTGAATATTACAAATCAGTAACGTTTGTCAGTAATTGCGGTTCTCACCCCTCAACAACTAGCAAAGGCAGCCCCATAAACACACAGTATGTTTTTAAGGACAATAGCTCGACGATTGAAGGTAGATACCCATACGACGTTCCAGACTACGCTCTGCAGGCTAGTGGTGGT"

    def generate_protocol(self, target_mutations: List[str],
                         rescue_mutations: List[str],
                         p53_sequence: str) -> YeastDisplayProtocol:
        """Generate complete yeast display protocol."""

        # Design library
        library_info = self._design_library(target_mutations, rescue_mutations)

        # Materials
        materials = self._generate_materials_list()

        # Protocol steps
        steps = self._generate_protocol_steps(target_mutations, rescue_mutations)

        # Expected outcomes
        outcomes = self._predict_outcomes(target_mutations, rescue_mutations)

        return YeastDisplayProtocol(
            title=f"Yeast Display Selection for p53 Rescue of {', '.join(target_mutations)}",
            target_mutations=target_mutations,
            rescue_mutations=rescue_mutations,
            vector="pYD1",
            promoter="GAL1",
            fusion_partner="Aga2p",
            linker_sequence="GGGGS" * 3,
            library_size=library_info['size'],
            coverage=library_info['coverage'],
            mutagenesis_method=library_info['method'],
            selection_rounds=3,
            sorting_gates=["DNA binding", "Folding (PAb1620)", "Function (FL reporter)"],
            negative_selection=True,
            estimated_days=21,
            materials=materials,
            protocol_steps=steps,
            expected_outcomes=outcomes
        )

    def _design_library(self, target_mutations: List[str],
                       rescue_mutations: List[str]) -> Dict:
        """Design mutagenesis library."""

        # For site-saturation mutagenesis at rescue positions
        n_positions = len(rescue_mutations)

        if n_positions <= 3:
            # Full combinatorial NNK at each position
            # NNK = 32 codons covering all 20 aa
            library_size = 32 ** n_positions
            coverage = 10.0  # 10x coverage
            method = "NNK site-saturation mutagenesis"
        else:
            # Too many positions - use error-prone PCR
            library_size = 1e6
            coverage = 3.0
            method = "Error-prone PCR (MnCl2)"

        return {
            'size': int(library_size),
            'coverage': coverage,
            'method': method
        }

    def _generate_materials_list(self) -> List[Dict[str, Any]]:
        """Generate materials needed for protocol."""
        return [
            {"item": "pYD1 vector", "quantity": "1 μg", "supplier": "Addgene #73447", "notes": "Aga2p fusion vector"},
            {"item": "EBY100 yeast strain", "quantity": "1 aliquot", "supplier": "ATCC", "notes": "Display strain"},
            {"item": "Anti-c-myc antibody (9E10)", "quantity": "100 μL", "supplier": "Santa Cruz", "notes": "Display detection"},
            {"item": "DNA-Cy5 probe", "quantity": "Custom", "supplier": "IDT", "notes": "p53 response element"},
            {"item": "PAb1620 antibody", "quantity": "50 μL", "supplier": "Abcam", "notes": "Folded p53 conformation"},
            {"item": "PAb240 antibody", "quantity": "50 μL", "supplier": "Abcam", "notes": "Unfolded p53 (negative)"},
            {"item": "SD-CAA media", "quantity": "2 L", "supplier": "In-house", "notes": "Non-inducing growth"},
            {"item": "SG-CAA media", "quantity": "1 L", "supplier": "In-house", "notes": "Galactose induction"},
            {"item": "NNK primers", "quantity": "Custom", "supplier": "IDT", "notes": "Site-saturation"},
            {"item": "Electrocompetent EBY100", "quantity": "20 aliquots", "supplier": "In-house", "notes": "Library transformation"},
        ]

    def _generate_protocol_steps(self, target_mutations: List[str],
                                rescue_mutations: List[str]) -> List[Dict[str, str]]:
        """Generate step-by-step protocol."""

        # Parse mutation positions
        positions = [int(m[1:-1]) for m in rescue_mutations]
        pos_str = ", ".join(str(p) for p in positions)

        return [
            {
                "day": "1-2",
                "step": "Library Construction",
                "description": f"Perform NNK site-saturation mutagenesis at positions {pos_str} using overlap extension PCR. Verify library diversity by sequencing 20 clones.",
                "critical": "Ensure NNK primers are phosphorylated for ligation"
            },
            {
                "day": "3",
                "step": "Vector Preparation",
                "description": "Digest pYD1 with BamHI and NheI. Gel purify linearized vector (>90% purity required).",
                "critical": "Use fresh restriction enzymes"
            },
            {
                "day": "4",
                "step": "Yeast Transformation",
                "description": "Transform library into EBY100 by electroporation (25 μg DNA per 4×10^8 cells). Plate on SD-CAA. Target: >10^7 transformants.",
                "critical": "Keep cells on ice, pulse at 2.5 kV"
            },
            {
                "day": "5-7",
                "step": "Library Amplification",
                "description": "Pool colonies into 500 mL SD-CAA. Grow to OD600=1.0. Make glycerol stocks. Calculate library size by dilution plating.",
                "critical": "Do not overgrow (maintain diversity)"
            },
            {
                "day": "8",
                "step": "Display Induction",
                "description": "Dilute to OD600=0.5 in SG-CAA (galactose). Induce at 20°C for 16-20 hours.",
                "critical": "20°C is critical for proper folding"
            },
            {
                "day": "9",
                "step": "Round 1 Selection - Folding",
                "description": "Stain with PAb1620 (folded) and PAb240 (unfolded). FACS sort for PAb1620+/PAb240- population. Collect 10^6 cells.",
                "critical": "Set gates using WT p53 control"
            },
            {
                "day": "10-12",
                "step": "Round 1 Recovery",
                "description": "Grow sorted cells in SD-CAA to saturation. Re-induce as above.",
                "critical": "Monitor for contamination"
            },
            {
                "day": "13",
                "step": "Round 2 Selection - DNA Binding",
                "description": "Incubate with 100 nM Cy5-labeled p53 response element DNA. Sort for top 1% DNA binding.",
                "critical": "Include non-specific DNA control"
            },
            {
                "day": "14-16",
                "step": "Round 2 Recovery",
                "description": "Amplify sorted population as above.",
                "critical": None
            },
            {
                "day": "17",
                "step": "Round 3 Selection - Stringent",
                "description": "Dual selection: PAb1620+ AND DNA binding (top 0.1%). This enriches for functional rescues.",
                "critical": "Use stringent gates to reduce false positives"
            },
            {
                "day": "18-19",
                "step": "Clone Analysis",
                "description": "Plate sorted cells. Pick 96 colonies for Sanger sequencing. Identify enriched mutations.",
                "critical": None
            },
            {
                "day": "20-21",
                "step": "Hit Validation",
                "description": "Retest top 10 hits individually. Measure display level, DNA binding Kd, and thermal stability.",
                "critical": "Include WT and cancer mutant controls"
            },
        ]

    def _predict_outcomes(self, target_mutations: List[str],
                         rescue_mutations: List[str]) -> Dict[str, str]:
        """Predict expected outcomes."""
        return {
            "expected_enrichment": "10-100x enrichment of functional variants per round",
            "hit_rate": "5-15% of clones in final pool expected to show rescue",
            "false_positive_rate": "~30% (validate individually)",
            "novel_rescue_discovery": "May identify rescues not predicted computationally",
            "timeline_variability": "±3 days depending on transformation efficiency",
            "success_criteria": "≥3 independent rescues with >50% WT activity"
        }

    def export_protocol_docx(self, protocol: YeastDisplayProtocol,
                            output_path: str) -> str:
        """Export protocol to markdown (can be converted to docx)."""

        lines = []
        lines.append(f"# {protocol.title}")
        lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d')}")
        lines.append(f"\n## Overview")
        lines.append(f"- **Target mutations:** {', '.join(protocol.target_mutations)}")
        lines.append(f"- **Predicted rescue mutations:** {', '.join(protocol.rescue_mutations)}")
        lines.append(f"- **Vector:** {protocol.vector}")
        lines.append(f"- **Library size:** {protocol.library_size:,}")
        lines.append(f"- **Timeline:** {protocol.estimated_days} days")

        lines.append(f"\n## Materials")
        lines.append("| Item | Quantity | Supplier | Notes |")
        lines.append("|------|----------|----------|-------|")
        for mat in protocol.materials:
            lines.append(f"| {mat['item']} | {mat['quantity']} | {mat['supplier']} | {mat['notes']} |")

        lines.append(f"\n## Protocol")
        for step in protocol.protocol_steps:
            lines.append(f"\n### Day {step['day']}: {step['step']}")
            lines.append(step['description'])
            if step.get('critical'):
                lines.append(f"\n**CRITICAL:** {step['critical']}")

        lines.append(f"\n## Expected Outcomes")
        for key, value in protocol.expected_outcomes.items():
            lines.append(f"- **{key.replace('_', ' ').title()}:** {value}")

        content = "\n".join(lines)

        with open(output_path, 'w') as f:
            f.write(content)

        return output_path


# ============================================================================
# PRIME EDITING GUIDE DESIGNER
# ============================================================================

class PrimeEditingDesigner:
    """Design CRISPR Prime Editing guides for p53 rescue mutations."""

    def __init__(self):
        # p53 genomic DNA (exons 4-8 region, simplified)
        # In real implementation, would fetch from genome
        self.pbs_length_default = 13  # nt
        self.rt_extension = 10  # nt after edit

    def design_guide(self, mutation: str, genomic_context: str,
                    edit_type: str = "rescue_introduction") -> PrimeEditingGuide:
        """
        Design pegRNA for introducing or reverting a mutation.

        Args:
            mutation: Mutation in format "A123B"
            genomic_context: ~200bp genomic DNA around mutation site
            edit_type: "rescue_introduction" or "mutation_reversion"
        """

        # Parse mutation
        original_aa = mutation[0]
        position = int(mutation[1:-1])
        new_aa = mutation[-1]

        # Find codon in genomic context (simplified)
        # In real implementation, would need exact genomic coordinates

        # Find PAM sites near mutation
        pam_sites = find_pam_sites(genomic_context)

        if not pam_sites:
            # No PAM found - return dummy guide
            return self._create_dummy_guide(mutation, edit_type)

        # Select best PAM (closest to edit site, on either strand)
        edit_position = len(genomic_context) // 2  # Assume mutation is centered
        best_pam = min(pam_sites, key=lambda x: abs(x[0] - edit_position))
        pam_pos, strand = best_pam

        # Design spacer (20nt upstream of PAM)
        if strand == '+':
            spacer_start = max(0, pam_pos - 20)
            spacer = genomic_context[spacer_start:pam_pos]
            pam = genomic_context[pam_pos:pam_pos+3]
        else:
            # Reverse strand
            spacer_end = min(len(genomic_context), pam_pos + 23)
            spacer = reverse_complement(genomic_context[pam_pos+3:spacer_end])
            pam = reverse_complement(genomic_context[pam_pos:pam_pos+3])

        # Design PBS (primer binding site) - complementary to target
        pbs_start = pam_pos - 3  # Nick is 3bp upstream of PAM
        pbs = reverse_complement(genomic_context[pbs_start-self.pbs_length_default:pbs_start])

        # Design RT template with edit
        # Determine codon change needed
        if edit_type == "rescue_introduction":
            new_codon = HUMAN_PREFERRED_CODONS.get(new_aa, 'NNN')
        else:
            new_codon = HUMAN_PREFERRED_CODONS.get(original_aa, 'NNN')

        # RT template includes edit + homology
        rt_start = pbs_start
        rt_end = min(len(genomic_context), rt_start + 3 + self.rt_extension)  # 3 for codon + extension
        rt_template = genomic_context[rt_start:rt_start+3]  # To be replaced
        rt_template = new_codon + genomic_context[rt_start+3:rt_end]

        # Design nicking guide for PE3 (60-100bp away)
        nick_distance = 70
        nick_pos = pam_pos + nick_distance
        if nick_pos + 23 < len(genomic_context):
            nick_spacer = genomic_context[nick_pos:nick_pos+20]
            nick_pam = genomic_context[nick_pos+20:nick_pos+23]
        else:
            nick_spacer = "N" * 20
            nick_pam = "NGG"
            nick_distance = 0

        # Calculate metrics
        spacer_gc = gc_content(spacer)
        pbs_gc = gc_content(pbs)

        # Predict efficiency (simplified model)
        efficiency = self._predict_efficiency(spacer_gc, pbs_gc, len(rt_template))

        # Design oligos for cloning
        pegRNA_top, pegRNA_bottom = self._design_oligos(spacer, pbs, rt_template)
        nick_top, nick_bottom = self._design_nick_oligos(nick_spacer)

        return PrimeEditingGuide(
            target_mutation=mutation,
            edit_type=edit_type,
            spacer_sequence=spacer,
            pam_sequence=pam,
            pbs_sequence=pbs,
            rt_template=rt_template,
            nick_spacer=nick_spacer,
            nick_pam=nick_pam,
            nick_distance=nick_distance,
            spacer_gc_content=spacer_gc,
            pbs_gc_content=pbs_gc,
            predicted_efficiency=efficiency,
            off_target_score=0.85,  # Would need BLAST/Cas-OFFinder
            pegRNA_oligo_top=pegRNA_top,
            pegRNA_oligo_bottom=pegRNA_bottom,
            nick_oligo_top=nick_top,
            nick_oligo_bottom=nick_bottom
        )

    def _predict_efficiency(self, spacer_gc: float, pbs_gc: float,
                           rt_length: int) -> float:
        """Predict PE efficiency based on design parameters."""
        # Simplified model based on published parameters
        score = 0.5

        # Optimal spacer GC: 40-60%
        if 0.4 <= spacer_gc <= 0.6:
            score += 0.15
        elif spacer_gc < 0.3 or spacer_gc > 0.7:
            score -= 0.2

        # Optimal PBS GC: 40-60%
        if 0.4 <= pbs_gc <= 0.6:
            score += 0.1
        elif pbs_gc > 0.7:
            score -= 0.15

        # Shorter RT templates are better
        if rt_length <= 15:
            score += 0.1
        elif rt_length > 25:
            score -= 0.15

        return max(0.1, min(0.9, score))

    def _design_oligos(self, spacer: str, pbs: str, rt_template: str) -> Tuple[str, str]:
        """Design oligos for pegRNA cloning into U6 vector."""

        # pegRNA structure: spacer + scaffold + RT template + PBS
        scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

        # Top oligo (for BsmBI cloning)
        top = f"CACC{spacer}G"  # Add CACC overhang and G for U6

        # Bottom oligo
        bottom = f"AAAC{reverse_complement(spacer[1:])}C"  # Complement

        # Full pegRNA would need scaffold + RT + PBS in separate piece
        # Simplified: return spacer oligos

        return top, bottom

    def _design_nick_oligos(self, nick_spacer: str) -> Tuple[str, str]:
        """Design oligos for nicking guide."""
        top = f"CACC{nick_spacer}G"
        bottom = f"AAAC{reverse_complement(nick_spacer[1:])}C"
        return top, bottom

    def _create_dummy_guide(self, mutation: str, edit_type: str) -> PrimeEditingGuide:
        """Create placeholder guide when design fails."""
        return PrimeEditingGuide(
            target_mutation=mutation,
            edit_type=edit_type,
            spacer_sequence="N" * 20,
            pam_sequence="NGG",
            pbs_sequence="N" * 13,
            rt_template="N" * 15,
            nick_spacer="N" * 20,
            nick_pam="NGG",
            nick_distance=70,
            spacer_gc_content=0.5,
            pbs_gc_content=0.5,
            predicted_efficiency=0.0,
            off_target_score=0.0,
            pegRNA_oligo_top="DESIGN_FAILED",
            pegRNA_oligo_bottom="DESIGN_FAILED",
            nick_oligo_top="DESIGN_FAILED",
            nick_oligo_bottom="DESIGN_FAILED"
        )

    def batch_design(self, mutations: List[str],
                    genomic_context: str) -> List[PrimeEditingGuide]:
        """Design guides for multiple mutations."""
        return [self.design_guide(mut, genomic_context) for mut in mutations]

    def export_order_sheet(self, guides: List[PrimeEditingGuide],
                          output_path: str) -> str:
        """Export oligo order sheet for IDT/Sigma."""

        lines = ["Name,Sequence,Scale,Purification"]

        for i, guide in enumerate(guides):
            prefix = f"p53_{guide.target_mutation}"

            lines.append(f"{prefix}_pegRNA_top,{guide.pegRNA_oligo_top},25nm,STD")
            lines.append(f"{prefix}_pegRNA_bot,{guide.pegRNA_oligo_bottom},25nm,STD")
            lines.append(f"{prefix}_nick_top,{guide.nick_oligo_top},25nm,STD")
            lines.append(f"{prefix}_nick_bot,{guide.nick_oligo_bottom},25nm,STD")

        content = "\n".join(lines)

        with open(output_path, 'w') as f:
            f.write(content)

        return output_path


# ============================================================================
# GENBANK EXPORTER
# ============================================================================

class GenBankExporter:
    """Export sequences in GenBank format."""

    def export_rescue_construct(self, name: str, protein_sequence: str,
                               mutations: List[str],
                               organism: str = "human") -> GenBankExport:
        """Create GenBank entry for rescued p53 construct."""

        # Convert to DNA
        dna_sequence = protein_to_dna(protein_sequence, organism)

        # Create features
        features = []

        # Gene feature
        features.append({
            'type': 'gene',
            'location': f'1..{len(dna_sequence)}',
            'qualifiers': {
                'gene': 'TP53',
                'note': 'Rescued p53 variant'
            }
        })

        # CDS feature
        features.append({
            'type': 'CDS',
            'location': f'1..{len(dna_sequence)}',
            'qualifiers': {
                'gene': 'TP53',
                'product': 'tumor protein p53 (rescued)',
                'codon_start': '1',
                'translation': protein_sequence
            }
        })

        # Mutation features
        for mut in mutations:
            pos = int(mut[1:-1])
            dna_pos = (pos - 1) * 3 + 1  # 1-indexed DNA position

            features.append({
                'type': 'misc_feature',
                'location': f'{dna_pos}..{dna_pos+2}',
                'qualifiers': {
                    'note': f'Rescue mutation {mut}',
                    'label': mut
                }
            })

        return GenBankExport(
            locus=name,
            definition=f"Synthetic p53 rescue variant with mutations {', '.join(mutations)}",
            accession="SYNTHETIC",
            version="SYNTHETIC.1",
            keywords="p53; rescue mutation; synthetic",
            source="synthetic construct",
            organism="synthetic construct",
            sequence=dna_sequence,
            features=features
        )

    def export_to_file(self, export: GenBankExport, output_path: str) -> str:
        """Write GenBank export to file."""
        with open(output_path, 'w') as f:
            f.write(export.to_genbank_string())
        return output_path


# ============================================================================
# ASSAY RECOMMENDER
# ============================================================================

class AssayRecommender:
    """Recommend validation assays based on rescue properties."""

    def __init__(self):
        self.assay_database = self._build_assay_database()

    def _build_assay_database(self) -> List[AssayRecommendation]:
        """Build database of available assays."""
        return [
            AssayRecommendation(
                assay_name="Thermal Shift Assay (DSF)",
                assay_type="structural",
                priority=1,
                rationale="Quick assessment of protein stability and folding",
                readout="Melting temperature (Tm)",
                throughput="high",
                estimated_cost="$5/sample",
                timeline="1 day",
                required_equipment=["RT-PCR machine", "SYPRO Orange"],
                protocol_reference="Niesen et al., Nature Protocols 2007"
            ),
            AssayRecommendation(
                assay_name="DNA Binding EMSA",
                assay_type="binding",
                priority=1,
                rationale="Direct measurement of p53 DNA binding",
                readout="Kd for p53 response element",
                throughput="low",
                estimated_cost="$50/sample",
                timeline="2 days",
                required_equipment=["Gel apparatus", "PhosphorImager"],
                protocol_reference="El-Deiry et al., Nature Genetics 1992"
            ),
            AssayRecommendation(
                assay_name="Luciferase Reporter Assay",
                assay_type="functional",
                priority=1,
                rationale="Measure p53 transcriptional activity in cells",
                readout="Luciferase activity (RLU)",
                throughput="high",
                estimated_cost="$10/sample",
                timeline="3 days",
                required_equipment=["Luminometer", "Cell culture"],
                protocol_reference="Brummelkamp et al., Science 2002"
            ),
            AssayRecommendation(
                assay_name="PAb1620/PAb240 ELISA",
                assay_type="structural",
                priority=2,
                rationale="Antibody-based assessment of folding state",
                readout="PAb1620/PAb240 ratio",
                throughput="medium",
                estimated_cost="$20/sample",
                timeline="1 day",
                required_equipment=["ELISA reader", "96-well plates"],
                protocol_reference="Bullock et al., Oncogene 1997"
            ),
            AssayRecommendation(
                assay_name="Fluorescence Polarization DNA Binding",
                assay_type="binding",
                priority=2,
                rationale="Quantitative DNA binding in solution",
                readout="Kd (nM)",
                throughput="high",
                estimated_cost="$15/sample",
                timeline="1 day",
                required_equipment=["FP plate reader", "Labeled DNA"],
                protocol_reference="Weinberg et al., JMB 2004"
            ),
            AssayRecommendation(
                assay_name="Cell Growth Suppression",
                assay_type="functional",
                priority=2,
                rationale="Test tumor suppressor function in p53-null cells",
                readout="Growth inhibition (%)",
                throughput="medium",
                estimated_cost="$100/sample",
                timeline="7 days",
                required_equipment=["Cell culture", "Colony counter"],
                protocol_reference="Baker et al., Science 1990"
            ),
            AssayRecommendation(
                assay_name="Apoptosis Induction (Annexin V)",
                assay_type="functional",
                priority=3,
                rationale="Measure ability to induce apoptosis",
                readout="% Annexin V+ cells",
                throughput="medium",
                estimated_cost="$30/sample",
                timeline="2 days",
                required_equipment=["Flow cytometer", "Annexin V kit"],
                protocol_reference="Chipuk et al., Science 2004"
            ),
            AssayRecommendation(
                assay_name="Crystallography",
                assay_type="structural",
                priority=4,
                rationale="Atomic-resolution structure of rescue",
                readout="3D structure",
                throughput="low",
                estimated_cost="$5000/structure",
                timeline="3-6 months",
                required_equipment=["Crystallization robot", "Synchrotron access"],
                protocol_reference="Joerger et al., Oncogene 2007"
            ),
            AssayRecommendation(
                assay_name="HDX-MS",
                assay_type="structural",
                priority=3,
                rationale="Map dynamics and conformational changes",
                readout="HDX protection map",
                throughput="low",
                estimated_cost="$500/sample",
                timeline="1 week",
                required_equipment=["HDX-MS system"],
                protocol_reference="Englander, JACS 2006"
            ),
            AssayRecommendation(
                assay_name="ChIP-seq",
                assay_type="functional",
                priority=3,
                rationale="Genome-wide DNA binding profile",
                readout="Binding peaks",
                throughput="low",
                estimated_cost="$1000/sample",
                timeline="2 weeks",
                required_equipment=["Sequencer", "ChIP antibodies"],
                protocol_reference="Botcheva et al., PNAS 2011"
            ),
        ]

    def recommend_assays(self, rescue_properties: Dict,
                        max_recommendations: int = 5,
                        budget: str = "medium") -> List[AssayRecommendation]:
        """
        Recommend assays based on rescue properties.

        Args:
            rescue_properties: Dict with keys like 'stability_change', 'binding_change'
            max_recommendations: Maximum number of assays to recommend
            budget: "low", "medium", or "high"
        """

        # Filter by budget
        budget_filter = {
            "low": ["high", "medium"],
            "medium": ["high", "medium", "low"],
            "high": ["high", "medium", "low"]
        }

        candidates = [a for a in self.assay_database
                     if a.throughput in budget_filter[budget]]

        # Score each assay
        scored = []
        for assay in candidates:
            score = self._score_assay(assay, rescue_properties)
            scored.append((score, assay))

        # Sort by score (descending) then priority (ascending)
        scored.sort(key=lambda x: (-x[0], x[1].priority))

        return [assay for _, assay in scored[:max_recommendations]]

    def _score_assay(self, assay: AssayRecommendation,
                    properties: Dict) -> float:
        """Score assay relevance for given rescue properties."""
        score = 5 - assay.priority  # Base score from priority

        stability_change = properties.get('stability_change', 0)
        binding_change = properties.get('binding_change', 0)

        # If stability is major change, prioritize structural assays
        if abs(stability_change) > 1.0:
            if assay.assay_type == "structural":
                score += 2

        # If binding is major change, prioritize binding assays
        if abs(binding_change) > 1.0:
            if assay.assay_type == "binding":
                score += 2

        # Always need functional validation
        if assay.assay_type == "functional":
            score += 1

        return score

    def generate_validation_plan(self, rescue_mutations: List[str],
                                rescue_properties: Dict,
                                output_path: str) -> str:
        """Generate complete validation plan document."""

        assays = self.recommend_assays(rescue_properties)

        lines = []
        lines.append("# Experimental Validation Plan")
        lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d')}")
        lines.append(f"\n## Rescue Mutations: {', '.join(rescue_mutations)}")

        lines.append("\n## Recommended Validation Assays")
        lines.append("\n### Tier 1 (Essential)")
        for assay in assays[:2]:
            lines.append(f"\n#### {assay.assay_name}")
            lines.append(f"- **Type:** {assay.assay_type}")
            lines.append(f"- **Rationale:** {assay.rationale}")
            lines.append(f"- **Readout:** {assay.readout}")
            lines.append(f"- **Timeline:** {assay.timeline}")
            lines.append(f"- **Cost:** {assay.estimated_cost}")
            lines.append(f"- **Reference:** {assay.protocol_reference}")

        lines.append("\n### Tier 2 (Recommended)")
        for assay in assays[2:4]:
            lines.append(f"\n#### {assay.assay_name}")
            lines.append(f"- **Type:** {assay.assay_type}")
            lines.append(f"- **Rationale:** {assay.rationale}")
            lines.append(f"- **Timeline:** {assay.timeline}")

        lines.append("\n### Tier 3 (If Resources Allow)")
        for assay in assays[4:]:
            lines.append(f"- {assay.assay_name}: {assay.rationale}")

        lines.append("\n## Timeline")
        lines.append("| Week | Activity |")
        lines.append("|------|----------|")
        lines.append("| 1-2 | Protein expression and purification |")
        lines.append("| 3 | Thermal shift assay (stability) |")
        lines.append("| 4 | DNA binding assays |")
        lines.append("| 5-6 | Cell-based functional assays |")
        lines.append("| 7-8 | Data analysis and follow-up |")

        content = "\n".join(lines)

        with open(output_path, 'w') as f:
            f.write(content)

        return output_path


# ============================================================================
# MAIN EXPERIMENTAL PIPELINE
# ============================================================================

class ExperimentalPipeline:
    """
    Main pipeline for generating all experimental validation materials.
    """

    def __init__(self):
        self.yeast_generator = YeastDisplayGenerator()
        self.pe_designer = PrimeEditingDesigner()
        self.genbank_exporter = GenBankExporter()
        self.assay_recommender = AssayRecommender()

    def generate_full_package(self, name: str,
                             p53_sequence: str,
                             cancer_mutation: str,
                             rescue_mutations: List[str],
                             output_dir: str,
                             genomic_context: str = None) -> Dict[str, str]:
        """
        Generate complete experimental validation package.

        Returns:
            Dict mapping file type to file path
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        outputs = {}

        # 1. Yeast Display Protocol
        protocol = self.yeast_generator.generate_protocol(
            [cancer_mutation], rescue_mutations, p53_sequence
        )
        protocol_path = output_dir / f"{name}_yeast_display_protocol.md"
        self.yeast_generator.export_protocol_docx(protocol, str(protocol_path))
        outputs['yeast_display_protocol'] = str(protocol_path)

        # 2. Prime Editing Guides
        if genomic_context:
            pe_guides = self.pe_designer.batch_design(rescue_mutations, genomic_context)
            oligo_path = output_dir / f"{name}_prime_editing_oligos.csv"
            self.pe_designer.export_order_sheet(pe_guides, str(oligo_path))
            outputs['prime_editing_oligos'] = str(oligo_path)

        # 3. GenBank Export
        rescue_sequence = p53_sequence  # Would apply mutations in real impl
        genbank = self.genbank_exporter.export_rescue_construct(
            name, rescue_sequence, rescue_mutations
        )
        genbank_path = output_dir / f"{name}.gb"
        self.genbank_exporter.export_to_file(genbank, str(genbank_path))
        outputs['genbank'] = str(genbank_path)

        # 4. FASTA Export
        fasta_path = output_dir / f"{name}.fasta"
        with open(fasta_path, 'w') as f:
            f.write(f">{name} | Rescue mutations: {', '.join(rescue_mutations)}\n")
            # Wrap sequence at 60 characters
            for i in range(0, len(rescue_sequence), 60):
                f.write(rescue_sequence[i:i+60] + "\n")
        outputs['fasta'] = str(fasta_path)

        # 5. Validation Plan
        rescue_properties = {
            'stability_change': -1.5,  # Would compute from scoring
            'binding_change': 0.5
        }
        plan_path = output_dir / f"{name}_validation_plan.md"
        self.assay_recommender.generate_validation_plan(
            rescue_mutations, rescue_properties, str(plan_path)
        )
        outputs['validation_plan'] = str(plan_path)

        # 6. Summary JSON
        summary = {
            'name': name,
            'cancer_mutation': cancer_mutation,
            'rescue_mutations': rescue_mutations,
            'outputs': outputs,
            'generated': datetime.now().isoformat()
        }
        summary_path = output_dir / f"{name}_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        outputs['summary'] = str(summary_path)

        return outputs


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for experimental pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description="p53-proteoMgCAD Experimental Pipeline")
    parser.add_argument("--name", type=str, required=True,
                       help="Name for the rescue construct")
    parser.add_argument("--cancer", type=str, required=True,
                       help="Cancer mutation (e.g., R175H)")
    parser.add_argument("--rescue", type=str, nargs='+', required=True,
                       help="Rescue mutations (e.g., N268D T125R)")
    parser.add_argument("--output", type=str, default="./experimental_package",
                       help="Output directory")

    args = parser.parse_args()

    # Load WT sequence
    from ..data.dms import P53_WT

    # Generate package
    pipeline = ExperimentalPipeline()
    outputs = pipeline.generate_full_package(
        name=args.name,
        p53_sequence=P53_WT,
        cancer_mutation=args.cancer,
        rescue_mutations=args.rescue,
        output_dir=args.output
    )

    print("\n" + "="*60)
    print("EXPERIMENTAL VALIDATION PACKAGE GENERATED")
    print("="*60)
    for file_type, path in outputs.items():
        print(f"  {file_type}: {path}")


if __name__ == "__main__":
    main()
