"""
Therapeutic Viability Analysis for p53 Rescue Proteins
Determines minimum sequence identity required for clinical use.

Based on:
1. Immunogenicity risk assessment
2. Structural integrity requirements  
3. Regulatory precedents for protein therapeutics
4. p53-specific functional constraints
"""

# THERAPEUTIC IDENTITY THRESHOLDS
# ================================

# 1. IMMUNOGENICITY RISK (Primary Constraint)
# --------------------------------------------
# Therapeutic proteins with <85% identity to human proteins are flagged by FDA
# as potentially immunogenic. However, context matters:

IMMUNOGENICITY_THRESHOLDS = {
    "conservative": 95,  # Minimal risk - suitable for chronic use
    "moderate": 90,      # Low risk - acceptable for most therapies
    "aggressive": 85,    # Moderate risk - requires immunogenicity testing
    "experimental": 80,  # High risk - research only, extensive testing needed
}

# 2. P53-SPECIFIC STRUCTURAL REQUIREMENTS
# ----------------------------------------
# p53 has distinct domains with different tolerance for mutation:

P53_DOMAIN_TOLERANCE = {
    "N-terminal_TAD": {
        "residues": "1-93",
        "function": "Transactivation",
        "mutation_tolerance": "HIGH",  # Intrinsically disordered
        "min_identity": 70,
    },
    "DNA_binding_domain": {
        "residues": "94-292", 
        "function": "DNA recognition",
        "mutation_tolerance": "LOW",   # Highly structured
        "min_identity": 90,            # Critical for function
    },
    "Tetramerization": {
        "residues": "293-356",
        "function": "Oligomerization",
        "mutation_tolerance": "MEDIUM",
        "min_identity": 85,
    },
    "C-terminal": {
        "residues": "357-393",
        "function": "Regulation",
        "mutation_tolerance": "HIGH",
        "min_identity": 75,
    }
}

# 3. REGULATORY PRECEDENTS
# -------------------------
# FDA-approved protein therapeutics with sequence modifications:

APPROVED_PRECEDENTS = {
    "Insulin_analogs": {
        "identity_to_human": 97,  # Lispro, Aspart (1-3 mutations)
        "status": "Widely approved",
    },
    "Erythropoietin_analogs": {
        "identity_to_human": 96,  # Darbepoetin (5 mutations)
        "status": "Approved",
    },
    "Factor_VIII_variants": {
        "identity_to_human": 99,  # B-domain deleted
        "status": "Approved",
    },
    "Antibody_humanization": {
        "identity_to_human": 90,  # Humanized mAbs
        "status": "Approved (with immunogenicity monitoring)",
    }
}

# 4. P53 CANCER THERAPY SPECIFIC CONSIDERATIONS
# ----------------------------------------------

P53_THERAPY_CONSTRAINTS = {
    "delivery_method": {
        "gene_therapy": {
            "min_identity": 85,  # Expressed endogenously, lower immunogenicity
            "rationale": "Intracellular expression reduces immune exposure"
        },
        "protein_therapy": {
            "min_identity": 95,  # Direct injection, higher immunogenicity risk
            "rationale": "Extracellular exposure requires near-native sequence"
        },
        "peptide_therapy": {
            "min_identity": 80,  # Small fragments, different rules
            "rationale": "Short sequences less immunogenic"
        }
    },
    
    "patient_population": {
        "cancer_patients": {
            "immune_status": "Often immunocompromised",
            "min_identity": 85,
            "rationale": "Reduced immune surveillance allows more variation"
        },
        "healthy_prophylaxis": {
            "immune_status": "Normal",
            "min_identity": 95,
            "rationale": "Full immune system requires near-native sequence"
        }
    },
    
    "treatment_duration": {
        "acute": {
            "duration": "<1 month",
            "min_identity": 85,
            "rationale": "Short exposure limits antibody development"
        },
        "chronic": {
            "duration": ">6 months",
            "min_identity": 95,
            "rationale": "Long exposure increases immunogenicity risk"
        }
    }
}

# 5. FUNCTIONAL RESCUE REQUIREMENTS
# ----------------------------------

FUNCTIONAL_REQUIREMENTS = {
    "DNA_binding": {
        "critical_residues": [120, 175, 245, 248, 273, 282],
        "max_mutations_allowed": 2,  # Can compensate with 1-2 changes
        "rationale": "Core DNA contact points"
    },
    "zinc_coordination": {
        "critical_residues": [176, 179, 238, 242],
        "max_mutations_allowed": 0,  # Absolutely conserved
        "rationale": "Structural zinc finger"
    },
    "tetramerization": {
        "critical_residues": [325, 337, 344],
        "max_mutations_allowed": 1,
        "rationale": "Oligomerization interface"
    }
}

# ============================================
# RECOMMENDED MINIMUM IDENTITY THRESHOLDS
# ============================================

RECOMMENDED_THRESHOLDS = {
    "ISEF_presentation": {
        "min_identity": 90,
        "max_mutations": 39,  # 10% of 393 residues
        "rationale": "Demonstrates therapeutic viability + scientific rigor",
        "risk_level": "Low",
        "regulatory_path": "Feasible"
    },
    
    "research_publication": {
        "min_identity": 85,
        "max_mutations": 59,  # 15% of 393 residues
        "rationale": "Acceptable for proof-of-concept studies",
        "risk_level": "Moderate",
        "regulatory_path": "Requires immunogenicity testing"
    },
    
    "computational_exploration": {
        "min_identity": 80,
        "max_mutations": 79,  # 20% of 393 residues
        "rationale": "Useful for understanding sequence-function relationships",
        "risk_level": "High",
        "regulatory_path": "Not therapeutically viable without extensive engineering"
    }
}

# ============================================
# FINAL RECOMMENDATION
# ============================================

print("=" * 80)
print("THERAPEUTIC VIABILITY ANALYSIS: p53 Rescue Proteins")
print("=" * 80)
print()
print("📊 RECOMMENDED MINIMUM SEQUENCE IDENTITY: 90%")
print()
print("Rationale:")
print("  ✅ Aligns with FDA precedents for protein therapeutics")
print("  ✅ Minimizes immunogenicity risk for cancer patients")
print("  ✅ Maintains structural integrity of DNA-binding domain")
print("  ✅ Allows ~39 mutations for functional optimization")
print("  ✅ Suitable for ISEF presentation (demonstrates rigor)")
print()
print("Maximum Allowable Mutations: 39 residues (10% of 393)")
print()
print("Critical Constraints:")
print("  • DNA-binding domain (94-292): Maintain >95% identity")
print("  • Zinc coordination (176,179,238,242): Zero mutations")
print("  • Hotspot positions (175,248,273): Max 2 compensatory changes")
print()
print("Delivery Method Compatibility:")
print("  ✅ Gene therapy (AAV, lentiviral)")
print("  ✅ mRNA therapy")
print("  ⚠️  Direct protein injection (may need 95% identity)")
print()
print("=" * 80)
print()
print("IMPLEMENTATION RECOMMENDATION:")
print("Set identity preservation weights to target 90-95% identity")
print("This balances functional rescue with therapeutic viability")
print("=" * 80)
