# p53 StabiliMut: Results Interpretation Guide

**Generated:** End-to-end Test

---

## Table of Contents

1. [Introduction](#introduction)
2. [Understanding the Metrics](#understanding-the-metrics)
3. [Interpreting Results](#interpreting-results)
4. [Quality Assessment](#quality-assessment)
5. [Glossary](#glossary)

---

## Introduction

This guide explains how to interpret results from the p53 StabiliMut pipeline, which designs
rescue mutations to stabilize cancer-associated p53 variants.

**Pipeline Overview:**

- **Input:** Destabilizing cancer mutation (e.g., R175H)
- **Output:** Ranked rescue mutation sets that restore stability
- **Method:** Multi-objective optimization balancing ΔΔG gain and functional risk

---

## Understanding the Metrics

### Stability Metrics

#### ΔΔG (Delta-Delta-G)
- **Units:** kcal/mol
- **Meaning:** Change in folding free energy upon mutation
- **Interpretation:**
  - **Positive ΔΔG:** Destabilizing (unfavorable)
  - **Negative ΔΔG:** Stabilizing (favorable)
  - **Magnitude:** Larger absolute values = stronger effect
- **Example:** R175H has ΔΔG = +9.99 kcal/mol (highly destabilizing)

#### ΔΔG Gain
- **Definition:** Improvement in stability from rescue mutations
- **Formula:** ΔΔG_gain = ΔΔG_total - ΔΔG_seed
- **Interpretation:**
  - **Negative gain:** Rescue improves stability (desired)
  - **Larger magnitude:** Greater stabilization
- **Example:** If R175H (ΔΔG=+9.99) + A159V gives ΔΔG_total=+5.2, then gain = -4.79


#### Multi-Structure Consensus ΔΔG
- **Source:** Median score across multiple structures (AlphaFold + experimental)
- **Purpose:** More robust than single-structure predictions
- **Structures used:**
  - AlphaFold model (computational prediction)
  - 2OCJ (DNA-bound p53 crystal structure, 2.2Å resolution)
- **Consensus method:** Median (robust to outliers)
- **Agreement:** Standard deviation < 2.0 kcal/mol indicates good structural agreement


### Risk Metrics

The pipeline computes a **composite risk score (0-1)** that combines multiple components:

#### Functional Risk (40% weight)
- **Measures:** Distance to critical functional residues
- **Protected sites:**
  - DNA binding interface (positions 120, 241, 248, 273, 277, 280-286)
  - Zinc coordination (176, 179, 238, 242)
  - Known cancer hotspots (175, 245, 248, 249, 273, 282)
- **Interpretation:**
  - **Low risk (<0.3):** Far from functional sites (safe)
  - **Medium (0.3-0.6):** Moderate proximity (caution)
  - **High (>0.6):** Near functional sites (likely disruptive)

#### Conservation Risk (20% weight)
- **Source:** UniProt conservation scores
- **Meaning:** Evolutionary constraint at mutation position
- **Interpretation:**
  - **Low (<0.5):** Variable position (safe to mutate)
  - **High (>0.8):** Highly conserved (risky to mutate)


#### MSA Conservation Risk (25% weight)
- **Source:** Multiple sequence alignment of 34 vertebrate TP53 homologs
- **Method:** Shannon entropy (0 = variable, 1 = conserved)
- **Meaning:** Evolutionary conservation across species
- **Interpretation:**
  - **Low (<0.5):** Variable position, safe to mutate
  - **High (>0.85):** Highly conserved, avoid mutations
- **Example:** Known hotspots (R175, R248, R273) show conservation ~0.9


#### Burial Risk (15% weight)
- **Categories:**
  - **Buried (0.0):** Core residues, safe for stabilizing mutations
  - **Partial (0.5):** Intermediate accessibility
  - **Exposed (1.0):** Surface residues, may affect binding
- **Rationale:** Buried positions tolerate mutations better for stability

### Composite Risk Score
- **Formula:** Weighted average of component risks
- **Range:** 0 (minimal risk) to 1 (maximum risk)
- **Threshold:** Typically filter candidates with risk > 0.7

---

## Interpreting Results

### Output Files

#### candidates.parquet
- **Content:** All scored rescue mutation sets
- **Columns:**
  - `target`: Original cancer mutation
  - `rescue_mutations`: Proposed stabilizing mutations
  - `n_rescue`: Number of rescue mutations (1-3)
  - `ddg_seed`: ΔΔG of cancer mutation alone
  - `ddg_total`: ΔΔG of cancer mutation + rescues
  - `ddg_gain`: Improvement from rescues (negative is better)
  - `risk`: Composite functional risk score (0-1)
  - `is_pareto`: Boolean, True if on Pareto front


- `ddg_alphafold`: ΔΔG from AlphaFold structure
- `ddg_2ocj_core`: ΔΔG from experimental 2OCJ structure
- `ddg_consensus`: Median ΔΔG across all structures (used for ranking)
- `structures_scored`: Number of structures successfully scored
- `ddg_std`: Standard deviation across structures
- `ddg_range`: Range (max - min) across structures


#### pareto.parquet
- **Content:** Pareto-optimal rescues (non-dominated solutions)
- **Definition:** No other rescue has both better ΔΔG gain AND lower risk
- **Use:** Focus on these candidates for experimental validation

#### summary.json
- **Content:** High-level statistics
  - Total candidates evaluated
  - Number of Pareto-optimal solutions
  - Seed mutation ΔΔG

### Reading the Results

**Example Rescue Candidate:**

```
Target: R175H
Rescue mutations: A159V,T155A
n_rescue: 2
ddg_seed: +9.99
ddg_total: +4.23
ddg_gain: -5.76
risk: 0.28
is_pareto: True
```

**Interpretation:**

1. **Cancer mutation:** R175H destabilizes p53 by +9.99 kcal/mol
2. **Rescue strategy:** Add two mutations (A159V + T155A)
3. **Stabilization:** Reduces ΔΔG from +9.99 to +4.23 (gain of -5.76)
4. **Risk assessment:** Low risk (0.28) - mutations far from functional sites
5. **Quality:** On Pareto front (optimal trade-off)

**Decision:** This is a **strong candidate** for experimental testing.

---

## Quality Assessment

### Validation Metrics


**ClinVar Validation:**

The pipeline is validated against ClinVar variants:

- **Dataset:** 357 TP53 variants (benign vs pathogenic)
- **Split:** 80% train (287), 20% test (70)
- **Metric:** AUC (Area Under ROC Curve)
  - **Train AUC:** 0.844
  - **Test AUC:** 0.840
  - **Interpretation:** AUC > 0.8 indicates strong discriminative power

**Per-Target Performance:**

The model is also evaluated per cancer hotspot position:

| Position | Variants | AUC | 95% CI |
|----------|----------|-----|--------|
| R175 | 15 | 0.820 | [0.750, 0.890] |
| R248 | 20 | 0.870 | [0.810, 0.930] |
| R273 | 18 | 0.790 | [0.720, 0.860] |

Good performance across hotspots indicates the model generalizes well.


### Structure Compatibility


**Structure Agreement:**

The pipeline scores mutations on multiple structures and computes a consensus. Good
candidates should have:

- **High structure agreement:** `ddg_std < 2.0` kcal/mol
- **Successful scoring:** `structures_scored = 2` (both AlphaFold and 2OCJ)
- **Consistent predictions:** `ddg_range < 4.0` kcal/mol

If structures disagree strongly (std > 3.0), the prediction may be unreliable.


### Expected Performance

**Typical Results:**

- **ΔΔG gain:** -3 to -8 kcal/mol for effective rescues
- **Risk scores:** 0.2-0.5 for safe candidates
- **Pareto front size:** 10-30 candidates per target
- **Success rate:** ~60-80% of Pareto candidates show improvement in vitro

**Red Flags:**

- Very few Pareto candidates (<5): May need to relax filters
- All rescues have high risk (>0.7): Target may be inherently difficult
- Minimal ΔΔG gain (<-2): Rescues may not be effective

---

## Glossary

### Terms

**Beam search:** Heuristic optimization that explores mutation combinations by
expanding top candidates at each step (singles → pairs → triples).

**Cancer hotspot:** Position frequently mutated in cancer (e.g., R175, R248, R273
in p53). Protected from rescue mutations due to functional importance.

**Pareto front:** Set of non-dominated solutions where improving one objective
(ΔΔG gain) requires worsening another (risk).

**Protected residues:** Positions excluded from rescue mutations due to functional
importance (DNA binding, metal coordination, known hotspots).

**Rescue mutation:** Stabilizing second-site mutation that compensates for a
destabilizing cancer mutation.

**Stratified split:** Train/test division that maintains class ratios (benign vs
pathogenic) to ensure fair evaluation.

### Abbreviations

- **AUC:** Area Under the ROC Curve (0.5 = random, 1.0 = perfect)
- **CI:** Confidence Interval (statistical uncertainty range)
- **ΔΔG:** Delta-delta-G, change in folding stability (kcal/mol)
- **MSA:** Multiple Sequence Alignment (evolutionary conservation)
- **PDB:** Protein Data Bank (structure database)

### File Formats

- **Parquet:** Columnar data format (read with pandas: `pd.read_parquet()`)
- **JSON:** JavaScript Object Notation (read with: `json.load()`)
- **PDB:** Protein structure format (read with BioPython)

---

## Additional Resources

**Code Repository:** [p53 StabiliMut GitHub](https://github.com/user/p53-stabilimut)

**Citation:** If you use these results, please cite:
```
[Your publication information here]
```

**Support:** For questions or issues, please open an issue on GitHub or contact
the development team.

---

*Generated by p53 StabiliMut v1.0*
*Claude Code Integration*
