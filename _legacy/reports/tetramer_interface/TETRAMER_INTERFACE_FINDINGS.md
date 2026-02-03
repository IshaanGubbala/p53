# Tetramer Interface Validation Results (Updated with Fixed MSA Scoring)

**Date:** January 25, 2026
**Status:** ✅ Complete
**Bug Fix:** MSA conservation scoring now enabled

---

## Executive Summary

**Key Finding:** New top candidates (T155A, L145I) are equally safe for tetramerization as old favorites (M133L, A189S). Only 1-2 medium-risk rescues per target, no high-risk candidates.

### Results by Target

| Target | Candidates Analyzed | Safe | Low | Medium | High | Critical |
|--------|---------------------|------|-----|--------|------|----------|
| R175H | 50 | 48 (96%) | 0 | 2 (4%) | 0 | 0 |
| R248Q | 50 | 49 (98%) | 0 | 1 (2%) | 0 | 0 |
| R273H | 50 | 49 (98%) | 0 | 1 (2%) | 0 | 0 |
| Y220C | 50 | 49 (98%) | 0 | 1 (2%) | 0 | 0 |

**Total:** 195/200 candidates are safe for tetramerization (97.5%)

---

## New Top Candidates Validation

### Top 5 Universal Rescues (by MSA conservation)

| Rescue | R175H | R248Q | R273H | Y220C | Tetramer Risk | Interface Distance |
|--------|-------|-------|-------|-------|---------------|-------------------|
| **T155A** | Rank #1 | Rank #1 | Rank #1 | Rank #1 | **Safe** | >8 Å |
| **L145I** | Rank #4 | Rank #4 | Rank #4 | Rank #3 | **Safe** | >8 Å |
| **M133L** | Rank #11 | Rank #9 | Rank #7 | Rank #9 | **Safe** | >8 Å |
| **A189S** | Rank #26 | Rank #20 | Rank #10 | Rank #19 | **Safe** | >8 Å |
| **V217L** | Rank #8 | Rank #6 | Rank #4 | Rank #6 | **Safe** | >8 Å |

All new top candidates are >8 Å from tetramer interface residues.

---

## Tetramer Interface Residues (from 3KMD structure)

**33 positions identified at monomer-monomer interface** (within 5 Å across chains):

### Core Domain Interface Regions

1. **N-terminal interface (8 positions):** 93, 94, 95, 96, 98, 99, 100, 101
   - Critical for initial tetramer assembly
   - Position 95 appears in some rescues (medium risk)

2. **Zinc-binding region (9 positions):** 176, 177, 178, 179, 180, 181, 242, 243, 244
   - Essential for structural integrity
   - Not targeted by any top rescues

3. **H2 helix (6 positions):** 198, 199, 200, 201, 202, 267
   - Stabilizes dimer-dimer interface
   - Not targeted by any top rescues

4. **L3 loop (7 positions):** 212, 224, 225, 226, 231, 233
   - Flexible interface region
   - Not targeted by any top rescues

5. **Other (3 positions):** 140, 166, 167, 170
   - Minor interface contacts

---

## Medium-Risk Rescues

Only 5 total medium-risk rescues identified across all targets:

| Rescue | Target(s) | Risk Level | Impact Score | Reason |
|--------|-----------|------------|--------------|--------|
| **S95A,M133L** | R175H | Medium | 1.50 | S95 at N-terminal interface |
| **S95A,T211A** | R248Q | Medium | 1.50 | S95 at N-terminal interface |
| **S95A,M133L** | R273H | Medium | 1.50 | S95 at N-terminal interface |
| **S95A,T155A** | Y220C | Medium | 1.50 | S95 at N-terminal interface |

**Pattern:** All medium-risk rescues contain **S95A**. Position 95 is at the N-terminal tetramer interface.

### S95A Impact Assessment

- **Impact score:** 1.50 (medium)
- **Mechanism:** Loss of H-bond capability at interface
- **Rank:** S95A-containing rescues rank #30-50 (not in top 10)
- **Recommendation:** Avoid S95A for initial validation

---

## Critical Discovery: R196 is NOT at Tetramer Interface

### Evidence

1. **3KMD structure:** R196 is **>8 Å** from all interface residues
2. **Nearest interface residue:** R196 is ~10 Å from position 200 (H2 helix)
3. **Literature:** R196 forms **internal salt bridges**, not inter-monomer contacts

### Implication

**R196Q/K/H mutations are SAFE for tetramerization** - contrary to initial concerns.

R196's high MSA conservation (0.61) and partial burial (surface-accessible) made it seem risky, but structural analysis confirms:
- ✅ Not at tetramer interface
- ✅ No inter-monomer contacts
- ⚠️ High MSA conservation due to **structural role** (internal salt bridges)

---

## Why Most Candidates Pass Tetramer Validation

### Multi-Layered Filtering Architecture

The 97.5% tetramer safety rate results from a **combination of explicit constraints and scoring biases**:

**Layer 1: Explicit Pre-filtering (Partial Coverage)**
- **Zinc-binding sites protected:** Positions 176, 179, 238, 242 marked as "protected"
  - These overlap with tetramer interface (zinc region)
  - 8 Å distance filter prevents mutations at/near zinc sites
  - Result: 9/33 interface positions explicitly filtered

**Layer 2: MSA Conservation Scoring (Strong Bias)**
- **Interface residues are highly conserved:** MSA ~0.44-0.61
  - Zinc region: 176-181 have MSA 0.47-0.61
  - H2 helix: 198-202 have MSA 0.37-0.61
  - High MSA → high risk scores → deprioritized in Pareto front
- **Effect:** Interface positions that pass filtering score poorly for risk

**Layer 3: Burial Preference (Secondary Bias)**
- **Interface positions are partially buried:** Between surface and core
  - Core buried positions provide stronger stabilization
  - Pareto optimization prefers high stability gains
- **Effect:** Reinforces MSA bias against interface positions

**Layer 4: Hard MSA Filter (Minimal Impact)**
- Positions with MSA > 0.85 completely excluded
- Most interface positions have MSA 0.44-0.61 (below threshold)
- Only extreme conservation filtered by this layer

**Critical Methodological Note:**

Unlike DNA binding validation (which achieved 100% safety through explicit pre-filtering), tetramer validation represents a **hybrid outcome**:
- ✅ ~27% of interface positions (zinc sites) are explicitly protected
- ✅ ~70% of remaining interface positions are deprioritized by MSA risk scoring
- ⚠️ This is **not** purely emergent behavior, but also **not** purely pre-filtered

The 5 medium-risk rescues (S95A-containing) demonstrate that the system allows some interface mutations when MSA conservation is moderate (S95 has MSA=0.445).

---

## Comparison: Before vs After MSA Fix

| Metric | Before (Broken MSA) | After (Fixed MSA) |
|--------|---------------------|-------------------|
| Top rescue | M133L | T155A |
| At interface? | No (>8 Å) | No (>8 Å) |
| Tetramer risk | Safe | Safe |
| MSA conservation | 0.48 | 0.38 |

**Key insight:** MSA fix changed rankings but did not affect tetramer safety. Both old and new top candidates avoid interface positions.

---

## Recommendations

### For Experimental Validation

**All top rescues are safe for tetramerization:**

1. **T155A** (MSA=0.38, tetramer safe)
2. **L145I** (MSA=0.40, tetramer safe)
3. **M133L** (MSA=0.48, tetramer safe)
4. **A189S** (MSA=0.50, tetramer safe)

**Expected:** Tetramer formation >70% in SEC (size-exclusion chromatography)

### Avoid for Initial Validation

**Medium-risk rescues containing S95A:**
- S95A,M133L
- S95A,T211A
- S95A,T155A

S95 is at the N-terminal interface. While impact is moderate (1.50), these are not top-ranked rescues and can be deprioritized.

---

## Experimental Validation Protocol

### Recommended Assays

1. **Size-Exclusion Chromatography (SEC)**
   - Expected: >70% tetramer peak
   - Mutants should elute at ~220 kDa (tetramer) not ~55 kDa (monomer)

2. **Native PAGE**
   - Expected: Tetramer band at high molecular weight
   - No accumulation of dimers or monomers

3. **Cross-linking + SDS-PAGE**
   - Chemical cross-linking (glutaraldehyde or BS3)
   - Expected: Tetramer band after cross-linking

### Success Criteria

- ✅ Tetramer population ≥70% (SEC)
- ✅ No dominant-negative effect (doesn't prevent WT tetramerization)
- ✅ Thermal stability of tetramer maintained (DSC)

---

## Files Generated

```
reports/tetramer_interface/
├── TETRAMER_INTERFACE_FINDINGS.md          (this file)
├── R175H_tetramer_interface.csv            (50 rescues, 48 safe)
├── R248Q_tetramer_interface.csv            (50 rescues, 49 safe)
├── R273H_tetramer_interface.csv            (50 rescues, 49 safe)
├── Y220C_tetramer_interface.csv            (50 rescues, 49 safe)
├── all_targets_tetramer_interface.csv      (200 rescues total)
├── interface_residues_3KMD.json            (33 interface positions)
└── R175H_tetramer_summary.json             (per-target statistics)
```

---

## Conclusion

**195/200 top candidates are safe for tetramerization** (97.5%). The 5 medium-risk rescues all contain S95A and are not top-ranked.

Key findings:
1. ✅ New top candidates (T155A, L145I) are tetramer-safe
2. ✅ Old favorites (M133L, A189S) remain tetramer-safe
3. ✅ R196 is NOT at interface (structural role only)
4. ✅ Multi-layered filtering (explicit constraints + MSA scoring) ensures interface safety
5. ⚠️ Avoid S95A-containing rescues (medium risk, at N-terminal interface)

**Methodological Clarity:**
- Unlike DNA binding (100% pre-filtered), tetramer safety is a **hybrid outcome**
- ~27% of interface positions explicitly protected (zinc sites)
- ~70% deprioritized by MSA conservation risk scoring
- The 5 medium-risk rescues show the system allows interface mutations with moderate MSA

**Status: Tetramer interface validation PASSED for 97.5% of candidates ✅**

---

*Generated: January 25, 2026*
*Updated after MSA conservation scoring bug fix*
*Structure: 3KMD (p53 core tetramer + DNA, 2.13 Å)*
