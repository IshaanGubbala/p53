# DNA Binding Validation Results (Updated with Fixed MSA Scoring)

**Date:** January 25, 2026
**Status:** ✅ Complete
**Bug Fix:** MSA conservation scoring now enabled

---

## Executive Summary

**Key Finding:** All top rescue candidates were confirmed to be 100% safe for DNA binding. This result validates the effectiveness of the pipeline's constraint-aware design, which successfully filters out all mutations that would pose a risk to DNA interaction before optimization begins.

### Results by Target

| Target | Candidates Analyzed | Safe | Low | Medium | High | Critical |
|--------|---------------------|------|-----|--------|------|----------|
| R175H | 50 | 50 (100%) | 0 | 0 | 0 | 0 |
| R248Q | 50 | 50 (100%) | 0 | 0 | 0 | 0 |
| R273H | 50 | 50 (100%) | 0 | 0 | 0 | 0 |
| Y220C | 50 | 50 (100%) | 0 | 0 | 0 | 0 |

**Total:** 200/200 candidates are safe for DNA binding (100%)

---

## New Top Candidates Validation

With proper MSA scoring, the top candidates changed, but all remain safe for DNA binding:

### Top 5 Universal Rescues (by MSA conservation)

| Rescue | R175H | R248Q | R273H | Y220C | DNA Risk | Min Distance |
|--------|-------|-------|-------|-------|----------|--------------|
| **T155A** | Rank #1 | Rank #1 | Rank #1 | Rank #1 | **Safe** | >15 Å |
| **L145I** | Rank #4 | Rank #4 | Rank #4 | Rank #3 | **Safe** | >15 Å |
| **M133L** | Rank #11 | Rank #9 | Rank #7 | Rank #9 | **Safe** | >15 Å |
| **A189S** | Rank #26 | Rank #20 | Rank #10 | Rank #19 | **Safe** | >15 Å |
| **V217L** | Rank #8 | Rank #6 | Rank #4 | Rank #6 | **Safe** | >15 Å |

All new top candidates maintain >15 Å separation from DNA-contacting residues.

---

## DNA-Contacting Residues (from 1TSR structure)

**17 positions identified as DNA contacts:**

### Strong DNA Contacts (≥5 contacts, <3.5 Å)
- **Q165** (6 contacts, min 2.89 Å) - H-bond, electrostatic
- **Q167** (5 contacts, min 3.05 Å) - H-bond
- **H168** (6 contacts, min 2.87 Å) - H-bond, electrostatic
- **R249** (9 contacts, min 2.74 Å) - Electrostatic, critical

### Medium DNA Contacts (2-4 contacts, <4.5 Å)
- R248, P250, S241, K120, R280, R282, K283, R273

### Weak DNA Contacts (<2 contacts)
- A159, C277, D281, E286

---

## Critical Discovery: R196 is NOT a DNA Contact

### Evidence

1. **1TSR structure:** R196 is **22.46 Å** from nearest DNA atom
2. **1TUP structure:** R196 is **18.20 Å** from nearest DNA atom
3. **Literature:** R196 forms **internal salt bridges** with D184, D186, S183

### Implication

**R196Q/K/H mutations are SAFE for DNA binding** - contrary to initial annotations in config file.

R196's high MSA conservation (0.61) reflects its **structural role** (stabilizing internal salt bridges), not DNA binding. These rescues are now correctly classified as:
- DNA binding risk: **Safe** (no DNA contact)
- MSA conservation risk: **High** (0.61, important for structure)
- Overall assessment: **Use with caution** (conserved for structural reasons)

---

## Why All Candidates Pass DNA Binding Validation

### Constraint-Aware Design Architecture

The 100% DNA binding safety rate is **by design, not by accident**. The pipeline enforces DNA safety through explicit constraints **before** Pareto optimization:

**Pipeline Order:**
1. **Candidate Generation:** Beam search generates all possible rescue mutations
2. **Constraint Filtering (CRITICAL STEP):**
   - 8 Å minimum distance filter applied to all candidates
   - Mutations at or near DNA-contacting residues are **rejected before scoring**
   - See `src/design/candidate_filters.py::filter_by_protected_distance()`
3. **Scoring:** Only pre-filtered candidates are scored for stability
4. **Pareto Optimization:** Runs on already-safe candidates

**Why This Works:**

1. **Pre-filtering guarantees safety:**
   - DNA-contacting positions (Q165, R249, R248, etc.) are marked as "protected"
   - Any mutation within 8 Å of protected sites is filtered out
   - This happens **before** Pareto optimization sees the candidates

2. **Surface vs buried bias reinforces constraint:**
   - DNA-contacting residues are surface-exposed (low burial)
   - Even if they passed filtering, they would score poorly for stability
   - This is a secondary effect, not the primary safety mechanism

3. **MSA conservation provides additional filtering:**
   - DNA contacts are highly conserved (MSA 0.5-0.6)
   - High-conservation positions are deprioritized by risk scoring
   - Again, this reinforces the constraint but doesn't replace it

**Critical Methodological Note:**

This validation **tests the constraint system**, not the optimization algorithm. The 100% pass rate confirms that:
- ✅ The 8 Å distance constraint is correctly implemented
- ✅ No bugs allow unsafe candidates to slip through filtering
- ✅ The constraint system successfully enforces DNA binding safety

It does **not** demonstrate that Pareto optimization "naturally" or "inherently" avoids DNA contacts. The optimization never sees unsafe candidates in the first place.

---

## Validation Details

### Structures Analyzed
- **1TSR:** p53-DNA complex, 1.8 Å resolution
- **1TUP:** Alternative p53-DNA complex, 2.2 Å resolution
- **Contact definition:** Protein atom within 5.0 Å of DNA atom

### Risk Classification
- **Safe:** No mutations at DNA-contacting positions, distance >8 Å
- **Low:** Distance 5-8 Å from DNA contact
- **Medium:** Distance 3-5 Å from DNA contact
- **High:** Mutations at weak DNA contact positions
- **Critical:** Mutations at strong DNA contact positions (Q165, R249, etc.)

---

## Comparison: Before vs After MSA Fix

| Metric | Before (Broken MSA) | After (Fixed MSA) |
|--------|---------------------|-------------------|
| Top rescue | M133L | T155A |
| MSA conservation | 0.48 | 0.38 |
| Risk score | 0.000 | 0.095 |
| DNA safe? | Yes (>15 Å) | Yes (>15 Å) |
| Stability | -5.6 kcal/mol | -1.5 kcal/mol |

**Key insight:** New #1 rescue (T155A) has **lower MSA conservation** but **weaker stability** than old #1 (M133L). DNA safety is unaffected by MSA fix.

---

## Recommendations

### For Experimental Validation

**Prioritize these DNA-safe rescues:**

1. **High stability, acceptable MSA:** M133L (-5.6 kcal/mol, MSA=0.48)
2. **Balanced:** L145I (-1.6 kcal/mol, MSA=0.40)
3. **Low MSA, weak stability:** T155A (-1.5 kcal/mol, MSA=0.38)

All three are:
- ✅ Safe for DNA binding (>15 Å from DNA)
- ✅ Not at functional sites
- ✅ Pass tetramer validation (see separate report)

### Avoid

**No rescues need to be avoided for DNA binding reasons.** All 200 top candidates are safe.

---

## Files Generated

```
reports/dna_binding/
├── DNA_BINDING_FINDINGS.md           (this file)
├── R175H_dna_binding.csv             (50 rescues, all safe)
├── R248Q_dna_binding.csv             (50 rescues, all safe)
├── R273H_dna_binding.csv             (50 rescues, all safe)
├── Y220C_dna_binding.csv             (50 rescues, all safe)
├── all_targets_dna_binding.csv       (200 rescues total)
├── dna_contacts.json                 (17 DNA-contacting residues)
└── R175H_dna_binding_summary.json    (per-target statistics)
```

---

## Conclusion

**All top rescue candidates preserve DNA binding capability** by avoiding DNA-contacting residues. This validation confirms that:

1. ✅ Constraint system successfully enforces DNA binding safety (8 Å filter works correctly)
2. ✅ R196 is NOT a DNA contact (structural role only, safe to mutate)
3. ✅ MSA scoring fix did not change DNA binding safety (constraint system unchanged)
4. ✅ New top candidates (T155A, L145I) are equally safe as old favorites (M133L)

**Methodological Clarity:**
- The 100% pass rate validates the **constraint-aware design architecture**
- This is **not** an emergent property of Pareto optimization
- Safety is guaranteed by pre-filtering, confirmed by post-hoc validation

**Status: DNA binding validation PASSED for all candidates ✅**

---

*Generated: January 25, 2026*
*Updated after MSA conservation scoring bug fix*
*Structures: 1TSR, 1TUP (p53-DNA complexes)*
