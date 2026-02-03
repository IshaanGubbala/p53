# Risk Scoring Bug Fix: Before vs After

**Date:** January 25, 2026
**Issue:** MSA conservation scoring was disabled due to path resolution bug
**Impact:** Risk scores artificially compressed to 0.000-0.075 range

---

## The Bug

### Path Resolution Error

```python
# Config had:
msa:
  precomputed_path: Data/processed/msa/P04637_conservation.json

# Code did:
processed_dir = Path("Data/processed")
msa_path = Path(config_path)  # "Data/processed/msa/..."
resolved = processed_dir / msa_path  # WRONG: Data/processed/Data/processed/msa/...

# Result: File not found, MSA map returned as empty {}
```

### Consequence

All 685 candidates had:
- `msa_conservation_risk = 0.0` (should be 0.16-0.61)
- `conservation_risk = 0.0` (no conservation map configured)
- `functional_risk = 0.0` (8 Å filter removed everything)
- Only `burial_risk` varied (0.0 buried, 0.5 partial)

**Effective risk formula:** `risk = 0.15 × burial_risk`
**Result:** Only 2 possible risk values: 0.000 (buried) or 0.075 (partial)

---

## Fix

Changed config path to be relative to `processed_dir`:

```yaml
msa:
  precomputed_path: msa/P04637_conservation.json  # Relative to Data/processed
```

---

## Impact on Rankings

### Before (Broken MSA Scoring)

| Target | Top 3 "Rescues" | Risk | Note |
|--------|----------------|------|------|
| R175H | M133L, A189S, S215A | 0.000 | All "perfect" |
| R248Q | M133L, S215A, L145I | 0.000 | All "perfect" |
| R273H | A189S, Y163F, M133L | 0.000 | All "perfect" |
| Y220C | M133L, A189S, A189G | 0.000 | All "perfect" |

- **142 candidates (20.7%) had risk = 0.000** - "perfect" scores
- **100% of candidates had risk < 0.1** - no discrimination
- Rankings were determined **solely by burial** (binary: buried or partial)

### After (Fixed MSA Scoring)

| Target | Top 3 Stabilizing Rescues (ΔΔG < -1.0) | Risk | ΔΔG |
|--------|----------------------------------------|------|-----|
| R175H | T155A, L145I, M133L,T155A | 0.095, 0.101, 0.108 | -1.53, -1.62, -7.13 |
| R248Q | T155A, L145I, M133L | 0.095, 0.101, 0.120 | -1.53, -1.62, -5.60 |
| R273H | T155A, L145I, M133L | 0.095, 0.101, 0.120 | -1.52, -1.61, -2.93 |
| Y220C | L145I, M133V, M133L | 0.101, 0.120, 0.120 | -1.46, -2.57, -5.60 |

- **0 candidates with risk = 0.000** - realistic discrimination
- **Risk range: 0.095 - 0.228** - 3x wider range
- **MSA conservation (0.38-0.61) now drives risk differences**

---

## Biological Meaning Restored

### MSA Conservation Now Matters

| Position | MSA Conservation | Interpretation | Example Rescues |
|----------|------------------|----------------|----------------|
| 155 | 0.38 | Evolutionarily flexible | T155A, T155S (risk=0.095) |
| 145 | 0.40 | Flexible | L145I, L145V (risk=0.101) |
| 133 | 0.48 | Moderate importance | M133L, M133V (risk=0.120) |
| 189 | 0.50 | Moderately conserved | A189S, A189G (risk=0.125) |
| 196 | 0.61 | Highly conserved | R196Q/K/H (risk=0.228) |
| 215 | 0.61 | Highly conserved | S215A (risk=0.153) |

**Key insight:** R196 has high MSA conservation (0.61) because it forms critical internal salt bridges, even though it doesn't contact DNA. This is biologically correct!

---

## Where Did the "Old Favorites" Go?

### M133L Performance

| Target | Old Rank | New Rank (All) | New Rank (Stabilizing) | ΔΔG | MSA |
|--------|----------|----------------|----------------------|-----|-----|
| R175H | #1 | #11 | **#6** | -5.60 | 0.48 |
| R248Q | #1 | #9 | **#3** | -5.60 | 0.48 |
| R273H | #3 | #7 | **#3** | -2.93 | 0.48 |
| Y220C | #1 | #9 | **#2** | -5.60 | 0.48 |

**M133L is still excellent** - ranks #2-6 among stabilizing rescues. It's no longer artificially "perfect" (risk=0.0), but it remains a top candidate due to strong stability gains.

### A189S Performance

| Target | Old Rank | New Rank (All) | New Rank (Stabilizing) | ΔΔG | MSA |
|--------|----------|----------------|----------------------|-----|-----|
| R175H | #2 | #26 | #14 | -5.21 | 0.50 |
| R248Q | #2 | #20 | #9 | -5.21 | 0.50 |
| R273H | #1 | #10 | **#5** | -4.51 | 0.50 |
| Y220C | #2 | #19 | #9 | -4.51 | 0.50 |

A189S also remains strong, with moderate MSA conservation (0.50) and good stability.

### A189G Performance

Similar to A189S (same position, same MSA conservation).

---

## Multi-Mutation Rescues

**Combinations with M133L dominate top rankings:**

| Target | Rescue | Risk | ΔΔG | Rank |
|--------|--------|------|-----|------|
| R175H | M133L,T155A | 0.108 | -7.13 | #3 |
| R175H | A189S,M133L | 0.123 | -10.81 | #10 |
| R248Q | A189S,M133L | 0.123 | -10.12 | #5 |
| R273H | A189S,M133L | 0.123 | -7.44 | #4 |
| Y220C | A189S,M133L | 0.123 | -10.12 | #6 |

Multi-mutation rescues achieve **stronger stabilization** (ΔΔG -7 to -11 kcal/mol) while maintaining moderate risk (0.11-0.12).

---

## Key Takeaways

### ✅ What's Better Now

1. **Realistic risk discrimination** - 3x wider range (0.095-0.228 vs 0.000-0.075)
2. **MSA conservation matters** - evolutionarily flexible positions favored
3. **No "perfect" candidates** - all have trade-offs (safety vs stability)
4. **Biologically meaningful** - R196 correctly identified as conserved (internal salt bridges)
5. **True Pareto optimization** - real stability/safety trade-off, not artifacts

### ⚠️ What Changed

1. **No more "zero risk" rescues** - everything has some conservation cost
2. **T155A is new universal #1** (lowest MSA) but weak stability (-1.5 kcal/mol)
3. **M133L dropped to #2-6** but remains top for stabilization
4. **Multi-mutations preferred** for strong stabilization with acceptable risk

### 🎯 Best Strategy

**For experimental validation, prioritize:**
1. **Single rescues:** M133L, L145I, A189S (balance of stability + safety)
2. **Multi-rescues:** M133L+T155A, A189S+M133L (strong stabilization)
3. **Avoid:** T155V, V217L (low risk but destabilizing or neutral)

---

## Technical Details

### Risk Component Breakdown (Example: M133L)

| Component | Weight | Before | After | Contribution |
|-----------|--------|--------|-------|--------------|
| Functional | 0.40 | 0.0 | 0.0 | 0.000 (≥8 Å filter) |
| Conservation | 0.20 | 0.0 | 0.0 | 0.000 (no map) |
| Burial | 0.15 | 0.0 | 0.0 | 0.000 (buried) |
| MSA | 0.25 | **0.0** | **0.48** | **0.120** |
| **Total** | 1.00 | **0.000** | **0.120** | |

### Theoretical Risk Range

With design filters:
- `min_distance_protected: 8.0 Å` → functional ≤ 0.0
- `max_conservation: 0.8` → conservation ≤ 0.8 (unused, no map)
- `max_msa_conservation: 0.85` → MSA ≤ 0.85
- `allow_exposed: false` → burial ≤ 0.5

**Maximum possible risk:** `0.40×0 + 0.20×0 + 0.15×0.5 + 0.25×0.85 = 0.288`
**Observed maximum:** 0.228 (positions have MSA ≤ 0.61 in practice)

---

## Validation Status

**Re-running validations with proper risk scoring:**
- ✅ Rescue design complete (all 4 targets)
- 🔄 DNA binding analysis (pending)
- 🔄 Tetramer interface analysis (pending)
- 🔄 Sensitivity analysis (pending)

**Expected changes:**
- Top candidates will have **higher** MSA conservation than before
- Validations may identify T155A/L145I as "safer" but weaker rescues
- M133L will still pass all validations (not at DNA contact or tetramer interface)

---

*Generated: January 25, 2026*
*Fixed bug identified through user questioning of suspicious "perfect" scores*
