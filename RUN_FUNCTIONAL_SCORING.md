# Running Functional Scoring Pipeline

## Quick Start

### Test on Single Target (R175H)

```bash
# Run R175H with functional scoring
python -m experiments.run_design_rescues --targets R175H --functional-scoring

# Expected runtime: ~1.5 hours (162 Pareto rescues)
# Output: Data/processed/rescues/R175H/pareto.parquet (with functional scores)
```

### Run All Targets

```bash
# Run all 4 targets with functional scoring
bash experiments/run_all_functional.sh

# Expected runtime: 3-5 hours total
# Targets: R175H, R248Q, R273H, Y220C
```

### Analyze Results

```bash
# Compare old vs new results across all targets
python experiments/analyze_functional_results.py

# Generates summary with:
# - Category distribution (excellent/good/acceptable/poor)
# - DNA binding pass rates
# - Tetramer interface pass rates
# - Top 10 rescues per target
# - Filtered-out rescues
```

---

## What's New

### Before (Stability-Only)

```yaml
Optimization objective: ΔΔG_folding (monomer stability)

Top rescue: M133L
  ΔΔG_folding = -5.6 kcal/mol
  Result: "Stabilizes p53 monomer"
```

### After (Functional Scoring)

```yaml
Optimization: Multi-dimensional functional restoration

Top rescue: M133L
  ΔΔG_folding   = -5.6 kcal/mol  (monomer stability)
  ΔΔG_binding   =  0.0 kcal/mol  (DNA binding preserved)
  ΔΔG_interface =  0.0 kcal/mol  (tetramer preserved)
  functional_score = 0.700
  Result: "Restores p53 function across all dimensions"
```

---

## Output Files

### Enhanced Pareto Front

**File:** `Data/processed/rescues/{TARGET}/pareto.parquet`

**New columns added:**
- `functional_score` - Composite 0-1 score (higher = better)
- `ddg_binding` - DNA binding affinity (kcal/mol)
- `ddg_interface` - Tetramer interface stability (kcal/mol)
- `folding_norm`, `binding_norm`, `interface_norm`, `risk_norm` - Normalized scores (0-1)
- `overall_category` - "excellent", "good", "acceptable", "poor"
- `binding_category`, `interface_category` - Per-dimension categories

**Existing columns preserved:**
- `rescue_mutations`, `ddg_gain`, `ddg_total`, `risk`, `rasp_gain`, etc.

---

## Understanding the Results

### Categories

**Excellent:** All 3 dimensions "good"
- ΔΔG_folding ≤ -3.0
- ΔΔG_binding ≤ 0.0
- ΔΔG_interface ≤ 0.0

**Good:** ≥2 dimensions "good", none "bad"

**Acceptable:** ≥1 dimension "good", or all "acceptable"

**Poor:** Any dimension "bad"
- Indicates rescue likely fails in experiments
- These are **false positives** from stability-only optimization

### Functional Score Formula

```python
functional_score = 0.35 × norm(ΔΔG_folding) +
                  0.35 × norm(ΔΔG_binding) +
                  0.20 × norm(ΔΔG_interface) +
                  0.10 × norm(risk)
```

**Weights (customizable in configs/functional_scoring.yaml):**
- 35% folding stability
- 35% DNA binding (equally important)
- 20% tetramer integrity
- 10% safety (MSA conservation)

---

## Expected Results

Based on R175H test (top 5 rescues):

### Pass Rates (estimated)

- **DNA binding:** 85-95% pass (ΔΔG ≤ 2.0 kcal/mol)
- **Interface:** 90-98% pass (ΔΔG ≤ 2.5 kcal/mol)
- **Overall:** 80-90% rated "excellent" or "good"

### Filtering

- **Top-ranked rescues:** Mostly "excellent" (all dimensions good)
- **Mid-ranked rescues:** Mix of "good" and "acceptable"
- **Lower-ranked rescues:** Some "poor" (false positives filtered out)

### Discoveries

- Some high-stability rescues **improve** DNA binding (e.g., S95A: -1.15 kcal/mol)
- Most rescues **preserve** DNA binding and interface (ΔΔG = 0.0)
- ~10-20% may fail one dimension (filtered as "poor")

---

## Computational Details

### Runtime

**Per rescue:**
- WT complex binding: ~10 sec
- Build mutant: ~10 sec
- Mutant complex binding: ~10 sec
- **Total:** ~30 sec per rescue

**Per target:**
- R175H: ~1.5 hours (162 Pareto rescues)
- R248Q: ~1 hour (~120 rescues)
- R273H: ~1 hour (~120 rescues)
- Y220C: ~1 hour (~120 rescues)

**Optimization:** Results are cached! Re-running is instant.

### Resources

- **CPU:** Uses existing EvoEF2 (no GPU needed)
- **RAM:** 8-16 GB recommended
- **Storage:** ~500 MB for cache
- **Parallelization:** One rescue at a time (EvoEF2 limitation)

---

## Troubleshooting

### Error: "EvoEF2 binary not found"

Check `configs/scoring.yaml`:
```yaml
evoef2:
  binary: /path/to/EvoEF2  # Update this path
```

### Error: "Structure not found"

Ensure structures exist:
```bash
ls Data/raw/experimental_pdbs/1TSR.pdb  # DNA-bound
ls Data/raw/experimental_pdbs/3KMD.pdb  # Tetramer
```

### Slow performance

**Normal:** 30 sec/rescue × 162 rescues = ~1.5 hours
**If slower:** Check EvoEF2 installation and system load

### Missing functional scores

Check that `--functional-scoring` flag was used:
```bash
python -m experiments.run_design_rescues --targets R175H --functional-scoring
```

---

## Validation

### Sanity Checks

After running, verify:

```python
import pandas as pd

# Load results
df = pd.read_parquet('Data/processed/rescues/R175H/pareto.parquet')

# Check functional scores present
assert 'functional_score' in df.columns
assert 'ddg_binding' in df.columns
assert 'ddg_interface' in df.columns

# Check values reasonable
assert df['functional_score'].between(0, 1).all()
assert df['ddg_binding'].between(-10, 10).all()  # Typical range
assert df['overall_category'].isin(['excellent', 'good', 'acceptable', 'poor']).all()

print("✅ Validation passed!")
```

### Compare with Test Results

Expected for R175H top 5:
- S95A: functional_score = 0.750, ddg_binding = -1.15
- M133L: functional_score = 0.700, ddg_binding = 0.00
- All rated "excellent"

---

## Next Steps

### 1. Run Full Pipeline

```bash
# Run all targets with functional scoring
bash experiments/run_all_functional.sh
```

### 2. Analyze Results

```bash
# Generate comprehensive report
python experiments/analyze_functional_results.py

# Output:
# - Category distribution across all targets
# - Pass rates for DNA binding and interface
# - Top 10 rescues per target
# - Filtered-out rescues (false positives)
```

### 3. Generate Validation Report

Create manuscript-ready report showing:
- Stability-only vs functional optimization comparison
- False positive identification and removal
- Top candidates meet all biological requirements

### 4. Prepare for Experimental Validation

Select top 3-5 "excellent" rescues per target for:
- ΔT_m measurements (validate ΔΔG_folding)
- K_d measurements (validate ΔΔG_binding)
- SEC/AUC (validate ΔΔG_interface)

---

## Configuration

### Adjust Weights

Edit `configs/functional_scoring.yaml`:

```yaml
functional_score:
  weights:
    folding: 0.35      # Increase if stability is most critical
    dna_binding: 0.35  # Increase if DNA binding is priority
    interface: 0.20    # Increase if avoiding dominant-negative
    risk: 0.10         # Increase if safety is paramount
    # Total must sum to 1.0
```

### Adjust Thresholds

```yaml
dna_binding:
  thresholds:
    good: 0.0        # ΔΔG ≤ 0 (neutral or improved)
    acceptable: 2.0  # ΔΔG ≤ 2.0 (modest penalty)
    bad: 4.0         # ΔΔG > 4.0 (significantly disrupted)
```

---

## References

- **Initiative 1 Report:** `INITIATIVE1_COMPLETE.md`
- **Test Results:** `experiments/test_evoef2_functional_scoring.py`
- **Source Code:** `src/scoring/functional/`

---

*Last Updated: January 26, 2026*
*Pipeline Version: 2.0 (with functional scoring)*
