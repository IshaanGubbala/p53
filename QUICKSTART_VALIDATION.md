# Quick Start: Running the Validation Suite

This guide shows you how to reproduce all three advanced validations on your rescue candidates.

## Prerequisites

```bash
# Required Python packages (already in environment)
pip install biopython pandas numpy matplotlib seaborn pyyaml

# Required PDB structures (already downloaded)
Data/raw/experimental_pdbs/
├── 1TSR.pdb  # p53-DNA complex
├── 1TUP.pdb  # p53-DNA complex (alternative)
├── 3KMD.pdb  # p53 core tetramer + DNA
└── 1OLG.pdb  # p53 tetramerization domain
```

---

## Step 1: DNA Binding Analysis (10 minutes)

**What it does:** Identifies which rescue mutations contact DNA (potential functional disruption)

```bash
python -m experiments.run_dna_binding_analysis
```

**Output:**
- `reports/dna_binding/*_dna_binding.csv` - Per-target results
- `reports/dna_binding/*_dna_binding_summary.json` - Summary statistics
- `reports/dna_binding/dna_contacts.json` - DNA-contacting residues from 1TSR

**Expected results:**
- ~0-5 rescues per target contact DNA
- Top candidates (M133L, A189S, A189G) show 0.00 DNA binding risk

**Interpretation:**
- `dna_binding_risk == "safe"` → No DNA contact, safe to proceed
- `dna_binding_risk == "critical"` → Disrupts DNA binding, exclude from consideration

---

## Step 2: Tetramer Interface Analysis (5 minutes)

**What it does:** Identifies which rescue mutations disrupt p53 tetramerization (dominant-negative risk)

```bash
python -m experiments.run_tetramer_interface_analysis
```

**Output:**
- `reports/tetramer_interface/*_tetramer_interface.csv` - Per-target results
- `reports/tetramer_interface/*_tetramer_summary.json` - Summary statistics

**Expected results:**
- ~1-5 rescues per target at interface positions
- Top candidates show 0.00 tetramer impact
- Only S95A shows medium risk (not in top candidates)

**Interpretation:**
- `tetramer_risk == "safe"` → No interface disruption
- `tetramer_risk == "critical"` → Disrupts tetramerization, exclude

---

## Step 3: Sensitivity Analysis (15 minutes)

**What it does:** Tests if rescue rankings are robust across different risk weighting schemes

```bash
python -m experiments.run_sensitivity_analysis
```

**Output:**
- `reports/sensitivity/*_sensitivity_rankings.csv` - Per-target rankings
- `reports/sensitivity/*_sensitivity_summary.json` - Robustness metrics

**Expected results:**
- Top 3 rescues IDENTICAL across all 4 weighting profiles
- 110-142 highly robust rescues per target
- 106 universal winners across all targets

**Interpretation:**
- `robustness_score == 1.0` → Highly robust (top 10 in all profiles)
- `robustness_score >= 0.75` → Robust (top 10 in ≥3/4 profiles)
- `robustness_score < 0.50` → Sensitive to weights (use with caution)

---

## Step 4: Generate Summary Figures (30 seconds)

**What it does:** Creates publication-ready visualizations

```bash
python visualization/create_validation_summary_figures.py
```

**Output:**
- `reports/figures/validation_summary.pdf` - Comprehensive validation summary
- `reports/figures/robustness_heatmap.pdf` - Rankings across profiles
- `reports/figures/validation_workflow.pdf` - Pipeline diagram

---

## Full Pipeline (30 minutes)

Run all three validations sequentially:

```bash
# DNA binding analysis
python -m experiments.run_dna_binding_analysis

# Tetramer interface analysis
python -m experiments.run_tetramer_interface_analysis

# Sensitivity analysis
python -m experiments.run_sensitivity_analysis

# Generate figures
python visualization/create_validation_summary_figures.py
```

**Check results:**
```bash
# DNA binding findings
cat reports/dna_binding/DNA_BINDING_FINDINGS.md

# Tetramer findings
cat reports/tetramer_interface/TETRAMER_INTERFACE_FINDINGS.md

# Sensitivity findings
cat reports/sensitivity/SENSITIVITY_FINDINGS.md

# Complete validation summary
cat reports/VALIDATION_COMPLETE.md
```

---

## Interpreting Results

### ✅ Safe Rescue Candidate Checklist

A rescue is considered **safe and robust** if:

1. **DNA binding:** `dna_binding_risk == "safe"` or `"low"`
   - Impact score < 1.0
   - Not at DNA-contacting positions

2. **Tetramerization:** `tetramer_risk == "safe"` or `"low"`
   - Impact score < 1.5
   - Not at interface positions

3. **Robustness:** `robustness_score >= 0.75`
   - Appears in top 10 for ≥3/4 weighting profiles
   - Rank variance < 10

4. **Pareto optimal:** `is_pareto == True`
   - On stability-safety Pareto frontier

### 🏆 Gold Standard Candidates

Our top 3 rescues (M133L, A189S, A189G) meet ALL criteria:
- ✅ DNA safe (risk: 0.00)
- ✅ Tetramer safe (risk: 0.00)
- ✅ Highly robust (score: 1.00)
- ✅ Pareto optimal (True)

---

## Customization

### Analyze Your Own Rescues

If you have rescue candidates in a different format:

```python
import pandas as pd
from src.eval.dna_binding_analysis import analyze_rescue_dna_binding
from src.eval.tetramer_interface_analysis import analyze_rescue_tetramer_impact
from src.eval.sensitivity_analysis import sensitivity_analysis

# Load your candidates
my_rescues = pd.read_csv("my_rescues.csv")

# DNA binding analysis
dna_results = analyze_rescue_dna_binding(
    my_rescues,
    dna_complex_pdb=Path("Data/raw/experimental_pdbs/1TSR.pdb"),
    protein_chains=["A"],
    dna_chains=["E", "F"],
)

# Tetramer analysis
tetramer_results = analyze_rescue_tetramer_impact(
    my_rescues,
    tetramer_pdb=Path("Data/raw/experimental_pdbs/3KMD.pdb"),
)

# Sensitivity analysis (requires risk component columns)
sensitivity_results = sensitivity_analysis(
    my_rescues,
    top_n=10,
)
```

### Change Validation Thresholds

Edit the source files to adjust risk thresholds:

**DNA binding:** `src/eval/dna_binding_analysis.py`
```python
# Line 232-248: Adjust impact score thresholds
if impact_score >= 4:      # Change from 4 to your threshold
    risk_level = "critical"
elif impact_score >= 2.5:  # Change from 2.5
    risk_level = "high"
```

**Tetramer:** `src/eval/tetramer_interface_analysis.py`
```python
# Line 193-206: Adjust interface risk thresholds
if impact_score >= 4:      # Change from 4
    risk_level = "critical"
```

---

## Troubleshooting

### Issue: "No candidates file for [target]"
**Solution:** Ensure you've run rescue design first:
```bash
python -m experiments.run_design_rescues --targets R175H R248Q R273H Y220C
```

### Issue: "PDB file not found: 1TSR.pdb"
**Solution:** Download structure manually:
```bash
curl -o Data/raw/experimental_pdbs/1TSR.pdb \
  "https://files.rcsb.org/download/1TSR.pdb"
```

### Issue: "Missing risk columns"
**Solution:** Candidates file needs `risk_components` column (JSON format):
```json
{
  "functional": 0.0,
  "conservation": 0.0,
  "burial": 0.0,
  "msa_conservation": 0.0
}
```

### Issue: BioPython import error
**Solution:** Install BioPython:
```bash
pip install biopython
```

---

## Advanced: Adding New Validations

### Template for New Validation Module

```python
# src/eval/my_validation.py

from pathlib import Path
import pandas as pd
from src.core.logging import get_logger

logger = get_logger(__name__)

def analyze_my_validation(
    rescues_df: pd.DataFrame,
    **kwargs
) -> pd.DataFrame:
    """
    Analyze rescues for [your validation].

    Args:
        rescues_df: DataFrame with rescue candidates
        **kwargs: Additional parameters

    Returns:
        DataFrame with validation columns added
    """
    results = []

    for idx, row in rescues_df.iterrows():
        rescue_mutations = row["rescue_mutations"]

        # Your validation logic here
        risk_score = compute_risk(rescue_mutations)

        results.append({
            "rescue_mutations": rescue_mutations,
            "my_validation_risk": risk_score,
        })

    # Merge with original
    results_df = pd.DataFrame(results)
    return rescues_df.merge(results_df, on="rescue_mutations", how="left")
```

### Add to Pipeline

```python
# experiments/run_my_validation.py

from src.eval.my_validation import analyze_my_validation

def run(args, configs):
    rescues = pd.read_parquet("Data/processed/rescues/R175H/candidates.parquet")

    validated = analyze_my_validation(rescues)

    validated.to_csv("reports/my_validation/R175H_results.csv", index=False)

    return 0
```

---

## Performance Tips

### Speed Up DNA Binding Analysis
Use smaller `top_n` in analysis:
```python
dna_binding_df = analyze_rescue_dna_binding(
    pareto_df.head(20),  # Analyze only top 20 instead of 50
    dna_complex_pdb=dna_complex_pdb,
)
```

### Parallel Processing
Run validations in parallel:
```bash
# Terminal 1
python -m experiments.run_dna_binding_analysis &

# Terminal 2
python -m experiments.run_tetramer_interface_analysis &

# Terminal 3
python -m experiments.run_sensitivity_analysis &

# Wait for all to complete
wait
```

### Cache Results
All analyses save CSV outputs. To avoid re-running:
```python
# Check if already run
if Path("reports/dna_binding/R175H_dna_binding.csv").exists():
    print("Already analyzed, loading cached results")
    df = pd.read_csv("reports/dna_binding/R175H_dna_binding.csv")
else:
    # Run analysis
    df = analyze_rescue_dna_binding(...)
```

---

## Expected Runtime

| Validation | Time | Bottleneck |
|-----------|------|------------|
| DNA binding | ~10 min | BioPython PDB parsing (1-2 sec per structure) |
| Tetramer | ~5 min | Interface distance calculations |
| Sensitivity | ~15 min | Re-ranking under multiple profiles |
| Figures | ~30 sec | Matplotlib rendering |
| **Total** | **~30 min** | - |

---

## Output File Reference

```
reports/
├── dna_binding/
│   ├── DNA_BINDING_FINDINGS.md           # Summary document
│   ├── R175H_dna_binding.csv             # Per-rescue DNA binding risk
│   ├── R248Q_dna_binding.csv
│   ├── R273H_dna_binding.csv
│   ├── Y220C_dna_binding.csv
│   ├── all_targets_dna_binding.csv       # Consolidated results
│   └── dna_contacts.json                 # DNA-contacting residues
│
├── tetramer_interface/
│   ├── TETRAMER_INTERFACE_FINDINGS.md    # Summary document
│   ├── R175H_tetramer_interface.csv      # Per-rescue tetramer risk
│   ├── R248Q_tetramer_interface.csv
│   ├── R273H_tetramer_interface.csv
│   ├── Y220C_tetramer_interface.csv
│   └── all_targets_tetramer_interface.csv
│
├── sensitivity/
│   ├── SENSITIVITY_FINDINGS.md           # Summary document
│   ├── R175H_sensitivity_rankings.csv    # Per-profile rankings
│   ├── R248Q_sensitivity_rankings.csv
│   ├── R273H_sensitivity_rankings.csv
│   ├── Y220C_sensitivity_rankings.csv
│   └── *_sensitivity_summary.json        # Robustness metrics
│
├── figures/
│   ├── validation_summary.pdf            # Comprehensive overview
│   ├── robustness_heatmap.pdf            # Ranking stability
│   └── validation_workflow.pdf           # Pipeline diagram
│
├── VALIDATION_COMPLETE.md                # Master summary
├── EXECUTIVE_SUMMARY.md                  # Stakeholder summary
└── EXPERIMENTAL_VALIDATION_PROTOCOL.md   # Experimental plan
```

---

## Success! What's Next?

After running all validations:

1. **Review findings:**
   ```bash
   cat reports/VALIDATION_COMPLETE.md
   ```

2. **Check top candidates:**
   ```bash
   cat reports/EXECUTIVE_SUMMARY.md
   ```

3. **Plan experiments:**
   ```bash
   cat reports/EXPERIMENTAL_VALIDATION_PROTOCOL.md
   ```

4. **Share figures:**
   ```bash
   open reports/figures/validation_summary.pdf
   ```

**Your rescue candidates are now comprehensively validated and ready for experimental testing!**

---

*Quickstart Guide v1.0*
*p53 Rescue Mutation Validation Suite*
*Last Updated: January 25, 2026*
