# Initiative 2: Quick Start Guide

**Purpose:** Step-by-step commands to run Initiative 2 from start to finish
**Prerequisites:** Tools installed (GROMACS, Vina, Open Babel)
**Duration:** ~2-3 weeks for full pipeline

---

## 🚀 One-Command Full Pipeline

```bash
# Run everything (top 10 rescues × 6 drugs = 60 pairs)
python experiments/run_initiative2_pipeline.py \
  --rescues top10 \
  --drugs all \
  --md-duration 50
```

**This will run all 4 phases automatically:**
1. MD simulations (50 ns × 10 rescues)
2. Trajectory analysis (extract 10 conformers per rescue)
3. Ensemble docking (10 rescues × 6 drugs × 10 conformers = 600 dockings)
4. Synergy scoring (rank all 60 rescue-drug pairs)

**Time estimate:**
- With GPU: ~2-3 days
- With CPU: ~2-3 weeks

---

## 📋 Step-by-Step (Manual Mode)

### Step 0: Download Drugs ✅

```bash
python src/md/download_drugs.py
```

**Output:** 6 drug structures in `Data/raw/drugs/structures/`

**Status:** ✅ Already done (you ran this successfully!)

---

### Step 1: MD Simulations (Per Rescue)

```bash
# Single rescue (for testing)
python src/md/setup_md_simulation.py \
  --rescue "A189S,M133L,S95T" \
  --target R175H \
  --duration 10

# Full 50 ns production run
python src/md/setup_md_simulation.py \
  --rescue "A189S,M133L,S95T" \
  --target R175H \
  --duration 50
```

**Time per rescue:**
- 10 ns pilot: ~30 min (GPU) or ~3 hours (CPU)
- 50 ns production: ~2-3 hours (GPU) or ~1-2 days (CPU)

**Output:**
- `Data/processed/md_simulations/A189S_M133L_S95T/md.xtc` (trajectory)
- `Data/processed/md_simulations/A189S_M133L_S95T/md.tpr` (topology)

---

### Step 2: Trajectory Analysis (Per Rescue)

```bash
python src/md/analyze_trajectory.py \
  --rescue "A189S_M133L_S95T" \
  --n-clusters 10
```

**Time:** ~5-10 minutes per rescue

**Output:**
- `Data/processed/md_simulations/A189S_M133L_S95T/analysis/rmsd.png`
- `Data/processed/md_simulations/A189S_M133L_S95T/analysis/rmsf.png`
- `Data/processed/md_simulations/A189S_M133L_S95T/analysis/clusters.png`
- `Data/processed/md_simulations/A189S_M133L_S95T/analysis/conformers/` (10 PDBs)

**Success check:**
- RMSD should converge < 3 Å
- RMSF should show expected flexible regions (loops high, core low)

---

### Step 3: Ensemble Docking (Per Rescue-Drug Pair)

```bash
python src/docking/run_ensemble_docking.py \
  --rescue "A189S_M133L_S95T" \
  --drug "COTI-2"
```

**Time:** ~30-60 minutes per pair (10 conformers × ~3-5 min each)

**Output:**
- `Data/processed/docking/results/A189S_M133L_S95T/COTI-2/docking_results.json`
- Contains: consensus_affinity, affinity_std, best_affinity

**Repeat for all drugs:**
```bash
for drug in COTI-2 APR-246 Nutlin-3a PHI-KAN-083 PRIMA-1 CP-31398; do
  python src/docking/run_ensemble_docking.py \
    --rescue "A189S_M133L_S95T" \
    --drug "$drug"
done
```

---

### Step 4: Synergy Scoring (All Pairs)

```bash
python src/docking/calculate_synergy.py
```

**Time:** ~1 minute

**Output:**
- `Data/processed/synergy/scores/all_pairs_ranked.csv` (all 60 pairs)
- `Data/processed/synergy/scores/superstar_pairs.csv` (synergy > 0.8)

**Columns:**
- `synergy_score`: Overall score (0-1, higher = better)
- `category`: superstar/promising/additive/poor
- `consensus_affinity`: Median binding affinity (kcal/mol)
- `functional_score`: From Initiative 1 (0-1)

---

## 🎯 Pilot Test (Fastest Way to Validate)

Run a quick pilot to test the full pipeline:

```bash
# Step 1: Run 10 ns MD on top rescue
python src/md/setup_md_simulation.py \
  --rescue "A189S,M133L,S95T" \
  --target R175H \
  --duration 10

# Step 2: Analyze trajectory
python src/md/analyze_trajectory.py \
  --rescue "A189S_M133L_S95T" \
  --n-clusters 10

# Step 3: Dock with one drug
python src/docking/run_ensemble_docking.py \
  --rescue "A189S_M133L_S95T" \
  --drug "COTI-2"

# Step 4: Calculate synergy (will only have 1 pair)
python src/docking/calculate_synergy.py
```

**Total time:** ~1-2 hours (GPU) or ~4-5 hours (CPU)

**If successful:** Scale to full 50 ns + all 10 rescues + all 6 drugs

---

## 📊 Monitoring Progress

### Check MD simulation progress
```bash
# View current RMSD (should converge)
cd Data/processed/md_simulations/A189S_M133L_S95T
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
# Select: 4 (Backbone)

# Plot RMSD
python -c "
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('rmsd.xvg', comments=['#', '@'])
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')
plt.savefig('rmsd_quick.png')
"
```

### Check trajectory analysis output
```bash
ls -lh Data/processed/md_simulations/A189S_M133L_S95T/analysis/
# Should see: rmsd.png, rmsf.png, clusters.png, conformers/
```

### Check docking results
```bash
cat Data/processed/docking/results/A189S_M133L_S95T/COTI-2/docking_results.json
# Should show: consensus_affinity, n_conformers, successful_dockings
```

### Check synergy scores
```bash
head -20 Data/processed/synergy/scores/all_pairs_ranked.csv
# Top rows = best rescue-drug pairs
```

---

## 🔧 Troubleshooting

### Issue: GROMACS not found
```bash
# Check installation
gmx --version

# If not found, install:
conda install -c conda-forge gromacs
```

### Issue: MD simulation too slow
**Solution:** Use shorter pilot (10 ns) or cloud GPU

```bash
# Launch AWS p3.2xlarge instance
# SSH into instance
# Run MD there (~10× faster)
```

### Issue: Docking fails (Vina not found)
```bash
# Check installation
vina --version

# If not found:
conda install -c conda-forge autodock-vina
```

### Issue: Can't convert SDF to PDBQT
**Solution:** Use RDKit to convert SDF to PDB first

```python
from rdkit import Chem
mol = Chem.SDMolSupplier('drug.sdf')[0]
Chem.MolToPDBFile(mol, 'drug.pdb')
```

### Issue: Trajectory analysis crashes
**Solution:** Check if trajectory exists and is readable

```bash
# Check file size (should be ~1-5 GB for 50 ns)
ls -lh Data/processed/md_simulations/*/md.xtc

# Try loading with MDAnalysis
python -c "
import MDAnalysis as mda
u = mda.Universe('md.tpr', 'md.xtc')
print(f'Frames: {len(u.trajectory)}')
"
```

---

## 🎯 Success Criteria

### After MD Simulations
- ✅ All 10 rescues have `md.xtc` files (~2-5 GB each)
- ✅ Mean RMSD < 3 Å (converged)
- ✅ No NaN or Inf values in trajectory

### After Trajectory Analysis
- ✅ All rescues have 10 conformer PDBs
- ✅ RMSD plots show convergence
- ✅ RMSF plots show realistic flexibility

### After Docking
- ✅ All 60 rescue-drug pairs have `docking_results.json`
- ✅ Consensus affinity < -5 kcal/mol (reasonable binding)
- ✅ Affinity std < 2 kcal/mol (consistent across ensemble)

### After Synergy Scoring
- ✅ `all_pairs_ranked.csv` has 60 rows
- ✅ ≥ 3 "superstar" pairs (synergy > 0.8)
- ✅ Top pairs are DNA-safe (ddg_binding ≤ 2.0)

---

## 📁 Expected Output Structure

```
Data/
├── raw/drugs/
│   └── structures/
│       ├── COTI-2.sdf ✅
│       ├── APR-246.sdf ✅
│       ├── Nutlin-3a.sdf ✅
│       ├── PHI-KAN-083.sdf ✅
│       ├── PRIMA-1.sdf ✅
│       └── CP-31398.sdf ✅
│
├── processed/
│   ├── md_simulations/
│   │   ├── A189S_M133L_S95T/
│   │   │   ├── md.xtc           # Trajectory
│   │   │   ├── md.tpr           # Topology
│   │   │   └── analysis/
│   │   │       ├── rmsd.png
│   │   │       ├── rmsf.png
│   │   │       ├── clusters.png
│   │   │       └── conformers/  # 10 PDBs
│   │   └── [9 more rescues...]
│   │
│   ├── docking/
│   │   └── results/
│   │       ├── A189S_M133L_S95T/
│   │       │   ├── COTI-2/
│   │       │   │   ├── docking_results.json
│   │       │   │   └── [10 docked PDBQTs]
│   │       │   └── [5 more drugs...]
│   │       └── [9 more rescues...]
│   │
│   └── synergy/
│       └── scores/
│           ├── all_pairs_ranked.csv     # ⭐ MAIN RESULT
│           └── superstar_pairs.csv      # ⭐ TOP CANDIDATES
```

---

## ⏱️ Time Estimates

### Pilot Run (1 rescue, 1 drug, 10 ns)
- With GPU: ~1-2 hours total
- With CPU: ~4-6 hours total

### Full Run (10 rescues, 6 drugs, 50 ns)
- With GPU: ~2-3 days
  - MD: ~1.5 days (parallel)
  - Analysis: ~2 hours
  - Docking: ~12 hours (parallel)
  - Synergy: ~1 minute

- With CPU: ~2-3 weeks
  - MD: ~2 weeks (sequential)
  - Analysis: ~2 hours
  - Docking: ~2 days (sequential)
  - Synergy: ~1 minute

### With 4 GPUs (Parallel)
- MD: ~12 hours (4 rescues in parallel)
- Docking: ~3 hours (multiple pairs in parallel)
- **Total: ~1 day**

---

## 🎉 Final Deliverable

After completion, you'll have:

1. **`all_pairs_ranked.csv`** - All 60 rescue-drug pairs ranked by synergy
2. **`superstar_pairs.csv`** - Top candidates for experimental validation
3. **Full docking results** - Binding modes and affinities
4. **MD trajectories** - Conformational ensembles for each rescue
5. **Analysis plots** - RMSD, RMSF, clustering for all rescues

**Next step:** Analyze top 5 superstar pairs, generate figures, write report!

---

*Quick start guide prepared: January 26, 2026*
*Ready to achieve 95%+ scientific rigor*
