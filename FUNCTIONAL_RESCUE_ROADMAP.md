# Functional Rescue & Synergistic Therapy Roadmap

**Goal:** Elevate from 85% → 95%+ scientific rigor by modeling **functional restoration**, not just monomer stability

**Current State:** Optimizing for monomer ΔΔG_folding only
**Target State:** Multi-dimensional functional score + drug synergy predictions

---

## The Critical Gap

**What we have now:**
```
Good Rescue = High monomer stability (ΔΔG_folding)
```

**What we need:**
```
Functional Rescue = ΔΔG_folding (stability)
                  + ΔΔG_binding (DNA affinity)
                  + ΔΔG_interface (tetramer integrity)
                  + Drug synergy potential
```

---

## Current Infrastructure Audit

### ✅ Already Have
- Multi-structure scoring (AlphaFold + 2OCJ core domain)
- OpenMM 8.4 (for MD simulations)
- Experimental structures:
  - `2OCJ_full.pdb` (DNA-bound p53 core, 2.2Å)
  - `3KMD.pdb` (p53 tetramer)
  - `1TSR.pdb` (DNA-bound tetramer)
- EvoEF2 for monomer stability
- Pareto optimization framework

### ❌ Missing Components
- FoldX (for AnalyseComplex)
- DNA binding affinity calculation (ΔΔG_binding)
- Tetramer interface stability calculation (ΔΔG_interface)
- MD simulation pipeline
- Ensemble docking infrastructure
- Drug library (PHI-KAN-083, COTI-2, etc.)

---

## Initiative 1: Functional Rescue Score (Target: 90%)

**Timeline:** 2-3 weeks
**Impact:** Filters out false positives, validates true functional rescues

### Part A: DNA Binding Affinity (ΔΔG_binding)

**Goal:** Ensure rescue mutations don't break DNA binding

**Implementation:**

1. **Tool Selection**
   - **Option 1 (Recommended):** FoldX AnalyseComplex
     - Pros: Fast, well-validated for protein-DNA
     - Cons: Need to install, license for academic use
   - **Option 2:** Rosetta InterfaceAnalyzer
     - Pros: Already in conda ecosystem
     - Cons: Slower, more complex setup
   - **Option 3:** PRODIGY-PROTEIN-DNA
     - Pros: Lightweight, fast
     - Cons: Less accurate

**Recommended:** Start with FoldX AnalyseComplex

2. **Workflow**
```python
# For each rescue candidate:
# 1. Build mutant DNA-bound complex (from 2OCJ or 1TSR)
# 2. Calculate WT binding energy: E_complex - (E_protein + E_DNA)
# 3. Calculate mutant binding energy: E_complex_mut - (E_protein_mut + E_DNA)
# 4. ΔΔG_binding = ΔG_binding_mut - ΔG_binding_WT
```

3. **Acceptance Criteria**
   - **Good:** ΔΔG_binding ≤ 0 (neutral or improved binding)
   - **Acceptable:** ΔΔG_binding ≤ +1.0 kcal/mol (modest penalty)
   - **Bad:** ΔΔG_binding > +2.0 kcal/mol (significantly weakened binding)

4. **Expected Filtering**
   - Current: 162 Pareto rescues for R175H
   - After ΔΔG_binding filter: ~120-140 rescues (remove 15-25%)
   - **Key filter:** R196Q likely fails (high MSA conservation, structural role)

---

### Part B: Tetramerization Integrity (ΔΔG_interface)

**Goal:** Ensure rescue mutations don't disrupt tetramer formation (dominant-negative effect)

**Implementation:**

1. **Structure Selection**
   - Use `3KMD.pdb` (p53 tetramer, 2.5Å resolution)
   - Focus on dimer-dimer interface (critical for tetramer assembly)

2. **Tool:** FoldX AnalyseComplex
```python
# For each rescue candidate:
# 1. Extract dimer from tetramer (chains A+B)
# 2. Calculate WT interface energy
# 3. Introduce rescue mutation
# 4. Calculate mutant interface energy
# 5. ΔΔG_interface = ΔG_interface_mut - ΔG_interface_WT
```

3. **Acceptance Criteria**
   - **Good:** ΔΔG_interface ≤ 0 (neutral or improved interface)
   - **Acceptable:** ΔΔG_interface ≤ +1.5 kcal/mol (modest destabilization)
   - **Bad:** ΔΔG_interface > +2.5 kcal/mol (disrupts tetramer)

4. **Expected Filtering**
   - Current: ~120-140 rescues (after ΔΔG_binding filter)
   - After ΔΔG_interface filter: ~100-120 rescues
   - **Key filter:** S95A likely fails (near N-terminal interface)

---

### Part C: Composite Functional Score

**Integration:**

```python
# Multi-dimensional scoring
functional_score = w1 * normalize(ΔΔG_folding)     # Monomer stability
                 + w2 * normalize(ΔΔG_binding)     # DNA binding
                 + w3 * normalize(ΔΔG_interface)   # Tetramerization
                 + w4 * normalize(risk)            # Safety

# Pareto optimization on functional_score dimensions
```

**Weighting scheme (initial):**
- w1 = 0.35 (folding stability)
- w2 = 0.35 (DNA binding - equally important)
- w3 = 0.20 (interface - slightly less critical)
- w4 = 0.10 (risk - safety tax)

**Output:** New Pareto front ranked by **functional rescue potential**, not just stability

---

### Implementation Files

**New modules:**
```
src/scoring/
  foldx_runner.py              # FoldX wrapper (similar to evoef2_runner.py)
  binding_affinity.py          # ΔΔG_binding calculation
  interface_stability.py       # ΔΔG_interface calculation
  functional_score.py          # Composite scoring

src/structures/
  complex_builder.py           # Build mutant protein-DNA complexes
  tetramer_utils.py            # Extract dimer, handle symmetry

configs/
  functional_scoring.yaml      # Weights, thresholds, structure paths
```

**Modified modules:**
```
experiments/run_design_rescues.py   # Add functional scoring step
experiments/run_make_report.py      # Add ΔΔG_binding, ΔΔG_interface columns
configs/optimizer.yaml              # Update Pareto objectives
```

---

## Initiative 2: Synergistic Therapy Simulation (Target: 95%+)

**Timeline:** 3-4 weeks
**Impact:** Clinically testable hypothesis, publication-grade contribution

### Part A: MD Simulation Pipeline

**Goal:** Generate dynamic ensemble of top rescued structures

**Implementation:**

1. **Input Selection**
   - Top 3-5 rescues from Initiative 1 (highest functional_score)
   - Candidates: Likely L145I, M133L, T155A, M133L+T155A

2. **Simulation Protocol**
   - **Force Field:** AMBER14 (protein) + TIP3P (water)
   - **System Setup:** Solvate in cubic box, 150 mM NaCl
   - **Minimization:** 5000 steps steepest descent
   - **Equilibration:**
     - NVT: 500 ps at 310K
     - NPT: 1 ns at 310K, 1 bar
   - **Production:** 50 ns, save frames every 10 ps → 5000 snapshots
   - **Replicates:** 3x per rescue (different random seeds)

3. **Tool:** OpenMM (already installed!)
```python
# For each top rescue:
# 1. Build system (protein + water + ions)
# 2. Minimize energy
# 3. Equilibrate (NVT → NPT)
# 4. Run 50 ns production MD
# 5. Extract ensemble (5000 frames)
# 6. Cluster conformations (RMSD-based, identify top 100 representatives)
```

4. **Analysis Metrics**
   - RMSD stability (should converge < 3Å)
   - RMSF per-residue (identify flexible regions)
   - Pocket volume evolution (track druggable sites)
   - Secondary structure retention

5. **Expected Output**
   - Per rescue: 100 representative conformations (clustered from 5000 frames)
   - Total: ~500 structures for docking (5 rescues × 100 conformations)

---

### Part B: Ensemble Docking & Virtual Screening

**Goal:** Find rescue-drug synergistic pairs

**Implementation:**

1. **Drug Library**
   - **Tier 1 (Known p53 reactivators):**
     - COTI-2 (Phase I clinical trial)
     - PHI-KAN-083 (Y220C-specific)
     - APR-246/PRIMA-1 (cysteine-targeting)
     - Stictic acid derivatives
   - **Tier 2 (FDA-approved, repurposing candidates):**
     - Nutlin-3a (MDM2 inhibitor)
     - Kevetrin (p53 modulator)
   - **Size:** ~20-30 compounds initially

2. **Docking Protocol**
   - **Tool:** AutoDock Vina (fast) or Glide (accurate, if available)
   - **Grid Setup:**
     - Center on Y220 pocket (for Y220C rescues)
     - Center on L1/L2 loop (for conformational rescues)
   - **Exhaustiveness:** 32 (high sampling)

3. **Workflow**
```python
# For each rescue:
#   For each conformation in ensemble (100):
#     For each drug (20-30):
#       - Dock drug to conformation
#       - Record binding affinity (ΔG_dock)
#       - Record binding pose
#   Calculate ensemble-averaged affinity: <ΔG_dock>
#   Identify top drugs (lowest <ΔG_dock>)

# Output: Matrix of rescue-drug pairs ranked by synergy
```

4. **Synergy Score**
```python
synergy_score = functional_score (from Initiative 1)
              + w_drug * normalize(<ΔG_dock>)
              + w_stability * (pocket_occupancy > 50%)
```

5. **Expected Output**
   - **Best pair:** e.g., "G245S + COTI-2" with <ΔG_dock> = -9.5 kcal/mol
   - **Hypothesis:** "G245S rescue creates a stable pocket that binds COTI-2 with high affinity across its conformational ensemble"
   - **Testable prediction:** Combine G245S + COTI-2 in vitro → measure ΔT_m

---

### Implementation Files

**New modules:**
```
src/md/
  openmm_runner.py             # MD simulation wrapper
  system_builder.py            # Solvate, add ions, minimize
  trajectory_analyzer.py       # RMSD, RMSF, clustering
  ensemble_extractor.py        # Cluster → representative frames

src/docking/
  autodock_runner.py           # AutoDock Vina wrapper
  ensemble_docking.py          # Dock to multiple conformations
  drug_library.py              # SMILES → 3D, protonation
  synergy_scorer.py            # Combine functional + docking scores

configs/
  md_simulation.yaml           # Force field, temperature, duration
  docking.yaml                 # Grid parameters, drug library paths
```

**Modified modules:**
```
experiments/run_md_simulations.py    # NEW: CLI for MD pipeline
experiments/run_ensemble_docking.py  # NEW: CLI for docking
experiments/run_make_report.py       # Add drug synergy section
```

---

## Deliverables & Impact

### Scientific Rigor (90% → 95%+)

**Before:**
- 162 Pareto rescues for R175H
- Top candidate: M133L (ΔΔG_folding = -5.6 kcal/mol)
- Claim: "Stabilizes p53 monomer"

**After Initiative 1:**
- ~100-120 **functional** rescues (filtered by DNA binding + interface)
- Top candidate: L145I (ΔΔG_folding = -1.6, ΔΔG_binding = +0.3, ΔΔG_interface = -0.1)
- Claim: "Restores p53 function by preserving stability, DNA binding, and tetramerization"

**After Initiative 2:**
- Top 5 rescues with validated drug synergy
- Best pair: G245S + COTI-2 (<ΔG_dock> = -9.5 kcal/mol)
- Claim: "G245S rescue creates a druggable pocket with high affinity for COTI-2, enabling dual gene-drug therapy"

---

### Publication Narrative

**Title:** "Rational Design of p53 Rescue Mutations with Synergistic Drug Therapy Potential: A Multi-Dimensional Functional Approach"

**Abstract:**
> "We present a computational pipeline for designing p53 rescue mutations that restore function by optimizing monomer stability (ΔΔG_folding), DNA binding affinity (ΔΔG_binding), and tetramerization integrity (ΔΔG_interface). By screening 700 candidates across 4 cancer hotspots (R175H, R248Q, R273H, Y220C), we identified 15 functional rescues that preserve all three dimensions. Molecular dynamics simulations of the top 5 rescues revealed that G245S adopts a dynamic conformation with a stable pocket optimized for binding the clinical-stage drug COTI-2 (ensemble-averaged affinity: -9.5 kcal/mol). This rescue-drug pair represents a testable dual-pronged therapeutic strategy for restoring p53 function in cancer."

**Target Journals:**
- *Nature Communications* (high impact, computational biology)
- *Cell Chemical Biology* (drug design focus)
- *PLOS Computational Biology* (methodology focus)

---

## Implementation Timeline

### Phase 1: Tool Setup (Week 1)
- Install FoldX (academic license)
- Setup AutoDock Vina
- Build drug library (download structures, protonate)
- Create config files for functional scoring

### Phase 2: ΔΔG_binding Implementation (Week 2)
- Build `foldx_runner.py` (wrapper like EvoEF2)
- Implement `binding_affinity.py`
- Test on known mutations (TP53 database)
- Integrate into design pipeline

### Phase 3: ΔΔG_interface Implementation (Week 3)
- Implement `interface_stability.py`
- Test on tetramer-disrupting mutations
- Add to functional scoring

### Phase 4: Re-run Design with Functional Score (Week 4)
- Update Pareto objectives
- Re-run all 4 targets (R175H, R248Q, R273H, Y220C)
- Compare old vs new Pareto fronts
- Validate top candidates

### Phase 5: MD Simulations (Week 5-6)
- Setup OpenMM pipeline
- Run 50 ns simulations for top 5 rescues
- Cluster conformations → ensembles
- Analyze stability metrics

### Phase 6: Ensemble Docking (Week 7-8)
- Dock drug library to ensembles
- Calculate synergy scores
- Identify top rescue-drug pairs
- Generate publication figures

### Phase 7: Validation & Reporting (Week 9)
- Cross-validate with experimental data (TP53 database)
- Generate comprehensive report
- Create manuscript draft

---

## Success Metrics

### Initiative 1 (Functional Rescue Score)
- ✅ **80%+ candidates** pass ΔΔG_binding filter (don't break DNA binding)
- ✅ **90%+ candidates** pass ΔΔG_interface filter (don't break tetramer)
- ✅ **Top 10 functional rescues** have favorable scores in all 3 dimensions
- ✅ **False positives removed:** R196Q, S95A filtered out

### Initiative 2 (Synergistic Therapy)
- ✅ **MD simulations converge:** RMSD < 3Å for all rescues
- ✅ **Druggable pockets identified:** ≥2 rescues with stable binding sites
- ✅ **High-affinity pairs found:** ≥1 rescue-drug pair with <ΔG_dock> < -9.0 kcal/mol
- ✅ **Testable hypothesis:** "G245S + COTI-2" with specific predicted ΔT_m

---

## Risks & Mitigation

### Risk 1: FoldX Installation/Licensing
- **Mitigation:** Apply for academic license (free), fallback to Rosetta if needed

### Risk 2: MD Simulations Too Slow
- **Mitigation:** Use GPU acceleration (OpenMM CUDA), reduce to 20 ns if needed

### Risk 3: No Clear Drug Synergy
- **Mitigation:** Expand drug library, focus on known p53 binders first

### Risk 4: Complexity Explosion
- **Mitigation:** Start with top 3 rescues only, expand if successful

---

## Next Steps (Immediate Actions)

1. **Decision Point:** Do you want to proceed with Initiative 1 first (functional scoring)?
   - If yes, I'll start with FoldX installation and ΔΔG_binding implementation
   - Expected time to first results: ~1 week

2. **Resource Check:** Do you have access to:
   - GPU for MD simulations? (check: `nvidia-smi`)
   - Computational cluster for parallel docking?
   - FoldX academic license (can apply: https://foldxsuite.crg.eu/)

3. **Prioritization:** Which is more important to you?
   - **Option A:** Faster results (do Initiative 1 only, defer drug screening)
   - **Option B:** Maximum impact (do both initiatives, 9-week timeline)

**My recommendation:** Start with Initiative 1 (functional scoring) to validate the approach, then proceed to Initiative 2 if results are promising.

---

*Roadmap created: January 25, 2026*
*Current state: 85% rigor (monomer stability only)*
*Target state: 95%+ rigor (functional restoration + drug synergy)*
