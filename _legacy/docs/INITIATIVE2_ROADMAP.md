# Initiative 2: Molecular Dynamics + Drug Synergy - Roadmap

**Goal:** Move from 90% → 95%+ scientific rigor via dynamic ensemble modeling and drug synergy prediction
**Status:** Planning Phase
**Timeline:** ~2-3 weeks for full implementation
**Date Started:** January 26, 2026

---

## Executive Summary

With Initiative 1 complete (functional scoring across stability, DNA binding, and interface), we now elevate the pipeline to model:

1. **Dynamic protein ensembles** via molecular dynamics (MD) simulations
2. **Drug binding** to rescue mutants (ensemble docking)
3. **Rescue-drug synergy** to identify optimal combination therapies

**Scientific upgrade:**
- **Before (Initiative 1):** Static structures, single conformer analysis
- **After (Initiative 2):** Dynamic ensembles, conformational heterogeneity, drug synergy

**Clinical impact:** Move from "Here are functional rescues" to "Here are rescue-drug pairs with predicted synergistic activity"

---

## Candidate Selection

### Top 28 Unique Rescues (from 749 total)

**Selection criteria:** Top 10 per target, deduplicated across targets

**Top 5 by Functional Score:**
1. **A189S,M133L,S95T** (R175H) - functional_score: 0.895
2. **M133L,R196Q,S95T** (R248Q) - functional_score: 0.883
3. **R196Q,S215A,S95T** (R273H) - functional_score: 0.880
4. **C229A,R196H,S95T** (R273H) - functional_score: 0.879
5. **R196H,S95T** (R273H/Y220C) - functional_score: 0.874

**Key observations:**
- All top 5 contain S95T (master rescue position)
- Mix of double and triple mutants
- Represent all 4 cancer hotspots

### Initial Focus: Top 10 for Pilot Study

**Recommendation:** Start with top 10 unique rescues for pilot MD/docking study
- Manageable computational cost (~100-200 GPU hours)
- Sufficient diversity to validate approach
- Can scale to all 28 if pilot succeeds

---

## Phase 1: MD Simulation Infrastructure (Week 1)

### Objective
Generate conformational ensembles for top 10 rescues via 50 ns MD simulations

### Tools & Requirements

**MD Engine:** GROMACS 2023+
- Fast, well-validated, GPU-accelerated
- Excellent for protein simulations
- Strong community support

**Force Field:** AMBER99SB-ILDN or CHARMM36m
- Recommendation: AMBER99SB-ILDN (well-tested for p53)
- Includes improvements for intrinsically disordered proteins

**Water Model:** TIP3P
- Standard explicit solvent model
- Good balance of accuracy and speed

**Computational Requirements:**
- GPU: NVIDIA RTX 3080 or better (or cloud GPU: AWS p3.2xlarge)
- Time: ~10-12 hours per rescue (50 ns @ ~4 ns/hour on RTX 3080)
- Storage: ~5 GB per trajectory (snapshots every 10 ps)
- Total: ~100-120 GPU hours for 10 rescues

### Workflow

#### 1.1 Structure Preparation
```bash
# For each rescue mutation:
# - Build mutant model from AlphaFold p53 structure
# - Add hydrogen atoms (pdb2gmx)
# - Define simulation box (10 Å padding)
# - Solvate with TIP3P water
# - Add ions (0.15 M NaCl) for charge neutralization
```

**Input:** Mutant PDB from EvoEF2 BuildMutant
**Output:** GROMACS topology (.top) + coordinate (.gro) files

#### 1.2 Energy Minimization
```bash
# Steepest descent (5000 steps) + Conjugate gradient
# Goal: Remove steric clashes, Fmax < 100 kJ/mol/nm
```

**Duration:** ~10 minutes per structure
**Output:** Minimized structure

#### 1.3 Equilibration
```bash
# NVT equilibration (100 ps, 300 K, Berendsen thermostat)
# NPT equilibration (100 ps, 1 bar, Parrinello-Rahman barostat)
# Restraints on protein heavy atoms (1000 kJ/mol/nm²)
```

**Duration:** ~30 minutes per structure
**Purpose:** Stabilize solvent and temperature/pressure

#### 1.4 Production MD
```bash
# 50 ns simulation, 300 K, 1 bar
# 2 fs timestep (LINCS constraints on bonds)
# Save frames every 10 ps (5000 snapshots)
# PME electrostatics, 10 Å cutoffs
```

**Duration:** ~10-12 hours per rescue
**Output:** Trajectory file (.xtc), 5000 conformers per rescue

#### 1.5 Trajectory Analysis
```bash
# RMSD, RMSF, Rg (radius of gyration)
# Secondary structure evolution (DSSP)
# Cluster analysis (gromos method, 2 Å cutoff)
# Extract representative conformers (top 10 clusters)
```

**Output:**
- MD metrics (convergence, stability)
- 10 representative conformers per rescue
- Total: 100 conformers for docking (10 rescues × 10 conformers)

### File Structure
```
Data/processed/md_simulations/
├── structures/                  # Input structures
│   ├── A189S_M133L_S95T.pdb
│   ├── M133L_R196Q_S95T.pdb
│   └── ...
├── topologies/                  # GROMACS topologies
│   ├── A189S_M133L_S95T.top
│   └── ...
├── trajectories/               # MD trajectories
│   ├── A189S_M133L_S95T.xtc
│   └── ...
├── analysis/                   # RMSD, RMSF, clusters
│   ├── A189S_M133L_S95T/
│   │   ├── rmsd.xvg
│   │   ├── rmsf.xvg
│   │   ├── clusters.pdb
│   │   └── representative_conformers.pdb  # 10 conformers
│   └── ...
└── ensembles/                  # Final ensembles for docking
    ├── A189S_M133L_S95T_ensemble.pdb  # Multi-model PDB
    └── ...
```

### Success Criteria (Phase 1)
- ✅ RMSD converges < 3 Å (stable trajectory)
- ✅ RMSF hotspots correlate with known flexible regions (loops)
- ✅ Secondary structure preserved (core β-sheets stable)
- ✅ 10 representative conformers per rescue

---

## Phase 2: Drug Library Preparation (Week 1)

### Objective
Compile and prepare known p53-targeting drugs for docking

### Drug Selection Criteria
1. **Clinical relevance:** FDA-approved or in clinical trials
2. **Mechanism:** Stabilize p53 mutants or enhance WT function
3. **Structural availability:** 3D structures available (PDB/ChEMBL/PubChem)

### Known p53 Drugs (Top Candidates)

#### 1. **COTI-2** (Cotinine derivative)
- **Mechanism:** Stabilizes mutant p53, restores WT function
- **Status:** Phase 1 clinical trials (gynecological cancers)
- **Target mutations:** Y220C, R175H, R273H
- **PubChem CID:** 71464531

#### 2. **PHI-KAN-083** (small molecule)
- **Mechanism:** Binds Y220C pocket, stabilizes core domain
- **Status:** Preclinical, high specificity for Y220C
- **Target mutation:** Y220C
- **Published structure:** Available in literature

#### 3. **APR-246 (Eprenetapopt)**
- **Mechanism:** Covalent modifier, restores DNA binding
- **Status:** Phase 3 clinical trials (AML, MDS)
- **Target mutations:** Multiple hotspots
- **PubChem CID:** 11388326

#### 4. **Nutlin-3a** (MDM2 inhibitor)
- **Mechanism:** Blocks p53-MDM2 interaction, stabilizes WT p53
- **Status:** Preclinical (many derivatives in trials)
- **Note:** Works on WT p53, but may synergize with rescues
- **PubChem CID:** 11433190

#### 5. **PRIMA-1 / PRIMA-1MET**
- **Mechanism:** Restores mutant p53 to WT conformation
- **Status:** Preclinical (PRIMA-1MET = APR-246)
- **PubChem CID:** 5719

#### 6. **CP-31398**
- **Mechanism:** Stabilizes p53 core domain
- **Status:** Preclinical
- **PubChem CID:** 2842

### Drug Preparation Workflow

#### 2.1 Structure Retrieval
```python
# Download 3D structures from PubChem
from pubchempy import get_compounds

drugs = {
    'COTI-2': 71464531,
    'APR-246': 11388326,
    'Nutlin-3a': 11433190,
    'PRIMA-1': 5719,
    'CP-31398': 2842,
}

for name, cid in drugs.items():
    compound = get_compounds(cid, 'cid')
    # Download SDF format
```

#### 2.2 3D Optimization
```bash
# Use Open Babel for 3D structure generation
obabel -ipubchem COTI-2.sdf -o pdb -O COTI-2.pdb --gen3d

# Energy minimize with MMFF94 force field
obminimize -ff MMFF94 -n 500 COTI-2.pdb > COTI-2_minimized.pdb
```

#### 2.3 Protonation State
```bash
# Add hydrogens at pH 7.4 (physiological)
obabel COTI-2_minimized.pdb -O COTI-2_pH7.4.pdb -p 7.4
```

#### 2.4 Format Conversion
```bash
# Convert to PDBQT format for AutoDock Vina
prepare_ligand4.py -l COTI-2_pH7.4.pdb -o COTI-2.pdbqt
```

### File Structure
```
Data/raw/drugs/
├── structures/
│   ├── COTI-2.pdb
│   ├── APR-246.pdb
│   ├── Nutlin-3a.pdb
│   ├── PRIMA-1.pdb
│   ├── CP-31398.pdb
│   └── PHI-KAN-083.pdb
├── prepared/                   # For docking (PDBQT format)
│   ├── COTI-2.pdbqt
│   └── ...
└── metadata.json              # Drug properties, PubChem IDs, mechanisms
```

### Success Criteria (Phase 2)
- ✅ All 6 drugs downloaded and prepared
- ✅ 3D structures optimized (MMFF94)
- ✅ Protonation states correct (pH 7.4)
- ✅ PDBQT format ready for docking

---

## Phase 3: Ensemble Docking (Week 2)

### Objective
Dock all drugs to all rescue ensemble conformers to identify binding modes and affinities

### Docking Strategy

**Tool:** AutoDock Vina (fast, accurate, well-validated)
- Scoring function: empirical + knowledge-based
- Flexible ligand, semi-flexible receptor (side chains)
- Exhaustiveness: 32 (high accuracy)

**Ensemble approach:**
- Dock each drug to 10 conformers per rescue
- Total: 10 rescues × 6 drugs × 10 conformers = **600 docking runs**
- Time: ~5 min per run = ~50 hours (parallelizable)

### Binding Site Definition

**Three key sites for p53:**

#### Site 1: DNA Binding Interface
- **Purpose:** Assess if drug disrupts DNA contact (bad) or enhances it (good)
- **Center:** Residues R248, R273, R280 (DNA-contacting loop)
- **Box size:** 20 × 20 × 20 Å

#### Site 2: Y220C Pocket (if applicable)
- **Purpose:** Known druggable pocket for Y220C rescues
- **Center:** Residue 220
- **Box size:** 15 × 15 × 15 Å

#### Site 3: Hydrophobic Core
- **Purpose:** Stabilization via core binding
- **Center:** Residues M133, L145, V157
- **Box size:** 20 × 20 × 20 Å

### Docking Workflow

#### 3.1 Receptor Preparation
```bash
# For each conformer:
# - Remove water molecules
# - Add hydrogens (reduce)
# - Convert to PDBQT (AutoDockTools)

prepare_receptor4.py -r conformer_1.pdb -o conformer_1.pdbqt
```

#### 3.2 Grid Box Definition
```python
# Define binding site coordinates from structure
config = """
receptor = conformer_1.pdbqt
ligand = COTI-2.pdbqt

center_x = 30.5
center_y = 25.0
center_z = 40.2

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = 32
num_modes = 10
"""
```

#### 3.3 Docking Execution
```bash
# Run AutoDock Vina
vina --config config.txt --out output.pdbqt --log log.txt

# Parse binding affinity from log
# Best binding affinity: -8.5 kcal/mol (mode 1)
```

#### 3.4 Consensus Scoring
```python
# For each rescue × drug pair:
# - Dock to 10 conformers
# - Extract best binding affinity per conformer
# - Calculate consensus: median ΔG_binding across ensemble
# - Standard deviation: measure of binding mode consistency

consensus_affinity = median([ΔG_conf1, ΔG_conf2, ..., ΔG_conf10])
consistency = stdev([ΔG_conf1, ΔG_conf2, ..., ΔG_conf10])
```

### File Structure
```
Data/processed/docking/
├── results/
│   ├── A189S_M133L_S95T/
│   │   ├── COTI-2/
│   │   │   ├── conformer_1_docked.pdbqt
│   │   │   ├── conformer_1_log.txt
│   │   │   ├── conformer_2_docked.pdbqt
│   │   │   └── ...
│   │   ├── APR-246/
│   │   └── ...
│   └── ...
├── consensus/                  # Consensus scores per rescue-drug pair
│   ├── A189S_M133L_S95T.csv   # All drugs for this rescue
│   └── ...
└── summary.csv                 # All rescue-drug pairs ranked
```

### Success Criteria (Phase 3)
- ✅ All 600 docking runs complete
- ✅ Consensus affinities calculated for all pairs
- ✅ Binding modes visualized (top 3 per rescue)
- ✅ Identification of "superstar" rescue-drug pairs (ΔG < -9 kcal/mol)

---

## Phase 4: Synergy Scoring (Week 2-3)

### Objective
Identify rescue-drug pairs with synergistic effects (rescue + drug > rescue alone + drug alone)

### Synergy Metrics

#### 4.1 Functional Enhancement Score
```python
# Does the drug enhance rescue function beyond additive effects?

FES = (ΔΔG_rescue + ΔΔG_drug + ΔG_binding) / ΔΔG_rescue

# Where:
# - ΔΔG_rescue: Functional score from Initiative 1 (stability + DNA + interface)
# - ΔΔG_drug: Drug binding affinity (docking)
# - ΔG_binding: Interaction energy (synergy term)

# Interpretation:
# FES > 1.2: Synergistic (rescue + drug better than sum)
# FES 0.8-1.2: Additive (no synergy)
# FES < 0.8: Antagonistic (drug interferes with rescue)
```

#### 4.2 DNA Binding Preservation
```python
# Does the drug bind near DNA interface?
# If yes, does it enhance or disrupt DNA binding?

DNA_safe = distance(drug_pose, DNA_interface) > 5.0  # Å

# If drug is near DNA interface:
#   - Calculate ΔΔG_DNA_with_drug using EvoEF2 ComputeBinding
#   - Compare to ΔΔG_DNA_without_drug
#   - Enhancement: ΔΔG improves (more negative)
#   - Disruption: ΔΔG worsens (less negative)
```

#### 4.3 Stability Synergy
```python
# Does the drug further stabilize the rescued protein?

# Run EvoEF2 ComputeStability on:
# 1. Rescue mutant alone
# 2. Rescue mutant + docked drug (as ligand)

ΔΔG_stability_synergy = ΔG_rescue_with_drug - ΔG_rescue_alone

# Synergy if ΔΔG_stability_synergy < -2 kcal/mol
```

### Composite Synergy Score
```python
synergy_score = (
    0.40 × normalize(FES) +
    0.30 × normalize(ΔG_binding) +
    0.20 × normalize(ΔΔG_stability_synergy) +
    0.10 × (1.0 if DNA_safe else 0.0)
)

# Range: 0.0 to 1.0 (higher = better synergy)
```

### Ranking & Prioritization

**Categories:**
- **Superstar (synergy_score > 0.8):** Strong synergy, top priority for validation
- **Promising (synergy_score 0.6-0.8):** Moderate synergy, worth testing
- **Additive (synergy_score 0.4-0.6):** No strong synergy, but both work
- **Poor (synergy_score < 0.4):** Weak binding or antagonistic

### File Structure
```
Data/processed/synergy/
├── scores/
│   ├── A189S_M133L_S95T_synergy.csv   # All drugs for this rescue
│   └── ...
├── ranked_pairs.csv                    # All rescue-drug pairs ranked by synergy
└── superstar_pairs.csv                 # Top candidates (synergy > 0.8)
```

### Success Criteria (Phase 4)
- ✅ Synergy scores calculated for all 60 rescue-drug pairs (10 rescues × 6 drugs)
- ✅ Identification of ≥5 "superstar" pairs
- ✅ Validation that top pairs are DNA-safe
- ✅ Mechanistic insights: why do certain pairs synergize?

---

## Phase 5: Validation & Reporting (Week 3)

### Objective
Generate comprehensive report with actionable rescue-drug combination recommendations

### Analyses

#### 5.1 Cross-Target Comparison
- Which drugs work best for R175H vs R248Q vs R273H vs Y220C?
- Are there universal rescue-drug pairs (work across all mutations)?
- Target-specific vs general strategies

#### 5.2 Mechanism Analysis
- Where do drugs bind? (DNA interface, core, Y220C pocket)
- How do they stabilize? (H-bonds, hydrophobic contacts)
- Structural basis for synergy (visualizations)

#### 5.3 Clinical Prioritization
```python
# Rank by clinical feasibility
clinical_score = (
    0.30 × synergy_score +
    0.25 × (1.0 if drug_in_trials else 0.5) +
    0.20 × functional_score +
    0.15 × (1.0 - drug_toxicity_score) +
    0.10 × synthesis_feasibility
)
```

#### 5.4 Experimental Validation Plan
**In vitro assays (immediate):**
1. Recombinant protein stability (DSF, CD spectroscopy)
2. DNA binding affinity (EMSA, SPR)
3. Drug binding (ITC, SPR)

**Cell-based assays (6-12 months):**
1. Transient transfection (H1299 cells, p53-null)
2. Transcriptional activity (p21, Bax luciferase reporters)
3. Apoptosis assays (Annexin V, caspase-3)
4. Drug synergy validation (Bliss independence, Loewe additivity)

### Report Contents

**Executive Summary:**
- Top 5 rescue-drug pairs (ranked by synergy)
- Key findings (universal enhancers, target-specific hits)
- Clinical recommendations

**Methods:**
- MD simulation protocol
- Docking parameters
- Synergy scoring formula

**Results:**
- Per-target performance
- Binding mode analysis
- Synergy mechanisms

**Discussion:**
- Comparison to literature (existing p53 drugs)
- Novelty of rescue-drug combination approach
- Limitations and future work

**Supplementary:**
- Full ranking of all 60 pairs
- MD trajectory metrics (RMSD, RMSF)
- Docking poses (top 3 per pair)

### File Structure
```
reports/
├── initiative2_results.md          # Main report
├── figures/
│   ├── synergy_heatmap.png        # Rescue × Drug synergy matrix
│   ├── top_binding_modes.png      # Structural visualizations
│   ├── md_convergence.png         # RMSD/RMSF plots
│   └── clinical_prioritization.png
└── supplementary/
    ├── all_pairs_ranked.csv
    ├── md_metrics_summary.csv
    └── docking_details.csv
```

### Success Criteria (Phase 5)
- ✅ Comprehensive report complete
- ✅ Top 5 rescue-drug pairs identified
- ✅ Mechanistic insights documented
- ✅ Experimental validation plan ready

---

## Implementation Timeline

### Week 1: Infrastructure Setup
- **Days 1-2:** GROMACS installation, force field setup, test MD on 1 rescue
- **Days 3-4:** Drug structure retrieval and preparation (all 6 drugs)
- **Day 5:** MD pipeline automation (scripts for all 10 rescues)

### Week 2: Production Simulations & Docking
- **Days 1-3:** Run all 10 MD simulations (parallel on GPU)
- **Days 4-5:** Trajectory analysis, extract representative conformers

### Week 2-3: Docking & Synergy
- **Days 6-8:** Run all 600 docking calculations (parallel)
- **Days 9-10:** Consensus scoring, identify top pairs

### Week 3: Synergy Scoring & Reporting
- **Days 11-12:** Calculate synergy scores for all pairs
- **Days 13-15:** Generate report, figures, prioritization

**Total:** ~3 weeks (can compress to 2 weeks with more parallelization)

---

## Resource Requirements

### Computational
**GPU:**
- Option 1: Local NVIDIA RTX 3080/3090 (~$1000, reusable)
- Option 2: Cloud GPU (AWS p3.2xlarge: $3.06/hour × 100 hours = $306)
- **Recommendation:** Cloud GPU for flexibility

**Storage:**
- MD trajectories: ~50 GB (10 rescues × 5 GB each)
- Docking results: ~5 GB
- Total: ~60 GB (easily fits on external drive)

### Software (All Free/Open-Source)
- GROMACS 2023+ (GPU-accelerated MD)
- AutoDock Vina (molecular docking)
- Open Babel (drug preparation)
- BioPython, MDAnalysis (analysis)
- PyMOL (visualization)

### External Data
- PubChem (drug structures) - free API
- PDB structures (already have 1TSR, 3KMD)
- Literature (binding site coordinates for Y220C drugs)

---

## Risk Mitigation

### Risk 1: MD Simulations Don't Converge
**Mitigation:**
- Extend simulation time to 100 ns if needed
- Check force field parameters
- Validate on known p53 WT structure first

### Risk 2: Poor Docking Scores (All Weak Binders)
**Mitigation:**
- Try multiple binding sites (DNA interface, core, Y220C pocket)
- Use flexible docking (allow side chain movement)
- Validate with known p53-drug complexes (if available)

### Risk 3: No Synergy Detected
**Mitigation:**
- Synergy may be subtle - look for additive effects
- Even additive combinations are clinically valuable
- Focus on DNA-safe drugs (avoid disruption)

### Risk 4: Computational Cost Higher Than Expected
**Mitigation:**
- Start with top 5 rescues (instead of 10)
- Reduce MD time to 25 ns (still useful)
- Use GPU efficiently (batch jobs)

---

## Success Criteria (Overall Initiative 2)

### Quantitative Goals
- ✅ 10 rescue ensembles generated (50 ns MD each)
- ✅ 600 docking runs completed (10 rescues × 6 drugs × 10 conformers)
- ✅ ≥5 "superstar" rescue-drug pairs identified (synergy > 0.8)
- ✅ 100% DNA-safe top pairs (no disruption of DNA binding)

### Qualitative Goals
- ✅ Mechanistic insights into rescue-drug synergy
- ✅ Clinical prioritization of top combinations
- ✅ Experimental validation plan ready
- ✅ Upgrade to 95%+ scientific rigor

### Publication Impact
**Updated narrative:**
> "We have identified rescue mutations that restore p53 function across stability, DNA binding, and tetramerization (Initiative 1), and predicted synergistic rescue-drug combinations using dynamic ensemble modeling and molecular docking (Initiative 2). Top combinations show strong binding affinity, preserve DNA interactions, and are ready for experimental validation."

**Suitable for:**
- *Nature Communications* (computational + experimental validation)
- *Cell Chemical Biology* (drug discovery focus)
- *Journal of Medicinal Chemistry* (drug-protein interactions)

---

## Next Steps After Initiative 2

### Option 1: Experimental Validation (Partner with Lab)
- Collaborate with structural biology lab
- Test top 3 rescue-drug pairs
- Validate synergy predictions

### Option 2: Expand Drug Library
- Screen FDA-approved drug library (>1000 compounds)
- Use virtual screening to find novel p53 binders
- Repurposing existing drugs

### Option 3: Patient-Specific Predictions
- Analyze TCGA tumor mutation data
- Match mutations to best rescue-drug pairs
- Precision medicine approach

---

## Status: Ready to Begin

**Current:** Roadmap complete, ready for implementation
**Next action:** Set up GROMACS environment + run pilot MD on top rescue
**Decision point:** Confirm computational resources (local GPU vs cloud)

**Let's elevate p53 rescue design to 95%+ rigor! 🚀**
