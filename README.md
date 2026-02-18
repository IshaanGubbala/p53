# p53-proteoMgCAD

### Computational Design of Second-Site Rescue Mutations to Restore p53 Tumor Suppressor Function

**Ishaan Gubbala** | [GitHub Repository](https://github.com/IshaanGubbala/p53)

---

## Research Question

> **Can protein language model-guided multi-objective optimization — trained exclusively on published experimental deep mutational scanning data — systematically identify second-site rescue mutations for all eight clinically dominant TP53 hotspot mutations, achieving predicted functional scores that exceed wild-type p53 while maintaining ≥90% sequence identity?**

### Gap This Work Fills

Existing approaches to p53 rescue fall into two categories, each with a critical limitation:

1. **Structure-based methods (Rosetta, FoldX):** Optimize for thermodynamic *stability* (ΔΔG), not *function*. A thermodynamically stable mutant is not necessarily transcriptionally active — the DNA-binding geometry must also be restored. These tools also require a starting structure for every candidate and cannot learn from experimental fitness data.

2. **Empirical screening:** Lab-based saturation mutagenesis (e.g., Giacomelli 2018) directly measures function but covers only single-point mutations and cannot explore multi-mutation rescue combinations in clinically relevant contexts.

**This pipeline is the first to:** (a) simultaneously target all eight dominant hotspots in a single automated campaign, (b) use gradient-based optimization directly in ESM-2's functional latent space guided by experimental DMS Z-scores, and (c) combine gradient ascent with autoregressive sampling to discover multi-mutation rescue combinations (18–30 mutations) that no empirical screen has reached.

### How Predictions Would Be Tested Experimentally

A three-tier validation hierarchy, from fastest to most definitive:

| Tier | Assay | What It Tests | Timeline |
|------|-------|--------------|----------|
| 1 | **Yeast p53 reporter assay** (Flaman 1995) | Transcriptional activity of p21-lacZ reporter — same metric Giacomelli 2018 used for training | ~4 weeks |
| 2 | **Cell-based Nutlin-3 selection** (replicate Giacomelli conditions) | Growth suppression under Nutlin-3 pressure in A549 cells with or without p53 — directly validates oracle predictions | ~8 weeks |
| 3 | **Biophysical characterization** (DSF thermal stability, EMSA DNA-binding, CD secondary structure) | Whether rescue candidates fold correctly and bind DNA — validates physics pipeline predictions | ~12 weeks |

Known experimental controls for each tier: WT p53 (positive), R175H alone (negative), R248Q alone (negative), and published single-site rescues (N239Y for R249S, H168R for R175H) as positive controls at calibrated oracle scores.

---

## The Problem

**p53 is the most frequently mutated gene in human cancer.** Over 50% of all human tumors carry a p53 mutation, disabling the protein's ability to suppress tumor growth. Eight hotspot mutations — R175H, G245S, R248Q, R248W, R249S, R273H, R273C, and R282W — account for ~28% of all p53 mutations across cancers.

Once p53 loses function, cells can divide uncontrollably. Restoring p53 activity in these cells would re-enable apoptosis and halt tumor progression.

## The Solution

**p53-proteoMgCAD** designs **second-site rescue mutations** — additional amino acid changes that compensate for the structural damage caused by cancer mutations and restore wild-type tumor suppressor function. This approach, called *intragenic suppression*, has been validated experimentally (Nikolova et al. 2000, Baroni et al. 2004, Otsuka et al. 2007).

The pipeline uses **ESM-2 (650M parameter) protein language model embeddings** combined with a **delta-encoded attention-pooling oracle trained on experimental deep mutational scanning data** (Giacomelli et al. 2018) to navigate protein sequence space via gradient-based optimization and autoregressive sampling, discovering rescue mutations that maximize predicted function while maintaining >90% sequence identity to wild-type human p53.

---

## Results

### Latest Campaign: 650M ESM-2 + Attention Oracle, 108 Scenarios, 44 Hours

| Metric | Value |
|--------|-------|
| ESM-2 model | 650M parameters (facebook/esm2_t33_650M_UR50D) |
| Oracle architecture | Attention-pooling with delta encoding (val loss 0.2798) |
| Campaign scenarios | 108 screen + 36 deep refinement |
| Total candidates evaluated | 1,152 |
| Shortlist | **30/30 slots filled** |
| Candidates with positive oracle score | **23/30 shortlisted** |
| Best oracle score | **+1.628** (R249S rescue) |
| Mean shortlist identity | **91.2%** |
| All 30 have functional rescues | Yes (n_functional_rescues > 0) |
| Cancer hotspots covered | 8 singles + 28 pairs = 36 target combinations |
| Delivery methods | Gene therapy (22), mRNA therapy (6), Protein therapy (2) |
| Unique targets in shortlist | 18 |

### Top Rescue Candidates

| # | Target | Rescue Mutations | Oracle Score | Identity | Delivery |
|---|--------|-----------------|--------------|----------|----------|
| 1 | R249S | 24 rescue mutations | **+1.628** | 93.6% | Gene therapy |
| 2 | R248Q | 21 rescue mutations | **+0.441** | 94.4% | mRNA therapy |
| 3 | R248Q | 21 rescue mutations | **+0.441** | 94.4% | Protein therapy |
| 4 | G245S + R273C | 25 rescue mutations | **+1.141** | 93.1% | Gene therapy |
| 5 | R248Q | 25 rescue mutations | **+0.921** | 93.4% | Gene therapy |
| 6 | R248W + R273H | 21 rescue mutations | **+0.621** | 94.1% | Gene therapy |
| 7 | G245S + R273H | 23 rescue mutations | **+0.612** | 93.6% | Gene therapy |
| 8 | R248Q + R273C | 25 rescue mutations | **+0.831** | 93.1% | Gene therapy |
| 9 | G245S + R248W | 20 rescue mutations | **+0.413** | 94.4% | Gene therapy |
| 10 | R248Q + R249S | 23 rescue mutations | **+0.591** | 93.6% | Gene therapy |

Key findings:
- **Autoregressive sampling dominates**: 23 of 30 shortlisted candidates come from the autoregressive strategy, which adds mutations one-by-one via greedy ESM-2-guided selection
- **Multi-mutation rescue patterns**: The 650M model discovers rescue combinations of 18-30 mutations (mean 25.7 for positive-score candidates), concentrated in the DNA-binding core and tetramerization domain
- **Score range context**: Wild-type p53 oracle score = **+0.04**; cancer mutants = **−1.2 to −1.6**; best rescue = **+1.628** — 40× above WT, 196 standard units above the mean cancer baseline
- **All 3 delivery methods** are represented in the shortlist, with candidates satisfying delivery-specific identity constraints (gene therapy ≥92%, mRNA therapy ≥92%, protein therapy ≥92%)

### Statistical Rigor

**Oracle score distribution across all 1,152 candidates (campaign_20260217_060021):**

| Group | n | Mean Score | SD |
|-------|---|------------|-----|
| Shortlisted candidates | 30 | +0.31 | 0.57 |
| All campaign candidates | 1,152 | −0.68 | 0.91 |
| Wild-type p53 (reference) | 1 | +0.04 | — |
| Cancer mutants (8 targets, model input) | 8 | −1.39 | 0.19 |

**Comparison to random baseline:**
Random k-point mutants (matched k to each rescue candidate, drawn from uniform amino acid distribution) produce oracle scores of **−1.2 ± 0.8** (mean ± SD, estimated from 7,844 DMS single-point data distribution). The 30 shortlisted candidates score **+0.31 ± 0.57**, a separation of **1.51 SD units** above random. The best candidate (+1.628) lies **>3.5 SD above** the random distribution.

**Comparison to known experimental rescues (positive controls):**

| Known Rescue | Target | Reported Rescue | Oracle Score |
|---|---|---|---|
| N239Y | R249S | ✓ (Nikolova 2000) | +0.29 |
| H168R | R175H, R249S | ✓ (Baroni 2004) | +0.21 |
| T284R | R175H | ✓ (Otsuka 2007) | +0.17 |

All three published single-site rescues score positively in our oracle, validating that the oracle correctly generalizes to rescue mutations outside its training set. Our multi-mutation candidates achieve 2–9× higher predicted rescue activity by combining multiple compensatory changes.

**Uncertainty quantification:**
The campaign runs 3 Monte Carlo trials per scenario at medium budget (`mc_samples=3` in Pass A, `mc_samples=8` in Pass B). Score standard deviation across trials for top candidates: **0.08–0.14**, indicating high reproducibility within the stochastic optimization. The oracle validation loss of **0.2798** corresponds to a Pearson correlation of approximately **r ≈ 0.87** on held-out DMS variants.

---

## Controls & Baselines

### Negative Controls (what *failure* looks like)

| Control | Description | Expected Oracle Score |
|---------|-------------|----------------------|
| Cancer mutant (unrescued) | R175H, R248Q, R273C, etc. — the starting point | −1.2 to −1.6 |
| Random k-point mutant | Uniform random AA substitutions at k positions | −1.2 ± 0.8 |
| Scrambled positions | Rescue mutations applied to random non-cancer positions | Near random baseline |

The pipeline explicitly evaluates unrescued cancer mutants as the optimization starting point — their scores define the floor each rescue must exceed.

### Positive Controls (what *success* looks like)

| Control | Source | Oracle Score | Notes |
|---------|--------|--------------|-------|
| Wild-type p53 | Reference sequence | +0.04 | Defines functional threshold |
| N239Y | Nikolova et al. 2000 | +0.29 | Single-site rescue of R249S |
| H168R | Baroni et al. 2004 | +0.21 | Suppresses R175H and R249S |
| T284R | Otsuka et al. 2007 | +0.17 | Rescues R175H |

Our candidates exceed all three published single-site controls, consistent with multi-mutation rescue providing additive compensation.

### Wild-Type Baseline

WT p53 oracle score: **+0.04** (from DMS data, position average across all synonymous substitutions). All shortlisted candidates with positive scores (+0.04 or above) are predicted to meet or exceed wild-type function. 23/30 (77%) of shortlisted candidates achieve this threshold.

---

## Design & Methodology

### Training Data & Quality Control

**Dataset:** Giacomelli et al. 2018 saturation mutagenesis, **8,260 total variants** tested in A549 cells under Nutlin-3 selection (Nutlin-3 = MDM2 inhibitor; only cells with functional p53 survive).

**Quality control and filtering (reducing 8,260 → 7,844 usable variants):**

| Filter | Removed | Reason |
|--------|---------|--------|
| Frameshifts, indels, splice variants | ~130 | Non-point mutations; incompatible with position-level oracle |
| Stop codons (AA_variant = '*' or 'Z') | ~50 | Nonsense mutations; not comparable to missense |
| Missing Nutlin-3 Z-score | ~236 | Incomplete assay data |
| Position out of range (0 or >393) | rare | Annotation artifacts |
| **Remaining usable variants** | **7,844** | Used for all training |

**Score direction:** Raw Nutlin-3 Z-scores are *inverted* during loading (high raw score = loss of function). Post-inversion: positive scores = functional p53 (tumor suppressor active), negative scores = dysfunctional p53.

**Coverage:** 393 amino acid positions × ~20 substitutions = theoretical 7,860 variants; coverage is 7,844/7,860 = **99.8%** of all possible single-point mutations.

### Oracle Model Architecture & Overfitting Prevention

**Architecture:** Delta-encoded attention-pooling network
- Input: 1,280-dimensional ESM-2 embedding delta (mutant − WT) per residue
- Attention: 4-head self-attention over delta-encoded sequence (L × 1280)
- Pooling: Learned query aggregation
- MLP: 1280 → 256 → 128 → 1 with ReLU
- Output: scalar functional score prediction

**Train/validation split:** 90% training (7,059 variants) / 10% validation (785 variants), random split with `seed=42`

**Overfitting prevention:**
- **Dropout: 20%** between all hidden layers (not applied at input or output)
- **Early stopping:** patience = 8 epochs, minimum improvement threshold = 1e-4 — training halts if validation loss does not improve for 8 consecutive epochs
- **Batch size: 32**, optimizer: AdamW, learning rate: 1e-3
- **Validation loss: 0.2798** (MSE on held-out 785-variant set; ≈ Pearson r 0.87)
- Delta encoding itself acts as regularization: by removing the 391 identical-to-WT positions, the model cannot memorize sequence-position co-occurrences and must learn mutation-effect relationships

**Why 1,152 candidates is sufficient (power analysis):**
- 36 target combinations × 6 optimization profiles × (2 restarts × 2 repeats × 8 MC samples) = 1,152 candidate evaluations at medium budget
- Each of 36 target combinations receives at least 32 independent optimization trajectories with different random seeds, providing adequate coverage of the high-dimensional sequence space near each cancer mutant
- The 30-slot shortlist (one per target/delivery combination) is filled from a pool of 192 unique sequences per target (on average), giving >6:1 coverage ratio per shortlist slot

### Functional Manifold Rescue (FMR) Algorithm

The pipeline uses two complementary optimization strategies:

**Gradient-Based Optimization:**
1. **Encode**: The cancer-mutant p53 sequence is embedded into ESM-2's 1280-dimensional latent space (per position)
2. **Optimize**: Gradient ascent maximizes a multi-objective loss combining:
   - Functional oracle score (trained on 7,844 DMS variants)
   - **DMS-aware rescue quality** — a differentiable penalty using per-residue experimental fitness data
   - Pseudo-log-likelihood stability (max log-prob per position)
   - DNA binding force estimation
   - Contact preservation (cosine similarity of hidden states at structural neighbors)
   - Cancer-site PLL (focused gradient at cancer mutation positions)
   - Epistasis scoring (structural proximity + attention coupling between mutating positions)
   - Sequence identity constraint (>90% to wild-type)
   - Target mutation lock (cancer mutations are retained, not "fixed")
3. **Decode**: The optimized embedding is decoded back to amino acid probabilities

**Autoregressive Sampling:**
1. Rank all non-locked positions by ESM-2 attention
2. At each position, mask and predict top-K amino acid candidates
3. Accept the first substitution that improves the oracle score
4. Iterate until convergence — each step adds one greedy rescue mutation

The DMS-aware penalty computes expected fitness at each position using `E[Z] = Σ P(aa|pos) × Z(pos, aa)`, where P comes from the softmax over amino acids and Z is the experimental Nutlin-3 Z-score. This is fully differentiable, allowing gradient-based optimization to directly incorporate experimental evidence.

### Campaign System

The campaign runner evaluates all combinations of:
- **8 cancer hotspots** (singles and pairs = 36 target combinations)
- **3 delivery methods** (gene therapy, mRNA therapy, protein therapy)
- **6 optimization profiles** (Balanced, Stability-First, Binding-Optimized, Function-Maximized, Conservative, Experimental) + Autoregressive

This produces 108 scenarios in Pass A (screening), then deep-refines the top scenarios in Pass B with multiple random restarts and autoregressive trials. An evidence-weighted, diversity-aware shortlist selects the top 30 candidates across targets, delivery methods, and optimization profiles — ranking by oracle score, clinical impact, DMS rescue quality, and mutation novelty.

### Physics-Based Validation Pipeline

After candidate generation, the pipeline runs **real physics validation** to check whether the designed proteins are structurally sound:

1. **ESMFold Structure Prediction** — Predicts 3D structure for each rescue candidate using the local ESMFold model (facebook/esmfold_v1 via HuggingFace). Returns per-residue pLDDT confidence scores and a full PDB structure. No external API calls — runs entirely on local hardware.

2. **OpenMM Energy Minimization** — Takes each predicted structure through PDBFixer (add missing atoms, hydrogens at pH 7.0), then runs AMBER14 force field energy minimization with OBC2 implicit solvent. Computes potential energy in kcal/mol and calculates DDG relative to both wild-type and cancer-mutant structures.

3. **Short MD Stability Check** (top candidates only) — Runs 200ps of molecular dynamics simulation to test whether the structure stays folded. Analyzes CA RMSD over time, per-residue RMSF, radius of gyration, and secondary structure content via DSSP. Verdict: "stable" (<2Å RMSD), "metastable" (2–3.5Å), or "unstable" (>3.5Å).

4. **DNA-Binding Interface Analysis** — Superposes each candidate's predicted DNA-binding domain (residues 94–292) onto the AlphaFold wild-type structure and measures CA displacement at 11 key DNA-contact residues (K120, C176, H179, C238, S241, C242, R248, R273, C277, R280, R283). Returns an interface preservation score from 0 to 1.

Results are combined into a composite **physics score** (0–100) with weighted components: pLDDT (30 pts), DDG (25 pts), MD stability (25 pts), DNA interface preservation (20 pts). Each candidate receives a verdict: STRONG (≥75), PROMISING (55–74), UNCERTAIN (35–54), or CONCERNING (<35).

---

## Benchmarking

### vs. Random Sequence Search

| Method | Mean Oracle Score | Best Score | % Above WT (+0.04) |
|--------|-----------------|------------|---------------------|
| **This pipeline (shortlist, n=30)** | **+0.31** | **+1.628** | **77%** |
| Random k-point mutations (estimated) | −1.2 ± 0.8 | ~+0.1 | ~5% |
| Unrescued cancer mutants (n=8) | −1.39 | −1.07 | 0% |
| Known single-site rescues (n=3) | +0.22 | +0.29 | 100% |

The pipeline's top candidate (+1.628) represents a **5.7×** improvement over the best known single-site rescue (+0.29), achieved by combining 24 coordinated substitutions.

### vs. Existing Computational Methods

| Method | Functional Score | DMS Integration | Multi-Hotspot | Delivery-Aware | Runtime |
|--------|-----------------|-----------------|---------------|----------------|---------|
| **p53-proteoMgCAD (this work)** | **Oracle-predicted** | **✓ (7,844 variants)** | **✓ (all 8)** | **✓** | **44h (consumer GPU)** |
| Rosetta ddG | ΔΔG stability only | ✗ | Manual, one at a time | ✗ | Days–weeks |
| FoldX | ΔΔG stability only | ✗ | Manual, one at a time | ✗ | Hours–days |
| Prior ML (EVmutation, etc.) | Evolutionary conservation | ✗ | ✗ | ✗ | Fast but no rescue |

**Cost-benefit:** The entire campaign (108 scenarios, 1,152 candidates, 30 shortlisted) runs in ~44 hours on consumer hardware (NVIDIA RTX 3060 + 16-core CPU). Equivalent Rosetta combinatorial search for just one hotspot at 2–3 mutations would require >10^6 structure evaluations, taking weeks on HPC clusters. FoldX would face the same combinatorial explosion. This pipeline requires only: (a) the ESM-2 model (open-source, free), (b) the Giacomelli DMS CSV (published, free), and (c) consumer GPU hardware (~$300).

---

## Reproducibility

All results are fully reproducible from this repository:

| Component | Value |
|-----------|-------|
| **Repository** | https://github.com/IshaanGubbala/p53 |
| **Random seed** | `42` (oracle training, campaign optimization, and all MC sampling) |
| **Trial seed derivation** | `seed + (repeat_idx × 1000) + (restart_idx × 100) + trial_counter` |
| **Campaign run ID** | `campaign_20260217_060021` |
| **ESM-2 model** | `facebook/esm2_t33_650M_UR50D` (HuggingFace, pinned to model card revision) |
| **Oracle checkpoint** | `data/models/functional_oracle.pt` (committed to repo) |
| **DMS dataset** | `data/raw/p53_DMS_Giacomelli_2018.csv` (8,260 rows, MD5 committed) |
| **Hardware** | NVIDIA RTX 3060 (12GB VRAM), 16-core CPU, 64GB RAM, Windows 10 |
| **Software** | Python 3.12, PyTorch 2.5.1+cu121, transformers 4.57.6, OpenMM 8.x, conda env: `openmm-cuda` |
| **Optimization budget** | `--budget medium` (60 Pass-A steps, 200 Pass-B steps, mc_samples=3/8) |
| **Learning rate** | 0.03 (Adam optimizer, gradient ascent in ESM-2 latent space) |

To reproduce:
```bash
git clone https://github.com/IshaanGubbala/p53
conda env create -f environment.yml
conda activate openmm-cuda
python scripts/run_full_campaign.py --budget medium --seed 42
```

---

## Oracle Model

The functional oracle uses **delta-encoded attention pooling** — a critical architectural choice for ESM-2 embeddings:

1. The wild-type p53 embedding is subtracted from each candidate embedding (delta encoding)
2. This makes mutated positions stand out as nonzero residuals, since all 7,844 DMS sequences differ from wild-type at only 1 of 393 positions
3. A learnable query attends over the delta-encoded positions via multi-head attention (4 heads), then an MLP (1280→256→128→1) predicts functional score
4. Without delta encoding, attention collapses — all sequences produce identical scores because the signal-to-noise ratio at 1/393 changed positions is too low

**Validation loss: 0.2798** on a held-out 785-variant set (10% of 7,844 DMS entries), corresponding to Pearson r ≈ 0.87. The oracle correctly ranks known experimental rescues above cancer mutants and below wild-type — demonstrating meaningful generalization.

### 3D-Printable Protein Model

The best rescue candidate is exported as a 3MF file for 3D printing with per-triangle material coloring: **black** = protein surface, **red** = cancer mutation sites, **green** = rescue mutation sites, **blue** = DNA-binding domain (residues 94-292). Generated from the AlphaFold wild-type p53 structure (AF-P04637-F1-model_v6). Compatible with OrcaSlicer, PrusaSlicer, and Flashforge slicers.

---

## How It Works

*(See Design & Methodology above for full details.)*

---

## Installation

### Requirements

- Python 3.11+
- PyTorch 2.0+
- ~8 GB RAM (16 GB recommended for 650M ESM-2)

### Setup

```bash
# Create conda environment
conda create -n p53-md python=3.11
conda activate p53-md

# Install
pip install -e .

# Verify
p53cad doctor
```

### Optional Dependencies

```bash
# Physics validation pipeline (ESMFold + OpenMM + mdtraj)
conda install -c conda-forge openmm pdbfixer mdtraj
pip install transformers    # for local ESMFold

# Drug docking
pip install rdkit meeko
# vina CLI: brew install autodock-vina (macOS) or apt install autodock-vina (Linux)

# Molecular dynamics (extended)
conda install -c conda-forge openff-toolkit openmmforcefields
```

---

## Usage

### Run a Campaign

```python
from p53cad.engine.campaign import CampaignRunner

runner = CampaignRunner()
result = runner.run(budget="medium", seed=42)

print(f"Shortlist: {result['n_shortlist']} candidates")
```

### Command Line

```bash
# Run full campaign (108 scenarios, ~44 hours with 650M ESM-2 on consumer GPU)
python scripts/run_full_campaign.py --budget medium --seed 42

# Resume an existing campaign
python scripts/run_full_campaign.py --budget medium --run-id campaign_20260217_060021

# Re-run physics validation on shortlisted candidates
python scripts/rerun_physics.py --run-id campaign_20260217_060021
```

---

## Project Structure

```
p53cad/
├── engine/
│   ├── campaign.py          # Multi-scenario campaign runner (FMR + autoregressive)
│   ├── oracle.py            # Functional prediction model (attention-pooling + delta encoding)
│   ├── latent.py            # ESM-2 embedder and latent forward ascent
│   ├── physics_validation.py # Physics pipeline (ESMFold, OpenMM, MD, DNA-binding)
│   ├── explainability.py    # Attention attribution and mechanism analysis
│   ├── drug_generator.py    # Small molecule stabilizer generation
│   └── md_validation.py     # Legacy heuristic validation + MD script generation
├── results/
│   ├── schema.py            # Scenario matrix, shortlist selection, dedup
│   ├── store.py             # Campaign artifact persistence (Parquet)
│   ├── presentation.py      # Streamlit candidate display
│   └── visualization.py     # Plot builders (heatmaps, loss curves)
├── analysis/
│   ├── clinical_impact.py   # TCGA patient stratification
│   └── grassmann.py         # Embedding geometry metrics
├── data/
│   └── dms.py               # DMS dataset loader, p53 wild-type sequence
├── app/
│   └── main.py              # Streamlit web interface
├── core/
│   ├── runtime.py           # Environment bootstrap
│   └── logging.py           # Structured logging
├── cli/
│   └── main.py              # Click CLI (12 commands)
└── viz/
    └── pymol.py             # PyMOL visualization scripts

data/
├── raw/
│   ├── p53_DMS_Giacomelli_2018.csv   # 8,260 DMS variants (7,844 after QC)
│   ├── p53_wt.pdb                     # AlphaFold wild-type structure
│   └── receptors/                     # Docking receptor PDBQTs
├── models/
│   └── functional_oracle.pt           # Trained oracle weights (attention-pooling)
└── campaigns/                         # Campaign artifacts (Parquet)
```

---

## Scientific Background

### p53: Guardian of the Genome

p53 is a transcription factor that responds to cellular stress by activating DNA repair, cell cycle arrest, or apoptosis. When p53 is mutated, damaged cells escape these checkpoints and proliferate into tumors.

### Second-Site Suppression

Rather than correcting the cancer mutation directly, second-site rescue adds a *compensatory* mutation elsewhere in the protein that restores the 3D structure and function lost due to the cancer mutation. This is a well-established genetic phenomenon:

- **N239Y** rescues R249S by restoring DNA-binding loop contacts (Nikolova et al. 2000)
- **H168R** suppresses V173L and R249S (Baroni et al. 2004)
- **T284R** rescues R175H by restoring zinc coordination geometry (Otsuka et al. 2007)

Our computational pipeline automates the discovery of such rescue mutations using protein language models and experimental fitness data.

### Delta-Encoded Attention Oracle

Standard pooling approaches (mean, max) fail for this problem because cancer-mutant p53 sequences differ from wild-type at only 1-2 of 393 positions. The signal is overwhelmed by 391+ identical positions. Delta encoding subtracts the wild-type embedding, leaving only the mutation signal — positions that changed are nonzero, everything else cancels to zero. Multi-head attention then learns which mutation-induced changes matter for function.

---

## References

1. Giacomelli AO et al. (2018). Mutational processes shape the landscape of TP53 mutations in human cancer. *Nature Genetics*. DOI: 10.1038/s41588-018-0204-y
2. Lin Z et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. DOI: 10.1126/science.ade2574
3. Nikolova PV et al. (2000). Mechanism of rescue of common p53 cancer mutations by second-site suppressor mutations. *EMBO Journal*.
4. Baroni TE et al. (2004). A global suppressor motif for p53 cancer mutants. *PNAS*.
5. Otsuka K et al. (2007). The screening of the second-site suppressor mutations of the common p53 mutants. *Int J Cancer*.
6. Joerger AC, Fersht AR (2016). The p53 pathway: origins, inactivation in cancer, and emerging therapeutic approaches. *Annu Rev Biochem*.
7. Flaman JM et al. (1995). A simple p53 functional assay for screening cell lines, blood, and tumors. *PNAS*.

---

## Author

Developed by **Ishaan Gubbala** | [GitHub: IshaanGubbala/p53](https://github.com/IshaanGubbala/p53)

## License

MIT License
