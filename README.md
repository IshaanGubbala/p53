# p53-proteoMgCAD

### Computational Design of Second-Site Rescue Mutations to Restore p53 Tumor Suppressor Function

**Ishaan Gubbala**

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
- **Score range**: Wild-type p53 scores +0.04; cancer mutants score -1.2 to -1.6; the best rescue reaches +1.628 — exceeding wild-type predicted function
- **All 3 delivery methods** are represented in the shortlist, with candidates satisfying delivery-specific identity constraints (gene therapy ≥92%, mRNA therapy ≥92%, protein therapy ≥92%)

### Oracle Model

The functional oracle uses **delta-encoded attention pooling** — a critical architectural choice for ESM-2 embeddings:

1. The wild-type p53 embedding is subtracted from each candidate embedding (delta encoding)
2. This makes mutated positions stand out as nonzero residuals, since all 7,844 DMS sequences differ from wild-type at only 1 of 393 positions
3. A learnable query attends over the delta-encoded positions via multi-head attention (4 heads), then an MLP (1280→256→128→1) predicts functional score
4. Without delta encoding, attention collapses — all sequences produce identical scores because the signal-to-noise ratio at 1/393 changed positions is too low

### 3D-Printable Protein Model

The best rescue candidate is exported as a 3MF file for 3D printing with per-triangle material coloring: **black** = protein surface, **red** = cancer mutation sites, **green** = rescue mutation sites, **blue** = DNA-binding domain (residues 94-292). Generated from the AlphaFold wild-type p53 structure (AF-P04637-F1-model_v6). Compatible with OrcaSlicer, PrusaSlicer, and Flashforge slicers.

---

## How It Works

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

### Experimental Validation

The functional oracle is trained on the **Giacomelli et al. 2018** saturation mutagenesis dataset — 7,844 p53 variants tested in a cell-based Nutlin-3 selection assay. This is real experimental data measuring whether each p53 variant retains tumor suppressor activity. The oracle achieves validation loss of 0.2798 and generalizes to multi-mutation combinations through delta encoding.

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

### Optional (for drug docking and MD simulations)

```bash
pip install rdkit meeko
conda install -c conda-forge openmm openff-toolkit
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
# Run full campaign (108 scenarios, ~44 hours with 650M ESM-2 on Apple MPS GPU)
p53cad campaign-run --budget medium --seed 42

# List past campaigns
p53cad campaign-list

# Generate report for a campaign
p53cad campaign-report --run-id <run_id>

# Launch web interface
streamlit run p53cad/app/main.py
```

### Generate 3D-Printable Model

```python
# After running a campaign, the best candidate can be exported as 3MF
# See data/campaigns/<run_id>/best_candidate_p53_rescue.3mf
```

---

## Project Structure

```
p53cad/
├── engine/
│   ├── campaign.py          # Multi-scenario campaign runner (FMR + autoregressive)
│   ├── oracle.py            # Functional prediction model (attention-pooling + delta encoding)
│   ├── latent.py            # ESM-2 embedder and latent forward ascent
│   ├── explainability.py    # Attention attribution and mechanism analysis
│   ├── drug_generator.py    # Small molecule stabilizer generation
│   └── md_validation.py     # Molecular dynamics simulation scripts
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
│   └── main.py              # Click CLI (10 commands)
└── viz/
    └── pymol.py             # PyMOL visualization scripts

data/
├── raw/
│   ├── p53_DMS_Giacomelli_2018.csv   # 8,258 DMS variants
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

---

## Author

Developed by **Ishaan Gubbala**

## License

MIT License
