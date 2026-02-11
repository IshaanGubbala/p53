# p53-proteoMgCAD

### Computational Design of Second-Site Rescue Mutations to Restore p53 Tumor Suppressor Function

**Ishaan Gubbala**

---

## The Problem

**p53 is the most frequently mutated gene in human cancer.** Over 50% of all human tumors carry a p53 mutation, disabling the protein's ability to suppress tumor growth. Eight hotspot mutations — R175H, G245S, R248Q, R248W, R249S, R273H, R273C, and R282W — account for ~28% of all p53 mutations across cancers.

Once p53 loses function, cells can divide uncontrollably. Restoring p53 activity in these cells would re-enable apoptosis and halt tumor progression.

## The Solution

**p53-proteoMgCAD** designs **second-site rescue mutations** — additional amino acid changes that compensate for the structural damage caused by cancer mutations and restore wild-type tumor suppressor function. This approach, called *intragenic suppression*, has been validated experimentally (Nikolova et al. 2000, Baroni et al. 2004, Otsuka et al. 2007).

The pipeline uses **ESM-2 protein language model embeddings** combined with a **functional oracle trained on experimental deep mutational scanning data** (Giacomelli et al. 2018) to navigate protein sequence space via gradient-based optimization, discovering rescue mutations that maximize predicted function while maintaining >90% sequence identity to wild-type human p53.

---

## Results

### Campaign: 108 Scenarios, 8 Cancer Hotspots, 3 Delivery Methods

| Metric | Value |
|--------|-------|
| Cancer hotspot targets | 10 unique target combinations |
| Rescue candidates designed | 23 shortlisted from 1,176 evaluated |
| Best oracle score | 1.714 (wild-type-like function) |
| Sequence identity to WT | >99% (1-5 mutations from wild-type) |
| Delivery methods covered | Gene therapy, mRNA therapy, Protein therapy |
| Target retention rate | 39.2% of deep candidates |

### Top Rescue Candidates

| Target Mutation | Rescue Mutation | Oracle Score | DMS Validation | Delivery |
|-----------------|-----------------|--------------|----------------|----------|
| R273H + R282W | G244C | 1.714 | Epistatic rescue | Gene therapy |
| R175H + R249S | S241R | 1.709 | Epistatic rescue | Gene therapy |
| R249S | M237R + S241R | 1.704 | Epistatic rescue | Gene therapy |
| R249S + R282W | G389E | 1.701 | F328I: Z = -0.31 (functional) | Gene therapy |
| R175H + R249S | S241H | 1.699 | Epistatic rescue | Gene therapy |
| G245S + R249S | V203D | 1.699 | Epistatic rescue | Gene therapy |
| R273H | F328I | 1.686 | **F328I: Z = -0.96 (strongly functional)** | Gene therapy |
| R273H + R282W | T118H | 1.667 | **T118H: Z = -0.10 (functional)** | mRNA therapy |
| R248Q | W91P | 1.598 | **W91P: Z = -0.57 (functional)** | Gene therapy |

Rescue mutations marked **bold** are independently functional in Giacomelli 2018 DMS data (negative Z-score = wild-type-like). Others work through **epistatic compensation** — they are harmful alone but restore function in combination with the cancer mutation.

### 3D-Printable Protein Model

The best rescue candidate (R273H+R282W rescued by G244C) is exported as a 3MF file for 3D printing, with cancer mutation sites colored red and rescue sites colored green. Generated from the AlphaFold wild-type p53 structure (AF-P04637-F1-model_v6).

---

## How It Works

### Functional Manifold Rescue (FMR) Algorithm

1. **Encode**: The cancer-mutant p53 sequence is embedded into ESM-2's 320-dimensional latent space (per position)
2. **Optimize**: Gradient ascent maximizes a multi-objective loss combining:
   - Functional oracle score (trained on 7,844 DMS variants)
   - Stability preservation
   - DNA binding maintenance
   - Sequence identity constraint (>90% to wild-type)
   - Target mutation lock (cancer mutations are retained, not "fixed")
3. **Decode**: The optimized embedding is decoded back to amino acid probabilities, yielding a rescued protein sequence

### Campaign System

The campaign runner evaluates all combinations of:
- **8 cancer hotspots** (singles and pairs = 36 target combinations)
- **3 delivery methods** (gene therapy, mRNA therapy, protein therapy)
- **3 optimization profiles** (Balanced, Stability-First, Binding-Optimized)

This produces 108 scenarios in Pass A (screening), then deep-refines the top 40% in Pass B with multiple random restarts and trials. A diversity-aware shortlist selects the top 30 candidates across targets, delivery methods, and optimization profiles.

### Experimental Validation

The functional oracle is trained on the **Giacomelli et al. 2018** saturation mutagenesis dataset — 7,844 p53 variants tested in a cell-based Nutlin-3 selection assay. This is real experimental data measuring whether each p53 variant retains tumor suppressor activity. The oracle achieves strong correlation with these experimental measurements and generalizes to multi-mutation combinations.

---

## Installation

### Requirements

- Python 3.11+
- PyTorch 2.0+
- ~8 GB RAM

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
# Run full campaign (108 scenarios, ~80 min on Apple MPS GPU)
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
│   ├── campaign.py          # Multi-scenario campaign runner (FMR optimizer)
│   ├── oracle.py            # Functional prediction model (ESM-2 + MLP)
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
│   └── functional_oracle.pt           # Trained oracle weights
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
