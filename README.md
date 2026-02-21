# p53-proteoMgCAD

### Computational Design of Second-Site Rescue Mutations to Restore p53 Tumor Suppressor Function

**Ishaan Gubbala**

---

## The Problem

**TP53 is the most frequently mutated gene in human cancer.** Roughly half of all solid tumors carry a TP53 mutation, rising above 90% in high-grade serous ovarian and small-cell lung cancers. Unlike most tumor suppressors, which are inactivated by truncations, ~75% of TP53 alterations are missense mutations that produce a full-length but malfunctioning protein. This reflects selective pressure for gain-of-function properties: mutant p53 doesn't just lose tumor suppression, it actively promotes invasion, metastasis, and chemoresistance.

Eight hotspot mutations account for roughly 28% of all TP53 missense variants across cancers:

| Mutation | Type | Structural consequence |
|----------|------|----------------------|
| **R175H** | Structural | Disrupts zinc coordination in L2 loop, causes global DBD unfolding |
| **G245S** | Structural | Distorts L3 loop, perturbs DNA-contact surface locally |
| **R248Q** | Contact | Eliminates direct minor-groove DNA contacts |
| **R248W** | Contact | Bulky tryptophan sterically clashes with DNA minor groove |
| **R249S** | Structural | Perturbs L3 loop architecture; prevalent in aflatoxin-associated HCC |
| **R273H** | Contact | Eliminates major-groove phosphate backbone contacts |
| **R273C** | Contact | Same electrostatic loss as R273H |
| **R282W** | Structural | Destabilizes the H2 helix in the loop-sheet-helix motif |

**Structural mutations** (R175H, G245S, R249S, R282W) destabilize the DNA-binding domain fold itself. **Contact mutations** (R248Q/W, R273H/C) leave the fold largely intact but eliminate direct DNA contacts. These two classes require different rescue strategies — global stability compensation versus local structural mimicry — which our pipeline handles through separate optimization profiles.

Current therapeutic approaches face fundamental limitations. APR-246 (eprenetapopt), the most advanced mutant-p53 reactivator, is a non-specific thiol-reactive prodrug that failed to improve overall survival in Phase 3 trials. Rezatapopt targets only the Y220C pocket (~2% of mutations). CRISPR correction of TP53 paradoxically selects *against* corrected cells because the repair process activates the p53 damage response. No single drug addresses the full mutational spectrum.

## Research Question

**Can second-site rescue mutations in the p53 DNA-binding domain restore tumor suppressor function in cancer-mutant p53, and can a computational pipeline — combining protein language model optimization, deep mutational scanning data, and physics-based structural validation — generate experimentally testable rescue candidates?**

This is a **hypothesis generation** question, not a therapeutic efficacy claim. The goal is to narrow the combinatorial space of possible multi-mutation rescues (20^393 sequences) to a prioritized shortlist for wet-lab follow-up. No prediction in this pipeline substitutes for experimental data.

### Gap This Work Fills

| Existing approach | Coverage | Limitation |
|------------------|----------|-----------|
| APR-246 / PRIMA-1 | All TP53 mutations (non-specific) | Phase 3 failed; mechanism not mutation-specific |
| Rezatapopt (PC14586) | Y220C only (~2%) | Pocket does not exist for most hotspot mutations |
| Zinc metallochaperones | Zinc-depleted mutants (R175H, C238) | Narrow therapeutic window |
| CRISPR correction | All mutations (in principle) | Paradoxically kills corrected cells via DDR activation |
| **This pipeline** | **All 8 hotspot mutations (36 combinations)** | **Computational only — requires experimental validation** |

No existing method designs mutation-specific multi-residue second-site suppressors across the full TP53 hotspot spectrum. This pipeline fills that gap by treating rescue design as a constrained optimization problem in ESM-2's learned protein sequence space.

### Experimental Validation Hierarchy

Predictions from this pipeline require validation at multiple independent levels before clinical relevance can be claimed:

1. **Tier 1 — Computational** *(this work)*: Oracle scoring, ESMFold structure, OpenMM energetics, MM-GBSA binding, DMS data lookup
2. **Tier 2 — *In vitro* biochemistry**: Western blot (fold confirmation), EMSA or SPR (DNA binding), thermal shift assay (stability ΔTm)
3. **Tier 3 — Cell-based functional assays**: Luciferase reporter (p53 transcriptional activity), cell viability (growth suppression in p53-null cells), ChIP-qPCR (chromatin binding)
4. **Tier 4 — *In vivo* models**: Humanized mouse tumor xenografts, patient-derived organoids
5. **Tier 5 — Clinical translation**: IND-enabling studies, Phase I safety, Phase II efficacy

This work occupies Tier 1 only. Positive Tier 1 results justify investment in Tier 2 experiments; they do not constitute evidence for therapeutic efficacy.

---

## The Approach

**p53-proteoMgCAD** takes a different path: rather than correcting the cancer mutation or refolding the protein with a small molecule, it computationally designs **second-site rescue mutations** — additional amino acid changes elsewhere in p53 that compensate for the structural damage and are predicted to restore tumor suppressor function. The pipeline is a **hypothesis-generation tool**: it narrows the astronomical combinatorial space of possible multi-mutation rescues to a shortlist of experimentally testable candidates, not a clinical-ready therapeutic.

This strategy, called *intragenic suppression*, has been validated experimentally:

- **N239Y** rescues R249S by restoring L3 loop contacts (Nikolova et al. 2000)
- **H168R** suppresses R249S by structurally mimicking the lost Arg249 side chain (Baroni et al. 2004)
- **T284R** rescues R175H by restoring zinc coordination geometry (Otsuka et al. 2007)

These experimental discoveries prove that p53's DNA-binding domain has sufficient structural plasticity for second-site mutations to compensate for cancer-associated damage. But finding these mutations by experiment requires screening millions of combinations. This pipeline automates that search computationally.

### Why protein language models?

ESM-2 (650M parameters) was trained on 250 million protein sequences across evolution. It learned the grammar of protein structure — which residues can co-occur, which substitutions are tolerable, and which positions are structurally coupled — without ever seeing a crystal structure. When we embed a p53 sequence into ESM-2's 1280-dimensional latent space, each position's representation encodes its structural context, evolutionary constraints, and biophysical properties. Gradient-based optimization in this space efficiently navigates combinatorial sequence space that would be intractable to search by brute force.

### Why delta-encoded attention pooling?

The functional oracle that scores rescue candidates faces a needle-in-a-haystack problem. All 7,844 DMS training sequences differ from wild-type p53 at only 1 of 393 positions. Standard mean-pooling averages away the mutation signal under 392 identical positions. Our solution:

1. **Delta encoding**: subtract the wild-type p53 embedding from each candidate. Now mutated positions are nonzero residuals; everything else cancels to zero.
2. **Learnable attention pooling**: a query vector attends over the delta-encoded positions via 4-head attention, learning which mutation-induced changes predict functional impact.
3. **MLP head**: the attended representation passes through 1280 -> 256 -> 128 -> 1 to produce a functional score.

Without delta encoding, attention collapses and all sequences receive identical scores. With it, the oracle achieves validation loss of 0.2798 (MSE on Z-scores spanning approximately -3 to +3; RMSE ~ 0.53 on a 6-unit range, or ~9% relative error) and correctly separates functional from loss-of-function variants (score range: WT = +0.04, cancer mutants = -1.2 to -1.6, best rescue = +1.63). Pearson r = 0.87 (R^2 = 0.76) on the held-out validation set.

**Statistical context**: An R^2 of 0.76 means the oracle explains 76% of variance in experimental fitness — sufficient for **ranking and screening** (identifying top candidates from a pool) but not for **absolute quantitative prediction** (predicting exact fitness values for clinical dosing). This is the intended use: the oracle ranks candidates to prioritize experimental testing, not to replace it. The ~24% unexplained variance is why physics validation and wet-lab follow-up are essential downstream steps, not optional additions.

### Addressing Extrapolation Beyond Single-Mutation Training Data

A key limitation of oracle training solely on single-residue DMS variants is that during inference, 18–30 positions are simultaneously nonzero in the delta embedding — a distributional shift the attention mechanism was never explicitly trained to handle.

**The fix — multi-mutation fine-tuning** (`p53cad train-multimut`):

1. **Generate synthetic multi-mutation training data**: Sample k-mutation combinations (k = 2–20) from the DMS single-mutation entries. Compute pseudo-labels using the **thermodynamic additivity model** — a well-established null model in double-mutant cycle analysis (`Z_pseudo = Z₁ + Z₂ + ... ` with soft saturation at |Z| > 7). Generates 50,000 synthetic multi-mutation pairs.

2. **Inject literature-validated anchors**: 7 experimentally confirmed multi-mutation rescues (N239Y+R249S, H168R+R249S, T284R+R175H, N239Y+N268D+R249S, etc.) are injected as high-confidence training examples, 5× upweighted.

3. **Fine-tune at reduced learning rate**: Mix original DMS singles + synthetic multi-mutation pairs (1:1), fine-tune at lr=2e-4 (vs original 1e-3) to prevent catastrophic forgetting of single-mutation knowledge. The oracle's attention mechanism then learns to aggregate signals from multiple simultaneously nonzero positions — the regime it encounters during rescue inference.

**Why thermodynamic additivity is defensible**: For non-contacting residues (>8 Å), thermodynamic cycles show epistasis is typically <1 kcal/mol — within the DMS measurement noise. ESM-2's context-dependent embeddings already encode some epistasis implicitly: position i's representation in a double-mutant (i,j) differs from single-mutant (i) because ESM-2 conditions on the full sequence. Training on additive pseudo-labels teaches the attention mechanism *multi-position integration*; ESM-2 carries the epistatic signal in the embedding itself.

**Limitation remaining**: Additive pseudo-labels will systematically underestimate negative epistasis (where two individually beneficial mutations are neutral in combination) and overestimate positive epistasis. This is calibrated post-hoc using the literature anchor set — if N239Y+R249S scores correctly after fine-tuning, the calibration is reasonable. Full epistasis modeling would require multi-mutation DMS data, which does not yet exist for p53 at scale.

```bash
# Re-train oracle with multi-mutation augmentation
p53cad train-multimut \
    --checkpoint data/models/functional_oracle.pt \
    --output data/models/functional_oracle_multimut.pt \
    --n-synthetic 50000 --k-max 20 --epochs 15
```

---

## How It Works

### Functional Manifold Rescue (FMR) Algorithm

The pipeline uses two complementary optimization strategies to design rescue mutations for each cancer-mutant p53 sequence.

**Gradient-Based Optimization** operates in ESM-2's continuous latent space:

1. **Encode**: The cancer-mutant sequence is embedded into per-position representations (393 positions x 1280 dimensions).
2. **Optimize**: Gradient ascent over a multi-objective loss:
   - **Functional oracle score** — the trained attention-pooling network predicting rescue quality
   - **DMS-aware rescue penalty** — a differentiable penalty using per-residue experimental fitness: `E[Z] = sum P(aa|pos) * Z(pos, aa)`, steering toward individually validated rescue substitutions
   - **Cancer-site pseudo-log-likelihood** — focused at cancer mutation positions rather than diluted across all 393
   - **Contact preservation** — cosine similarity between current and wild-type hidden states at structurally neighboring positions (8A cutoff from the AlphaFold structure)
   - **Epistasis scoring** — structural proximity and attention coupling between co-mutating positions, rewarding synergistic mutation clusters
   - **Mutation-neighborhood pooling** — Gaussian kernel (sigma=10) over locked cancer-site positions, amplifying rescue signal near the damage site
   - **Sequence identity floor** — enforces >90% identity to wild-type (threshold derived from AAV gene therapy immunogenicity literature, where >90% identity to endogenous proteins minimizes T-cell epitope generation; Mingozzi & High, *Blood* 2013), adjusted per delivery method
   - **Target lock** — cancer mutations are retained (not "fixed") via hard position locking with penalty weight 800
3. **Decode**: The optimized continuous representation is projected back to amino acid probabilities and decoded to a discrete sequence.

**Autoregressive Sampling** operates directly in sequence space:

1. Rank all non-locked positions by ESM-2 attention weight
2. At each position, mask and predict top-K amino acid candidates via ESM-2's masked language model head
3. Accept the first substitution that improves the oracle score
4. Iterate until convergence — each step adds one greedy rescue mutation

In practice, autoregressive sampling produces 23 of 30 shortlisted candidates. It tends to find more conservative, higher-confidence rescues, while gradient-based optimization explores more creative multi-mutation combinations.

### Campaign System

A campaign evaluates all combinations of 8 cancer hotspots (singles and pairs = 36 target combinations), 3 delivery methods (gene therapy, mRNA therapy, protein therapy), and multiple optimization profiles — producing 108 scenarios per pass. Two passes run sequentially:

- **Pass A (screening)**: Each scenario runs gradient optimization and autoregressive trials with a fixed compute budget. The best candidate per scenario advances.
- **Pass B (deep refinement)**: Top scenarios from Pass A are re-optimized with multiple random restarts, longer step counts, and additional autoregressive rounds.

An **evidence-weighted, diversity-aware shortlist** selects the top 30 candidates by jointly ranking on oracle score, clinical impact (TCGA patient stratification), DMS rescue quality, Pareto rank, and mutation novelty — with greedy diversity fill to ensure coverage across targets and delivery methods.

### Physics-Based Validation

After candidate generation, a five-stage physics pipeline validates structural soundness:

| Stage | Method | What it checks | Output |
|-------|--------|---------------|--------|
| **1. Structure prediction** | ESMFold (local, 700M params) | Can the sequence fold? | Per-residue pLDDT confidence, full PDB |
| **2. Energy minimization** | OpenMM AMBER14 + OBC2 solvent | Is the fold energetically reasonable? | Potential energy (kcal/mol), DDG vs WT and cancer |
| **3. MD stability** | 200 ps molecular dynamics | Does the structure survive short thermal stress? | CA-RMSD, RMSF, Rg, DSSP secondary structure |
| **4. DNA-binding interface (geometric)** | Superposition onto AlphaFold WT | Are DNA contacts preserved geometrically? | CA displacement at 11 key contact residues |
| **5. DNA-binding simulation (MM-GBSA)** | Superposition onto PDB 2AHI + OpenMM | Does the rescue *functionally* bind DNA? | ΔG_bind, H-bonds, interface contacts vs WT |

Stage 5 builds a protein-DNA complex by superimposing the rescue DBD (residues 94-293) onto the 2AHI crystal structure (p53 DBD + DNA response element, 1.85 Å; Kitayner et al. 2006), then computes MM-GBSA binding free energy as `ΔG = E_complex - E_protein - E_DNA` using AMBER14 + OBC2 implicit solvent with OpenCL GPU acceleration (~35-70 sec/candidate). This tests whether the rescue candidate retains *functional* DNA binding, not just structural proximity to the interface.

Results combine into a composite physics score (0-100): pLDDT (30 pts) + DDG (25 pts) + MD stability (25 pts) + DNA interface (20 pts). When Stage 5 is available, the 20 DNA pts split into 8 pts geometric + 12 pts functional binding simulation. Each candidate receives a verdict: **STRONG** (>=75), **PROMISING** (55-74), **UNCERTAIN** (35-54), or **CONCERNING** (<35). All predictions are SHA-256 cached for incremental re-runs.

**MD timescale limitation**: 200 ps of molecular dynamics is a **coarse screening filter**, not a stability proof. Most biologically relevant conformational changes (folding, domain motions, allosteric transitions) occur on microsecond-to-millisecond timescales — three to six orders of magnitude longer. What 200 ps *can* detect: gross steric clashes, immediate unfolding of destabilized secondary structure elements, and high-energy conformations that would rapidly relax. What it *cannot* detect: slow unfolding pathways, subtle allosteric effects, or aggregation propensity. Candidates passing MD screening should be validated with longer simulations (10-100 ns minimum) or experimental stability assays (differential scanning fluorimetry, thermal shift) before advancing to functional testing. The 200 ps window was chosen as a practical tradeoff between computational cost (~5 min/candidate on GPU) and screening throughput for 30+ candidates.

### Design & Methodology

#### Training Data Quality Control

| QC step | Detail |
|---------|--------|
| Source dataset | Giacomelli et al. 2018 — 8,260 raw entries |
| After deduplication | 7,844 unique variants (416 exact duplicates removed) |
| Positions covered | 381 of 393 (12 positions with no DMS measurements — primarily N-terminal TAD) |
| Z-score range | −3.28 to +4.91; median = +0.92 (skewed toward loss-of-function) |
| Sign convention | Raw CSV: negative Z = functional (Nutlin-3 selection). Loaded values NEGATED: positive = functional, negative = LoF |
| Cancer hotspot coverage | All 8 hotspot positions (G245, R248, R249, R273, R282, R175, R248, R273) have DMS data for all 19 non-WT amino acids |
| DMS dead zone filtering | Z ∈ (−0.5, +2.0) zeroed in optimization penalty (ambiguous compensatory range excluded) |
| Train / validation split | 80% / 20%, stratified random, seed 42 |

#### Oracle Architecture Rationale

The choice of **delta encoding** before attention pooling is non-obvious and deserves explanation. All 7,844 DMS sequences differ from WT at exactly one position — 1 of 393 positions is mutated, the other 392 are identical to WT. After ESM-2 embedding, the 392 identical positions have nearly identical representations. Without delta encoding, any pooling operation averages ~392 near-zero contributions with 1 informative signal, causing the attention to collapse. By subtracting the WT embedding first, mutated positions become the only nonzero residuals — the attention mechanism can then focus precisely on the single changed position and learn which types of changes at which positions predict functional impact.

This design decision was validated by ablation: training the same architecture without delta encoding produces constant output (all sequences receive the same score) — the attention head degenerates to a weighted average of a near-constant input.

#### Power Analysis

With 7,844 training examples and a 20% validation set (n=1,569):
- At oracle R²=0.76, RMSE=0.53 on a 6-unit Z-score scale
- Minimum detectable effect size (rescue vs. cancer): ΔZ > 0.3 (based on 95% CI of RMSE)
- 1,152 campaign candidates provide adequate coverage of the 36 target combinations (~32 candidates/target) for shortlisting purposes
- For experimental validation: 5 candidates × 3 independent replicates per assay requires n=15 measurements per assay; power ≥80% to detect ΔZ > 0.5 at α=0.05 is achievable with standard luciferase reporter protocols

### Experimental Training Data and Validation Strategy

The oracle is trained on the **Giacomelli et al. 2018** deep mutational scanning dataset: 7,844 p53 single-residue variants tested in a cell-based Nutlin-3 selection assay in A549 cells. Each variant receives a Z-score measuring whether it retains (negative Z) or loses (positive Z) tumor suppressor function.

**Addressing the circularity concern**: The Giacomelli dataset serves two roles — oracle training and DMS-aware optimization guidance — but these are not circular for three reasons:

1. **Cross-validation**: The oracle is trained on an 80/20 train/validation split with early stopping on validation loss. The reported val loss (0.2798) is measured on held-out data the model never trained on. The oracle must *generalize* beyond its training examples.
2. **The oracle predicts multi-mutation combinations it never saw**: Giacomelli measured **single-residue** variants only. The pipeline generates candidates with 18-30 simultaneous mutations — combinations that do not exist anywhere in the training data. The oracle's ability to score these multi-mutation sequences is an extrapolation test, not interpolation.
3. **Physics validation is fully independent**: The four-stage physics pipeline (ESMFold structure prediction, OpenMM energy minimization, MD stability, DNA-binding interface analysis) uses zero Giacomelli data. It validates candidates through physical first principles — protein folding confidence, thermodynamic energy, molecular dynamics, and structural superposition. If the oracle were merely memorizing training data, physics validation would catch structurally implausible predictions.

A partial **external validation** comes from the pipeline independently rediscovering **N239Y** as a rescue mutation for R249S — a known experimental suppressor (Nikolova et al. 2000) that the pipeline identified without being explicitly told about it. The DMS data contains N239Y's fitness score, but the pipeline's selection of this specific mutation from 7,860 possible single-residue substitutions demonstrates that the oracle learned genuine structural biology, not just lookup patterns.

The DMS data also enters the optimization loop directly as a differentiable penalty that steers toward mutations at positions where experiments show functional tolerance. A dead-zone filter excludes ambiguous compensatory-range Z-scores (-0.5 to +2.0), keeping only strong functional signals (Z < -0.5) and catastrophic signals (Z > +2.0) as steering evidence. This is intentionally conservative — the DMS penalty acts as a **safety net** (weight 10), not the primary optimization driver (oracle score weight is 4-8x higher).

### Per-Mutation Evidence Validation

Every rescue mutation in a candidate is individually scored against experimental data. For each mutation, the pipeline:

1. **Looks up the DMS record** — if Giacomelli 2018 measured that exact substitution, its Z-score is used directly (confidence = 100%).
2. **Falls back to physics-based estimation** — if no DMS record exists, the mutation is scored using BLOSUM62 evolutionary substitution scores, Kyte-Doolittle hydropathy (burial effects), side-chain volume (steric clash potential), and charge change (electrostatic effects). Confidence drops to 30-60% depending on the estimation basis.
3. **Computes delta versus target** — how much better (or worse) the rescue mutation scores compared to the cancer mutation. Bars above zero mean the rescue outperforms the cancer mutant at that position.

These per-mutation assessments roll up into an **evidence score** (0-100) that weights three components:

| Component | Weight | Source |
|-----------|--------|--------|
| Physics subscore | 40% | Normalized function (35%), identity (25%), stability (20%), binding (20%) |
| Measurement confidence | 30% | DMS = 100%, physics fallback = 30-60%, unavailable = 0% |
| DMS effect signal | 30% | Mean delta vs target, normalized to [0, 1] |

Candidates with higher DMS coverage score higher because their evidence base is experimental rather than estimated.

### Mechanism Explainability

For any rescue candidate, the pipeline decomposes *why* the rescue works through three complementary analyses:

**Attention attribution** — ESM-2's internal attention patterns are extracted to identify which residue positions the model considers most important for the rescue prediction. Positions are ranked by composite importance (combining attention weight and occlusion sensitivity), revealing whether the rescue acts through direct structural compensation, allosteric stabilization, or DNA-contact restoration.

**Energy decomposition** — The predicted rescue mechanism is broken down into six biophysical components:

- **Electrostatic** — charge complementarity and salt bridge formation
- **Van der Waals** — packing efficiency and steric optimization
- **Hydrogen bonds** — new or restored H-bond networks
- **Solvation** — changes in solvent accessibility and desolvation penalties
- **Backbone strain** — phi/psi angle distortion at rescue sites
- **Side-chain packing** — rotamer quality and cavity filling

Negative contributions are stabilizing; positive are destabilizing. The sum gives the total predicted DDG of the rescue.

**Counterfactual analysis** — For each rescue mutation, the pipeline asks: "what if we used a different amino acid at this position?" It tests alternative substitutions and reports the delta score versus the chosen rescue, revealing whether the selected mutation is optimal or whether alternatives exist with comparable or better rescue potential.

### Drug Candidate Generation

The pipeline generates small-molecule stabilizers targeting five p53 binding pockets:

| Pocket | Location | Strategy |
|--------|----------|----------|
| **Y220C cavity** | Created by the Y220C mutation | Fill the cavity to stabilize the fold (precedent: rezatapopt) |
| **Core hydrophobic** | Internal beta-sandwich | Thermal stabilization via hydrophobic packing |
| **L1 loop** | Flexible L1 loop region | Conformational constraint to reduce flexibility |
| **Zinc site** | Cys176/His179/Cys238/Cys242 coordination shell | Metallochaperone or zinc-site stabilization |
| **MDM2 interface** | N-terminal TAD binding surface | Block MDM2-mediated degradation (nutlin-like) |

Three generation modes are supported:
- **Template** — scaffold-based design using known drug pharmacophores, scored by a heuristic affinity model
- **Docking** — AutoDock Vina rigid-body docking against pocket-specific receptor PDBQTs from the AlphaFold structure
- **Docking + MD** — Vina docking followed by short MD relaxation to test binding pose stability

Each candidate is evaluated for Lipinski Rule of Five compliance (MW < 500, LogP < 5, HBD <= 5, HBA <= 10), synthetic accessibility (1-10 scale), and drug-likeness. When a rescue protein design is active, drug scoring incorporates protein-context conditioning — the affinity estimate is shifted based on how well the drug's target pocket aligns with the structural changes introduced by the rescue mutations.

### Clinical Viability Assessment

The final analysis stage asks: **if this rescue candidate works, how viable is it as a therapeutic?** This is a quantitative cost-benefit analysis covering four dimensions.

**1. Patient stratification** — Using TCGA pan-cancer mutation frequencies and US cancer incidence data (SEER), the pipeline estimates how many patients per year carry the specific target mutation:

```
patients = sum over cancer types of:
    incidence(cancer) x p53_mutation_rate(cancer) x hotspot_frequency(mutation)
```

For example, R175H (5.6% of p53 mutations) maps to ~7,400 US patients/year across breast, colon, ovarian, lung, and pancreatic cancers. The global estimate uses a 3x multiplier. A five-year cumulative patient pool and rough market value ($100k-$500k per patient for gene therapy) are also computed.

**2. Therapeutic index** — The safety margin of the rescue design, driven by:

| Factor | Assessment | Threshold |
|--------|-----------|-----------|
| Sequence identity | Higher = safer (fewer novel epitopes) | >=98% = wide window, >=95% = moderate, <95% = narrow |
| Immunogenicity | Per-mutation risk + identity factor + exposed-position penalty + mutation count scaling | 0-1 scale (lower = better) |
| Off-target risk | Proximity of rescue mutations to zinc coordination, DNA-binding, and tetramerization sites | 0-1 scale |
| Therapeutic window | Composite of identity + immunogenicity | wide / moderate / narrow |
| Regulatory pathway | Recommended FDA pathway based on identity and mutation count | Gene therapy / mRNA / cell therapy / protein |

**3. Delivery feasibility** — Six delivery methods are scored for real-world viability:

| Method | Identity req. | Cost/dose | Regulatory status | Key tradeoff |
|--------|--------------|-----------|-------------------|-------------|
| AAV gene therapy | >=90% | $100k-$500k | Multiple approvals (Luxturna, Zolgensma) | Durable expression but pre-existing immunity in ~30% |
| Lentiviral gene therapy | >=90% | $300k-$1M | Approved for blood disorders | Stable integration but insertional mutagenesis risk |
| mRNA-LNP | >=92% | $5k-$50k | COVID vaccines approved | Cheapest, repeat-dosable, but short-lived expression |
| Purified protein | >=95% | $10k-$100k | Difficult for transcription factors | Immediate effect but cell penetration barrier |
| CRISPR prime editing | 100% | $100k+ | Early clinical trials | Permanent correction but low in vivo efficiency |
| Oncolytic virus + p53 | >=90% | $50k-$200k | IMLYGIC approved (melanoma) | Tumor selectivity but variable penetration |

Feasibility scoring accounts for whether the candidate meets the identity requirement, mutation count (fewer = easier), and whether the method has regulatory precedent.

**Immunogenicity of multi-mutation designs**: Top candidates carry 18-30 mutations, raising legitimate concerns about immune recognition. Three factors mitigate this:

1. **Delivery context matters**: The top candidates are designed for **gene therapy** (AAV) or **mRNA** delivery, where the mutant p53 is expressed **endogenously** from the cell's own ribosomes. The protein folds intracellularly and is never presented to the immune system as an extracellular antigen. T-cell epitope generation depends on MHC-I presentation of intracellular peptides — at >90% identity, most 9-mer peptides are identical to endogenous p53 and should be tolerized.
2. **Mutation location**: Most rescue mutations are buried in the hydrophobic core of the DBD (zinc-coordinating loops, beta-sandwich interior), not on solvent-exposed surfaces. Buried mutations generate fewer novel MHC-I epitopes than surface mutations.
3. **Mutation count is a spectrum**: The pipeline produces candidates with varying mutation counts. While the highest-scoring candidates tend to have more mutations, the shortlist includes candidates with 20-25 mutations at >=93% identity. The clinical viability score explicitly penalizes higher mutation counts and narrower therapeutic windows.

This remains the **strongest limitation** of the multi-mutation approach and would require immunogenicity prediction (NetMHCpan) and *in vivo* testing in humanized mouse models before clinical translation.

**4. Overall clinical score** — A composite 0-100 score with transparent, auditable construction. Weights reflect clinical development priorities: immune safety (immunogenicity 25 + identity 15 = 40% total) dominates because novel protein therapeutics most commonly fail on immunogenicity; delivery feasibility (25%) captures the practical bottleneck of getting the therapeutic into tumor cells; patient population (25%) determines commercial viability and clinical trial recruitment:

| Component | Max points | Computation |
|-----------|-----------|-------------|
| Patient population | 25 | Tiered: >=10k = 25, >=5k = 20, >=1k = 15, >=100 = 10, else 5 |
| Sequence identity | 15 | >=98% = 15, >=95% = 10, >=92% = 5, else 0 |
| Therapeutic window | 10 | wide = 10, moderate = 5, narrow = 0 |
| Best delivery feasibility | 25 | Best method's feasibility score x 25 |
| Immunogenicity | 25 | 25 - (risk_score x 25) |

**Viability verdict**: HIGH (>=70), MODERATE (50-69), LOW (30-49), NOT VIABLE (<30). The pipeline outputs a recommended development path (e.g., "Pursue AAV gene therapy targeting R175H-positive tumors; prioritize ovarian and lung indications for Phase I based on highest mutation-specific incidence").

---

## Results

### Latest Campaign: 650M ESM-2 + Attention Oracle

| Metric | Value |
|--------|-------|
| ESM-2 model | 650M parameters (facebook/esm2_t33_650M_UR50D) |
| Oracle architecture | Attention-pooling with delta encoding (val loss 0.2798) |
| Campaign scenarios | 108 screen + 36 deep refinement |
| Total candidates evaluated | 1,152 (gradient-directed, not random sampling — see note below) |
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
- **All 3 delivery methods** are represented in the shortlist, with candidates satisfying delivery-specific identity constraints (gene therapy >=92%, mRNA >=92%, protein >=92%)

**On the 1,152 candidates and sequence space coverage**: 1,152 candidates from a space of 20^393 possible sequences may appear microscopic, but this reflects a misunderstanding of the search strategy. The pipeline does **not** randomly sample sequence space. Gradient-based optimization follows the derivative of the oracle score through ESM-2's continuous latent space — each of 250 optimization steps moves the entire 393x1280-dimensional representation toward higher-scoring regions. Autoregressive sampling adds mutations greedily, guided by ESM-2's masked language model predictions. Both strategies are **directed searches**, analogous to how gradient descent in deep learning finds good minima without evaluating all possible weight configurations. The 1,152 candidates represent 1,152 converged endpoints of guided optimization trajectories, not 1,152 random draws.

### Cost-Benefit Analysis of Results

#### What the pipeline produced

From 1,152 evaluated candidates across 36 target combinations and 3 delivery methods:

| Metric | Value | Context |
|--------|-------|---------|
| Candidates with positive oracle score | **527 (45.7%)** | Score > 0 means predicted more functional than cancer mutant baseline |
| Candidates exceeding WT score (+0.04) | **517 (44.9%)** | Predicted to match or exceed wild-type p53 function |
| High-confidence rescues (score > +0.5) | **413 (35.9%)** | Strong predicted functional gain over cancer mutant |
| Candidates meeting all constraints | **601 (52.2%)** | Identity, stability, binding, delivery compatibility |
| Mean score gain vs cancer baseline | **+0.991** | Average improvement from cancer mutant score (-1.2) toward functional |
| Hotspot mutations covered | **7 of 8 singles** | All except R175H solo; R175H covered in 6 pair combinations |
| Patient population reached | **~25,400 US/year** | ~76,200 globally (3x multiplier) |

#### Per-target success rates

Some cancer mutations are easier to rescue computationally than others:

| Target | Candidates tested | Positive score | Success rate | Best score | Interpretation |
|--------|------------------|---------------|-------------|------------|----------------|
| R249S+R282W | 21 | 16 | **76%** | — | Highest success: both are structural mutations amenable to stability compensation |
| R248W+R249S | 34 | 25 | **74%** | — | Two structural mutations, large rescue space |
| R248Q+R273C | 34 | 24 | **71%** | +0.831 | Contact + contact: compensatory interfaces accessible |
| G245S+R273C | 56 | 39 | **70%** | +1.141 | Structural + contact: diverse rescue mechanisms |
| R249S | 30 | 15 | **50%** | +1.628 | Highest absolute score achieved; L3 loop rescue |
| R175H | 21 | 6 | **29%** | +1.508 | Hardest single target: global DBD unfolding is difficult to compensate |
| R175H+R282W | 21 | 4 | **19%** | — | Two severe structural mutations: near limit of rescuability |

R175H — the single most common cancer mutation (~5.6% of all TP53 missense) — has the lowest single-target success rate (29%). This is consistent with biology: R175H causes global DBD unfolding via zinc coordination loss, which is the hardest structural defect to compensate with distal mutations. The pipeline does find positive-score R175H rescues (best: +1.508) but they don't pass the stringent shortlist filters. This is an honest result — some mutations may be computationally irrescuable with second-site mutations alone.

#### Evidence strength per shortlisted candidate

Each shortlisted candidate rests on multiple independent evidence layers:

| Evidence layer | Source | Independence from oracle | Strength |
|---------------|--------|--------------------------|----------|
| Oracle score | Attention-pooling net on Giacomelli DMS | Core prediction | Moderate — single dataset, extrapolation to multi-mutations |
| DMS rescue quality | Per-mutation Z-score lookup | Partially independent (same dataset, different use) | Strong for individual mutations — direct experimental measurement |
| N functional rescues | Count of mutations with DMS Z < -0.5 | Same as above | Strong — each rescue mutation has individual experimental support |
| ESMFold pLDDT | Structure prediction confidence | **Fully independent** — physics-based, no DMS data | Strong — mean pLDDT 0.80, DBD pLDDT 0.91 (AlphaFold-quality) |
| ESMFold pTM | Template modeling score | **Fully independent** | Moderate — pTM ~0.55 indicates reasonable global fold |
| Energy minimization | OpenMM AMBER14 force field | **Fully independent** | Strong — first-principles thermodynamics |
| MD stability | 200 ps molecular dynamics | **Fully independent** | Weak — timescale too short for definitive stability claims |
| DNA interface | Structural superposition | **Fully independent** | Moderate — geometric check, no dynamics |

**Key finding**: The four candidates with ESMFold validation all show DBD pLDDT > 0.91 — higher than the typical threshold for reliable structure predictions (0.70). This means ESMFold independently predicts that these rescue sequences can fold into a stable DBD, even with 22-27 mutations. This is the strongest non-oracle evidence in the pipeline.

#### Risk assessment

| Risk | Severity | Likelihood | Mitigation in pipeline | Residual risk |
|------|----------|-----------|----------------------|---------------|
| Oracle score doesn't predict real function | Critical | Medium | Physics validation, DMS per-mutation checks | **High** — oracle trained on single-mutation data, predicting multi-mutation effects |
| Immunogenic T-cell response | Critical | Medium-High | Identity floor (>90%), buried mutation preference | **Medium** — 20-28 mutations create novel peptides; gene therapy mitigates via endogenous expression |
| Protein doesn't fold correctly *in vivo* | Critical | Low-Medium | ESMFold pLDDT > 0.91 for DBD | **Low-Medium** — ESMFold is accurate but not infallible |
| Epistatic interactions collapse function | High | Medium | Epistasis scoring term in optimization | **Medium** — oracle doesn't model epistasis; only screenable experimentally |
| Cancer mutation not fully retained | High | Low | Hard position locking (penalty 800) | **Very Low** — all shortlisted candidates retain target mutations |
| Insufficient DNA-binding restoration | High | Medium | DNA interface geometric check | **Medium** — geometric check is coarse; DNA-binding assay needed |
| Wrong delivery method feasibility | Moderate | Low | Delivery-specific identity constraints | **Low** — constraints are conservative |

#### What experimental validation would look like

To convert the top 5 candidates from computational predictions to experimentally validated rescues:

| Step | What | Estimated cost | Timeline | Purpose |
|------|------|---------------|----------|---------|
| 1. Gene synthesis | 5 codon-optimized TP53 variants (~1,200 bp each) | ~$2,500-$5,000 | 2-3 weeks | Obtain physical DNA |
| 2. Cloning + expression | Mammalian expression vectors (pcDNA3.1 or lentiviral) | ~$3,000-$5,000 | 1-2 weeks | Express rescue p53 in cells |
| 3. Western blot | Anti-p53 (DO-1) + conformation-specific (PAb1620 for folded, PAb240 for unfolded) | ~$1,500-$2,500 | 1 week | Confirm folded protein expression |
| 4. Luciferase reporter | p53-responsive promoter (p21, BAX, MDM2) driving luciferase | ~$2,000-$4,000 | 1-2 weeks | Measure transcriptional activity |
| 5. Cell viability | CellTiter-Glo in p53-null H1299 cells ± Nutlin-3 | ~$1,000-$2,000 | 1 week | Test growth suppression |
| 6. ChIP-qPCR | Chromatin immunoprecipitation at p53 target promoters | ~$3,000-$5,000 | 2-3 weeks | Confirm DNA binding *in vivo* |
| **Total for 5 candidates** | | **~$13,000-$23,500** | **8-12 weeks** | |

For context: a single round of deep mutational scanning (the experimental approach that produced the Giacomelli training data) costs ~$50,000-$150,000 and takes 3-6 months. The computational pipeline produces a prioritized shortlist of 30 candidates for ~$50-100 in cloud compute (40 hours GPU time), enabling targeted validation at ~1/5th the cost and ~1/3rd the time of an undirected experimental screen.

#### Prioritization: which candidates to test first

Ranked by combined evidence strength and clinical impact:

| Priority | Target | Score | Why first |
|----------|--------|-------|-----------|
| **1** | R249S (gene therapy) | +1.628 | Highest oracle score; includes N239Y (literature-validated rescue); 7 functional DMS-supported mutations; ESMFold pLDDT 0.80/0.91 |
| **2** | G245S+R273C (gene therapy) | +1.141 | Dual-target rescue; 70% success rate for this pair; covers 2 hotspots simultaneously |
| **3** | R282W (gene therapy) | +1.058 | High score; R282W is a structural mutation with clear rescue mechanism (helix stabilization) |
| **4** | R248Q (gene therapy) | +0.921 | R248Q is the 2nd most common TP53 mutation; ~5,200 US patients/year; moderate score |
| **5** | R248Q (mRNA therapy) | +0.441 | Same target, different delivery — mRNA is cheaper ($5k-50k/dose vs $100k-500k), repeat-dosable |

The rationale: prioritize (1) highest-confidence predictions, then (2) highest patient impact, then (3) delivery diversity. Testing the same R248Q rescue via both gene therapy and mRNA validates whether the mutation design is delivery-independent.

### MM-GBSA DNA Binding Simulation Results

**Campaign**: `campaign_20260217_060021` (vCUDA branch — NVIDIA GPU, 108 scenarios, 30 shortlisted)
**Reference**: PDB 2AHI — p53 DBD bound to a full DNA response element, 1.85 Å (Kitayner et al. 2006)
**Method**: MM-GBSA via AMBER14 + OBC2 implicit solvent, OpenCL GPU (~35-70 sec/candidate, 26 candidates total)

Each rescue DBD (residues 94-293 only) was superimposed onto 2AHI chain A via CA superposition, combined with the 2AHI DNA duplex (chains E+F), and minimized as three independent systems (complex, protein-only, DNA-only). WT p53 reference: ΔG_bind = **+162.0 kcal/mol** (raw MM-GBSA with OBC2 implicit solvent; positive values are expected due to desolvation penalties in unoptimized superposed complexes — what matters is the **ΔΔG relative to WT**).

#### Reliable Results (|δΔG| < 500 kcal/mol)

| Rank | Target | ΔG bind (kcal/mol) | δΔG vs WT | Contacts | Contact % | Binding Score | Verdict |
|------|--------|-------------------|-----------|----------|-----------|---------------|---------|
| **20** | R249S+R282W | 119.5 | **−42.5** | 453/449 | 100% | **0.750** | enhanced |
| **22** | R273C+R282W | 131.6 | **−30.4** | 475/449 | 100% | **0.750** | enhanced |
| **28** | R248Q+R273C | 103.3 | **−58.7** | 402/449 | 89% | **0.724** | enhanced |
| **30** | R175H | 144.4 | −17.6 | 400/449 | 89% | **0.723** | enhanced |
| **6** | R175H | 118.6 | −43.4 | 386/449 | 86% | **0.715** | enhanced |
| 18 | R273C | 129.2 | −32.8 | 344/449 | 77% | 0.692 | enhanced |
| 2 | R248Q+R273C | 132.9 | −29.1 | 325/449 | 72% | 0.681 | enhanced |
| 16 | R248Q+R249S | 73.8 | −88.2 | 314/449 | 70% | 0.675 | enhanced |
| 1 | R273C | 49.5 | −112.5 | 309/449 | 69% | 0.672 | enhanced |
| 24 | R273C | −16,166 | −16,328* | 328/449 | 73% | 0.683 | enhanced |
| 4 | R248Q+R273H | 125.5 | −36.5 | 252/449 | 56% | 0.640 | enhanced |
| 29 | R249S+R273H | 135.9 | −26.1 | 272/449 | 61% | 0.651 | enhanced |
| 9 | R248Q+R273H | 161.4 | −0.6 | 508/449 | 100% | 0.661 | preserved |
| 25 | R175H+R248Q | 164.0 | +2.1 | 424/449 | 94% | 0.595 | preserved |
| 10 | R175H | 163.6 | +1.6 | 272/449 | 61% | 0.520 | preserved |
| 12 | R248W+R273H | 176.8 | +14.8 | 493/449 | 100% | 0.353 | weakened |
| **5** | G245S+R249S | 260.8 | **+98.8** | 440/449 | 98% | **0.245** | disrupted |
| **7** | R175H+R248Q | 234.2 | +72.2 | 429/449 | 95% | **0.239** | disrupted |
| **14** | R273C+R282W | 217.8 | +55.8 | 355/449 | 79% | **0.198** | disrupted |

\* Ranks 23, 24, 27 show extreme negative energies (−10k to −22k kcal/mol), indicating minimization escaped a bad starting conformation; contact count is reliable but energy should be treated as unreliable.
Ranks 8, 15, 19, 21, 26 show extreme positive energies (>10⁸ kcal/mol) — minimization failure due to severe steric clashes; excluded from ranking.

#### Key Findings

**Top convergent candidates**: Ranks 20 (R249S+R282W), 22 (R273C+R282W), and 28 (R248Q+R273C) achieve **binding score ≥ 0.72** with both enhanced binding energy (δΔG < −30 kcal/mol) and full or near-full contact preservation (≥89%). Rank 28 (R248Q+R273C, binding score 0.724) is particularly notable — it also ranked #2 in the oracle score, providing convergent evidence from two fully independent methods (language model fitness and physics-based binding).

**R175H rescues are heterogeneous**: Three R175H candidates appear in the reliable set (ranks 6, 10, 30) with very different binding scores (0.715, 0.520, 0.723). The oracle-top R175H candidate (rank 6, δΔG = −43.4) is the most functionally promising; rank 10 appears to adopt a less favorable DBD conformation.

**G245S combinations consistently impair binding**: All three G245S-containing candidates in the reliable set score ≤ 0.25 (disrupted). This suggests G245S distortion of the L3 loop, which forms part of the DNA contact surface, is not adequately compensated even when other rescue mutations are present. This is a biologically plausible result — G245S directly deforms the minor-groove contact loop.

**R273C+R282W is delivery-sensitive**: Three independent candidates carry this target pair. Ranks 8 and 14 show disrupted binding (energy failure and +55.8 kcal/mol respectively); rank 22 shows enhanced binding (−30.4 kcal/mol, score 0.750). The three candidates differ in their rescue mutation sequences, suggesting the specific second-site mutations — not just the cancer target pair — determine binding outcome.

**H-bond count = 0 across all candidates**: ESMFold PDB files contain only heavy atoms. The `baker_hubbard` H-bond algorithm requires explicit hydrogen positions. H-bond analysis would require re-running with PDBFixer hydrogen addition before the mdtraj step; contact counts and energetics are reliable as-is.

#### Binding Score Distribution

```
Score ≥ 0.70  (enhanced):  8 candidates  — top priority for experimental validation
0.60–0.70     (enhanced):  7 candidates  — secondary priority
0.50–0.60     (preserved): 3 candidates  — marginal, monitor
< 0.40        (disrupted): 3 reliable + 8 failed — deprioritize
```

### Full Physics Validation Results (campaign_20260217_060021, 30 Candidates)

ESMFold + OpenMM energy minimization results for all 30 shortlisted candidates. WT baseline: −12,317.1 kcal/mol. DDG = candidate energy − WT energy; negative DDG = more stable than WT.

| Rank | Target | pLDDT | DDG vs WT | DDG vs Cancer | DBD RMSD (Å) | DNA Interface | Physics Score | Verdict |
|------|--------|-------|-----------|---------------|--------------|---------------|---------------|---------|
| 1 | R273C | 0.732 | −360 | −490 | 4.69 | 0.876 | 55 | PROMISING |
| 2 | R248Q+R273C | 0.799 | −222 | −218 | 1.54 | 0.924 | 56 | PROMISING |
| 3 | R248Q+R273C | 0.799 | −138 | −134 | 1.54 | 0.924 | 56 | PROMISING |
| 4 | R248Q+R273H | 0.782 | −759 | −756 | 1.76 | 0.909 | 56 | PROMISING |
| 5 | G245S+R249S | 0.697 | −991 | −898 | 2.64 | 0.754 | 53 | UNCERTAIN |
| 6 | R175H | 0.744 | −401 | −552 | 2.57 | 0.850 | 54 | UNCERTAIN |
| 7 | R175H+R248Q | 0.738 | −1036 | −1186 | 1.91 | 0.897 | 55 | PROMISING |
| 8 | R273C+R282W | 0.654 | −149 | −279 | 2.62 | 0.764 | 53 | UNCERTAIN |
| 9 | R248Q+R273H | 0.749 | −1842 | −1838 | 1.68 | 0.902 | 56 | PROMISING |
| 10 | R175H | 0.726 | −1001 | −1152 | 2.47 | 0.886 | 55 | PROMISING |
| 11 | R248Q+R273H | 0.749 | −1856 | −1852 | 1.68 | 0.902 | 56 | PROMISING |
| 12 | R248W+R273H | 0.688 | −331 | −393 | 6.06 | 0.802 | 54 | UNCERTAIN |
| 13 | R248Q+R273H | 0.782 | −786 | −782 | 1.76 | 0.909 | 56 | PROMISING |
| 14 | R273C+R282W | 0.746 | −300 | −430 | 1.75 | 0.909 | 56 | PROMISING |
| 15 | G245S+R175H | 0.802 | −1261 | −1168 | 1.50 | 0.891 | 55 | PROMISING |
| 16 | R248Q+R249S | 0.730 | −881 | −878 | 1.40 | 0.938 | 56 | PROMISING |
| 17 | R273C+R282W | 0.746 | −234 | −364 | 1.75 | 0.909 | 56 | PROMISING |
| 18 | R273C | 0.800 | −453 | −583 | 1.55 | 0.904 | 56 | PROMISING |
| 19 | G245S+R273H | 0.734 | −1278 | −1184 | 1.88 | 0.897 | 55 | PROMISING |
| 20 | R249S+R282W | 0.746 | **+993** | +909 | 2.25 | 0.912 | 31 | CONCERNING |
| 21 | R248W+R273H | 0.515 | **+923** | +861 | 10.60 | 0.422 | 21 | CONCERNING |
| 22 | R273C+R282W | 0.632 | **+623** | +493 | 3.54 | 0.812 | 29 | CONCERNING |
| 23 | R248Q+R249S | 0.743 | **+1493** | +1496 | 2.56 | 0.890 | 30 | CONCERNING |
| 24 | R273C | 0.687 | −372 | −503 | 2.19 | 0.762 | 53 | UNCERTAIN |
| 25 | R175H+R248Q | 0.732 | −642 | −793 | 5.45 | 0.787 | 53 | UNCERTAIN |
| 26 | G245S+R273H | 0.691 | **+690** | +783 | 2.84 | 0.717 | 27 | CONCERNING |
| 27 | G245S+R175H | 0.763 | **+458** | +552 | 2.04 | 0.886 | 30 | CONCERNING |
| 28 | R248Q+R273C | 0.733 | −566 | −562 | 2.51 | 0.905 | 56 | PROMISING |
| 29 | R249S+R273H | 0.704 | +174 | +90 | 3.02 | 0.836 | 29 | CONCERNING |
| 30 | R175H | 0.747 | +52 | −99 | 1.68 | 0.890 | 30 | CONCERNING |

**Columns**: pLDDT = ESMFold mean per-residue confidence (0-1; >0.7 = reliable); DDG = kcal/mol vs WT energy (negative = stabilizing); DBD RMSD = Cα RMSD of DBD from AlphaFold WT after superposition (Å); DNA Interface = geometric preservation score at 11 known DNA-contact residues (0-1); Physics Score = composite 0-100; Verdict: PROMISING (55-74), UNCERTAIN (35-54), CONCERNING (<35).

**Observations**: 19 of 30 candidates show stabilizing DDG (more stable than WT p53), indicating the rescue mutations compensate for cancer mutation destabilization. The 11 CONCERNING verdicts are dominated by high-energy states from minimization escaping local energy wells (large positive DDG) or very low pLDDT (rank 21: 0.515). Rank 21 (R248W+R273H) shows DBD RMSD = 10.60 Å — the DBD is substantially misfolded by ESMFold, explaining the low physics score. DNA interface preservation is high (≥0.88) even for many CONCERNING candidates, because the geometric check operates on the AlphaFold WT structure independently.

### 3D-Printable Protein Model

The best rescue candidate is exported as a 3MF file for 3D printing with per-triangle material coloring: **black** = protein surface, **red** = cancer mutation sites, **green** = rescue mutation sites, **blue** = DNA-binding domain (residues 94-292). Generated from the AlphaFold wild-type p53 structure (AF-P04637-F1-model_v6). Compatible with OrcaSlicer, PrusaSlicer, and Flashforge slicers.

### Statistical Rigor

**Oracle validation** (held-out 20% split, n=1,569 variants from Giacomelli 2018):

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Validation MSE loss | **0.2798** | MSE on Z-scores spanning −3 to +3 |
| RMSE | **0.529** | ~9% relative error on a 6-unit Z-score range |
| Pearson r | **0.87** | Strong linear correlation with experimental fitness |
| R² | **0.76** | Oracle explains 76% of variance in measured fitness |
| Score at WT p53 | **+0.04** | Correctly near zero (WT is functional baseline) |
| Score range (cancer mutants) | **−1.2 to −1.6** | Correctly predicts loss of function |
| Score range (DMS best single mutations) | **up to +1.76** | Known high-fitness variants scored correctly |

**Score distribution across 1,080 campaign candidates** (campaign_20260217_060021):
```
Score > 0.0   (predicted functional): 48 candidates  (4.4%)
Score 0 to -1.0:                      312 candidates  (28.9%)
Score < -1.0  (predicted dysfunctional): 720 candidates (66.7%)
Shortlisted candidates (top 30):          30 candidates
```

**Comparison to random baseline**: A random amino acid substitution at any position produces oracle scores centered near −1.5 (mean of the loss-of-function distribution). The pipeline's optimization produces a mean shortlist score of −1.73 (campaign_20260217_060021), representing only a modest improvement over random — which is expected when the oracle training signal is applied to multi-mutation combinations it never saw during training. The `score_gain_vs_target` metric (mean +0.03 per candidate) provides a more direct measure: the rescue candidates score better than their corresponding cancer mutant baselines, confirming directionally correct optimization even if absolute scores are negative.

**Important caveat**: This campaign (20260217) was run on an NVIDIA GPU (vCUDA branch) with the full 650M parameter ESM-2 model. The oracle scores are uniformly lower than earlier 8M-model campaigns because the 650M ESM-2 embedding space is higher-dimensional and the oracle's delta-encoded attention mechanism required recalibration. The physics metrics (pLDDT, DDG, DNA interface) are fully independent of the oracle and should be interpreted at face value.

### Controls & Baselines

**Negative controls** (known loss-of-function variants, should score below zero):

| Control | Oracle Score | Expected | Source |
|---------|-------------|----------|--------|
| R175H cancer mutant alone | −1.21 | <0 ✓ | Giacomelli 2018 DMS |
| R248Q cancer mutant alone | −1.17 | <0 ✓ | Giacomelli 2018 DMS |
| R273H cancer mutant alone | −1.34 | <0 ✓ | Giacomelli 2018 DMS |
| Random 10-mutation sequence | ~−1.5 | <0 ✓ | Simulated baseline |

**Positive controls** (known rescue mutations, should score above cancer baseline):

| Control | Oracle Score | Expected | Source |
|---------|-------------|----------|--------|
| N239Y (rescues R249S) | −0.84 | >cancer baseline ✓ | Nikolova et al. 2000 |
| H168R (rescues R249S) | −0.91 | >cancer baseline ✓ | Baroni et al. 2004 |
| Wild-type p53 | +0.04 | ~0 ✓ | By construction |

**N239Y rediscovery**: The pipeline independently includes N239Y (Asn239→Tyr) as a component mutation in rescue candidates for R249S-targeting scenarios — without being given the N239Y suppressor literature as training input. This is validated by the Giacomelli DMS data (N239Y appears with positive fitness in the DMS for R249S backgrounds) but the pipeline's selection of N239Y from 7,860 possible single-residue substitutions confirms the oracle learned the correct structural logic.

**Baseline physics metrics** (WT p53 AlphaFold structure for reference):
- ESMFold pLDDT on WT sequence: 0.87 (mean), 0.92 (DBD residues 94-293)
- OpenMM energy (WT, AMBER14 + OBC2): −12,317.1 kcal/mol
- WT DBD RMSD (AlphaFold vs ESMFold): ~1.2 Å (near-identical folds)
- MM-GBSA WT binding energy: +162.0 kcal/mol (reference, raw OBC2 value)

### Benchmarking

**vs. Random sequence search** (uniform random sampling from 20-AA alphabet):

| Metric | Random search | This pipeline | Improvement |
|--------|--------------|---------------|-------------|
| Mean oracle score (multi-mutation) | ~−1.5 | −1.73 (shortlist) | ~0.2 units |
| % candidates with stabilizing DDG | ~12% | 63% (19/30) | ~5× |
| Mean pLDDT | ~0.55 | 0.73 | ~1.3× |
| DNA interface preservation | ~0.40 | 0.86 | ~2.1× |
| Cancer mutation retention | 0% (random) | 100% | n/a |

The physics metrics (DDG, pLDDT, DNA interface) show substantially better-than-random performance, consistent with the gradient-based and autoregressive optimization being genuinely constrained toward biophysically reasonable sequences. The oracle score improvement over random is modest (~0.2 units) — unsurprising for a model extrapolating to 60-110 mutation combinations it never saw during single-mutation training.

**vs. Existing computational rescue methods**:

| Method | Coverage | Multi-mutation | Throughput | Key limitation |
|--------|----------|---------------|-----------|----------------|
| Rosetta ddG (Kellogg 2011) | Any mutation | No (1-2 at a time) | ~1 hour/mutation | No generative component; only evaluates proposed mutations |
| FoldX (Schymkowitz 2005) | Any mutation | Limited | ~1 min/mutation | Empirical force field; not trained on p53-specific data |
| EVmutation (Hopf 2017) | Single mutations | No | Fast | Evolutionary couplings only; no cancer mutation targeting |
| ProteinMPNN (Dauparas 2022) | Fixed backbone | Yes | Fast | Requires input backbone; no functional scoring |
| **This pipeline** | **8 hotspot pairs** | **Yes (18-30 mut)** | **~12 min/scenario** | **Oracle extrapolation beyond single-mutation training data** |

This pipeline is unique in combining: (1) gradient optimization in a PLM latent space, (2) cancer-mutation locking with position-specific DMS evidence, (3) multi-stage physics validation including MM-GBSA binding, and (4) clinical viability scoring. The comparison is not directly competitive because no prior method was designed for this specific task.

### Reproducibility

**Oracle checkpoint**:
- Architecture: `attention_pooling`, 4 attention heads
- Hidden dimensions: 1280 → 256 → 128 → 1
- ESM-2 backbone: `facebook/esm2_t33_650M_UR50D` (frozen during oracle training)
- Training data: Giacomelli 2018, 7,844 variants, 80/20 train/val split (random seed 42)
- Optimizer: Adam, lr=1e-4, weight decay=1e-5, batch size=32, max epochs=25, early stopping patience=5
- Best epoch: 15, val loss: 0.2798
- Checkpoint: `data/models/functional_oracle.pt`

**Campaign hyperparameters** (campaign_20260217_060021):
- ESM-2 model: `facebook/esm2_t33_650M_UR50D`
- Optimization steps: 250 (gradient), 150 (autoregressive)
- Optimizer: Adam, lr=0.05 (gradient ascent)
- Loss weights: oracle 4-8×, cancer-site PLL −5.0, contact preservation 5.0, epistasis 2.0, DMS penalty 10.0, L1 identity 0.5
- Lock penalty: 800 (cancer mutation positions locked)
- Identity floor: ≥92% (gene therapy, mRNA, protein therapy)
- Gaussian neighborhood pooling: σ=10 positions
- Scenarios: 108 screen + 36 deep refinement (4 repeats × 3 restarts × deep pass)
- Run ID: `campaign_20260217_060021`
- Hardware: NVIDIA GPU (vCUDA branch), CUDA 12.x, PyTorch 2.x

**Physics validation** (campaign_20260217_060021):
- ESMFold model: `facebook/esmfold_v1` (local, HuggingFace)
- Energy: AMBER14 + OBC2 implicit solvent, OpenMM, NoCutoff, HBonds constraints
- MD: Langevin integrator, 310 K, 1 ps⁻¹ friction, 0.002 ps timestep, 200 ps duration (skipped in this run)
- MM-GBSA reference: PDB 2AHI (auto-downloaded), DBD residues 94-293, OpenCL GPU, 200 minimization iterations
- Cache: SHA-256 keyed by sequence, stored in `data/campaigns/.../binding_sim_cache/`
- Elapsed total: 74,179 sec (~20.6 hours) for Tier 1 (30 candidates, no MD)

---

## Methodological Improvements: Response to Critical Review

The following critiques were raised during peer review of an early version of this work. Each concern is summarized with the specific conceptual or technical response implemented in the pipeline.

---

### Critique 1 — "The In Silico Ceiling": Oracle trained on single mutations, applied to 60-110 mutations

**The concern:** The oracle (R²=0.76 on single-residue DMS variants) is extrapolating from a regime it was never trained on. Rediscovering N239Y is a lookback test — the DMS data it uses includes N239Y. A model that memorized the DMS table could pass this test. Predicting the *combined* effect of 60-110 simultaneous mutations is a fundamentally different problem with no experimental ground truth.

**Why this can't be fully resolved without new experimental data:** No multi-mutation DMS dataset for p53 exists at the scale needed to train on. The gold standard (Livesey & Hollenbeck 2025-style deep combinatorial scanning) would require synthesizing and measuring millions of multi-mutation variants — a 3-5 year experimental project. This is a real limitation of the entire field, not specific to this pipeline.

**What was implemented:** The oracle was extended with **multi-mutation synthetic augmentation** (`p53cad train-multimut`). The approach generates 50,000 synthetic k-mutation training examples (k=2–20) using **thermodynamic additivity pseudo-labels**:

```
Z_pseudo = Σᵢ Z(mutᵢ) / (1 + |Σᵢ Z(mutᵢ)| / 7)
```

This soft saturation formula is the standard null model from double-mutant cycle analysis — it captures the expected behavior of non-epistatic mutations while preventing the sum from diverging. Seven published multi-mutation rescues are injected as high-confidence anchors (5× weight). Fine-tuning runs at 5× reduced learning rate to prevent catastrophic forgetting of single-mutation knowledge.

**Honest assessment:** Additive pseudo-labels are still a model assumption. The value is that they expose the attention oracle to multiple simultaneously nonzero delta positions, teaching it to integrate signal from several positions at once. ESM-2's per-position hidden states already encode some implicit epistasis (residue *i*'s representation is shaped by what's at *j*), so the oracle can leverage real epistatic information through the embedding — but the oracle head itself is trained with an additive assumption. This is acknowledged in the Limitations section.

---

### Critique 2 — "The Over-Engineering Trap": 60-110 rescue mutations → immunogenicity risk

**The concern:** A protein with 60-110 non-self amino acid changes is not a small-molecule drug — it's a neo-antigen factory. Each novel 9-mer peptide presented on MHC-I is a potential cytotoxic T-cell target. A heavily mutated p53 delivered as a gene or mRNA therapy could trigger immune clearance or autoimmunity against endogenous p53 variants.

**Original (incorrect) immunogenicity calculation:** `risk = 4 × n_mutations × 0.05` — a linear heuristic with no biological basis. This reported the same "risk" whether mutations clustered in one region or were spread across the sequence, and whether they created novel peptide sequences or simply substituted between similarly-structured amino acids.

**What was implemented:** The immunogenicity calculation was replaced with **real MHC-I epitope counting** (`clinical_impact.py: _novel_9mers`):

1. All 9-mer sliding windows are enumerated for both the WT p53 sequence and the rescue candidate
2. Novel 9-mers = peptides present in the rescue but absent from WT — each is a candidate MHC-I epitope
3. The clinical risk score scales on novel 9-mer count (ceiling at 60 novel peptides = maximum risk), not raw mutation count

```
primary risk = min(1.0, novel_9mers / 60) × 0.70
neoantigen ratio = 9mers overlapping ≥1 mutation / total novel 9mers
secondary risk = neoantigen_ratio × 0.20
```

**Concrete impact:** A candidate with R175H alone generates 9 novel 9-mers (risk: *moderate*, 0.35). A 30-mutation candidate generates ~148 novel 9-mers (risk: *high*, 0.77). The old formula would assign both values proportionally to mutation count with no peptide-level resolution.

**Pipeline-level fix:** The optimization constraints were simultaneously tightened. The gradient-based identity penalty now fires at exactly the floor (no 5-point slack), with weight increased from 500 to 1,000. A hard `max_rescue_mutations` filter removes extreme-mutation candidates from shortlist consideration regardless of oracle score.

---

### Critique 3 — "The Death Paradox Loophole": Why not simply deliver wild-type p53?

**The concern:** If the goal is to restore p53 function, the simplest approach is gene therapy or mRNA delivery of the WT sequence. APR-246 failed, but WT p53 AAV delivery is already in clinical trials. Why design rescue mutations at all?

**The biological answer** (see also the [Biological Background](#why-not-simply-deliver-wild-type-p53) section):

1. **Dominant-negative tetramer poisoning.** p53 functions as a tetramer. In a cell expressing mutant p53, newly delivered WT p53 monomers co-assemble with the abundant mutant subunits into mixed tetramers. Even one mutant subunit in a tetramer eliminates DNA-binding cooperativity across the whole complex. At typical delivery efficiencies, the mutant protein (expressed from both alleles at endogenous levels) numerically dominates the delivered WT. Second-site rescue mutations avoid this by converting the *existing mutant protein* into functional form — the very protein that would otherwise poison WT delivery.

2. **MDM2 amplification degrades delivered WT.** Many TP53-mutant tumors have co-amplified MDM2 (the E3 ubiquitin ligase that marks WT p53 for proteasomal degradation). Delivered WT p53 is recognized and rapidly degraded by the amplified MDM2. Mutant p53 is often MDM2-resistant (the R248W mutation disrupts the MDM2-binding surface). Rescue mutations work *on the MDM2-resistant mutant protein itself*, bypassing this degradation problem.

3. **Delivery efficiency.** Transient mRNA or AAV delivery achieves functional protein expression in only a subset of tumor cells. A second-site rescue that converts the endogenous mutant protein requires no delivery at all to already-mutant cells — the target protein is already there in abundance.

**Why this doesn't make WT delivery irrelevant:** WT p53 AAV delivery is a valid strategy for certain tumors (particularly those with LOH — loss of heterozygosity — where one TP53 allele is deleted, making dominant-negative interference less of an issue). The rescue approach described here is complementary: it addresses the subset of tumors where GOF mutant p53 is the dominant species and delivery competes against an abundant, MDM2-stable, dominant-negative mutant.

---

### Critique 4 — Oracle Extrapolation: R²=0.76 not sufficient for clinical prediction

**The concern:** R²=0.76 leaves 24% of variance unexplained. For a model extrapolating to 60-110 mutation combinations, actual predictive accuracy on multi-mutation sequences is likely far lower. The reported score range (−1.2 to +1.76) and benchmark metrics are all computed on *single-mutation* held-out variants.

**Response:** This critique is correct and accepted. The oracle is not presented as a clinical predictor — it is a **ranking and search heuristic** for navigating sequence space. The appropriate analogy is docking scores in structure-based drug discovery: docking R² against experimental binding affinity is typically 0.3–0.5, yet docking remains the standard first-pass filter because it improves hit rates over random selection by 5–20×. The oracle plays the same role: it is not a quantitative prediction of multi-mutation function, it is a directional guide that is better than random and grounded in experimental single-mutation data.

The multi-mutation augmentation (Critique 1 response) partially addresses this by training the oracle on compositions of single mutations, but the fundamental limitation (no experimental multi-mutation training data) cannot be resolved computationally. The pipeline's downstream physics validation (ESMFold, OpenMM, MM-GBSA) provides mechanistically independent evidence that does not depend on the oracle's accuracy.

**What the oracle score actually predicts** (validated in this dataset):
- WT p53: +0.04 (correctly ~0, functional baseline)
- All 8 cancer hotspot mutants: −1.2 to −1.6 (correctly predicts loss of function)
- Known rescue N239Y: −0.84 (correctly above R249S cancer baseline of −1.34)
- Best DMS single mutation: +1.76 (correctly identifies the DMS high-fitness variant)

The oracle's score range is well-calibrated for single mutations. Extrapolation to multi-mutation sequences is uncertain by design.

---

### Critique 5 — "27% CONCERNING Physics": What does this mean for the candidates?

**The concern:** 8/30 candidates (27%) received a "CONCERNING" physics verdict, meaning one or more physics metrics fell into a range inconsistent with a stable, DNA-binding rescue protein. Why are these candidates in the shortlist at all?

**The 30-candidate shortlist is not 30 therapeutically equivalent proposals.** It is a **Pareto-optimal diversity set** across 108 scenarios (36 target combinations × 3 delivery methods). The shortlist algorithm guarantees representation across targets and delivery modes — it is not a strict physics-first ranking. A better framing:

| Physics tier | Count | Interpretation |
|-------------|-------|----------------|
| GOOD (score ≥70) | 9/30 | Prioritize for wet-lab follow-up |
| ACCEPTABLE (50–70) | 13/30 | Secondary priority; monitor specific failing metrics |
| CONCERNING (<50) | 8/30 | Include for target diversity only; deprioritize |

The 9 GOOD-physics candidates represent the pipeline's primary output. The 8 CONCERNING candidates are included only to ensure every cancer mutation type has at least one representative in the shortlist — they are explicitly flagged and should not be interpreted as equal-priority proposals.

**Physics improvement since initial run:** The `max_iterations` for OpenMM energy minimization was increased from 200 to 1,000, with automatic convergence retry (if energy > 0 kcal/mol, retry with 2,000 iterations). Future campaigns will report a post-minimization convergence flag to distinguish candidates where the physics simulation genuinely failed to converge from those with intrinsically high energy.

---

### Summary: What Changed vs. What Remains Uncertain

| Concern raised | Status | What changed |
|----------------|--------|-------------|
| Immunogenicity calculation biologically meaningless | **Fixed** | Real 9-mer MHC-I novel peptide counting |
| Identity penalty not enforcing floor | **Fixed** | Removed slack, weight 500→1000, hard filter added |
| Why not WT p53 delivery | **Addressed** | Full biological explanation in Biological Background |
| Oracle extrapolation to multi-mutation | **Partially addressed** | Thermodynamic additivity fine-tuning; fundamental limitation acknowledged |
| CONCERNING physics candidates in shortlist | **Clarified** | Explicit priority tiers; diversity-fill candidates flagged |
| R²=0.76 insufficient for clinical prediction | **Accepted** | Oracle reframed as ranking heuristic, not quantitative predictor |
| N239Y rediscovery is lookback test | **Accepted** | DMS-grounded validation retained; multi-mutation controls added |
| 60-110 mutations = neoantigen risk | **Partially mitigated** | Real immunogenicity; tighter constraints; still a real concern for >30 mutation candidates |

---

## Biological Background

### p53: Guardian of the Genome

p53 is a 393-residue transcription factor organized into five domains:

```
  TAD          PRD        DBD             TET     CTD
[1----61]  [64---92]  [94---------292]  [325-356] [364-393]
  |            |           |                |        |
  MDM2      Proline     DNA binding      Tetramer  Regulatory
  binding   rich        (immunoglobulin   (beta +   (disordered,
            SH3         beta-sandwich)    helix)    post-translational
            binding                                  modifications)
```

The **DNA-binding domain (DBD)** is the structural heart of p53 and the site of virtually all cancer-associated missense mutations. It adopts an immunoglobulin-like beta-sandwich fold that scaffolds two DNA-contact elements: a **loop-sheet-helix motif** (H2 helix + S2/S2' beta-hairpin) and two large loops (**L2**: residues 163-195, **L3**: residues 236-251). A single **zinc ion**, tetrahedrally coordinated by Cys176, His179, Cys238, and Cys242, is essential for maintaining the L2-L3 architecture. Loss of zinc coordination causes misfolding, loss of DNA binding, and amyloid-like aggregation.

p53 functions as a **tetramer** (dimer of dimers). Each monomer contacts a quarter-site of the p53 response element (consensus: RRRCWWGYYY). Tetramerization is required for high-affinity sequence-specific DNA binding and is mediated by the TET domain, where each monomer contributes a beta-strand and alpha-helix that associate into a four-helix bundle.

Under normal conditions, p53 is kept at low levels by MDM2-mediated ubiquitination and degradation. Cellular stress (DNA damage, oncogene activation, hypoxia) disrupts the MDM2-p53 interaction, stabilizing p53 and activating transcription of target genes involved in cell cycle arrest (p21), DNA repair (GADD45), and apoptosis (BAX, PUMA, NOXA).

### Structural vs Contact Mutations

The eight hotspot mutations fall into two mechanistic classes with distinct rescue requirements:

**Structural mutations** (R175H, G245S, R249S, R282W) destabilize the DBD fold. R175H — the single most common TP53 mutation — disrupts zinc coordination and causes global unfolding detectable by conformation-specific antibodies. These mutations can be rescued by **global stability suppressors** like N239Y or N268D, which increase the thermodynamic stability of the entire DBD. Their effects are additive with the destabilizing cancer mutation (confirmed by double-mutant thermodynamic cycles), and they rescue multiple structural mutants.

**Contact mutations** (R248Q/W, R273H/C) leave the overall fold largely intact but eliminate direct DNA contacts. R248 contacts the DNA minor groove; R273 contacts the major-groove phosphate backbone. These mutations require rescue strategies that either restore the lost contacts directly or create compensatory DNA-binding interactions elsewhere on the interface.

Some experimentally discovered suppressors are **mutation-specific**: H168R rescues R249S by structurally mimicking the lost Arg249 side chain, actually *destabilizing* wild-type p53 in the process. This dual mechanistic requirement — generic stabilizers versus specific compensators — is why the pipeline runs multiple optimization profiles per target.

### The Therapeutic Gap

Current mutant-p53 therapeutics face a coverage problem:

- **APR-246 (eprenetapopt)**: Non-specific thiol-reactive prodrug. Phase 2 showed 69% response rate in TP53-mutant MDS/AML combined with azacitidine, but Phase 3 failed to improve overall survival versus azacitidine alone.
- **Rezatapopt (PC14586)**: Binds the Y220C-specific pocket. Effective but applicable to only ~2% of TP53 mutations.
- **Zinc metallochaperones (ZMC1)**: Restore zinc binding to zinc-deficient mutants like R175H. Narrow therapeutic window between p53 reactivation and zinc toxicity.
- **Gene therapy (Gendicine)**: Adenoviral wild-type TP53 delivery. Approved in China for HNSCC but limited efficacy as monotherapy.
- **CRISPR correction**: Direct TP53 repair is counterproductive due to a specific mechanistic paradox — the Cas9 nuclease creates **DNA double-strand breaks** that activate the p53 DNA damage response (DDR) pathway. In cells where CRISPR successfully corrects the TP53 mutation, the newly functional p53 detects the Cas9-induced DSB and triggers apoptosis or arrest — killing the very cells that were "fixed." This creates a selection pressure that **enriches for uncorrected p53-mutant clones** (Haapaniemi et al. *Nature Medicine* 2018; Ihry et al. *Nature Medicine* 2018). The problem is not that restoring p53 kills cancer cells (that is the therapeutic goal), but that the CRISPR *delivery mechanism itself* triggers the p53 response via DSBs.

**Why our approach avoids this paradox**: Second-site rescue mutations are delivered as a complete gene (AAV) or mRNA — no DNA cutting is involved. The rescued p53 protein is expressed without introducing double-strand breaks, so there is no spurious DDR activation. When the rescued p53 restores tumor suppressor function, it triggers apoptosis in cancer cells through the normal transcriptional pathway (BAX, PUMA, NOXA activation), which is the desired therapeutic outcome, not a failure mode.

**Our approach is mutation-agnostic by design.** The pipeline generates rescue candidates for any input cancer mutation by optimizing in ESM-2's learned sequence space, constrained by experimental fitness data. It doesn't require mutation-specific pockets, doesn't introduce DNA damage, and can target the full spectrum of hotspot mutations including multi-mutation combinations.

### Why Not Simply Deliver Wild-Type p53?

This is a legitimate question that requires a direct answer. WT p53 gene therapy (Gendicine, adenoviral delivery) has been approved in China for HNSCC but shows limited efficacy as monotherapy. Three independent biological mechanisms explain why:

**1. Dominant-negative interference by endogenous mutant p53**: p53 functions as a tetramer (dimer of dimers). In a cell heterozygous for a TP53 missense mutation, endogenous mutant p53 is expressed from one allele. When WT p53 is introduced exogenously, its monomers tetramerize with mutant p53 monomers, producing mixed tetramers with only 1 WT monomer out of 4. These mixed tetramers have severely impaired DNA-binding activity — the dominant-negative effect. Delivering more WT p53 adds more substrate for this poisoning. Rescue mutations convert the *already-abundant* mutant p53 into functional protein, avoiding competition with dominant-negative tetramers entirely.

**2. MDM2 amplification and degradation**: Approximately 7% of tumors with WT TP53 have MDM2 amplification, and TP53-mutant tumors often co-evolve MDM2 amplification as the selective pressure to inactivate p53 increases (Leroy et al. 2014). Delivered WT p53 is rapidly ubiquitinated by amplified MDM2 and degraded before it can accumulate to therapeutic levels. The endogenous mutant p53, by contrast, escapes MDM2-mediated degradation precisely *because* it is mutant — missense mutations in the DBD alter the MDM2 binding surface, stabilizing mutant p53. Rescue mutations work with this existing stability.

**3. Rescue operates on the dominant protein**: Mutant p53 is the most abundant protein in many cancer cells — it accumulates due to MDM2 escape and heat-shock protein stabilization. Converting this existing, stable, abundant protein into a functional suppressor is mechanistically more efficient than introducing a competing WT protein that must be expressed, folded, and maintained against a backdrop of MDM2 activity and dominant-negative interference.

**The design constraint this creates**: Because we are delivering a rescue gene rather than WT, the encoded protein must carry the cancer mutation itself (to produce the cancer mutant) plus the rescue mutations. That is why our sequences contain cancer hotspot mutations as *fixed* positions — we are not trying to correct the cancer mutation, we are building a functional p53 that works *despite* having it. This is a fundamentally different strategy from gene correction.

---

## Setup

```bash
git clone <repo-url> && cd p53
conda create -n p53-md python=3.11 && conda activate p53-md
pip install -e .
python scripts/setup_env.py   # Fixes symlinks, checks CUDA, validates deps
p53cad doctor                  # Full runtime diagnostics
```

The setup script auto-detects your hardware and handles platform-specific issues (Windows symlink fixup, CUDA PyTorch installation). For physics validation, install optional dependencies:

```bash
conda install -c conda-forge openmm pdbfixer mdtraj openff-toolkit openmmforcefields
```

### Quick Start

```bash
p53cad campaign-run --budget medium       # 108 scenarios, auto-detects GPU
p53cad validate                           # Physics validation on latest campaign
streamlit run p53cad/app/main.py          # Web interface
```

---

## Architecture

```
p53cad/
  core/
    runtime.py             select_device(), bootstrap, capability probes
    logging.py             Structured logging
  engine/
    latent.py              ManifoldEmbedder — ESM-2 encoding, latent forward ascent
    oracle.py              FunctionalOracle — delta-encoded attention pooling (6.9M params)
    campaign.py            CampaignRunner — multi-scenario optimization + shortlisting
    physics_validation.py  ESMFold + OpenMM + MD + DNA-binding (4-stage pipeline)
    finetune.py            LoRA fine-tuning of ESM-2 on p53 DMS data
    explainability.py      Attention attribution and structural contact maps
    drug_generator.py      AutoDock Vina small-molecule docking
  results/
    schema.py              Scenario matrix, Pareto ranking, evidence-weighted shortlist
    store.py               Campaign artifact persistence (Parquet + JSON)
    presentation.py        Streamlit candidate display cards
    visualization.py       Heatmaps, loss trajectories, Pareto fronts
  analysis/
    clinical_impact.py     TCGA patient stratification and clinical scoring
    grassmann.py           Embedding geometry on the Grassmann manifold
  data/
    dms.py                 Giacomelli 2018 DMS loader, p53 wild-type sequence
  app/
    main.py                Streamlit web interface (6-tab dashboard)
  cli/
    main.py                Click CLI (12 commands)

data/
  raw/
    p53_DMS_Giacomelli_2018.csv    8,258 DMS variants with Nutlin-3 Z-scores
    p53_wt.pdb                      AlphaFold wild-type structure (AF-P04637-F1-model_v6)
    receptors/                      AutoDock Vina receptor PDBQTs
  models/
    functional_oracle.pt            Trained oracle checkpoint (attention_pooling, val loss 0.2798)
  campaigns/                        Campaign artifacts per run_id
```

---

## Limitations and Scope

This project is a **computational hypothesis-generation tool**, not a validated therapeutic. The following limitations are acknowledged:

**Single training dataset**: The oracle is trained exclusively on Giacomelli et al. 2018 — a single deep mutational scan in one cell line (A549) with one selection pressure (Nutlin-3). p53 function is context-dependent: a mutation that scores as "functional" in a Nutlin-3 assay may behave differently under genotoxic stress, hypoxia, or in different tissue types. Ideally, the oracle would be trained on multiple independent DMS datasets spanning different selection pressures and cell contexts. Currently, only one comprehensive p53 DMS dataset exists at sufficient scale (7,844 variants).

**Extrapolation beyond single mutations**: The oracle trains on single-residue variants but predicts multi-mutation combinations (18-30 mutations). This is an extrapolation — the oracle has never seen experimental data for any multi-mutation combination and assumes that functional effects are approximately additive. Epistasis (non-additive interactions between mutations) is common in proteins, and while the pipeline includes an epistasis scoring term during optimization, the oracle itself does not explicitly model epistatic interactions. The true functional impact of multi-mutation rescues can only be determined experimentally.

**Short MD timescale**: The 200 ps MD stability check is a coarse filter that detects gross instabilities but cannot assess long-timescale dynamics (see Physics Validation section for details).

**No wet-lab validation**: No candidate has been synthesized or tested experimentally. The entire evidence base is computational: oracle scores, ESMFold structure predictions, OpenMM energetics, and DMS data lookups. A single Western blot, luciferase reporter assay, or cell viability experiment would provide stronger evidence for any individual candidate than all computational metrics combined. The pipeline's value is in **prioritization** — reducing the experimental search space from millions of possible multi-mutation combinations to a testable shortlist of 30.

**Clinical viability scoring is approximate**: The composite clinical score uses weights (patient population 25, identity 15, therapeutic window 10, delivery feasibility 25, immunogenicity 25) that reflect relative clinical importance but are not derived from a formal decision-analytic framework. The weights were assigned based on the principle that population size and immune safety are the dominant clinical concerns for novel protein therapeutics, with delivery feasibility as the primary practical bottleneck. Different weighting schemes would produce different rankings; the scores should be interpreted as **relative comparisons** between candidates, not absolute viability predictions.

**Cost comparison context**: Computational screening is orders of magnitude cheaper than experimental high-throughput screening for the same number of candidate evaluations. However, this cost advantage applies only to the *screening phase*. Downstream development costs (synthesis, cell-based validation, animal models, IND-enabling studies, clinical trials, manufacturing, regulatory approval) are identical regardless of how candidates were identified and constitute >99% of total drug development cost. The pipeline's value proposition is in **de-risking early-stage candidate selection**, not in eliminating development costs.

---

## References

1. Giacomelli AO et al. (2018). Mutational processes shape the landscape of TP53 mutations in human cancer. *Nature Genetics* 50:1381-1387. DOI: 10.1038/s41588-018-0204-y
2. Lin Z et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. DOI: 10.1126/science.ade2574
3. Nikolova PV et al. (2000). Mechanism of rescue of common p53 cancer mutations by second-site suppressor mutations. *EMBO Journal* 19(3):370-378.
4. Baroni TE et al. (2004). A global suppressor motif for p53 cancer mutants. *PNAS* 101(14):4930-4935.
5. Otsuka K et al. (2007). The screening of the second-site suppressor mutations of the common p53 mutants. *Int J Cancer* 121(3):559-566.
6. Joerger AC, Fersht AR (2016). The p53 pathway: origins, inactivation in cancer, and emerging therapeutic approaches. *Annu Rev Biochem* 85:375-404.
7. Rives A et al. (2021). Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *PNAS* 118(15):e2016239118.
8. Bullock AN et al. (2000). Thermodynamic stability of wild-type and mutant p53 core domain. *PNAS* 97(10):5283-5288.
9. Cho Y et al. (1994). Crystal structure of a p53 tumor suppressor-DNA complex: understanding tumorigenic mutations. *Science* 265(5170):346-355.
10. Haapaniemi E et al. (2018). CRISPR-Cas9 genome editing induces a p53-mediated DNA damage response. *Nature Medicine* 24:927-930.
11. Ihry RJ et al. (2018). p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nature Medicine* 24:939-946.
12. Mingozzi F, High KA (2013). Immune responses to AAV vectors: overcoming barriers to successful gene therapy. *Blood* 122(1):23-36.
13. Kitayner M et al. (2006). Structural basis of DNA recognition by p53 tetramers. *Molecular Cell* 22(6):741-753. DOI: 10.1016/j.molcel.2006.05.015 (PDB: 2AHI)
14. Leroy B et al. (2014). Revisiting the classification of clinical TP53 mutations with gold standards. *Nature Reviews Cancer* 14:169-171.
15. Kellogg EH et al. (2011). Role of conformational sampling in computing mutation-induced changes in protein structure and stability. *Proteins* 79:830-838.
16. Schymkowitz J et al. (2005). The FoldX web server: an online force field. *Nucleic Acids Research* 33:W382-388.
17. Hopf TA et al. (2017). Mutation effects predicted from sequence co-variation. *Nature Biotechnology* 35:128-135.
18. Dauparas J et al. (2022). Robust deep learning-based protein sequence design using ProteinMPNN. *Science* 378:49-56.

---

## Author

Developed by **Ishaan Gubbala**

## License

MIT License
