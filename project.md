# ProteinForge-p53 (PMG-p53)
## Constraint-aware generative design of stability-rescuing suppressor mutations for cancer-destabilized p53

**Project name (recommended):** ProteinForge-p53  
**Acronym:** PMG-p53 (Mutation Protein Generator for p53)  
**Category:** Computational Biology and Bioinformatics (CBIO)  
**Project type:** Fully computational pipeline. No wet lab. No patient recruitment.

---

## 1) One-paragraph overview

TP53 is one of the most frequently altered genes in cancer, and many TP53 missense mutations destabilize the p53 DNA-binding core domain, increasing misfolding and reducing functional protein. ProteinForge-p53 reframes “p53 rescue” as a constrained generative design problem: given a destabilizing p53 mutant, generate and rank small second-site suppressor mutation sets (1–3 additional edits) predicted to improve folding stability while avoiding key functional and structural residues, including the zinc coordination site (Cys176, His179, Cys238, Cys242) and DNA-contact regions. The pipeline integrates public sequence and annotation data (UniProt), a structure model (AlphaFold or PDB), ClinVar-labeled variants for validation, and stability scoring via ΔΔG estimation (FoldX as the default engine, optionally with a second scorer for consensus). Outputs include a variant stability benchmark, Pareto-optimal rescue designs for cancer hotspot mutants, ablation studies showing why constraints are necessary, and poster-ready figures and tables.

---

## 2) What you ship (end product)

### 2.1 Primary deliverables
1. **A reproducible software pipeline** (CLI + configs) that:
   - downloads and caches public inputs (TP53 sequence, structure, ClinVar variants, hotspot lists)
   - builds clean label sets for validation (benign vs pathogenic where available)
   - generates rescue candidates under constraints
   - scores stability and risk
   - runs multi-objective selection (Pareto front)
   - exports figures and tables to a `reports/` directory

2. **Rescue libraries** for selected p53 hotspot mutants:
   - CSV tables ranked by predicted stability rescue and filtered by risk
   - Pareto plots for stability–risk tradeoff per mutant
   - 3D structure visualizations highlighting protected regions and proposed rescue sites

3. **Validation report**:
   - benign vs pathogenic separation on ClinVar-labeled TP53 missense variants using ΔΔG as a stability signal
   - effect sizes, AUROC/AUPRC, and confidence intervals via bootstrapping
   - robustness checks across scoring replicates and parameter stress
   - ablation results: what happens when you remove constraints

### 2.2 Optional deliverables (nice for science fairs)
- A simple web viewer or notebook that:
  - lets you click a target mutant
  - shows the top rescue sets
  - shows the Pareto front and structure highlights

---

## 3) Research question and hypothesis

### 3.1 Research question
Can a constraint-aware generative design system propose plausible second-site suppressor mutation sets that improve predicted stability of common cancer p53 mutants, while staying away from functional and highly conserved regions?

### 3.2 Hypotheses
- **H1 (validation):** ClinVar-labeled pathogenic TP53 missense variants will, on average, be more destabilizing (higher ΔΔG) than benign variants.
- **H2 (design):** For multiple destabilizing hotspot mutants, ProteinForge-p53 will produce non-trivial Pareto fronts with candidates that improve predicted stability without violating constraints.
- **H3 (constraints):** Removing biological constraints will produce designs with “better” stability scores that are substantially less plausible (more conserved sites, closer to functional residues, more surface-risk), demonstrating that constraints are essential.

### 3.3 Claim boundaries (say this clearly on your board)
- This project does **not** claim a therapy, clinical benefit, or patient outcome.
- This project does **not** claim restored p53 transcriptional function in cells.
- This project outputs **computationally ranked candidate designs** and evaluates them with transparent, reproducible metrics.

---

## 4) Scientific motivation

### 4.1 Why p53 and TP53 mutations matter
- p53 is a key tumor suppressor controlling stress response pathways.
- TP53 is among the most frequently mutated genes in many cancers.
- Many mutations occur in the DNA-binding core domain and include recurrent hotspots.

### 4.2 Why stability rescue is a meaningful computational objective
- A large subset of p53 mutants reduce folding stability, shifting the ensemble toward misfolded states.
- Second-site suppressor mutations can sometimes restore stability by improving core packing or local structural support, while leaving functional surfaces untouched.
- This is a natural fit for constrained optimization: maximize stability gain while minimizing functional risk.

---

## 5) Data plan and how to obtain it

Everything below is public, downloadable, and suitable for a fully computational project.

### 5.1 TP53 sequence and annotations
- **Source:** UniProt for TP53 (canonical protein)
- **What you use:**
  - sequence for reference checks
  - domain boundaries and functional annotations when available

**Files created**
- `data/raw/uniprot/TP53.*`
- `data/interim/tp53.fasta`
- `data/interim/uniprot_features.json`

### 5.2 Structure model
Pick one consistent structural starting point for scoring:
- **AlphaFold TP53 model** (convenient and consistent), or
- **PDB core domain structures** (if you want experimental structures)

**Files created**
- `data/raw/structures/TP53_model.pdb`
- `data/interim/structures/tp53_core_clean.pdb`

### 5.3 Variant labels for validation
- **Source:** ClinVar bulk downloads (XML)  
- **What you use:**
  - TP53 missense variants
  - clinical significance labels (benign/likely benign vs pathogenic/likely pathogenic)
  - review status as a quality filter (optional)

**Files created**
- `data/raw/clinvar/ClinVarVariationRelease.xml.gz`
- `data/interim/tp53_clinvar_missense.parquet`
- `data/processed/label_sets/benign.parquet`, `pathogenic.parquet`

### 5.4 Hotspot selection for case studies
- **Sources:** cancer mutation resources (for example, NCI TP53 database tables or curated hotspot lists)
- **What you use:**
  - shortlist 3–7 target mutants for deep design case studies

**Files created**
- `data/interim/hotspots_tp53.csv`

---

## 6) Methods and algorithms

### 6.1 Pipeline overview
1. **Ingest and normalize variants** (ClinVar)
2. **Define the p53 design target region** (core domain)
3. **Build constraints** (protected residues, conservation, geometry filters)
4. **Generate candidates** (1–3 edits)
5. **Score stability and risk** (ΔΔG + proxies)
6. **Select Pareto-optimal candidates**
7. **Validate and report** (variant separation, ablations, robustness)

### 6.2 Domain definition
- Use a consistent “design domain” such as the p53 DNA-binding core domain.
- Restrict scoring and candidate generation to this region to keep interpretation clean.

### 6.3 Constraint system (biological safety rules)
**Hard constraints (reject candidates):**
- Never mutate zinc-binding residues: **Cys176, His179, Cys238, Cys242**
- Never mutate residues annotated as essential for DNA binding or other core function
- Reject any candidate edit within a minimum distance (example: 6 Å) of protected residues

**Soft constraints (penalties in risk score):**
- conservation penalty from multiple sequence alignment (MSA)
- surface hydrophobic patch penalty (aggregation proxy)
- exposure penalty (avoid introducing hydrophobics on exposed surface)

Why this matters:
- Without constraints, the optimizer tends to “cheat” by editing functional residues that can improve computed stability but would likely break biology.

### 6.4 Candidate generation (generative design)
ProteinForge-p53 generates candidates in a bounded way to avoid combinatorial explosion.

**Step A: select designable positions**
- exclude protected residues
- exclude highly conserved positions above a threshold
- prefer buried or partially buried residues (good for packing stabilization)

**Step B: propose substitutions**
- use a controlled substitution palette (conservative swaps) for buried sites
- optionally use protein language model suggestions to prioritize “natural” substitutions

**Step C: grow mutation sets**
- generate single edits
- expand to pairs and triples using beam search:
  - keep only the top K candidates by preliminary score at each depth

### 6.5 Stability scoring (ΔΔG)
Default: FoldX BuildModel to estimate ΔΔG for each candidate mutation set.

Key implementation details:
- normalize scoring outputs consistently
- run multiple replicates per candidate if needed
- cache by `(structure_hash, mutation_set)` to avoid recomputation
- parallelize scoring across CPU cores

### 6.6 Functional risk scoring (proxy)
Risk is a composite score that keeps the design biologically plausible:
- distance-to-protected residues
- conservation penalty
- exposure and hydrophobic patch penalty

### 6.7 Multi-objective selection (Pareto)
Each candidate has at least two competing objectives:
- maximize stability rescue (more negative ΔΔG gain is better)
- minimize risk score (closer to 0 is safer)

Pareto selection yields a set of “best tradeoff” designs rather than one arbitrary winner.

---

## 7) Implementation plan and codebase structure

### 7.1 Repository layout
```
proteinforge-p53/
  configs/
  data/
    raw/
    interim/
    processed/
  src/
    cli.py
    core/
    data_fetch/
    variants/
    target/
    structure/
    design/
    scoring/
    optimize/
    eval/
    viz/
  experiments/
  reports/
    figures/
    tables/
    logs/
```

### 7.2 What each module does (high level)
- `data_fetch/`: downloads and caches UniProt, structures, ClinVar, hotspot lists
- `variants/`: parses ClinVar, normalizes HGVS, maps to protein positions, builds label sets
- `target/`: defines core domain, protected residues, conservation map
- `structure/`: cleans structure, computes burial and residue distances
- `design/`: generates mutation candidates, applies filters
- `scoring/`: runs ΔΔG scoring and computes risk proxies
- `optimize/`: beam search or evolutionary search; Pareto selection
- `eval/`: benchmark on labeled variants; rescue quality metrics; ablations
- `viz/`: poster-ready plots and structure highlighting outputs

### 7.3 Performance and reproducibility built-in from day 1
- caching of all expensive computations (ΔΔG runs, distance maps, conservation maps)
- parallel scoring with joblib or multiprocessing
- Parquet tables for large intermediate datasets
- deterministic seeds for candidate sampling and search steps
- run metadata: tool versions, config hashes, timestamps in `reports/logs/`

---

## 8) Metrics and evaluation

### 8.1 Variant stability benchmark
Purpose: demonstrate the ΔΔG signal aligns with real human variant labels.

**Outputs**
- distribution plots of ΔΔG for benign vs pathogenic
- AUROC and AUPRC
- effect sizes and confidence intervals

**Interpretation**
- strong separation supports stability as a major contributor for a subset of pathogenic variants
- overlap is expected because many variants act via other mechanisms

### 8.2 Rescue design evaluation
For each target mutant:
- number of feasible candidates generated
- number of candidates rejected by constraints (shows constraints matter)
- Pareto front shape and size
- top candidates with score breakdown and structure highlights

### 8.3 Ablation studies (critical for judge trust)
Run the same design pipeline with:
- no conservation penalty
- no distance-to-protected filter
- no surface-risk penalty

Expected: stability scores improve but risk skyrockets, proving constraints prevent biologically implausible solutions.

### 8.4 Robustness
- replicate scoring runs for a subset of candidates
- rank correlation across replicates
- parameter stress testing (beam width, penalty weights)

---

## 9) Current status and preliminary output interpretation

### 9.1 What exists in the current report outputs
- Rescue tables:
  - `reports/tables/rescues_R175H.csv`
  - `reports/tables/rescues_R248Q.csv`
  - `reports/tables/rescues_R273H.csv`
- Pareto plots:
  - `reports/figures/pareto_R175H.png`
  - `reports/figures/pareto_R248Q.png`
  - `reports/figures/pareto_R273H.png`
- Variant separation benchmark is currently missing because the label sets were empty and the benchmark step was skipped.

### 9.2 Summary of top rescues (from the existing CSV outputs)
- **R175H** best stability gain candidate:
  - `A189S, M133L, S95A` with `ddg_gain ≈ -17.04`, `risk ≈ 0.033`
- **R248Q** best stability gain candidate:
  - `M133L, R196Q, S95A` with `ddg_gain ≈ -22.30`, `risk ≈ 0.0667`
- **R273H** best stability gain candidates include:
  - `C229A, R196Q, S95A` with `ddg_gain ≈ -21.42`, `risk ≈ 0.10`
  - `R196Q, S215A, S95A` with `ddg_gain ≈ -20.29`, `risk ≈ 0.0667`

### 9.3 Patterns across targets (observed)
- Recurring edits appear across top-ranked designs:
  - `S95A`, `M133L`, and `R196Q` show up frequently.
  - `C229A` or `C229S` appears prominently for the R273H target.
- All top-20 rows are triple mutation rescues under current settings. No singles or doubles appear in the top-20.

### 9.4 How to interpret these patterns
- The recurrence suggests the pipeline may be discovering “global stabilizers,” edits that improve stability broadly across multiple mutant contexts.
- It can also indicate search-space bias or insufficient complexity penalties.
- The dominance of triple-mutation sets indicates the objective currently rewards cumulative stability improvements more than it penalizes design complexity.

### 9.5 Immediate quality checks to run (high priority)
If `ddg_gain` magnitudes appear unusually large, run these checks:
1. Score several single mutations alone (for example, `S95A`, `M133L`, `R196Q`) and confirm their individual ΔΔG effects are reasonable.
2. Confirm FoldX output parsing extracts the correct ΔΔG field and averages replicates rather than summing.
3. Confirm `ddg_gain` is computed against the intended baseline consistently.

### 9.6 Fixes to improve scientific realism and presentation strength
- Add an explicit penalty for mutation count:
  - treat “number of edits” as a third objective or add a regularization term to the score
- Enforce that single and double mutation candidates remain in the Pareto set by designing a 3D Pareto selection:
  - stability rescue
  - risk
  - mutation count
- Restore and validate the ClinVar label benchmark by fixing why label sets are empty:
  - common causes include overly strict filters, transcript mismatches, HGVS parsing, or mapping errors

---

## 10) Expected results (realistic for a strong CBIO project)

This section describes what a well-functioning pipeline should produce when end-to-end data ingestion and scoring are working.

### 10.1 Variant benchmark expected outcomes
- Pathogenic TP53 missense variants should skew more destabilizing than benign variants, but overlap is expected.
- AUROC and AUPRC should be meaningfully above chance if stability is capturing a real component of pathogenicity.

### 10.2 Rescue design expected outcomes
For each target mutant:
- a non-trivial Pareto front (dozens of candidates)
- multiple design strategies:
  - conservative, low-risk rescues with moderate stability improvement
  - higher-gain rescues that trade off some risk or complexity
- recurring stabilizing themes:
  - improved core packing at buried residues
  - avoidance of functional surfaces and conserved residues

### 10.3 Ablation expected outcomes
- Removing constraints increases apparent stability improvements but yields designs closer to protected residues or at highly conserved sites.
- The constrained version produces slightly smaller stability gains but much more plausible designs.

### 10.4 Robustness expected outcomes
- candidate ranking should remain stable across replicate scoring runs for top candidates
- Pareto fronts should remain broadly similar under small parameter changes

---

## 11) Science fair presentation plan

### 11.1 Poster sections
1. Motivation: why p53 stability rescue matters
2. Pipeline diagram: end-to-end overview
3. Constraint map: protected residues and designable region
4. Benchmark: benign vs pathogenic ΔΔG separation (ClinVar)
5. Case studies: R175H, R248Q, R273H
   - Pareto plots
   - top rescue tables
   - structure highlights
6. Ablations: “why constraints matter”
7. Limitations and next steps

### 11.2 Judge-proof talking points
- This is generative design under constraints, not a black-box classifier.
- Validation uses real variant labels to test whether stability scoring behaves sensibly.
- The strongest evidence is not one number, it is the combination of:
  - variant benchmark
  - Pareto tradeoffs
  - constraint ablations
  - robustness checks

---

## 12) Limitations and responsible framing

### 12.1 What the project cannot prove
- It cannot prove that designed rescues restore transcriptional activity in cells.
- It cannot prove therapeutic benefit.
- It cannot fully capture complex p53 biology (dominant-negative effects, interaction networks, context dependence).

### 12.2 Why it is still scientifically valuable
- It produces a reproducible, explainable design engine.
- It generates testable hypotheses prioritized by stability and safety constraints.
- It provides a rigorous computational basis for experimental follow-up.

---

## 13) Next build milestones (concrete checklist)

### Milestone A: Fix label sets and run benchmark
- Make sure benign and pathogenic label sets are non-empty.
- Produce:
  - ΔΔG distribution figure
  - AUROC/AUPRC table with CIs

### Milestone B: Add mutation-count objective
- Ensure singles and doubles compete on the Pareto front.
- Export 3D Pareto or separate “best per mutation count” tables.

### Milestone C: Add second stability scorer (optional)
- Implement consensus ranking to reduce dependence on one engine.

### Milestone D: Produce final report bundle
- `reports/figures/*` poster-quality
- `reports/tables/*` final tables for case studies
- `reports/reports.md` autogenerated writeup

---

## 14) References (starter set)

You should maintain a dedicated `references.bib` or `sources.md` with:
- p53 biology and hotspot context
- ClinVar documentation and related NAR paper
- AlphaFold DB paper
- FoldX paper and documentation
- p53 suppressor mutation literature (mechanism of rescue)

Add your final 6–10 core papers to the sources page in the report.
