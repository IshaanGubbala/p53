# p53 Rescue Mutations: Experimental Validation Protocol

## Overview

This protocol outlines a systematic experimental validation plan for p53 rescue mutations, organized into three phases from proof-of-concept to clinical development.

**Priority Candidates for Validation:**
1. **M133L** (Universal #1 rescue)
2. **A189S** (Universal #2 rescue)
3. **A189G** (Universal #3 rescue)

**Target Cancer Mutations:**
- R175H (DNA-binding domain hotspot)
- R248Q (DNA-binding domain hotspot)
- R273H (DNA-binding domain hotspot)
- Y220C (Structural hotspot)

---

## Phase 1: Biophysical Validation (3-6 months)

### Objective
Confirm that rescue mutations stabilize mutant p53 while preserving DNA binding and tetramerization.

---

### 1.1 Protein Expression and Purification

**Constructs to Generate:**
- Wild-type p53 (core domain, residues 94-312)
- Single mutants: R175H, R248Q, R273H, Y220C
- Rescued mutants:
  - R175H + M133L
  - R175H + A189S
  - R175H + A189G
  - R248Q + M133L
  - R248Q + A189S
  - R273H + M133L
  - R273H + A189S
  - Y220C + M133L
  - Y220C + A189S

**Expression System:** E. coli BL21(DE3)
- Vector: pET-based (N-terminal His6-tag)
- Growth: 37°C to OD₆₀₀ = 0.6, induce with 0.5 mM IPTG
- Expression: 18°C overnight

**Purification Protocol:**
1. Ni-NTA affinity chromatography
2. TEV protease cleavage (remove His-tag)
3. Size exclusion chromatography (SEC)
4. Concentrate to 5-10 mg/mL

**Expected Yield:** 5-10 mg per liter culture

**Buffer:** 20 mM Tris-HCl pH 7.5, 150 mM NaCl, 1 mM DTT, 5% glycerol

---

### 1.2 Thermal Stability Analysis

**Assay:** Differential Scanning Fluorimetry (DSF) / ThermoFluor

**Setup:**
- Protein concentration: 5 μM
- SYPRO Orange dye: 5X
- Buffer: PBS pH 7.4
- Temperature ramp: 20-95°C, 1°C/min
- Instrument: Real-time PCR machine

**Measurements:**
- Melting temperature (T_m): Temperature at 50% unfolding
- ΔT_m: Difference between mutant and wild-type
- ΔΔT_m: Difference between rescued and unrescued mutant

**Expected Results:**

| Construct | Expected T_m (°C) | Expected ΔT_m vs WT (°C) | Expected ΔΔT_m vs Mutant (°C) |
|-----------|-------------------|--------------------------|-------------------------------|
| Wild-type p53 | 42-44 | 0 | - |
| R175H | 35-37 | -7 to -9 | 0 |
| R175H + M133L | 40-43 | -2 to -4 | +5 to +7 ✅ |
| R175H + A189S | 39-42 | -3 to -5 | +4 to +6 ✅ |
| R248Q | 36-38 | -6 to -8 | 0 |
| R248Q + M133L | 41-43 | -1 to -3 | +5 to +7 ✅ |

**Success Criteria:**
- ✅ Rescued mutants show ΔΔT_m > +4°C compared to unrescued mutant
- ✅ Rescued mutants approach wild-type T_m (within 2-3°C)

**Replicates:** n=3 biological replicates, n=3 technical replicates

---

### 1.3 Protein Folding and Secondary Structure

**Assay:** Circular Dichroism (CD) Spectroscopy

**Setup:**
- Protein concentration: 0.2 mg/mL (5-10 μM)
- Buffer: 10 mM sodium phosphate pH 7.4, 50 mM NaF
- Cuvette: 1 mm path length quartz
- Wavelength range: 190-260 nm
- Temperature: 20°C

**Measurements:**
1. **Far-UV CD (190-260 nm):** Secondary structure content
   - α-helix: Negative peaks at 208 and 222 nm
   - β-sheet: Negative peak at 218 nm
   - Random coil: Negative peak at 195 nm

2. **Thermal denaturation (CD melting curve):**
   - Monitor ellipticity at 222 nm
   - Temperature ramp: 10-90°C, 1°C/min
   - Calculate T_m from inflection point

**Expected Results:**
- Wild-type and rescued mutants: ~30% α-helix, ~20% β-sheet (native fold)
- Unrescued mutants (R175H, R248Q): Reduced secondary structure (partially unfolded)
- Rescued mutants: Secondary structure restored to near wild-type levels

**Success Criteria:**
- ✅ Rescued mutants show CD spectrum similar to wild-type
- ✅ Secondary structure content within 10% of wild-type
- ✅ CD-derived T_m confirms DSF results

**Replicates:** n=2 biological replicates

---

### 1.4 DNA Binding Affinity

**Assay:** Electrophoretic Mobility Shift Assay (EMSA)

**DNA Probe:**
- p53 consensus sequence: 5'-GGACATGCCCGGGCATGTCC-3' (20-mer)
- 5' end-labeled with Cy5 or ³²P
- Final concentration: 10 nM

**Setup:**
- Protein concentrations: 0, 10, 25, 50, 100, 250, 500, 1000 nM
- Buffer: 20 mM HEPES pH 7.5, 50 mM KCl, 5 mM MgCl₂, 0.1 mg/mL BSA, 5% glycerol
- Incubation: 30 min at room temperature
- Gel: 6% native polyacrylamide, 0.5X TBE, run at 4°C
- Detection: Fluorescence imaging or autoradiography

**Measurements:**
- K_d (dissociation constant): Fit binding curve to one-site binding model
- Hill coefficient: Cooperativity of binding

**Expected Results:**

| Construct | Expected K_d (nM) | Expected Hill Coefficient |
|-----------|-------------------|---------------------------|
| Wild-type p53 | 20-50 | 2-3 (cooperative) |
| R175H | >500 | <1.5 (impaired) |
| R175H + M133L | 30-80 | 2-3 ✅ |
| R175H + A189S | 40-100 | 2-3 ✅ |
| R248Q | >1000 | <1.5 (very impaired) |
| R248Q + M133L | 50-150 | 2-3 ✅ |

**Success Criteria:**
- ✅ Rescued mutants show K_d < 200 nM (functional DNA binding)
- ✅ Rescued mutants show Hill coefficient > 1.5 (cooperative binding maintained)
- ✅ Rescued mutants bind within 2-5X of wild-type affinity

**Replicates:** n=3 biological replicates

---

### 1.5 Tetramerization Analysis

**Assay 1:** Size Exclusion Chromatography (SEC)

**Setup:**
- Column: Superdex 200 Increase 10/300 GL
- Buffer: 20 mM Tris-HCl pH 7.5, 150 mM NaCl, 1 mM DTT
- Flow rate: 0.5 mL/min
- Protein load: 100 μg at 2 mg/mL
- Detection: UV absorbance at 280 nm

**Measurements:**
- Elution volume (V_e): Compare to molecular weight standards
- Peak shape: Symmetry indicates homogeneous oligomeric state
- Calculate apparent molecular weight:
  - Monomer: ~25 kDa (core domain)
  - Dimer: ~50 kDa
  - Tetramer: ~100 kDa

**Expected Results:**
- Wild-type: Major peak at tetramer size (~100 kDa)
- Unrescued mutants: Mixed populations (monomer/dimer/tetramer)
- Rescued mutants: Restored tetramer peak similar to wild-type

**Assay 2:** Native PAGE

**Setup:**
- Gel: 6% native polyacrylamide, Tris-glycine buffer pH 8.3
- Protein load: 10 μg per lane at 1 mg/mL
- Run: 150V for 2 hours at 4°C
- Staining: Coomassie Blue

**Expected Results:**
- Wild-type: Single band at tetramer position
- Unrescued mutants: Multiple bands (monomer, dimer, tetramer)
- Rescued mutants: Predominant tetramer band

**Success Criteria:**
- ✅ Rescued mutants show >70% tetramer population (SEC)
- ✅ Rescued mutants show clear tetramer band (Native PAGE)
- ✅ Tetramer stability comparable to wild-type

**Replicates:** n=2 biological replicates

---

### 1.6 Zinc Binding (for R175H, R248Q, R273H)

**Assay:** Inductively Coupled Plasma Mass Spectrometry (ICP-MS)

**Setup:**
- Protein preparation: Extensive dialysis against Chelex-treated buffer
- Protein concentration: 50 μM
- Digestion: Dilute in 2% HNO₃
- ICP-MS measurement of Zn²⁺ content

**Expected Results:**
- Wild-type: 1.0 ± 0.1 Zn per protein molecule
- Unrescued mutants: 0.3-0.7 Zn (impaired coordination)
- Rescued mutants: 0.8-1.0 Zn ✅ (restored coordination)

**Alternative Assay:** 4-(2-pyridylazo)resorcinol (PAR) spectroscopic assay
- Zn-PAR complex: λ_max = 500 nm, ε = 66,000 M⁻¹cm⁻¹

**Success Criteria:**
- ✅ Rescued mutants show ≥0.8 Zn per protein
- ✅ Zinc coordination restored to near wild-type levels

---

## Phase 2: Cellular Functional Validation (6-12 months)

### Objective
Demonstrate that rescued p53 restores tumor suppressor function in cells.

---

### 2.1 Cell Line Selection

**Primary Cell Line:** H1299 (human non-small cell lung carcinoma)
- **Genotype:** p53-null (homozygous deletion)
- **Advantage:** Clean background, no endogenous p53 interference

**Secondary Cell Lines (optional):**
- SAOS-2 (osteosarcoma, p53-null)
- MEF p53⁻/⁻ (mouse embryonic fibroblasts)

---

### 2.2 Stable Cell Line Generation

**Lentiviral Expression System:**
- Vector: pLenti-CMV-p53-IRES-GFP (or similar)
- Promoter: CMV or EF1α (constitutive expression)
- Selection: Puromycin resistance or GFP sorting

**Constructs:**
- Empty vector (negative control)
- Wild-type p53 (positive control)
- R175H mutant (negative control)
- R175H + M133L (test)
- R175H + A189S (test)
- R248Q mutant (negative control)
- R248Q + M133L (test)

**Protocol:**
1. Transfect HEK293T cells with lentiviral packaging plasmids
2. Harvest viral supernatant at 48-72 hours
3. Transduce H1299 cells (MOI = 5)
4. Select with puromycin (2 μg/mL) for 7 days
5. Verify expression by Western blot

**Quality Control:**
- p53 expression level: Similar across all constructs (1-2X endogenous)
- Localization: Nuclear (immunofluorescence)
- Stability: No aggregation (cell fractionation)

---

### 2.3 Transcriptional Activity

**Reporter Assay:** p53-responsive luciferase reporter

**Setup:**
- Transfect p53 reporter construct:
  - pGL4.38[luc2P/p53 RE/Hygro] (Promega) or similar
  - Contains 4 copies of p53 response element (RE)
- Co-transfect Renilla luciferase (normalization control)
- Measure after 24-48 hours

**Measurements:**
- Luciferase activity (Firefly/Renilla ratio)
- Fold induction over empty vector

**Expected Results:**

| Construct | Expected Fold Induction | Expected Activity vs WT (%) |
|-----------|-------------------------|------------------------------|
| Empty vector | 1.0 | - |
| Wild-type p53 | 50-100 | 100% |
| R175H | 1-5 | <10% (loss of function) |
| R175H + M133L | 30-70 | 60-70% ✅ |
| R175H + A189S | 25-60 | 50-70% ✅ |

**Alternative Assay:** qRT-PCR of endogenous p53 target genes
- Measure mRNA levels of:
  - CDKN1A (p21) - cell cycle arrest
  - MDM2 - negative feedback
  - BAX - apoptosis
  - PUMA - apoptosis

**Success Criteria:**
- ✅ Rescued mutants show ≥50% of wild-type transcriptional activity
- ✅ Endogenous target genes upregulated (≥5-fold vs empty vector)

**Replicates:** n=3 biological replicates, n=3 technical replicates

---

### 2.4 Cell Cycle Arrest

**Assay:** Flow Cytometry (PI staining for cell cycle distribution)

**Setup:**
1. Induce p53 activity:
   - Option A: Nutlin-3 treatment (10 μM, 24h) - disrupts MDM2-p53 interaction
   - Option B: Doxorubicin treatment (0.5 μg/mL, 24h) - DNA damage
2. Harvest cells and fix (70% ethanol)
3. RNase treatment and propidium iodide (PI) staining
4. Analyze by flow cytometry

**Measurements:**
- % cells in G1, S, G2/M phases
- G1 arrest: Increase in G1 population after stress

**Expected Results:**

| Construct | % G1 Phase (Baseline) | % G1 Phase (Nutlin-3) | G1 Arrest (%) |
|-----------|------------------------|------------------------|---------------|
| Empty vector | 45-50% | 50-55% | +5% (no arrest) |
| Wild-type p53 | 45-50% | 70-75% | +25% (strong arrest) |
| R175H | 45-50% | 50-55% | +5% (no arrest) |
| R175H + M133L | 45-50% | 60-70% | +15-20% ✅ |
| R175H + A189S | 45-50% | 60-70% | +15-20% ✅ |

**Success Criteria:**
- ✅ Rescued mutants show ≥15% G1 arrest upon stress
- ✅ At least 60-75% of wild-type arrest activity

**Replicates:** n=3 biological replicates

---

### 2.5 Apoptosis Induction

**Assay:** Annexin V / PI Flow Cytometry

**Setup:**
1. Induce apoptosis:
   - Doxorubicin (1 μg/mL, 48h) or
   - Etoposide (50 μM, 48h)
2. Harvest cells (include floating cells)
3. Stain with Annexin V-FITC and PI
4. Analyze by flow cytometry

**Measurements:**
- Early apoptosis: Annexin V⁺ / PI⁻
- Late apoptosis: Annexin V⁺ / PI⁺
- Total apoptosis: Sum of early + late

**Expected Results:**

| Construct | % Total Apoptosis (Baseline) | % Total Apoptosis (Doxorubicin) |
|-----------|------------------------------|----------------------------------|
| Empty vector | 5-10% | 15-20% |
| Wild-type p53 | 5-10% | 40-50% (strong apoptosis) |
| R175H | 5-10% | 15-25% (impaired) |
| R175H + M133L | 5-10% | 30-45% ✅ |
| R175H + A189S | 5-10% | 28-40% ✅ |

**Alternative Assay:** Western blot for apoptotic markers
- Cleaved Caspase-3
- Cleaved PARP
- Cytochrome c release

**Success Criteria:**
- ✅ Rescued mutants show ≥30% apoptosis upon stress
- ✅ At least 60-75% of wild-type apoptotic activity

**Replicates:** n=3 biological replicates

---

### 2.6 Colony Formation Assay

**Assay:** Clonogenic Survival

**Setup:**
1. Seed cells at low density (500-1000 cells per 6-well dish)
2. Culture for 10-14 days in selection medium
3. Fix and stain with crystal violet
4. Count colonies (≥50 cells = 1 colony)

**Measurements:**
- Colony number
- Colony size (area)
- Plating efficiency: (# colonies / # seeded) × 100%

**Expected Results:**

| Construct | Expected Colony Number | Expected Plating Efficiency (%) |
|-----------|------------------------|----------------------------------|
| Empty vector | 200-300 | 40-60% (high proliferation) |
| Wild-type p53 | 50-100 | 10-20% (growth suppression) |
| R175H | 200-300 | 40-60% (no suppression) |
| R175H + M133L | 80-150 | 15-30% ✅ (restored suppression) |
| R175H + A189S | 90-160 | 18-32% ✅ |

**Success Criteria:**
- ✅ Rescued mutants show ≥50% reduction in colonies vs unrescued mutant
- ✅ Colony suppression ≥50% of wild-type activity

**Replicates:** n=3 biological replicates, n=3 technical replicates (wells)

---

## Phase 3: In Vivo Validation (12-36 months)

### Objective
Demonstrate tumor suppression in animal models.

---

### 3.1 Xenograft Tumor Growth

**Model:** Nude mice (athymic, immunodeficient)

**Setup:**
1. Inject H1299 cells expressing p53 constructs
   - Cell number: 5 × 10⁶ cells per mouse
   - Injection site: Subcutaneous (right flank)
   - Vehicle: Matrigel (50% v/v in PBS)
2. Monitor tumor growth:
   - Measure tumor dimensions (caliper) 2-3 times per week
   - Calculate volume: V = (L × W²) / 2
3. Endpoint: Tumor volume ≥1500 mm³ or 6-8 weeks

**Groups (n=8-10 mice per group):**
- Empty vector (negative control)
- Wild-type p53 (positive control)
- R175H mutant
- R175H + M133L
- R175H + A189S

**Expected Results:**

| Construct | Expected Tumor Volume at 6 weeks (mm³) | Expected Tumor Growth Delay (days) |
|-----------|----------------------------------------|-------------------------------------|
| Empty vector | 1200-1500 | 0 (baseline) |
| Wild-type p53 | 300-600 | +15-20 days ✅ |
| R175H | 1100-1400 | +1-3 days (no effect) |
| R175H + M133L | 500-900 | +10-15 days ✅ |
| R175H + A189S | 550-950 | +8-14 days ✅ |

**Survival Analysis:**
- Kaplan-Meier curves: Time to tumor volume ≥1000 mm³
- Log-rank test: Statistical significance (p < 0.05)

**Success Criteria:**
- ✅ Rescued mutants show ≥50% tumor growth reduction vs unrescued mutant
- ✅ Tumor growth delay ≥10 days vs unrescued mutant
- ✅ Survival significantly improved (p < 0.05)

**Histology:**
- H&E staining: Tumor architecture
- Ki-67 staining: Proliferation index (expected: 30-50% reduction)
- TUNEL assay: Apoptotic cells (expected: 2-3X increase)
- IHC for p53: Confirm expression

---

### 3.2 Orthotopic Tumor Model (Optional)

**Model:** Lung cancer orthotopic injection (H1299 cells into lung)

**Advantages:**
- Mimics natural tumor microenvironment
- Tests metastatic potential
- More clinically relevant

**Endpoint:**
- Lung tumor burden (bioluminescence imaging if cells are luciferase⁺)
- Metastasis (distant organ colonization)
- Survival time

---

### 3.3 Patient-Derived Xenograft (PDX) Model (Future)

**Setup:**
- Obtain tumor samples from cancer patients with R175H, R248Q, or R273H mutations
- Implant into NSG mice (NOD-SCID-gamma)
- Treat with AAV-delivered rescue mutations or small molecule stabilizers

**Goal:** Validate rescue strategy in human tumor context

---

## Timeline and Resources

### Phase 1: Biophysical (3-6 months)
**Personnel:**
- 1 postdoc or senior grad student
- 1 undergrad assistant

**Equipment Required:**
- Real-time PCR machine (DSF)
- CD spectrometer
- Gel imaging system (EMSA)
- FPLC system (SEC)

**Estimated Cost:** $15,000-20,000
- Reagents and consumables: $10,000
- DNA synthesis and cloning: $3,000
- Protein purification columns: $2,000

### Phase 2: Cellular (6-12 months)
**Personnel:**
- 1 postdoc or grad student
- 1 undergrad assistant

**Equipment Required:**
- Tissue culture facility
- Flow cytometer
- Plate reader (luciferase assays)
- Fluorescence microscope

**Estimated Cost:** $25,000-35,000
- Cell culture reagents: $8,000
- Lentivirus production: $5,000
- Flow cytometry: $4,000
- Assay kits: $8,000

### Phase 3: In Vivo (12-36 months)
**Personnel:**
- 1 postdoc
- Animal care technician

**Equipment Required:**
- Animal facility (IACUC-approved)
- Calipers for tumor measurement
- Imaging system (optional)

**Estimated Cost:** $80,000-120,000
- Mouse purchase and housing: $40,000
- Tumor monitoring: $10,000
- Histology and IHC: $15,000
- Labor: $35,000

**Total Estimated Cost:** $120,000-175,000 over 2-3 years

---

## Success Milestones

### Milestone 1 (Month 6): Biophysical Validation Complete
- ✅ At least 2/3 top rescues show ΔΔT_m > +4°C
- ✅ DNA binding restored (K_d < 200 nM)
- ✅ Tetramer formation confirmed

**Decision Point:** Proceed to cellular validation

### Milestone 2 (Month 12): Cellular Validation Complete
- ✅ Transcriptional activity ≥50% of wild-type
- ✅ G1 arrest ≥15% upon stress
- ✅ Apoptosis ≥30% upon stress

**Decision Point:** Proceed to in vivo validation

### Milestone 3 (Month 24): Xenograft Proof-of-Concept
- ✅ Tumor growth reduced ≥50% vs unrescued mutant
- ✅ Survival significantly improved (p < 0.05)

**Decision Point:** Consider therapeutic development or PDX models

---

## Publication Strategy

### Manuscript 1: Computational Design (Submit immediately)
**Title:** "Systematic Design of Rescue Mutations for Cancer-Associated p53 Mutants Using Pareto Optimization and Multi-Dimensional Validation"

**Highlights:**
- Comprehensive computational pipeline
- 106 universal rescue candidates across 4 hotspots
- DNA binding, tetramer interface, and sensitivity analysis
- Novel discovery: R196 is structural, not DNA contact

**Target Journal:** Nature Communications, eLife, or PNAS

### Manuscript 2: Experimental Validation (Submit after Phase 1-2)
**Title:** "M133L and A189S Rescue Mutations Restore Tumor Suppressor Function to Cancer-Associated p53 Mutants"

**Highlights:**
- Biophysical validation (T_m, DNA binding, tetramer)
- Cellular functional validation (transcription, cell cycle, apoptosis)
- Mechanism: Stabilization without functional disruption

**Target Journal:** Science Translational Medicine, Cancer Cell, or Molecular Cell

### Manuscript 3: Therapeutic Translation (Submit after Phase 3)
**Title:** "In Vivo Tumor Suppression by Rescued p53 Mutants: Therapeutic Implications for Cancer Treatment"

**Highlights:**
- Xenograft tumor growth suppression
- Survival benefit
- Rationale for small molecule stabilizer development

**Target Journal:** Nature Medicine, Science Translational Medicine, or Clinical Cancer Research

---

## Alternative Approaches (If Primary Validation Fails)

### If DNA Binding Is Impaired:
- **Hypothesis:** Rescue stabilizes protein but reduces DNA affinity
- **Solution:** Test double rescues (e.g., M133L + A189S together)
- **Alternative:** Focus on rescues in different regions (test Y163F, S215A)

### If Tetramerization Is Disrupted:
- **Hypothesis:** Rescue affects dimer-dimer interface
- **Solution:** Use full-length p53 (includes tetramerization domain)
- **Alternative:** Test rescues that don't touch interface (e.g., buried core only)

### If Cellular Function Is Weak:
- **Hypothesis:** Context-dependent effects (post-translational modifications)
- **Solution:** Use isogenic cell lines (CRISPR knock-in of rescues)
- **Alternative:** Test in multiple cell lines (breast, colon, lung cancer)

---

## Conclusion

This protocol provides a comprehensive roadmap for experimental validation of p53 rescue mutations, from biophysical characterization to in vivo tumor suppression.

**Key Advantages:**
1. Phased approach allows early decision-making
2. Multiple orthogonal assays reduce false positives
3. Quantitative success criteria at each stage
4. Flexible alternatives if primary strategy fails

**Timeline:** 2-3 years from initiation to publication-ready results

**Expected Outcome:** Demonstration that M133L and A189S are universal rescues that restore tumor suppressor function to multiple cancer-associated p53 mutants.

---

*Protocol Version 1.0*
*Last Updated: January 25, 2026*
*Contact: [Your Lab/Institution]*
