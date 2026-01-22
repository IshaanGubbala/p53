Locked. Here’s a **p53 StabiliMut codebase plan at function level**, plus **exact data acquisition**, plus **how to build it optimized from day 1** (caching, parallel scoring, reproducibility, minimal recompute).

I’m designing this so you can actually finish it and so a judge can trace every result back to a file and a function.

---

## Project end product

A **software pipeline** that takes:

* a destabilizing TP53 missense mutant (ex: R175H, R248Q, R273H)
* the p53 core-domain structure

…and outputs:

1. top “rescue” mutation sets (1–3 second-site edits) predicted to **increase stability**
2. strict constraint checks: avoids DNA-binding residues, Zn-binding residues, highly conserved sites
3. Pareto front: stability gain vs functional-risk proxies
4. validation report: destabilizing vs benign variants separation using ClinVar labels

Everything is computational.

Key biological anchor: the p53 DNA-binding domain (core) is stabilized by a zinc ion coordinated by **Cys176, His179, Cys238, Cys242**. ([PMC][1])

---

## Data acquisition plan

### A) TP53 sequence + annotations

Source: UniProt P04637 (reviewed) ([UniProt][2])
You need:

* full sequence
* domain boundaries (core domain ~94–292 is common in literature)
* annotated functional regions (DNA-binding, etc.)

**How you’ll use it**

* to map variants to protein positions
* to define protected residues and motifs

### B) Structure model(s)

Primary: AlphaFold entry for P04637 ([alphafold.ebi.ac.uk][3])
Optional: also use experimental PDBs for the core domain if you want (not required for v1).

**How you’ll use it**

* compute solvent accessibility + burial
* compute distances to DNA-contact and Zn-binding residues
* run ΔΔG tools on a stable starting structure

### C) Variant labels for validation

ClinVar bulk downloads via FTP (XML/VCF) ([NCBI][4])
You’ll extract:

* missense variants in TP53
* clinical significance (benign/likely benign vs pathogenic/likely pathogenic)
* review status stars (optional for quality filtering)

Supplement: NCI TP53 Database somatic and germline datasets (tab-delimited) ([DCEG][5])
You’ll use this for:

* common tumor hotspots to demonstrate relevance
* “design targets” set (R175H, G245S, R248Q, R249S, R273H, R282W, etc.)

### D) ΔΔG scoring engine

FoldX BuildModel docs ([foldxsuite.crg.eu][6])
This is the default stability predictor in the plan because it’s widely used and scriptable.

You’ll implement the scoring layer so you can swap engines later.

---

## Optimized-from-start design rules

If you follow these, the project stays fast and finishable:

1. **Cache everything** by content hash

   * structure hash + mutation string → score
   * never recompute a mutation you’ve already scored

2. Use **Parquet** for big tables

   * variants, candidates, scores
   * fast load, fast groupby

3. Every heavy step supports **parallelism**

   * candidate scoring across mutation sets
   * use joblib multiprocessing

4. Everything is **config-driven**

   * domains, protected residues, tool paths, thresholds
   * so you can run ablations by changing YAML only

5. Reproducibility

   * fixed seeds
   * deterministic sampling
   * logged environment and tool versions

---

## Repository layout (detailed)

```
p53-stabilimut/
├── configs/
│   ├── p53.yaml
│   ├── paths.yaml
│   ├── scoring.yaml
│   ├── optimizer.yaml
│   └── report.yaml
│
├── data/
│   ├── raw/              # downloaded as-is
│   ├── interim/          # parsed + normalized
│   └── processed/        # candidates + scores
│
├── src/
│   ├── cli.py
│   ├── core/
│   │   ├── types.py
│   │   ├── hashing.py
│   │   ├── logging.py
│   │   └── config.py
│   │
│   ├── data_fetch/
│   │   ├── fetch_uniprot.py
│   │   ├── fetch_alphafold.py
│   │   ├── fetch_clinvar.py
│   │   └── fetch_nci_tp53.py
│   │
│   ├── variants/
│   │   ├── parse_clinvar_xml.py
│   │   ├── normalize_variants.py
│   │   ├── map_to_protein.py
│   │   └── label_sets.py
│   │
│   ├── target/
│   │   ├── p53_domain.py
│   │   ├── p53_protected.py
│   │   └── conservation.py
│   │
│   ├── structure/
│   │   ├── clean_structure.py
│   │   ├── residue_geometry.py
│   │   └── sasa_burial.py
│   │
│   ├── design/
│   │   ├── candidate_generator.py
│   │   ├── candidate_filters.py
│   │   └── mutation_sets.py
│   │
│   ├── scoring/
│   │   ├── foldx_runner.py
│   │   ├── stability_consensus.py
│   │   └── risk_scores.py
│   │
│   ├── optimize/
│   │   ├── beam_search.py
│   │   ├── evolutionary.py
│   │   └── pareto.py
│   │
│   ├── eval/
│   │   ├── variant_separation.py
│   │   ├── rescue_quality.py
│   │   ├── ablations.py
│   │   └── robustness.py
│   │
│   └── viz/
│       ├── plots_variants.py
│       ├── plots_pareto.py
│       ├── plots_structures.py
│       └── export_figures.py
│
└── experiments/
    ├── run_fetch.py
    ├── run_build_dataset.py
    ├── run_score_variants.py
    ├── run_design_rescues.py
    └── run_make_report.py
```

---

## Script-by-script and function-level plan

### `src/cli.py`

Purpose: single entrypoint with subcommands.

Functions:

* `build_parser() -> argparse.ArgumentParser`
* `main(argv: list[str]) -> int`
* Subcommands:

  * `fetch` (downloads)
  * `build-dataset` (parses + maps)
  * `score-variants` (ΔΔG on ClinVar sets)
  * `design` (generate rescue candidates for chosen mutants)
  * `report` (figures + tables)

How it does it:

* loads YAML configs once
* routes to pipeline functions in `experiments/`

---

# Data fetching

### `src/data_fetch/fetch_uniprot.py`

Goal: download TP53 sequence + feature annotations.

Functions:

* `download_uniprot_entry(uniprot_id: str, out_path: Path) -> None`
* `extract_sequence(entry_html_or_txt: str) -> str`
* `extract_features(entry_html_or_txt: str) -> dict`
* `save_fasta(seq: str, out_fasta: Path) -> None`
* `save_annotations(features: dict, out_json: Path) -> None`

Optimization:

* store raw response in `data/raw/uniprot/`
* store parsed JSON in `data/interim/annotations/`

Source: UniProt entry P04637 ([UniProt][2])

---

### `src/data_fetch/fetch_alphafold.py`

Goal: download AlphaFold structure for P04637.

Functions:

* `download_alphafold_pdb(uniprot_id: str, out_pdb: Path) -> None`
* `download_alphafold_pae(uniprot_id: str, out_json: Path) -> None` (optional)
* `verify_model(out_pdb: Path) -> dict` (length, missing residues, pLDDT summary)

Source: AlphaFold DB entry ([alphafold.ebi.ac.uk][3])

---

### `src/data_fetch/fetch_clinvar.py`

Goal: download ClinVar XML release.

Functions:

* `download_clinvar_xml_release(out_dir: Path) -> Path`
* `download_clinvar_vcf(out_dir: Path, assembly: str="GRCh38") -> Path` (optional)
* `record_release_metadata(out_dir: Path) -> None`

Source: ClinVar download documentation ([NCBI][4])

---

### `src/data_fetch/fetch_nci_tp53.py`

Goal: download NCI TP53 database tables for hotspot list and somatic variant frequency.

Functions:

* `download_nci_tp53_tables(out_dir: Path) -> list[Path]`
* `parse_hotspots(table_path: Path) -> list[dict]`
* `save_hotspots(hotspots: list[dict], out_path: Path) -> None`

Source: NCI TP53 Database download availability ([DCEG][5])

---

# Variant parsing and mapping

### `src/variants/parse_clinvar_xml.py`

Goal: extract TP53 missense variants + labels.

Functions:

* `iter_clinvar_variants(xml_path: Path) -> Iterator[dict]`
* `is_tp53(record: dict) -> bool`
* `is_missense(record: dict) -> bool`
* `extract_hgvs(record: dict) -> dict`  (c., p., transcript)
* `extract_significance(record: dict) -> dict`
* `to_flat_row(record: dict) -> dict`

Performance:

* streaming parse with lxml iterparse
* write Parquet in chunks

Output:

* `data/interim/parsed_variants.parquet`

---

### `src/variants/normalize_variants.py`

Goal: clean inconsistent variant representations.

Functions:

* `normalize_protein_change(p_hgvs: str) -> tuple[int, str, str]`

  * returns (pos, refAA, altAA)
* `filter_to_canonical_transcript(df) -> df`
* `drop_conflicted_labels(df) -> df`
* `dedupe_variants(df) -> df`

Output:

* `data/interim/tp53_missense_normalized.parquet`

---

### `src/variants/map_to_protein.py`

Goal: ensure variant positions match UniProt numbering.

Functions:

* `load_uniprot_sequence(fasta_path: Path) -> str`
* `validate_variant_ref(seq: str, pos: int, refAA: str) -> bool`
* `map_and_filter(df, seq) -> df`  (drop mismatched refs)
* `restrict_to_domain(df, start: int, end: int) -> df`

This avoids silent off-by-one errors.

---

### `src/variants/label_sets.py`

Goal: create clean validation cohorts.

Functions:

* `make_benign_set(df) -> df`
* `make_pathogenic_set(df) -> df`
* `quality_filter(df, min_review_stars: int=1) -> df` (optional)
* `export_label_sets(benign_df, path_df, out_dir) -> None`

---

# Target definition and constraints

### `src/target/p53_domain.py`

Goal: define domain boundaries and design region.

Functions:

* `get_core_domain() -> tuple[int,int]`  (ex: 94–312, configurable)
* `get_designable_region() -> list[int]`

  * core domain minus protected residues minus low-confidence residues (optional)

---

### `src/target/p53_protected.py`

Goal: define “do not mutate” residues.

Functions:

* `get_zinc_binding_residues() -> set[int]`

  * includes {176,179,238,242} ([PMC][1])
* `get_dna_contact_residues(annotations: dict) -> set[int]`

  * from UniProt features + optional literature list
* `get_hotspot_residues(nci_hotspots: list[dict]) -> set[int]`
* `combine_protected_sets(...) -> set[int]`

---

### `src/target/conservation.py`

Goal: conservation penalty via MSA.

Functions:

* `fetch_homologs(uniprot_id: str, out_fasta: Path) -> None`
* `run_mafft(in_fasta: Path, out_aln: Path) -> None`
* `compute_conservation(aln_path: Path) -> dict[int, float]`

  * returns per-position conservation score [0..1]
* `conservation_penalty(pos: int, cons_map: dict, gamma: float) -> float`

Optimization:

* cache alignment outputs
* compute conservation once

---

# Structure geometry

### `src/structure/clean_structure.py`

Goal: prepare a stable structure file.

Functions:

* `extract_domain_pdb(in_pdb: Path, start: int, end: int, out_pdb: Path) -> None`
* `remove_altloc(out_pdb: Path) -> None`
* `standardize_residue_numbering(out_pdb: Path) -> None`
* `validate_structure(out_pdb: Path) -> dict`

---

### `src/structure/sasa_burial.py`

Goal: compute burial to prioritize core packing mutations.

Functions:

* `compute_sasa(pdb: Path) -> dict[int, float]`
* `classify_burial(sasa_map, thresholds) -> dict[int, str]`  (buried/partial/exposed)

---

### `src/structure/residue_geometry.py`

Goal: distance-based functional risk.

Functions:

* `compute_distance_matrix(pdb: Path) -> dict[tuple[int,int], float]`
* `min_distance_to_set(pos: int, protected: set[int], dist_map) -> float`
* `risk_from_distance(d: float, cutoff: float) -> float`

---

# Candidate generation and filtering

### `src/design/mutation_sets.py`

Goal: represent mutation sets with deterministic IDs.

Functions:

* `format_mutation(pos:int, ref:str, alt:str) -> str`  (“R175H”)
* `canonicalize_set(muts: list[str]) -> tuple[str,...]`
* `mutation_set_id(muts: tuple[str,...]) -> str` (hash)

---

### `src/design/candidate_generator.py`

Goal: propose second-site stabilizers.

Functions:

* `select_design_positions(burial_map, protected, cons_map, top_n:int) -> list[int]`

  * prefers buried / partially buried
  * excludes protected
  * excludes high-conservation
* `propose_substitutions(pos:int, refAA:str) -> list[str]`

  * conservative palette (hydrophobic↔hydrophobic etc.)
* `generate_single_mutants(positions) -> list[MutationSet]`
* `expand_to_pairs(best_singles, max_pairs:int) -> list[MutationSet]`
* `expand_to_triples(best_pairs, max_triples:int) -> list[MutationSet]`

Optimization:

* “beam width” controls growth
* never explode combinatorics

---

### `src/design/candidate_filters.py`

Goal: reject risky designs before expensive scoring.

Functions:

* `filter_by_protected_distance(mset, dist_map, min_angstrom: float) -> bool`
* `filter_by_conservation(mset, cons_map, max_cons: float) -> bool`
* `filter_by_surface_patch(mset, burial_map) -> bool` (proxy)
* `apply_all_filters(candidates) -> list[candidates_passed]`

---

# Scoring layer (FoldX)

### `src/scoring/foldx_runner.py`

Goal: compute ΔΔG for mutation sets using FoldX BuildModel.

Functions:

* `prepare_foldx_workdir(base_pdb: Path, workdir: Path) -> None`
* `write_mutant_file(mset: MutationSet, out_file: Path) -> None` ([foldxsuite.crg.eu][7])
* `run_buildmodel(workdir: Path, mutant_file: Path, n_runs:int) -> Path` ([foldxsuite.crg.eu][6])
* `parse_foldx_ddg(output_files: list[Path]) -> float`
* `score_mutation_set(mset, pdb, cache) -> float`

Performance:

* hashed cache key = (pdb_hash, mutation_set_id)
* run in parallel across candidates
* reuse the same prepared structure and FoldX repaired model

---

### `src/scoring/risk_scores.py`

Goal: compute non-stability risk metrics.

Functions:

* `functional_risk(mset, dist_map, protected) -> float`
* `conservation_risk(mset, cons_map) -> float`
* `burial_risk(mset, burial_map) -> float` (avoid too-exposed changes)
* `aggregate_risk(weights, components) -> float`

---

### `src/scoring/stability_consensus.py`

Goal: combine multiple scorers if you add another tool later.

Functions:

* `consensus_ddg(ddg_foldx: float, ddg_other: Optional[float]) -> float`
* `rank_stability(ddg: float) -> float`  (lower is better)

---

# Optimization

### `src/optimize/pareto.py`

Goal: multi-objective ranking.

Functions:

* `dominates(a, b) -> bool`
* `pareto_front(items, objectives) -> list`
* `crowding_distance(front) -> list`

---

### `src/optimize/beam_search.py`

Goal: efficient search that never blows up.

Functions:

* `beam_step(current_sets, generator, scorer, beam_width:int) -> list`
* `run_beam_search(seed_mutant, constraints, beam_width, depth) -> list[MutationSet]`

The seed_mutant is your destabilizing cancer mutant; you search for second-site sets that rescue stability.

---

### `src/optimize/evolutionary.py`

Goal: optional richer search.

Functions:

* `init_population(candidates, k:int) -> list`
* `mutate(mset) -> mset`
* `crossover(a,b) -> child`
* `select(pop, objectives) -> next_pop`
* `run_ea(...) -> results`

You likely only need beam search for v1; EA is v2.

---

# Evaluation

### `src/eval/variant_separation.py`

Goal: prove ΔΔG correlates with known clinical labels.

Functions:

* `score_variant_set(df_variants, scorer) -> df_scored`
* `compute_auc(df_scored) -> dict`
* `plot_distributions(df_scored) -> fig`
* `bootstrap_ci(metric_fn, n=2000) -> (lo,hi)`

ClinVar FTP provides labels; you only use them for evaluation. ([NCBI][4])

---

### `src/eval/rescue_quality.py`

Goal: show rescues are meaningful and safe.

Functions:

* `select_target_mutants(nci_hotspots) -> list[str]` ([DCEG][5])
* `run_rescue_for_target(target_mut: str) -> results`
* `summarize_pareto(results) -> table`
* `sanity_checks(results) -> dict` (no protected violations, etc.)

---

### `src/eval/ablations.py`

Goal: judge-proof: show constraints matter.

Functions:

* `run_without_conservation_constraint()`
* `run_without_distance_constraint()`
* `compare_risk_profiles()`

You’ll show that removing constraints “improves” stability but increases risk, which justifies your design.

---

# Visualization

### `src/viz/plots_pareto.py`

Functions:

* `plot_pareto_front(results, x="risk", y="ddg")`
* `plot_tradeoff_matrix(results)`

### `src/viz/plots_structures.py`

Functions:

* `render_residue_highlights(pdb, residues, out_png)`
* `render_mutation_set(pdb, mset, out_png)`

### `src/viz/plots_variants.py`

Functions:

* `plot_ddg_by_label(df_scored)`
* `plot_hotspot_heatmap(df_scored)`

---

## Execution plan (what you run)

### 1) Fetch

`python -m src.cli fetch --all`

Outputs:

* UniProt sequence + features ([UniProt][2])
* AlphaFold PDB ([alphafold.ebi.ac.uk][3])
* ClinVar XML ([NCBI][4])
* NCI TP53 tables ([DCEG][5])

### 2) Build dataset

`python -m src.cli build-dataset`

Outputs:

* normalized TP53 missense variants table
* benign/pathogenic cohorts

### 3) Score variants

`python -m src.cli score-variants --engine foldx --parallel 8`

Outputs:

* ΔΔG distributions + AUC report

### 4) Design rescues for hotspot mutants

`python -m src.cli design --targets R175H R248Q R273H --max-muts 3`

Outputs:

* Pareto fronts
* Top rescue sets per target

### 5) Report

`python -m src.cli report`

Outputs:

* poster-grade figures and summary tables

---

## What we start coding first (optimized path)

To start building code “optimized already,” implement in this order:

1. `core/config.py`, `core/hashing.py`, `core/logging.py`
2. `fetch_uniprot.py`, `fetch_alphafold.py`
3. `clean_structure.py`, `sasa_burial.py`, `residue_geometry.py`
4. `fetch_clinvar.py` + streaming `parse_clinvar_xml.py`
5. `map_to_protein.py` (prevents silent mapping errors)
6. `foldx_runner.py` with caching and parallel scoring
7. candidate generator + filters
8. beam search + Pareto
9. reports + plots

That sequence prevents rewriting and keeps everything measurable early.

---


[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4922747/?utm_source=chatgpt.com "Reactivating mutant p53 using small molecules as zinc ..."
[2]: https://www.uniprot.org/uniprotkb/P04637/entry?utm_source=chatgpt.com "TP53 - Cellular tumor antigen p53 - Homo sapiens (Human)"
[3]: https://alphafold.ebi.ac.uk/entry/P04637?utm_source=chatgpt.com "P04637 - AlphaFold Protein Structure Database"
[4]: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/?utm_source=chatgpt.com "Accessing and using data in ClinVar - NCBI - NIH"
[5]: https://dceg.cancer.gov/tools/public-data/tp53-database?utm_source=chatgpt.com "The TP53 Database - NCI - DCEG - National Cancer Institute"
[6]: https://foldxsuite.crg.eu/command/BuildModel?utm_source=chatgpt.com "BuildModel"
[7]: https://foldxsuite.crg.eu/parameter/mutant-file?utm_source=chatgpt.com "mutant-file"
