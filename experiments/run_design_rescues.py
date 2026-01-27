from __future__ import annotations

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from tqdm import tqdm
import sys

# Add project root to path
sys.path.append(str(Path.cwd()))

from src.core.hashing import dict_sha256, file_sha256
from src.core.logging import get_logger
from src.design.candidate_filters import apply_all_filters
from src.design.candidate_generator import generate_single_mutants, select_design_positions
from src.design.checkpointing import (
    compute_config_hash,
    delete_checkpoint,
    load_checkpoint,
    save_checkpoint,
    should_resume,
)
from src.design.mutation_sets import canonicalize_set, mutation_set_id
from src.scoring.evoef2_runner import compute_stability, score_mutation_set
from src.scoring.multi_structure import score_mutation_set_multi
# [RaSP Integration]
from src.scoring.rasp_scorer import score_mutation_set_rasp
from src.scoring.risk_scores import (
    aggregate_risk,
    burial_risk,
    conservation_risk,
    functional_risk,
    msa_conservation_risk,
)
from src.target.p53_protected import (
    combine_protected_sets,
    get_dna_contact_residues,
    get_hotspot_residues,
    get_zinc_binding_residues,
)
from src.scoring.allosteric_scorer import allosteric_distance_score
from src.scoring.docking_scorer import run_docking
from src.scoring.rasp_scorer import score_mutation_set_rasp


MutationSet = tuple[str, ...]
MUTATION_RE = re.compile(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$")


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_root = Path(paths_cfg.get("data_root", "data"))
    data_cfg = paths_cfg.get("data", {})
    
    raw_dir = Path(data_cfg.get("raw", "raw"))
    if not raw_dir.is_absolute(): raw_dir = data_root / raw_dir
        
    interim_dir = Path(data_cfg.get("interim", "interim"))
    if not interim_dir.is_absolute(): interim_dir = data_root / interim_dir
        
    processed_dir = Path(data_cfg.get("processed", "processed"))
    if not processed_dir.is_absolute(): processed_dir = data_root / processed_dir
        
    cache_dir = Path(paths_cfg.get("cache_dir", "processed/cache"))
    if not cache_dir.is_absolute(): cache_dir = data_root / cache_dir
        
    project_root = Path(paths_cfg.get("project_root", Path.cwd()))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
        "processed": processed_dir,
        "cache": cache_dir,
        "project_root": project_root,
    }


def _structure_basename(uniprot_id: str, start: int, end: int) -> str:
    return f"{uniprot_id}_core_{start}_{end}"


def _load_json_map(path: Path) -> dict[int, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {int(key): value for key, value in payload.items()}


def _load_distance_map(path: Path) -> dict[tuple[int, int], float]:
    df = pd.read_parquet(path)
    dist_map: dict[tuple[int, int], float] = {}
    for row in df.itertuples(index=False):
        i = int(row.pos_i)
        j = int(row.pos_j)
        dist = float(row.distance)
        dist_map[(i, j)] = dist
        dist_map[(j, i)] = dist
    return dist_map


def _parse_mutation(token: str) -> tuple[int, str, str]:
    match = MUTATION_RE.match(token.strip())
    if not match:
        raise ValueError(f"Invalid mutation token: {token}")
    ref, pos_text, alt = match.groups()
    return int(pos_text), ref.upper(), alt.upper()


def _load_sequence(fasta_path: Path) -> str:
    lines = []
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    return "".join(lines)


def _load_conservation_map(p53_cfg: dict[str, Any], interim_dir: Path) -> dict[int, float]:
    cons_cfg = p53_cfg.get("conservation", {})
    path_value = cons_cfg.get("path")
    if not path_value:
        return {}
    cons_path = Path(path_value)
    if not cons_path.is_absolute():
        cons_path = (interim_dir / cons_path).resolve()
    if not cons_path.exists():
        return {}
    payload = json.loads(cons_path.read_text(encoding="utf-8"))
    return {int(k): float(v) for k, v in payload.items()}


def _load_msa_conservation_map(p53_cfg: dict[str, Any], processed_dir: Path) -> dict[int, float]:
    """Load MSA-based conservation scores if enabled."""
    msa_cfg = p53_cfg.get("msa", {})
    if not msa_cfg.get("enabled", False):
        return {}
    path_value = msa_cfg.get("precomputed_path")
    if not path_value:
        return {}
    msa_path = Path(path_value)
    if not msa_path.is_absolute():
        msa_path = (processed_dir / msa_path).resolve()
    if not msa_path.exists():
        return {}
    payload = json.loads(msa_path.read_text(encoding="utf-8"))
    return {int(k): float(v) for k, v in payload.items()}


def _filter_targets(targets: Iterable[str]) -> list[str]:
    cleaned = []
    for token in targets:
        token = token.strip().upper()
        if not token:
            continue
        if not MUTATION_RE.match(token):
            raise ValueError(f"Invalid target mutation: {token}")
        cleaned.append(token)
    if not cleaned:
        raise ValueError("No valid target mutations provided")
    return cleaned


def _seed_ddg_from_scores(scores_path: Path) -> dict[str, float]:
    if not scores_path.exists():
        return {}
    df = pd.read_parquet(scores_path)
    if not {"pos", "ref", "alt", "ddg"}.issubset(df.columns):
        return {}
    df = df.copy()
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    df["mutation"] = df["ref"] + df["pos"].astype(str) + df["alt"]
    return dict(zip(df["mutation"], df["ddg"]))


def _validate_seed(sequence: str, seed: str) -> None:
    pos, ref, _ = _parse_mutation(seed)
    if pos < 1 or pos > len(sequence):
        raise ValueError(f"Seed position out of range: {seed}")
    if sequence[pos - 1].upper() != ref.upper():
        raise ValueError(f"Seed reference mismatch vs UniProt sequence: {seed}")


def _rank_candidates(
    candidates: list[MutationSet],
    scores: dict[MutationSet, dict[str, Any]],
    limit: int,
) -> list[MutationSet]:
    ranked = sorted(candidates, key=lambda mset: scores[mset]["ddg_gain"])
    return ranked[:limit] if limit else ranked


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))

    p53_cfg = configs.get("p53", {})
    optimizer_cfg = configs.get("optimizer", {})
    design_cfg = optimizer_cfg.get("design", {})
    beam_cfg = optimizer_cfg.get("beam_search", {})
    pareto_cfg = optimizer_cfg.get("pareto", {})

    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    core_domain = p53_cfg.get("core_domain", {})
    domain_start = int(core_domain.get("start", 1))
    domain_end = int(core_domain.get("end", 9999))

    cli_targets = getattr(args, "targets", None)
    if cli_targets:
        targets = _filter_targets(cli_targets)
    else:
        targets = p53_cfg.get("design", {}).get("targets") or p53_cfg.get("design_targets")
        if not targets:
            raise RuntimeError("No design targets provided; pass --targets")
        targets = _filter_targets(targets)

    max_mutations = int(getattr(args, "max_muts", None) or p53_cfg.get("design", {}).get("max_mutations", 3))
    max_mutations = max(1, min(max_mutations, 3))

    top_positions = int(design_cfg.get("top_positions", 80))
    min_distance = float(design_cfg.get("min_distance_protected", 8.0))
    max_cons = float(design_cfg.get("max_conservation", 0.8))
    max_msa_cons = float(design_cfg.get("max_msa_conservation", 0.85))
    allow_exposed = bool(design_cfg.get("allow_exposed", False))
    max_singles = int(design_cfg.get("max_singles", 200))
    max_pairs = int(design_cfg.get("max_pairs", 300))
    max_triples = int(design_cfg.get("max_triples", 200))

    beam_width = int(getattr(args, "beam_width", None) or beam_cfg.get("beam_width", 50))
    depth = int(getattr(args, "depth", None) or beam_cfg.get("depth", max_mutations))
    depth = min(depth, max_mutations)

    # Checkpointing configuration
    checkpoint_cfg = beam_cfg.get("checkpointing", {})
    checkpointing_enabled = checkpoint_cfg.get("enabled", True)
    verify_config_hash = checkpoint_cfg.get("verify_config_hash", True)
    recompute = getattr(args, "recompute", False)

    scoring_cfg = configs.get("scoring", {})
    evoef2_cfg = dict(scoring_cfg.get("evoef2", {}))
    if not evoef2_cfg:
        raise RuntimeError("Missing evoef2 configuration in configs/scoring.yaml")

    parallel_cfg = scoring_cfg.get("parallel", {})
    parallel = int(getattr(args, "parallel", None) or parallel_cfg.get("n_jobs", 2))
    parallel = max(1, parallel)

    # [RaSP Integration]
    rasp_cfg = scoring_cfg.get("rasp", {})
    rasp_enabled = rasp_cfg.get("enabled", False)
    if rasp_enabled:
        logger.info("RaSP scoring enabled.")

    # [Allosteric & Docking Integration]
    allosteric_cfg = scoring_cfg.get("allosteric", {})
    allosteric_enabled = allosteric_cfg.get("enabled", False)
    
    docking_cfg = scoring_cfg.get("docking", {})
    docking_enabled = docking_cfg.get("enabled", False)
    ligand_name = docking_cfg.get("ligand", "PHI-KAN-083")

    cache_cfg = scoring_cfg.get("cache", {})
    cache_dir = Path(cache_cfg.get("dir", paths["cache"]))
    cache_dir.mkdir(parents=True, exist_ok=True)
    work_root = cache_dir / "evoef2_design"
    work_root.mkdir(parents=True, exist_ok=True)

    fasta_path = paths["interim"] / "uniprot" / f"{uniprot_id}.fasta"
    sequence = _load_sequence(fasta_path)

    structure_dir = paths["interim"] / "structure"
    base_name = _structure_basename(uniprot_id, domain_start, domain_end)
    burial_path = structure_dir / f"{base_name}_burial.json"
    dist_path = structure_dir / f"{base_name}_distances.parquet"

    burial_map = _load_json_map(burial_path)
    dist_map = _load_distance_map(dist_path)

    annotations_path = paths["interim"] / "annotations" / f"{uniprot_id}.json"
    hotspots_path = paths["interim"] / "nci_tp53" / "hotspots.json"
    protected_manual = set(p53_cfg.get("design", {}).get("protected_residues", []))

    protected = combine_protected_sets(
        get_zinc_binding_residues(p53_cfg),
        get_dna_contact_residues(annotations_path),
        get_hotspot_residues(hotspots_path),
        protected_manual,
    )

    cons_map = _load_conservation_map(p53_cfg, paths["interim"])
    msa_cons_map = _load_msa_conservation_map(p53_cfg, paths["processed"])
    if msa_cons_map:
        logger.info("MSA conservation scores loaded for %d positions", len(msa_cons_map))

    # Check if multi-structure scoring is enabled
    use_multi_structure = "structures" in evoef2_cfg and len(evoef2_cfg.get("structures", [])) > 1
    if use_multi_structure:
        logger.info("Multi-structure scoring enabled with %d structures", len(evoef2_cfg["structures"]))
        # Validate structure paths
        for struct_cfg in evoef2_cfg["structures"]:
            struct_pdb = Path(struct_cfg["pdb"])
            if not struct_pdb.is_absolute():
                struct_pdb = (paths["project_root"] / struct_pdb).resolve()
                struct_cfg["pdb"] = str(struct_pdb)
            if not struct_pdb.exists():
                logger.warning("Structure %s not found at %s", struct_cfg["id"], struct_pdb)

    base_pdb = Path(evoef2_cfg.get("repaired_pdb") or paths["raw"] / "alphafold").expanduser()
    if not base_pdb.is_absolute():
        base_pdb = (paths["project_root"] / base_pdb).resolve()
    if base_pdb.is_dir():
        pdb_candidates = list(base_pdb.glob("*.pdb"))
        if not pdb_candidates:
            raise RuntimeError(f"No PDB files found in {base_pdb}")
        base_pdb = max(pdb_candidates, key=lambda path: path.stat().st_mtime)
    if not base_pdb.exists():
        raise RuntimeError(f"Base PDB not found: {base_pdb}")

    if rasp_enabled:
        # Pre-initialize RaSP for base_pdb to avoid race conditions in parallel loop
        try:
            logger.info("Pre-initializing RaSP for base structural model...")
            score_mutation_set_rasp([], base_pdb, cache_dir, rasp_cfg)
            logger.info("RaSP initialization complete.")
        except Exception as e:
            logger.warning(f"RaSP pre-initialization failed: {e}")

    base_energy_dir = work_root / "base_energy"
    base_energy_dir.mkdir(parents=True, exist_ok=True)
    base_energy = compute_stability(base_pdb, evoef2_cfg, base_energy_dir)
    seed_ddg_lookup = _seed_ddg_from_scores(paths["processed"] / "variant_scores.parquet")

    design_out_root = paths["processed"] / "rescues"
    design_out_root.mkdir(parents=True, exist_ok=True)

    risk_weights = design_cfg.get("risk_weights", {"functional": 0.5, "conservation": 0.3, "burial": 0.2})

    for seed in targets:
        _validate_seed(sequence, seed)
        seed_pos, _, _ = _parse_mutation(seed)
        protected_for_design = set(protected)
        protected_for_design.add(seed_pos)

        positions = select_design_positions(
            burial_map,
            protected_for_design,
            cons_map=cons_map,
            top_n=top_positions,
            max_cons=max_cons,
            allow_exposed=allow_exposed,
        )

        single_sets = generate_single_mutants(positions, sequence)
        single_sets = apply_all_filters(
            single_sets,
            dist_map=dist_map,
            protected=protected,
            burial_map=burial_map,
            min_angstrom=min_distance,
            cons_map=cons_map,
            max_cons=max_cons,
            msa_cons_map=msa_cons_map,
            max_msa_cons=max_msa_cons,
            allow_exposed=allow_exposed,
        )

        candidate_tokens = sorted({mset[0] for mset in single_sets})
        if not candidate_tokens:
            raise RuntimeError(f"No candidate positions found for target {seed}")

        seed_ddg = seed_ddg_lookup.get(seed)
        if seed_ddg is None:
            seed_ddg = score_mutation_set([seed], base_pdb, cache_dir, evoef2_cfg, work_root, base_energy=base_energy)

        # Compute RaSP seed (cancer mutation alone) for gain calculation
        rasp_seed = 0.0
        if rasp_enabled:
            try:
                rasp_seed = score_mutation_set_rasp([seed], base_pdb, cache_dir, rasp_cfg)
                logger.info(f"RaSP seed for {seed}: {rasp_seed:+.3f} kcal/mol")
            except Exception as e:
                logger.warning(f"RaSP seed scoring failed for {seed}: {e}")

        ddg_cache: dict[MutationSet, dict[str, Any]] = {}

        def _score_candidate(mset: MutationSet) -> dict[str, Any]:
            if mset in ddg_cache:
                return ddg_cache[mset]
            full_set = canonicalize_set([seed, *mset])

            # Choose scoring method based on configuration
            if use_multi_structure:
                # Multi-structure scoring
                multi_result = score_mutation_set_multi(
                    list(full_set),
                    evoef2_cfg["structures"],
                    evoef2_cfg,
                    cache_dir,
                    work_root,
                    consensus_method=evoef2_cfg.get("consensus_method", "median"),
                    require_all=evoef2_cfg.get("require_all", False),
                )
                ddg_total = multi_result["ddg_consensus"]
                # Store per-structure scores for later output
                structure_scores = {k: v for k, v in multi_result.items() if k.startswith("ddg_")}
            else:
                # Single-structure scoring
                ddg_total = score_mutation_set(
                    list(full_set),
                    base_pdb,
                    cache_dir,
                    evoef2_cfg,
                    work_root,
                    base_energy=base_energy,
                )
                structure_scores = {}
            
            # [RaSP Integration]
            rasp_score = 0.0
            rasp_gain = 0.0
            if rasp_enabled:
                try:
                    rasp_score = score_mutation_set_rasp(
                        list(full_set),
                        base_pdb, # RaSP usually uses the base PDB
                        cache_dir,
                        rasp_cfg
                    )
                    # Compute RaSP gain (comparable to EvoEF2 gain)
                    rasp_gain = rasp_seed - rasp_score
                except Exception as e:
                    logger.warning(f"RaSP scoring failed for {mset}: {e}")
                    rasp_score = 0.0 # Fallback
                    rasp_gain = 0.0

            # [Allosteric & Docking Integration]
            allo_dist = 999.0
            if allosteric_enabled:
                allo_dist = allosteric_distance_score(
                    list(full_set), 
                    dist_map, 
                    sites=allosteric_cfg.get("sites", ["Y220_POCKET"])
                )
            
            docking_affinity = 0.0
            if docking_enabled:
                # Docker usually runs on the structure. 
                # Since we don't necessarily have the mutant PDB yet (unless cached/available), 
                # we use the base PDB + mutations context or just a mock affinity for now.
                # Real implementation might need the repaired PDB of the mutant.
                docking_affinity = run_docking(base_pdb, ligand_name, mock=True)

            ddg_gain = ddg_total - seed_ddg

            # Compute risk components
            func_risk = functional_risk(mset, dist_map, protected, cutoff=min_distance)
            cons_risk = conservation_risk(mset, cons_map)
            bur_risk = burial_risk(mset, burial_map)
            msa_risk = msa_conservation_risk(mset, msa_cons_map)

            risk_components = {
                "functional": func_risk,
                "conservation": cons_risk,
                "burial": bur_risk,
                "msa_conservation": msa_risk,
            }
            risk = aggregate_risk(risk_weights, risk_components)

            payload = {
                "ddg_total": ddg_total,
                "ddg_gain": ddg_gain,
                "rasp_ddg": rasp_score,
                "rasp_gain": rasp_gain,  # Comparable to ddg_gain
                "allosteric_dist": allo_dist,
                "docking_affinity": docking_affinity,
                "risk": risk,
                "risk_components": risk_components,
            }

            # Add multi-structure metadata if available
            if use_multi_structure:
                payload.update(structure_scores)
                payload["structures_scored"] = multi_result.get("structures_scored", 0)
                payload["ddg_std"] = multi_result.get("ddg_std", 0.0)
                payload["ddg_range"] = multi_result.get("ddg_range", 0.0)

            ddg_cache[mset] = payload
            return payload

        def _score_candidates(candidates: list[MutationSet], desc: str) -> list[MutationSet]:
            if not candidates:
                return []
            results: dict[MutationSet, dict[str, Any]] = {}
            errors: list[str] = []
            with ThreadPoolExecutor(max_workers=parallel) as executor:
                futures = {executor.submit(_score_candidate, mset): mset for mset in candidates}
                for future in tqdm(as_completed(futures), total=len(futures), desc=desc, unit="set"):
                    mset = futures[future]
                    try:
                        results[mset] = future.result()
                    except Exception as exc:
                        errors.append(f"{mset}: {exc}")
            if errors:
                sample = "\n".join(errors[:5])
                raise RuntimeError(f"Scoring failed for {len(errors)} candidates:\n{sample}")
            for mset, payload in results.items():
                ddg_cache[mset] = payload
            return list(results.keys())

        candidate_positions = {token: _parse_mutation(token)[0] for token in candidate_tokens}

        def _expand_candidates(current: list[MutationSet]) -> list[MutationSet]:
            expanded: list[MutationSet] = []
            for mset in current:
                used_positions = {_parse_mutation(mut)[0] for mut in mset}
                for token in candidate_tokens:
                    pos = candidate_positions[token]
                    if pos in used_positions:
                        continue
                    expanded_set = canonicalize_set([*mset, token])
                    expanded.append(expanded_set)
            return expanded

        seed_out_dir = design_out_root / seed
        seed_out_dir.mkdir(parents=True, exist_ok=True)

        # Setup checkpointing
        checkpoint_dir = seed_out_dir / "checkpoints"
        checkpoint_config = {
            "seed": seed,
            "beam_width": beam_width,
            "depth": depth,
            "design_cfg": design_cfg,
            "beam_cfg": beam_cfg,
        }
        config_hash = compute_config_hash(checkpoint_config)

        # Check for existing checkpoint
        checkpoint = None
        start_step = 1
        if checkpointing_enabled and not recompute:
            checkpoint = load_checkpoint(
                checkpoint_dir,
                expected_config_hash=config_hash if verify_config_hash else None,
                verify_config=verify_config_hash,
            )

        if checkpoint:
            logger.info(f"Resuming from checkpoint at step {checkpoint['step']} for {seed}")
            start_step = checkpoint["step"] + 1
            results_list = checkpoint["candidates"].to_dict(orient="records")
            current = checkpoint["beam_state"]

            # Reconstruct ddg_cache from checkpoint
            for row in results_list:
                mset_str = row.get("rescue_mutations", "")
                if mset_str:
                    mset = tuple(mset_str.split(","))
                else:
                    mset = tuple()
                ddg_cache[mset] = {
                    "ddg_total": row["ddg_total"],
                    "ddg_gain": row["ddg_gain"],
                    "rasp_ddg": row.get("rasp_ddg", 0.0),
                    "allosteric_dist": row.get("allosteric_dist", 999.0),
                    "docking_affinity": row.get("docking_affinity", 0.0),
                    "risk": row["risk"],
                    "risk_components": json.loads(row["risk_components"]),
                }
        else:
            if recompute and checkpoint_dir.exists():
                logger.info(f"Deleting existing checkpoint (--recompute flag)")
                delete_checkpoint(checkpoint_dir)
            results_list = []
            current = [tuple()]

        results: list[dict[str, Any]] = results_list
        step_limits = {1: max_singles, 2: max_pairs, 3: max_triples}

        for step in range(start_step, depth + 1):
            candidates = _expand_candidates(current)
            # Deduplicate expanded sets
            candidates = sorted(list(set(candidates)))
            candidates = apply_all_filters(
                candidates,
                dist_map=dist_map,
                protected=protected,
                burial_map=burial_map,
                min_angstrom=min_distance,
                cons_map=cons_map,
                max_cons=max_cons,
                msa_cons_map=msa_cons_map,
                max_msa_cons=max_msa_cons,
                allow_exposed=allow_exposed,
            )

            scored = _score_candidates(candidates, desc=f"{seed} step {step}")
            ranked = _rank_candidates(scored, ddg_cache, step_limits.get(step, 0))
            current = ranked[:beam_width]

            if step == 1:
                candidate_tokens = sorted({mset[0] for mset in ranked})
                candidate_positions = {token: _parse_mutation(token)[0] for token in candidate_tokens}

            for mset in ranked:
                info = ddg_cache[mset]
                result_dict = {
                    "target": seed,
                    "rescue_mutations": ",".join(mset),
                    "full_mutations": ",".join(canonicalize_set([seed, *mset])),
                    "n_rescue": len(mset),
                    "ddg_seed": seed_ddg,
                    "ddg_total": info["ddg_total"],
                    "ddg_gain": info["ddg_gain"],
                    "rasp_seed": rasp_seed,  # RaSP score of cancer alone
                    "rasp_ddg": info.get("rasp_ddg", 0.0),  # RaSP score of cancer+rescue
                    "rasp_gain": info.get("rasp_gain", 0.0),  # RaSP rescue benefit
                    "allosteric_dist": info.get("allosteric_dist", 999.0),
                    "docking_affinity": info.get("docking_affinity", 0.0),
                    "risk": info["risk"],
                    "risk_components": json.dumps(info["risk_components"], sort_keys=True),
                    "set_id": mutation_set_id(mset),
                }

                # Add multi-structure columns if available
                if use_multi_structure:
                    # Add per-structure scores
                    for key in ["ddg_alphafold", "ddg_2ocj_core", "ddg_consensus"]:
                        if key in info:
                            result_dict[key] = info[key]
                    # Add metadata
                    result_dict["structures_scored"] = info.get("structures_scored", 0)
                    result_dict["ddg_std"] = info.get("ddg_std", 0.0)
                    result_dict["ddg_range"] = info.get("ddg_range", 0.0)

                results.append(result_dict)

            # Save checkpoint after each step
            if checkpointing_enabled:
                results_df_checkpoint = pd.DataFrame(results)
                save_checkpoint(
                    checkpoint_dir=checkpoint_dir,
                    step=step,
                    candidates=results_df_checkpoint,
                    beam_state=current,
                    config_hash=config_hash,
                    metadata={"seed": seed, "beam_width": beam_width, "depth": depth},
                )

        if not results:
            raise RuntimeError(f"No rescue candidates generated for {seed}")

        results_df = pd.DataFrame(results)
        pareto_objectives = pareto_cfg.get("objectives", ["ddg_gain", "risk"])
        pareto_items = results_df.to_dict(orient="records")
        front_items = []
        for item in pareto_items:
            dominated = False
            for other in pareto_items:
                if other is item:
                    continue
                if all(other[obj] <= item[obj] for obj in pareto_objectives) and any(
                    other[obj] < item[obj] for obj in pareto_objectives
                ):
                    dominated = True
                    break
            if not dominated:
                front_items.append(item)

        pareto_ids = {item["set_id"] for item in front_items}
        results_df["is_pareto"] = results_df["set_id"].isin(pareto_ids)

        candidates_path = seed_out_dir / "candidates.parquet"
        pareto_path = seed_out_dir / "pareto.parquet"
        summary_path = seed_out_dir / "summary.json"

        results_df.to_parquet(candidates_path, index=False)

        # Get Pareto dataframe
        pareto_df = results_df[results_df["is_pareto"]].copy()

        # Apply functional scoring if requested
        if args.functional_scoring:
            logger.info("Applying functional scoring to %d Pareto rescues for %s", len(pareto_df), seed)
            try:
                from src.scoring.functional.functional_score_evoef2 import score_rescue_candidates_evoef2
                from src.scoring.functional.evoef2_binding import load_evoef2_config

                # Load EvoEF2 config
                evoef2_cfg = load_evoef2_config(args.config_score)

                # Score Pareto rescues with functional scoring
                pareto_scored = score_rescue_candidates_evoef2(
                    pareto_df,
                    evoef2_cfg=evoef2_cfg,
                    cache_dir=str(paths["cache"] / "functional_scoring"),
                    verbose=True
                )

                # Use scored version as the Pareto dataframe
                pareto_df = pareto_scored
                logger.info("Functional scoring complete for %s", seed)

            except Exception as e:
                logger.error("Functional scoring failed for %s: %s", seed, e)
                logger.warning("Continuing without functional scores")

        # Save Pareto dataframe (with or without functional scores)
        pareto_df.to_parquet(pareto_path, index=False)

        summary = {
            "target": seed,
            "seed_ddg": seed_ddg,
            "candidates": int(len(results_df)),
            "pareto": int(results_df["is_pareto"].sum()),
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

        signature = {
            "target": seed,
            "sequence_sha256": file_sha256(fasta_path),
            "burial_sha256": file_sha256(burial_path),
            "dist_sha256": file_sha256(dist_path),
            "base_pdb_sha256": file_sha256(base_pdb),
            "evoef2_cfg": evoef2_cfg,
            "design_cfg": design_cfg,
            "beam_cfg": beam_cfg,
        }
        signature_path = seed_out_dir / "signature.json"
        signature_path.write_text(
            json.dumps({"signature_sha256": dict_sha256(signature), "inputs": signature}, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        # Clean up checkpoints on successful completion
        if checkpointing_enabled and checkpoint_dir.exists():
            delete_checkpoint(checkpoint_dir)
            logger.info("Cleaned up checkpoints for %s", seed)

        logger.info("Rescue design complete for %s (%d candidates)", seed, len(results_df))

    return 0


def main():
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Run beam search design for p53 rescue")
    parser.add_argument("--targets", type=str, nargs="+", help="Target mutations (e.g. R175H)")
    parser.add_argument("--max_muts", type=int, help="Max rescue mutations (1-3)")
    parser.add_argument("--beam_width", type=int, help="Beam width for search")
    parser.add_argument("--depth", type=int, help="Search depth")
    parser.add_argument("--parallel", type=int, help="Number of parallel jobs")
    parser.add_argument("--recompute", action="store_true", help="Recompute from scratch")
    parser.add_argument("--config_p53", type=Path, default="configs/p53.yaml")
    parser.add_argument("--config_opt", type=Path, default="configs/optimizer.yaml")
    parser.add_argument("--config_score", type=Path, default="configs/scoring.yaml")
    parser.add_argument("--config_paths", type=Path, default="configs/paths.yaml")
    parser.add_argument("--functional-scoring", action="store_true",
                        help="Apply EvoEF2-based functional scoring (DNA binding + interface) to Pareto rescues")

    args = parser.parse_args()

    # Load All Configs
    configs = {}
    for key, path in [
        ("p53", args.config_p53),
        ("optimizer", args.config_opt),
        ("scoring", args.config_score),
        ("paths", args.config_paths),
    ]:
        if path.exists():
            with path.open("r") as f:
                configs[key] = yaml.safe_load(f)
        else:
            configs[key] = {}

    return run(args, configs)


if __name__ == "__main__":
    sys.exit(main())
