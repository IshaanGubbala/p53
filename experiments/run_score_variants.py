from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm
import sys

# Add project root to path
sys.path.append(str(Path.cwd()))

from src.core.hashing import dict_sha256, file_sha256
from src.core.logging import get_logger
from src.scoring.evoef2_runner import (
    canonicalize_mutation_set,
    compute_stability,
    mutation_set_id,
    score_mutation_set,
)


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    processed_dir = Path(data_cfg.get("processed", "data/processed"))
    cache_dir = Path(paths_cfg.get("cache_dir", processed_dir / "cache"))
    project_root = Path(paths_cfg.get("project_root", Path.cwd()))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
        "processed": processed_dir,
        "cache": cache_dir,
        "project_root": project_root,
    }


def _find_latest_pdb(raw_dir: Path) -> Path:
    candidates = list(raw_dir.glob("*.pdb"))
    if not candidates:
        raise RuntimeError(f"No PDB files found in {raw_dir}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _normalize_cfg(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _normalize_cfg(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_normalize_cfg(item) for item in value]
    return value


def _load_variants(variant_path: Path) -> pd.DataFrame:
    df = pd.read_parquet(variant_path)
    
    # If we already have a 'mutation' column (like search/generative results), use it
    if "mutation" in df.columns:
        df = df.dropna(subset=["mutation"]).copy()
        df["mutation"] = df["mutation"].astype(str)
        return df

    # Otherwise, require and build from pos/ref/alt
    if not {"pos", "ref", "alt"}.issubset(df.columns):
        raise RuntimeError("Variants missing required columns: 'mutation' OR 'pos/ref/alt'")
    
    df = df.dropna(subset=["pos", "ref", "alt"]).copy()
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    df["mutation"] = df["ref"] + df["pos"].astype(str) + df["alt"]
    return df


def _resolve_repaired_pdb(evoef2_cfg: dict[str, Any], project_root: Path) -> Path | None:
    repaired = evoef2_cfg.get("repaired_pdb")
    if not repaired:
        return None
    repaired_path = Path(repaired).expanduser()
    if not repaired_path.is_absolute():
        repaired_path = (project_root / repaired_path).resolve()
    evoef2_cfg["repaired_pdb"] = str(repaired_path)
    return repaired_path


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))

    scoring_cfg = configs.get("scoring", {})
    engine = (getattr(args, "engine", None) or scoring_cfg.get("engine") or "evoef2").lower()
    if engine != "evoef2":
        raise RuntimeError(f"Unsupported scoring engine: {engine}")

    evoef2_cfg = dict(scoring_cfg.get("evoef2", {}))
    if not evoef2_cfg:
        raise RuntimeError("Missing evoef2 configuration in configs/scoring.yaml")

    parallel_cfg = scoring_cfg.get("parallel", {})
    parallel = getattr(args, "parallel", None)
    if parallel is None:
        parallel = int(parallel_cfg.get("n_jobs", 1))
    parallel = max(1, int(parallel))

    cache_cfg = scoring_cfg.get("cache", {})
    cache_enabled = bool(cache_cfg.get("enabled", True))
    cache_dir = Path(cache_cfg.get("dir", paths["cache"]))
    recompute = bool(scoring_cfg.get("recompute", False))
    if not cache_enabled:
        recompute = True

    variant_path_str = getattr(args, "variants", None)
    if variant_path_str:
        variant_path = Path(variant_path_str)
    else:
        variant_path = paths["interim"] / "tp53_missense_normalized.parquet"
        
    if not variant_path.exists():
        raise RuntimeError(f"Missing variants path: {variant_path}")

    df = _load_variants(variant_path)
    mutations = sorted(df["mutation"].unique())
    
    # [RaSP Integration]
    rasp_cfg = scoring_cfg.get("rasp", {})
    rasp_enabled = rasp_cfg.get("enabled", False)
    
    logger.info("Scoring %d unique mutations with %s %s", 
                len(mutations), engine, "(+ RaSP)" if rasp_enabled else "")

    repaired_pdb = _resolve_repaired_pdb(evoef2_cfg, paths["project_root"])
    if repaired_pdb:
        base_pdb = repaired_pdb
    else:
        base_pdb = _find_latest_pdb(paths["raw"] / "alphafold")

    if not base_pdb.exists():
        raise RuntimeError(f"Base PDB not found: {base_pdb}")

    scores_path = paths["processed"] / "variant_scores.parquet"
    if variant_path_str:
        scores_path = variant_path.with_name(f"{variant_path.stem}_scored.parquet")
    
    signature_path = scores_path.with_suffix(".signature.json")

    signature = {
        "variants_sha256": file_sha256(variant_path),
        "pdb_sha256": file_sha256(base_pdb),
        "engine": engine,
        "engine_cfg": _normalize_cfg(evoef2_cfg),
    }
    signature_hash = dict_sha256(signature)
    if not recompute and scores_path.exists() and signature_path.exists():
        existing = json.loads(signature_path.read_text(encoding="utf-8"))
        if existing.get("signature_sha256") == signature_hash:
            logger.info("Variant scores up to date; skipping scoring")
            return 0

    work_root = cache_dir / "evoef2_work"
    work_root.mkdir(parents=True, exist_ok=True)

    base_hash = file_sha256(base_pdb)
    base_energy_cache_key = dict_sha256(
        {"pdb_sha256": base_hash, "engine_cfg": _normalize_cfg(evoef2_cfg)}
    )
    base_energy_cache = cache_dir / f"evoef2_base_{base_energy_cache_key}.json"
    base_energy = None
    if base_energy_cache.exists() and not recompute:
        payload = json.loads(base_energy_cache.read_text(encoding="utf-8"))
        base_energy = float(payload["base_energy"])

    if base_energy is None:
        logger.info("Computing EvoEF2 base stability energy")
        base_workdir = work_root / "base_energy"
        base_workdir.mkdir(parents=True, exist_ok=True)
        with tqdm(total=1, desc="Computing base stability", unit="step") as pbar:
            base_energy = compute_stability(base_pdb, evoef2_cfg, base_workdir)
            pbar.update(1)
        base_energy_cache.write_text(
            json.dumps(
                {
                    "base_energy": base_energy,
                    "pdb_sha256": base_hash,
                    "engine_cfg": _normalize_cfg(evoef2_cfg),
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )

    scores: dict[str, float] = {}
    pending: list[str] = []

    for mut in mutations:
        # Support multi-mutation strings like "M1A,V2L"
        mut_list = [m.strip() for m in mut.split(",")]
        canonical = canonicalize_mutation_set(mut_list)
        key = mutation_set_id(base_hash, canonical, evoef2_cfg)
        cache_path = cache_dir / f"evoef2_{key}.json"
        if cache_path.exists() and not recompute:
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            if "ddg" in payload:
                scores[mut] = {"ddg": float(payload["ddg"])}
                # Try to get RaSP too if enabled
                if rasp_enabled:
                    rasp_cache = cache_dir / f"rasp_full_{key}.json"
                    if rasp_cache.exists():
                        r_payload = json.loads(rasp_cache.read_text(encoding="utf-8"))
                        scores[mut]["rasp_ddg"] = float(r_payload.get("rasp_ddg", 0))
                continue
        pending.append(mut)

    errors: list[str] = []
    pbar = tqdm(total=len(mutations), desc="Scoring mutations", unit="mut")
    pbar.update(len(scores))

    def _score(mut: str) -> tuple[str, dict[str, float], float]:
        start = time.perf_counter()
        mut_list = [m.strip() for m in mut.split(",")]
        res = {"ddg": score_mutation_set(
            mut_list,
            base_pdb,
            cache_dir,
            evoef2_cfg,
            work_root,
            base_energy=base_energy,
            recompute=recompute,
        )}
        
        if rasp_enabled:
            # Import here to avoid forcing dependency if not used
            from src.scoring.rasp_scorer import score_mutation_set_rasp
            res["rasp_ddg"] = score_mutation_set_rasp(mut_list, base_pdb, cache_dir, rasp_cfg)
            
        elapsed = time.perf_counter() - start
        return mut, res, elapsed

    if pending:
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = {executor.submit(_score, mut): mut for mut in pending}
            for future in as_completed(futures):
                mut = futures[future]
                try:
                    key, results, _elapsed = future.result()
                    scores[key] = results
                except Exception as exc:
                    errors.append(f"{mut}: {exc}")
                pbar.update(1)
    pbar.close()

    if errors:
        sample = "\n".join(errors[:5])
        raise RuntimeError(f"Scoring failed for {len(errors)} mutations:\n{sample}")

    scored_df = df.copy()
    scored_df["ddg"] = scored_df["mutation"].apply(lambda x: scores.get(x, {}).get("ddg"))
    if rasp_enabled:
        scored_df["rasp_ddg"] = scored_df["mutation"].apply(lambda x: scores.get(x, {}).get("rasp_ddg"))
    
    scored_df = scored_df.dropna(subset=["ddg"])
    scored_df["engine"] = engine

    scores_path.parent.mkdir(parents=True, exist_ok=True)
    scored_df.to_parquet(scores_path, index=False)
    logger.info("Saved variant scores to %s", scores_path)

    signature_path.write_text(
        json.dumps({"signature_sha256": signature_hash, "inputs": signature}, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    meta_path = paths["processed"] / "variant_scores_meta.json"
    meta_path.write_text(
        json.dumps(
            {
                "engine": engine,
                "base_pdb": str(base_pdb),
                "base_energy": base_energy,
                "variants": len(df),
                "unique_mutations": len(mutations),
                "scored_mutations": len(scores),
                "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    logger.info("Variant scoring complete")
    return 0


if __name__ == "__main__":
    import argparse
    import yaml

    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", help="Path to variants parquet/csv file")
    parser.add_argument("--engine", help="Scoring engine to use")
    parser.add_argument("--parallel", type=int, help="Number of parallel jobs")
    parser.add_argument("--config", default="configs/scoring.yaml", help="Path to scoring config")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        scoring_cfg = yaml.safe_load(f)
    
    # Load default paths
    with open("configs/paths.yaml", "r") as f:
        paths_cfg = yaml.safe_load(f)
        
    configs = {"scoring": scoring_cfg, "paths": paths_cfg}
    sys.exit(run(args, configs))
