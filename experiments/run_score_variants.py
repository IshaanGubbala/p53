from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

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
    if not {"pos", "ref", "alt"}.issubset(df.columns):
        raise RuntimeError("Normalized variants missing required columns: pos/ref/alt")
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

    variant_path = paths["interim"] / "tp53_missense_normalized.parquet"
    if not variant_path.exists():
        raise RuntimeError(f"Missing normalized variants: {variant_path}")

    df = _load_variants(variant_path)
    mutations = sorted(df["mutation"].unique())
    logger.info("Scoring %d unique mutations with EvoEF2", len(mutations))

    repaired_pdb = _resolve_repaired_pdb(evoef2_cfg, paths["project_root"])
    if repaired_pdb:
        base_pdb = repaired_pdb
    else:
        base_pdb = _find_latest_pdb(paths["raw"] / "alphafold")

    if not base_pdb.exists():
        raise RuntimeError(f"Base PDB not found: {base_pdb}")

    scores_path = paths["processed"] / "variant_scores.parquet"
    signature_path = paths["processed"] / "variant_scores_signature.json"

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
        canonical = canonicalize_mutation_set([mut])
        key = mutation_set_id(base_hash, canonical, evoef2_cfg)
        cache_path = cache_dir / f"evoef2_{key}.json"
        if cache_path.exists() and not recompute:
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            if "ddg" in payload:
                scores[mut] = float(payload["ddg"])
                continue
        pending.append(mut)

    errors: list[str] = []
    pbar = tqdm(total=len(mutations), desc="Scoring mutations", unit="mut")
    pbar.update(len(scores))

    def _score(mut: str) -> tuple[str, float, float]:
        start = time.perf_counter()
        ddg = score_mutation_set(
            [mut],
            base_pdb,
            cache_dir,
            evoef2_cfg,
            work_root,
            base_energy=base_energy,
            recompute=recompute,
        )
        elapsed = time.perf_counter() - start
        return mut, ddg, elapsed

    if pending:
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = {executor.submit(_score, mut): mut for mut in pending}
            for future in as_completed(futures):
                mut = futures[future]
                try:
                    key, ddg, _elapsed = future.result()
                    scores[key] = ddg
                except Exception as exc:
                    errors.append(f"{mut}: {exc}")
                pbar.update(1)
    pbar.close()

    if errors:
        sample = "\n".join(errors[:5])
        raise RuntimeError(f"EvoEF2 scoring failed for {len(errors)} mutations:\n{sample}")

    scored_df = df.copy()
    scored_df["ddg"] = scored_df["mutation"].map(scores)
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
