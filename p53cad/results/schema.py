from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations
import json
from typing import Any, Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd

SCHEMA_VERSION = "1.0.0"

BIG8_HOTSPOTS = [
    "R175H",
    "R248Q",
    "R248W",
    "R273H",
    "R273C",
    "G245S",
    "R249S",
    "R282W",
]

DEFAULT_DELIVERY_METHODS = ["gene_therapy", "mrna_therapy", "protein_therapy"]


@dataclass(frozen=True)
class ScenarioSpec:
    scenario_id: str
    target_label: str
    targets: tuple[str, ...]
    delivery_method: str


def build_run_id(prefix: str = "campaign") -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"{prefix}_{ts}"


def _normalize_mutation_label(muts: Sequence[str]) -> str:
    return "+".join(sorted(str(m).strip().upper() for m in muts if str(m).strip()))


def _normalize_mutation_token(value: Any) -> str:
    token = str(value).strip().upper()
    return "".join(ch for ch in token if not ch.isspace())


def _parse_mutation_list(value: Any, *, plus_separator: bool = False) -> List[str]:
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []

    raw_items: List[Any]
    if isinstance(value, (list, tuple, set)):
        raw_items = list(value)
    else:
        text = str(value).strip()
        if not text:
            return []
        try:
            parsed = json.loads(text)
            if isinstance(parsed, (list, tuple, set)):
                raw_items = list(parsed)
            else:
                raw_items = [text]
        except Exception:
            split_char = "+" if plus_separator and ("+" in text) else ","
            text_norm = text.replace(";", ",") if split_char == "," else text
            raw_items = [part for part in text_norm.split(split_char)]

    out = []
    for item in raw_items:
        token = _normalize_mutation_token(item)
        if token:
            out.append(token)
    return out


def _targets_for_row(row: pd.Series) -> List[str]:
    for key in ("targets_json", "targets"):
        if key in row and pd.notna(row.get(key)):
            vals = _parse_mutation_list(row.get(key))
            if vals:
                return vals

    label = str(row.get("target_label", "")).strip().upper()
    if label:
        vals = _parse_mutation_list(label, plus_separator=True)
        if vals:
            return vals

    scenario_id = str(row.get("scenario_id", ""))
    if "__" in scenario_id:
        target_part = scenario_id.split("__", 1)[0]
        vals = _parse_mutation_list(target_part, plus_separator=True)
        if vals:
            return vals

    return []


def _mutations_for_row(row: pd.Series) -> List[str]:
    for key in ("mutations_json", "mutations"):
        if key in row and pd.notna(row.get(key)):
            vals = _parse_mutation_list(row.get(key))
            if vals or str(row.get(key)).strip() == "[]":
                return vals
    return []


def _jaccard_overlap(set_a: set[str], set_b: set[str]) -> float:
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 0.0
    return float(len(set_a & set_b) / len(union))


def build_scenario_matrix(
    hotspots: Sequence[str] | None = None,
    delivery_methods: Sequence[str] | None = None,
    include_pairs: bool = True,
) -> List[ScenarioSpec]:
    base = [str(m).strip().upper() for m in (hotspots or BIG8_HOTSPOTS) if str(m).strip()]
    modes = [str(m).strip().lower() for m in (delivery_methods or DEFAULT_DELIVERY_METHODS) if str(m).strip()]

    singles = [(m,) for m in base]
    pairs = list(combinations(base, 2)) if include_pairs else []
    targets = singles + pairs

    scenarios: List[ScenarioSpec] = []
    for mut_group in targets:
        target_label = _normalize_mutation_label(mut_group)
        for delivery in modes:
            scenario_id = f"{target_label}__{delivery}"
            scenarios.append(
                ScenarioSpec(
                    scenario_id=scenario_id,
                    target_label=target_label,
                    targets=tuple(mut_group),
                    delivery_method=delivery,
                )
            )
    return scenarios


def scenarios_to_frame(scenarios: Iterable[ScenarioSpec]) -> pd.DataFrame:
    rows = []
    for s in scenarios:
        rows.append(
            {
                "scenario_id": s.scenario_id,
                "target_label": s.target_label,
                "targets": list(s.targets),
                "delivery_method": s.delivery_method,
                "n_targets": len(s.targets),
            }
        )
    return pd.DataFrame(rows)


def select_presentation_shortlist(
    candidates_df: pd.DataFrame,
    top_n: int = 30,
    max_per_target: int = 4,
    delivery_methods: Sequence[str] | None = None,
    min_per_delivery: int = 2,
    max_per_base_mutation: int = 12,
) -> pd.DataFrame:
    """
    Rank and select presentation candidates with strict rescue semantics.

    Strict policy:
      - retain all target mutation(s)
      - include at least one additional second-site mutation
      - dedupe by sequence and canonical mutation set

    Ranking base:
      - trust-adjusted score (desc)
      - clinical_score (desc)
      - uncertainty (asc)
      - ood_distance (asc)
      - n_mutations (asc)

    Diversity stage:
      - quota seeds across delivery, profile, and target hotspots
      - delivery floor ensures min_per_delivery per method
      - base-mutation cap prevents any single hotspot from dominating via combos
      - greedy diversity fill using normalized score + novelty - overlap penalty
    """
    if candidates_df is None or candidates_df.empty:
        cols = list(candidates_df.columns) if candidates_df is not None else []
        for needed in ["selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason", "presentation_rank"]:
            if needed not in cols:
                cols.append(needed)
        return pd.DataFrame(columns=cols)

    work = candidates_df.copy().reset_index(drop=True)

    for col, default in [
        ("score", -np.inf),
        ("clinical_score", -np.inf),
        ("uncertainty", np.inf),
        ("ood_distance", np.inf),
        ("n_mutations", np.inf),
        ("target_label", "unknown"),
        ("delivery_method", "unknown"),
        ("profile", "Unknown"),
        ("candidate_uid", ""),
        ("sequence", ""),
    ]:
        if col not in work.columns:
            work[col] = default

    work["score"] = pd.to_numeric(work["score"], errors="coerce").fillna(-np.inf)
    work["clinical_score"] = pd.to_numeric(work["clinical_score"], errors="coerce").fillna(-np.inf)
    work["uncertainty"] = pd.to_numeric(work["uncertainty"], errors="coerce").fillna(np.inf)
    work["ood_distance"] = pd.to_numeric(work["ood_distance"], errors="coerce").fillna(np.inf)
    work["n_mutations"] = pd.to_numeric(work["n_mutations"], errors="coerce").fillna(np.inf)

    work["_target_list"] = work.apply(_targets_for_row, axis=1)
    work["_mutation_list"] = work.apply(_mutations_for_row, axis=1)
    work["_target_set"] = work["_target_list"].apply(lambda xs: set(xs))
    work["_mutation_set"] = work["_mutation_list"].apply(lambda xs: set(xs))
    work["_rescue_set"] = work.apply(lambda r: set(r["_mutation_set"]) - set(r["_target_set"]), axis=1)
    work["_retains_targets"] = work.apply(
        lambda r: bool(r["_target_set"]) and set(r["_target_set"]).issubset(set(r["_mutation_set"])),
        axis=1,
    )
    work["_has_second_site"] = work["_rescue_set"].apply(lambda xs: len(xs) >= 1)

    # Strict retain+add filter.
    work = work[work["_retains_targets"] & work["_has_second_site"]].copy()
    if work.empty:
        cols = list(candidates_df.columns)
        for needed in ["selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason", "presentation_rank"]:
            if needed not in cols:
                cols.append(needed)
        return pd.DataFrame(columns=cols)

    ranking_cols = ["score", "clinical_score", "uncertainty", "ood_distance", "n_mutations"]
    ranking_asc = [False, False, True, True, True]

    work = work.sort_values(by=ranking_cols, ascending=ranking_asc, kind="mergesort").reset_index(drop=True)
    work["_candidate_uid"] = work["candidate_uid"].astype(str)
    work["_sequence_key"] = work.apply(
        lambda r: str(r.get("sequence", "")).strip().upper() or f"uid:{r['_candidate_uid']}",
        axis=1,
    )
    work["_mutation_key"] = work.apply(
        lambda r: "|".join(sorted(r["_mutation_set"])) or f"uid:{r['_candidate_uid']}",
        axis=1,
    )
    work["_delivery_key"] = work["delivery_method"].astype(str).str.strip().str.lower()

    # Dedupe: keep highest-ranked row per (sequence, delivery), then per (mutation set, delivery).
    # Including delivery_method ensures each delivery method retains its own representative
    # even when the underlying protein sequence is identical (e.g. mrna vs protein therapy).
    work = work.drop_duplicates(subset=["_sequence_key", "_delivery_key"], keep="first")
    work = work.drop_duplicates(subset=["_mutation_key", "_delivery_key"], keep="first")
    work = work.sort_values(by=ranking_cols, ascending=ranking_asc, kind="mergesort").reset_index(drop=True)

    # Normalize scores to [0, 1] so diversity terms are on a comparable scale.
    finite_scores = work["score"][np.isfinite(work["score"])]
    if len(finite_scores) > 1:
        score_min = float(finite_scores.min())
        score_max = float(finite_scores.max())
    else:
        score_min, score_max = 0.0, 1.0
    score_range = max(score_max - score_min, 1e-6)
    work["_score_norm"] = ((work["score"] - score_min) / score_range).clip(0.0, 1.0)

    desired_delivery = [str(d).strip().lower() for d in (delivery_methods or DEFAULT_DELIVERY_METHODS)]
    selected_indices: List[int] = []
    selected_index_set: set[int] = set()
    stage_reason: Dict[int, str] = {}
    target_counts: Dict[str, int] = {}
    base_mut_counts: Dict[str, int] = {}

    def _can_take(idx: int) -> bool:
        if idx in selected_index_set:
            return False
        target = str(work.at[idx, "target_label"])
        if target_counts.get(target, 0) >= int(max_per_target):
            return False
        # Base-mutation cap: prevent any single hotspot from dominating via combos.
        target_set = work.at[idx, "_target_set"]
        for mut in target_set:
            if base_mut_counts.get(mut, 0) >= int(max_per_base_mutation):
                return False
        return True

    def _take(idx: int, reason: str) -> bool:
        if not _can_take(idx):
            return False
        selected_indices.append(idx)
        selected_index_set.add(idx)
        target = str(work.at[idx, "target_label"])
        target_counts[target] = target_counts.get(target, 0) + 1
        for mut in work.at[idx, "_target_set"]:
            base_mut_counts[mut] = base_mut_counts.get(mut, 0) + 1
        stage_reason[idx] = reason
        return True

    # Quota seed 1: at least one per delivery if available.
    delivery_series = work["delivery_method"].astype(str).str.lower()
    for delivery in desired_delivery:
        matches = work.index[delivery_series == delivery].tolist()
        for idx in matches:
            if _take(int(idx), f"quota_delivery_{delivery}"):
                break

    # Quota seed 2: at least one per profile if available.
    for profile in sorted(work["profile"].astype(str).unique()):
        matches = work.index[work["profile"].astype(str) == profile].tolist()
        for idx in matches:
            if _take(int(idx), f"quota_profile_{profile}"):
                break

    # Quota seed 3: target diversity — at least one per Big-8 hotspot.
    target_label_upper = work["target_label"].astype(str).str.upper()
    for hotspot in BIG8_HOTSPOTS:
        # Match any target_label containing this hotspot (single or combo).
        matches = work.index[target_label_upper.str.contains(hotspot, regex=False)].tolist()
        for idx in matches:
            if _take(int(idx), f"quota_target_{hotspot}"):
                break

    # Quota seed 4: delivery floor — ensure min_per_delivery per method.
    delivery_selected_counts: Dict[str, int] = {}
    for idx in selected_indices:
        d = str(work.at[idx, "delivery_method"]).strip().lower()
        delivery_selected_counts[d] = delivery_selected_counts.get(d, 0) + 1
    for delivery in desired_delivery:
        current = delivery_selected_counts.get(delivery, 0)
        if current >= int(min_per_delivery):
            continue
        matches = work.index[delivery_series == delivery].tolist()
        for idx in matches:
            if delivery_selected_counts.get(delivery, 0) >= int(min_per_delivery):
                break
            if _take(int(idx), f"quota_delivery_floor_{delivery}"):
                delivery_selected_counts[delivery] = delivery_selected_counts.get(delivery, 0) + 1

    novelty_weight = 0.35
    overlap_penalty = 0.20

    def _greedy_utility(idx: int, selected: Sequence[int]) -> tuple[float, float, float]:
        rescue_set = set(work.at[idx, "_rescue_set"])
        if not selected:
            novelty = 1.0
            max_overlap = 0.0
        else:
            overlaps = [_jaccard_overlap(rescue_set, set(work.at[s, "_rescue_set"])) for s in selected]
            max_overlap = float(max(overlaps)) if overlaps else 0.0
            novelty = float(1.0 - np.mean(overlaps)) if overlaps else 1.0
        utility = float(work.at[idx, "_score_norm"]) + (novelty_weight * novelty) - (overlap_penalty * max_overlap)
        return utility, novelty, max_overlap

    # Fill remaining slots with greedy diversity utility.
    while len(selected_indices) < int(top_n):
        remaining = [int(i) for i in work.index if i not in selected_index_set and _can_take(int(i))]
        if not remaining:
            break

        best_idx = -1
        best_utility = -float("inf")
        best_novelty = 0.0
        best_overlap = 0.0
        for idx in remaining:
            utility, novelty, max_overlap = _greedy_utility(idx, selected_indices)
            if utility > best_utility + 1e-12:
                best_idx = idx
                best_utility = utility
                best_novelty = novelty
                best_overlap = max_overlap
            elif abs(utility - best_utility) <= 1e-12:
                # Stable deterministic tie-break by rank order (lower index is better).
                if best_idx < 0 or idx < best_idx:
                    best_idx = idx
                    best_utility = utility
                    best_novelty = novelty
                    best_overlap = max_overlap

        if best_idx < 0:
            break
        if _take(best_idx, "diversity_fill"):
            stage_reason[best_idx] = (
                f"diversity_fill_utility={best_utility:.3f}_novelty={best_novelty:.3f}_overlap={best_overlap:.3f}"
            )

    if not selected_indices:
        cols = list(candidates_df.columns)
        for needed in ["selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason", "presentation_rank"]:
            if needed not in cols:
                cols.append(needed)
        return pd.DataFrame(columns=cols)

    # Recompute selection metrics in selected order for deterministic persisted metadata.
    selected_rows = []
    prior_selected: List[int] = []
    for rank, idx in enumerate(selected_indices[: int(top_n)], start=1):
        utility, novelty, max_overlap = _greedy_utility(idx, prior_selected)
        row = work.loc[idx].copy()
        row["selection_score"] = float(utility)
        row["diversity_novelty"] = float(novelty)
        row["max_mutation_overlap"] = float(max_overlap)
        row["selection_reason"] = str(stage_reason.get(idx, "trust_adjusted_rank"))
        row["presentation_rank"] = int(rank)
        selected_rows.append(row)
        prior_selected.append(idx)

    result = pd.DataFrame(selected_rows)
    drop_cols = [c for c in result.columns if c.startswith("_")]
    result = result.drop(columns=drop_cols, errors="ignore")

    # Keep predictable column order with new fields near the end.
    base_cols = [c for c in candidates_df.columns if c in result.columns]
    extra_cols = [c for c in ["selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason", "presentation_rank"] if c in result.columns]
    ordered_cols = base_cols + [c for c in extra_cols if c not in base_cols]
    result = result.loc[:, ordered_cols]
    return result.reset_index(drop=True)
