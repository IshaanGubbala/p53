from __future__ import annotations

import json
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

PROFILE_COLORS = {
    "Balanced": "#2563EB",
    "Stability-First": "#0EA5E9",
    "Binding-Optimized": "#10B981",
    "Function-Maximized": "#14B8A6",
    "Conservative": "#1D4ED8",
    "Experimental": "#0891B2",
}


def _to_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                return parsed
            return []
        except json.JSONDecodeError:
            return []
    return []


def _to_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    return default


def _normalize_trajectory_frame(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    for col in [
        "step",
        "score",
        "score_raw",
        "score_calibrated",
        "stability",
        "binding",
        "identity",
        "n_mutations",
        "uncertainty",
        "ood_distance",
        "lx",
        "ly",
        "lz",
        "loss_total",
    ]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    if "step" in out.columns:
        out = out.sort_values("step")
    return out


def _mutation_set_from_row(row: pd.Series) -> set[str]:
    muts = [
        str(m).strip().upper()
        for m in _to_list(row.get("mutations_json", row.get("mutations", [])))
        if str(m).strip()
    ]
    return set(muts)


def _jaccard_overlap(set_a: set[str], set_b: set[str]) -> float:
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 0.0
    return float(len(set_a & set_b) / len(union))


def hydrate_display_candidates(
    candidates_df: pd.DataFrame,
    trajectories_df: Optional[pd.DataFrame] = None,
    top_df: Optional[pd.DataFrame] = None,
    max_candidates: int = 30,
) -> List[Dict[str, Any]]:
    """
    Convert campaign artifact tables into Streamlit candidate objects.
    """
    if candidates_df is None or candidates_df.empty:
        return []

    candidates = candidates_df.copy()
    for col in [
        "score",
        "score_raw",
        "score_calibrated",
        "stability",
        "binding",
        "identity",
        "n_mutations",
        "uncertainty",
        "ood_distance",
        "score_gain_vs_target",
        "selection_score",
        "diversity_novelty",
        "max_mutation_overlap",
    ]:
        if col in candidates.columns:
            candidates[col] = pd.to_numeric(candidates[col], errors="coerce")

    if top_df is not None and not top_df.empty and "candidate_uid" in top_df.columns:
        ranked = top_df.copy()
        if "presentation_rank" in ranked.columns:
            ranked["presentation_rank"] = pd.to_numeric(ranked["presentation_rank"], errors="coerce")
            ranked = ranked.sort_values("presentation_rank")
        ranked_uids = ranked["candidate_uid"].astype(str).tolist()
        candidates = candidates[candidates["candidate_uid"].astype(str).isin(ranked_uids)].copy()
        candidates["_rank_order"] = pd.Categorical(candidates["candidate_uid"].astype(str), categories=ranked_uids, ordered=True)
        candidates = candidates.sort_values("_rank_order")
    else:
        sort_cols = [c for c in ["selection_score", "score", "clinical_score"] if c in candidates.columns]
        if sort_cols:
            asc = [False for _ in sort_cols]
            candidates = candidates.sort_values(sort_cols, ascending=asc)

    if "pass_name" in candidates.columns:
        deep_only = candidates[candidates["pass_name"].astype(str).str.lower() == "deep"]
        if not deep_only.empty:
            candidates = deep_only

    # Backfill missing selection metadata for older artifacts deterministically.
    for col in ["selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason"]:
        if col not in candidates.columns:
            candidates[col] = np.nan if col != "selection_reason" else ""

    prior_sets: List[set[str]] = []
    for idx in candidates.index:
        row = candidates.loc[idx]
        mut_set = _mutation_set_from_row(row)
        if prior_sets:
            overlaps = [_jaccard_overlap(mut_set, prev) for prev in prior_sets]
            novelty = float(1.0 - np.mean(overlaps))
            max_overlap = float(max(overlaps))
        else:
            novelty = 1.0
            max_overlap = 0.0

        if pd.isna(candidates.at[idx, "diversity_novelty"]):
            candidates.at[idx, "diversity_novelty"] = novelty
        if pd.isna(candidates.at[idx, "max_mutation_overlap"]):
            candidates.at[idx, "max_mutation_overlap"] = max_overlap
        if pd.isna(candidates.at[idx, "selection_score"]):
            base_score = float(pd.to_numeric(row.get("score", np.nan), errors="coerce"))
            if np.isnan(base_score):
                base_score = -np.inf
            candidates.at[idx, "selection_score"] = base_score + (0.35 * novelty) - (0.20 * max_overlap)
        reason = str(candidates.at[idx, "selection_reason"]).strip()
        if not reason or reason.lower() == "nan":
            candidates.at[idx, "selection_reason"] = f"artifact_load_fallback_rank_{int(len(prior_sets) + 1)}"

        prior_sets.append(mut_set)

    trajectories = _normalize_trajectory_frame(
        trajectories_df if trajectories_df is not None else pd.DataFrame()
    )
    traj_by_uid: Dict[str, List[Dict[str, Any]]] = {}
    if not trajectories.empty and "candidate_uid" in trajectories.columns:
        for uid, sub in trajectories.groupby(trajectories["candidate_uid"].astype(str)):
            rows: List[Dict[str, Any]] = []
            for _, row in sub.iterrows():
                payload = row.to_dict()
                payload["mutations"] = _to_list(payload.get("mutations_json", payload.get("mutations", [])))
                payload["mut_positions"] = [
                    int(p)
                    for p in _to_list(payload.get("mut_positions_json", payload.get("mut_positions", [])))
                    if str(p).strip()
                ]
                rows.append(payload)
            traj_by_uid[str(uid)] = rows

    display: List[Dict[str, Any]] = []
    for idx, (_, row) in enumerate(candidates.head(int(max_candidates)).iterrows(), start=1):
        uid = str(row.get("candidate_uid", f"cand_{idx}"))
        profile = str(row.get("profile", "Balanced"))
        muts = [str(m).strip().upper() for m in _to_list(row.get("mutations_json", row.get("mutations", []))) if str(m).strip()]
        mut_positions = [
            int(p)
            for p in _to_list(row.get("mut_positions_json", row.get("mut_positions", [])))
            if str(p).strip()
        ]
        trajectory = traj_by_uid.get(uid, [])

        display.append(
            {
                "candidate_id": idx,
                "candidate_uid": uid,
                "profile": profile,
                "color": PROFILE_COLORS.get(profile, "#2563EB"),
                "restart": int(row.get("restart_idx", row.get("restart", 1)) or 1),
                "trial_idx": int(row.get("trial_idx", idx) or idx),
                "sequence": str(row.get("sequence", "")),
                "score": float(row.get("score", np.nan)),
                "score_raw": float(row.get("score_raw", row.get("score", np.nan))),
                "score_calibrated": float(row.get("score_calibrated", row.get("score", np.nan))),
                "score_gain_vs_target": float(row.get("score_gain_vs_target", np.nan)),
                "stability": float(row.get("stability", np.nan)),
                "binding": float(row.get("binding", np.nan)),
                "identity": float(row.get("identity", np.nan)),
                "n_mutations": int(row.get("n_mutations", len(muts)) or len(muts)),
                "mutations": muts,
                "mut_positions": mut_positions,
                "uncertainty": float(row.get("uncertainty", np.nan)),
                "ood_distance": float(row.get("ood_distance", np.nan)),
                "selection_score": float(row.get("selection_score", np.nan)),
                "diversity_novelty": float(row.get("diversity_novelty", np.nan)),
                "max_mutation_overlap": float(row.get("max_mutation_overlap", np.nan)),
                "selection_reason": str(row.get("selection_reason", "artifact_load")),
                "delivery_method": str(row.get("delivery_method", "unknown")),
                "target_label": str(row.get("target_label", "unknown")),
                "trajectory": trajectory,
                "meets_constraints": _to_bool(row.get("meets_constraints", True), default=True),
            }
        )

    return display


def bundle_tables(bundle: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    return {
        "manifest": bundle.get("manifest", {}),
        "scenarios": bundle.get("scenarios", pd.DataFrame()),
        "candidates": bundle.get("candidates", pd.DataFrame()),
        "trajectories": bundle.get("trajectories", pd.DataFrame()),
        "clinical": bundle.get("clinical", pd.DataFrame()),
        "top30": bundle.get("top30", pd.DataFrame()),
    }
