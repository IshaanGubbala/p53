from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import griddata


def _parse_mutation_position(mut: str) -> int | None:
    text = str(mut).strip().upper()
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        return None
    try:
        return int(digits)
    except ValueError:
        return None



def build_candidate_position_heatmap(
    candidates: Sequence[Dict[str, Any]],
    max_positions: int = 40,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      - matrix_df: rows candidates, cols informative positions, values 0/1 mutation presence
      - freq_df: position-level mutation frequency summary
    """
    rows: List[Dict[str, Any]] = []
    position_counts: Dict[int, int] = {}

    for cand in candidates:
        mut_positions = cand.get("mut_positions") or []
        if not mut_positions:
            mut_positions = [
                p for p in (_parse_mutation_position(m) for m in (cand.get("mutations") or [])) if p is not None
            ]
        mut_set = {int(p) for p in mut_positions}
        for pos in mut_set:
            position_counts[pos] = position_counts.get(pos, 0) + 1
        rows.append(
            {
                "candidate_label": f"#{cand.get('candidate_id', '?')} {cand.get('profile', 'Unknown')}",
                "score": float(cand.get("score", 0.0)),
                "positions": mut_set,
            }
        )

    if not rows:
        return pd.DataFrame(), pd.DataFrame(columns=["position", "count", "frequency"]) 

    n_candidates = len(rows)
    freq_df = pd.DataFrame(
        [
            {
                "position": int(pos),
                "count": int(count),
                "frequency": float(count / max(n_candidates, 1)),
            }
            for pos, count in sorted(position_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        ]
    )

    if freq_df.empty:
        return pd.DataFrame(), freq_df

    selected_positions = freq_df["position"].tolist()[: max(1, int(max_positions))]
    matrix = []
    for row in rows:
        matrix.append(
            {
                "candidate_label": row["candidate_label"],
                **{f"pos_{pos}": (1 if pos in row["positions"] else 0) for pos in selected_positions},
            }
        )
    matrix_df = pd.DataFrame(matrix)
    return matrix_df, freq_df



def build_multi_metric_frame(
    traj_df: pd.DataFrame,
    *,
    normalize: bool,
    min_function: float,
    min_stability: float,
    min_binding: float,
    min_identity: float,
) -> pd.DataFrame:
    out = traj_df.copy()
    if out.empty:
        return out

    for col in ["score", "stability", "binding", "identity"]:
        if col not in out.columns:
            out[col] = np.nan
        out[col] = pd.to_numeric(out[col], errors="coerce")

    if not normalize:
        return out

    # Threshold-anchored normalization (not per-trace min/max)
    score_hi = max(float(min_function) + 1.0, float(np.nanpercentile(out["score"], 95)) if out["score"].notna().any() else float(min_function) + 1.0)
    out["score"] = np.clip((out["score"] - float(min_function)) / max(score_hi - float(min_function), 1e-6), 0.0, 1.0)

    out["stability"] = np.clip((out["stability"] - float(min_stability)) / max(0.0 - float(min_stability), 1e-6), 0.0, 1.0)
    out["binding"] = np.clip((out["binding"] - float(min_binding)) / max(10.0 - float(min_binding), 1e-6), 0.0, 1.0)
    out["identity"] = np.clip((out["identity"] - float(min_identity)) / max(100.0 - float(min_identity), 1e-6), 0.0, 1.0)
    return out



def build_loss_mesh(
    traj_df: pd.DataFrame,
    grid_size: int = 35,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """
    Build interpolated loss mesh from trajectory points.

    Returns (grid_x, grid_y, grid_z) or None if insufficient points.
    """
    needed = {"lx", "ly", "loss_total"}
    if traj_df is None or traj_df.empty or not needed.issubset(set(traj_df.columns)):
        return None

    work = traj_df[list(needed)].copy()
    work = work.apply(pd.to_numeric, errors="coerce").dropna()
    if len(work) < 12:
        return None

    # Need enough unique points for interpolation
    if work[["lx", "ly"]].drop_duplicates().shape[0] < 8:
        return None

    x = work["lx"].to_numpy()
    y = work["ly"].to_numpy()
    z = work["loss_total"].to_numpy()

    gx = np.linspace(float(np.min(x)), float(np.max(x)), int(max(grid_size, 12)))
    gy = np.linspace(float(np.min(y)), float(np.max(y)), int(max(grid_size, 12)))
    grid_x, grid_y = np.meshgrid(gx, gy)

    grid_z = griddata((x, y), z, (grid_x, grid_y), method="linear")
    if np.isnan(grid_z).all():
        return None

    # Fill sparse holes with nearest so surface is viewable.
    mask = np.isnan(grid_z)
    if np.any(mask):
        nearest = griddata((x, y), z, (grid_x, grid_y), method="nearest")
        grid_z[mask] = nearest[mask]

    return grid_x, grid_y, grid_z
