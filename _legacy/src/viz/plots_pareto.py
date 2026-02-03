from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_pareto_front(
    results: pd.DataFrame,
    out_path: Path,
    x: str = "risk",
    y: str = "ddg_gain",
    pareto_flag: str = "is_pareto",
) -> None:
    if x not in results.columns or y not in results.columns:
        raise ValueError(f"Results missing required columns: {x}, {y}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    if pareto_flag in results.columns:
        front = results[results[pareto_flag] == True]
        rest = results[results[pareto_flag] == False]
        if not rest.empty:
            ax.scatter(rest[x], rest[y], s=18, alpha=0.4, label="Candidates", color="#457b9d")
        if not front.empty:
            ax.scatter(front[x], front[y], s=28, alpha=0.9, label="Pareto front", color="#e63946")
    else:
        ax.scatter(results[x], results[y], s=20, alpha=0.6, color="#457b9d")

    ax.set_xlabel("Risk score (lower is safer)")
    ax.set_ylabel("DDG gain vs seed (negative is better)")
    ax.set_title("Rescue candidates: stability vs risk")
    ax.axhline(0.0, color="#888888", linewidth=0.8, linestyle="--")
    ax.margins(x=0.07, y=0.1)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
