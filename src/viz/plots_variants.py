from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_ddg_by_label(df: pd.DataFrame, out_path: Path, bins: int = 50) -> None:
    if not {"ddg", "label"}.issubset(df.columns):
        raise ValueError("Dataframe must include ddg and label columns")

    benign = df.loc[df["label"] == 0, "ddg"].astype(float)
    pathogenic = df.loc[df["label"] == 1, "ddg"].astype(float)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(benign, bins=bins, alpha=0.6, label=f"Benign (n={len(benign)})", color="#2a9d8f")
    ax.hist(pathogenic, bins=bins, alpha=0.6, label=f"Pathogenic (n={len(pathogenic)})", color="#e76f51")
    ax.set_xlabel("Delta delta G (EvoEF2)")
    ax.set_ylabel("Count")
    ax.set_title("TP53 missense stability by ClinVar label")
    ax.margins(x=0.02, y=0.08)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
