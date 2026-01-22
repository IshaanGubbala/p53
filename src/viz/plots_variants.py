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

    plt.figure(figsize=(8, 5))
    plt.hist(benign, bins=bins, alpha=0.6, label=f"Benign (n={len(benign)})", color="#2a9d8f")
    plt.hist(pathogenic, bins=bins, alpha=0.6, label=f"Pathogenic (n={len(pathogenic)})", color="#e76f51")
    plt.xlabel("Delta delta G (EvoEF2)")
    plt.ylabel("Count")
    plt.title("TP53 missense stability by ClinVar label")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
