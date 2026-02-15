#!/usr/bin/env python3
"""Generate a self-contained Colab notebook from the p53cad source tree.

Reads all necessary package source files and embeds them as %%writefile
cells so the notebook can run on Google Colab with just the data files
uploaded alongside it.
"""

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ── Notebook cell helpers ────────────────────────────────────────────────

cells = []


def md(source: str):
    """Add a markdown cell."""
    lines = source.split("\n")
    src = [line + "\n" for line in lines[:-1]]
    if lines[-1]:
        src.append(lines[-1])
    cells.append({"cell_type": "markdown", "metadata": {}, "source": src})


def code(source: str):
    """Add a code cell."""
    lines = source.split("\n")
    src = [line + "\n" for line in lines[:-1]]
    if lines[-1]:
        src.append(lines[-1])
    cells.append(
        {
            "cell_type": "code",
            "metadata": {},
            "source": src,
            "outputs": [],
            "execution_count": None,
        }
    )


def writefile_cell(rel_path: str, content: str):
    """Add a %%writefile cell that writes a source file to disk."""
    full = f"%%writefile {rel_path}\n{content}"
    code(full)


def read_src(rel_path: str) -> str:
    """Read a source file from the project tree."""
    return (PROJECT_ROOT / rel_path).read_text(encoding="utf-8")


# ── Build the notebook ───────────────────────────────────────────────────

# ── Title ──
md(
    """\
# p53 Rescue Mutation Discovery — Full Simulation (Google Colab)

This notebook runs the **p53-proteoMgCAD** computational protein design pipeline to discover
second-site rescue mutations that restore tumor suppressor function in mutant p53.

**Fully self-contained** — upload this notebook (plus data files) to Google Colab and run.

**Pipeline overview:**
1. Install dependencies & write package source code
2. Upload data files (DMS CSV + oracle weights)
3. Load models (ESM-2 protein language model + functional oracle)
4. Build scenario matrix (8 cancer hotspots × 3 delivery methods)
5. Run campaign (Pass A screening → Pass B deep refinement)
6. Analyze & visualize results (Top-30 shortlist, clinical impact, heatmaps)

**Required data files** (upload when prompted):
- `p53_DMS_Giacomelli_2018.csv` — Deep mutational scanning data (~1.9 MB)
- `functional_oracle.pt` — Trained oracle model weights (~29 MB)
- `p53_wt.pdb` — Wild-type structure (optional, for contact penalties)"""
)

# ── Section 0: Install Dependencies ──
md("## 0. Install Dependencies & Build Package")

code(
    """\
# Install required packages (takes ~2 min on Colab)
!pip install -q torch transformers pandas numpy matplotlib seaborn scipy pyarrow

import os, sys

# Runtime guards
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

print("Dependencies installed successfully")"""
)

# ── Section 0a: Create package directory structure ──
md(
    """\
### 0a. Write Package Source Code

The following cells write the `p53cad` package source files to disk using `%%writefile`.
You can collapse these cells — they just recreate the package structure."""
)

code(
    """\
import os

# Create package directory structure
for d in [
    'p53cad', 'p53cad/core', 'p53cad/data', 'p53cad/engine',
    'p53cad/analysis', 'p53cad/results',
    'data/raw', 'data/models', 'data/campaigns',
]:
    os.makedirs(d, exist_ok=True)
    init_path = os.path.join(d, '__init__.py')
    if d.startswith('p53cad') and not os.path.exists(init_path):
        open(init_path, 'w').close()

print("Package directory structure created")"""
)

# ── Write all source files ──
source_files = [
    ("p53cad/core/logging.py", "p53cad/core/logging.py"),
    ("p53cad/core/runtime.py", "p53cad/core/runtime.py"),
    ("p53cad/data/dms.py", "p53cad/data/dms.py"),
    ("p53cad/engine/latent.py", "p53cad/engine/latent.py"),
    ("p53cad/engine/oracle.py", "p53cad/engine/oracle.py"),
    ("p53cad/engine/pareto.py", "p53cad/engine/pareto.py"),
    ("p53cad/engine/campaign.py", "p53cad/engine/campaign.py"),
    ("p53cad/analysis/clinical_impact.py", "p53cad/analysis/clinical_impact.py"),
    ("p53cad/results/schema.py", "p53cad/results/schema.py"),
    ("p53cad/results/store.py", "p53cad/results/store.py"),
    ("p53cad/results/visualization.py", "p53cad/results/visualization.py"),
]

for write_path, src_path in source_files:
    content = read_src(src_path)
    writefile_cell(write_path, content)

# ── Section 0b: Upload Data Files ──
md(
    """\
### 0b. Upload Data Files

Upload the required data files when prompted:
1. **`p53_DMS_Giacomelli_2018.csv`** — DMS functional scores (~1.9 MB)
2. **`functional_oracle.pt`** — Trained oracle weights (~29 MB)
3. **`p53_wt.pdb`** — Wild-type PDB structure (optional, ~243 KB)"""
)

code(
    """\
import shutil

# Check if running in Colab
IN_COLAB = 'google.colab' in str(get_ipython()) if 'get_ipython' in dir() else False

if IN_COLAB:
    from google.colab import files
    print("Upload data files (DMS CSV, oracle weights, optionally PDB)...")
    print("You can select multiple files at once.\\n")
    uploaded = files.upload()

    # Move files to expected locations
    for fname, content in uploaded.items():
        if fname.endswith('.csv'):
            dest = f'data/raw/{fname}'
            with open(dest, 'wb') as f:
                f.write(content)
            print(f"  -> {dest}")
        elif fname.endswith('.pt'):
            dest = f'data/models/{fname}'
            with open(dest, 'wb') as f:
                f.write(content)
            print(f"  -> {dest}")
        elif fname.endswith('.pdb'):
            dest = f'data/raw/{fname}'
            with open(dest, 'wb') as f:
                f.write(content)
            print(f"  -> {dest}")
        else:
            print(f"  Skipped unknown file: {fname}")
else:
    # Running locally — check if data files exist
    import os
    data_files = {
        'DMS CSV': 'data/raw/p53_DMS_Giacomelli_2018.csv',
        'Oracle weights': 'data/models/functional_oracle.pt',
        'WT PDB': 'data/raw/p53_wt.pdb',
    }
    for label, path in data_files.items():
        exists = os.path.exists(path)
        status = "found" if exists else "MISSING"
        print(f"  {label}: {path} [{status}]")

print("\\nData setup complete!")"""
)

# ── Section 1: Environment Setup ──
md("## 1. Environment Setup")

code(
    """\
from p53cad.core.runtime import bootstrap_runtime, get_runtime_capabilities

bootstrap_runtime(seed=42)
caps = get_runtime_capabilities()

print("Runtime Capabilities")
print("=" * 50)
for key, val in caps.items():
    print(f"  {key:<35} {val}")"""
)

# ── Section 2: Inspect Data ──
md("## 2. Inspect Data: Wild-Type p53 & DMS Dataset")

code(
    """\
from p53cad.data.dms import P53_WT, get_dms_data

print(f"p53 Wild-Type Sequence ({len(P53_WT)} amino acids):")
# Print in blocks of 60 with position markers
for i in range(0, len(P53_WT), 60):
    chunk = P53_WT[i:i+60]
    print(f"  {i+1:>4}  {chunk}")

print(f"\\nLoading DMS data (Giacomelli 2018)...")
dms_df = get_dms_data()
print(f"  Variants: {len(dms_df):,}")
print(f"  Columns: {list(dms_df.columns)}")
dms_df.head()"""
)

code(
    """\
import matplotlib.pyplot as plt
import numpy as np

# Visualize DMS score distribution
score_col = [c for c in dms_df.columns if "Nutlin" in c and "Z" in c]
if score_col:
    scores = dms_df[score_col[0]].dropna()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.hist(scores, bins=80, color="#2563EB", alpha=0.7, edgecolor="white")
    ax.axvline(0, color="red", linestyle="--", label="Neutral (Z=0)")
    ax.axvline(-0.5, color="green", linestyle="--", alpha=0.7, label="Functional threshold (Z=-0.5)")
    ax.set_xlabel("Nutlin-3 Z-score")
    ax.set_ylabel("Count")
    ax.set_title("DMS Functional Score Distribution (Giacomelli 2018)")
    ax.legend()
    plt.tight_layout()
    plt.show()
    print(f"  Mean Z-score: {scores.mean():.3f}")
    print(f"  Functional (Z<0): {(scores < 0).sum():,} / {len(scores):,} ({100*(scores < 0).mean():.1f}%)")
else:
    # Fallback for legacy format
    if 'score' in dms_df.columns:
        scores = dms_df['score'].dropna()
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        ax.hist(scores, bins=80, color="#2563EB", alpha=0.7, edgecolor="white")
        ax.axvline(0, color="red", linestyle="--", label="Neutral (score=0)")
        ax.set_xlabel("Functional Score")
        ax.set_ylabel("Count")
        ax.set_title("DMS Functional Score Distribution")
        ax.legend()
        plt.tight_layout()
        plt.show()"""
)

# ── Section 3: Scenario Matrix ──
md(
    """\
## 3. Scenario Matrix

The campaign explores rescue mutations for 8 major p53 cancer hotspots (the "BIG8"),
plus pairwise combinations, across 3 delivery methods."""
)

code(
    """\
from p53cad.results.schema import BIG8_HOTSPOTS, DEFAULT_DELIVERY_METHODS, build_scenario_matrix

print("BIG8 Cancer Hotspots:")
for i, h in enumerate(BIG8_HOTSPOTS, 1):
    print(f"  {i}. {h}")

print(f"\\nDelivery Methods: {DEFAULT_DELIVERY_METHODS}")

# Build the full scenario matrix
scenarios = build_scenario_matrix(
    hotspots=BIG8_HOTSPOTS,
    delivery_methods=DEFAULT_DELIVERY_METHODS,
    include_pairs=True,
)
print(f"\\nTotal scenarios: {len(scenarios)}")
print(f"  Single-hotspot: {sum(1 for s in scenarios if '+' not in s.target_label)}")
print(f"  Pair combos:    {sum(1 for s in scenarios if '+' in s.target_label)}")

# Show first few
print(f"\\nFirst 10 scenarios:")
for s in scenarios[:10]:
    print(f"  {s.scenario_id:<40} targets={s.targets}  delivery={s.delivery_method}")"""
)

# ── Section 4: Run the Campaign ──
md(
    """\
## 4. Run the Campaign

**Budget options:**
| Budget | Pass A Steps | Pass B Steps | Approx. Time |
|--------|-------------|-------------|---------------|
| `fast` | 40 steps × 1 restart | 120 steps × 1 restart | ~5 min (GPU) |
| `medium` | 60 steps × 1 restart | 200 steps × 2 restarts | ~30 min (GPU) |
| `high` | 80 steps × 1 restart | 280 steps × 2 restarts | ~2 hrs (GPU) |

Change `BUDGET` below to control run time. On Colab free tier (T4 GPU), `fast` is recommended."""
)

code(
    """\
# ===== CONFIGURATION =====
BUDGET = "fast"          # "fast", "medium", or "high"
SEED = 42
INCLUDE_PAIRS = True     # Include pairwise hotspot combos
SHORTLIST_N = 30         # Number of top candidates to shortlist
# ========================="""
)

code(
    """\
import time
from p53cad.engine.campaign import CampaignRunner

print(f"Initializing CampaignRunner...")
runner = CampaignRunner()

# Report loaded model info
if runner.embedder is not None:
    print(f"  ESM-2 model: {runner.embedder.model_name}")
    print(f"  Hidden dim:  {runner.embedder.hidden_size}")
if runner.oracle is not None:
    print(f"  Oracle input_dim: {runner.oracle.input_dim}")
    arch = "attention_pooling" if hasattr(runner.oracle.model, 'attn') else "legacy_mlp"
    print(f"  Oracle architecture: {arch}")
print(f"  Pairwise DMS entries: {len(runner._pairwise_dms)}")
print(f"  Contact map entries:  {len(runner._wt_contacts)}")
print(f"\\nReady to run.")"""
)

code(
    """\
print(f"Starting campaign (budget={BUDGET}, seed={SEED})...")
print(f"This may take a while depending on budget.\\n")

t0 = time.time()

result = runner.run(
    budget=BUDGET,
    seed=SEED,
    include_pairs=INCLUDE_PAIRS,
    shortlist_n=SHORTLIST_N,
    with_clinical=True,
)

elapsed = time.time() - t0

print(f"\\n{'=' * 60}")
print(f"  CAMPAIGN COMPLETE \\u2014 {elapsed/60:.1f} min")
print(f"{'=' * 60}")
print(f"  Run ID:      {result['run_id']}")
print(f"  Scenarios:   {result['n_scenarios']}")
print(f"  Candidates:  {result['n_candidates']}")
print(f"  Shortlist:   {result['n_shortlist']}")
print(f"  Run dir:     {result['run_dir']}")"""
)

# ── Section 5: Results Analysis ──
md("## 5. Results Analysis")

code(
    """\
import json
import os
import pandas as pd
import numpy as np

run_dir = result["run_dir"]

# Load all candidates
cand = pd.read_parquet(os.path.join(run_dir, "candidates.parquet"))
print(f"Total candidates: {len(cand)}")
print(f"Columns: {list(cand.columns)}")

# Separate deep-refined candidates
if "pass_name" in cand.columns:
    deep = cand[cand["pass_name"] == "deep"].copy()
    print(f"Deep-refined candidates: {len(deep)}")
else:
    deep = cand.copy()
    print(f"All candidates (no pass separation): {len(deep)}")

cand.describe()"""
)

code(
    """\
# Target retention analysis
retains = 0
for _, row in deep.iterrows():
    targets = json.loads(str(row.get("targets_json", "[]")))
    muts = json.loads(str(row.get("mutations_json", "[]")))
    if set(targets).issubset(set(muts)):
        retains += 1

print(f"Target Retention: {retains}/{len(deep)} ({100*retains/max(len(deep),1):.1f}%)")

# Pareto ranking
if "pareto_rank" in deep.columns:
    rank1 = (deep["pareto_rank"] == 1).sum()
    max_rank = int(deep["pareto_rank"].replace([np.inf], np.nan).dropna().max())
    print(f"Pareto rank-1 (non-dominated): {rank1}")
    print(f"Total Pareto fronts: {max_rank}")

# DMS quality
if "rescue_dms_mean" in deep.columns:
    valid_dms = deep[deep["rescue_dms_mean"].notna() & (deep["rescue_dms_mean"] != 0)]
    print(f"\\nCandidates with DMS data: {len(valid_dms)}")
    if len(valid_dms):
        print(f"  Mean rescue DMS Z: {valid_dms['rescue_dms_mean'].mean():.3f}")
        func = (valid_dms["rescue_dms_mean"] < 0).sum()
        print(f"  Functional rescues (Z<0): {func}/{len(valid_dms)} ({100*func/len(valid_dms):.0f}%)")

# Score stats
if "score" in deep.columns:
    print(f"\\nOracle Score Stats:")
    print(f"  Mean: {deep['score'].mean():.4f}")
    print(f"  Max:  {deep['score'].max():.4f}")
    print(f"  Min:  {deep['score'].min():.4f}")"""
)

# ── 5a: Top-30 Shortlist ──
md("### 5a. Top-30 Shortlist")

code(
    """\
from p53cad.results.schema import select_presentation_shortlist

top = select_presentation_shortlist(deep, top_n=SHORTLIST_N)
print(f"Shortlist: {len(top)} candidates\\n")

# Build display table
display_rows = []
for _, row in top.iterrows():
    targets = json.loads(str(row.get("targets_json", "[]")))
    muts = json.loads(str(row.get("mutations_json", "[]")))
    rescue = [m for m in muts if m not in targets]
    dms_z = row.get("rescue_dms_mean", None)
    dms_str = f"{dms_z:+.2f}" if dms_z and dms_z != 0 else "N/A"
    pr = row.get("pareto_rank", None)
    pr_str = f"{int(pr)}" if pr and pr != np.inf else "?"
    display_rows.append({
        "Rank": int(row.get("presentation_rank", 0)),
        "Target": row["target_label"],
        "Rescue Mutations": "+".join(rescue),
        "Oracle Score": f"{float(row['score']):.3f}",
        "DMS Z-score": dms_str,
        "Pareto Rank": pr_str,
        "Delivery": row["delivery_method"],
    })

shortlist_df = pd.DataFrame(display_rows)
shortlist_df"""
)

code(
    """\
# Delivery method distribution
if len(top):
    delivery_counts = top["delivery_method"].value_counts()
    print("Delivery Distribution:")
    for d, c in delivery_counts.items():
        print(f"  {d}: {c}")

    unique_targets = top["target_label"].nunique()
    print(f"\\nUnique target combos: {unique_targets}")"""
)

# ── 5b: Score Distribution Plots ──
md("### 5b. Score Distribution Plots")

code(
    """\
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 3, figsize=(16, 4))

# 1. Oracle score distribution
if "score" in deep.columns:
    axes[0].hist(deep["score"].dropna(), bins=40, color="#2563EB", alpha=0.7, edgecolor="white")
    axes[0].set_title("Oracle Score Distribution")
    axes[0].set_xlabel("Oracle Score")
    axes[0].set_ylabel("Count")

# 2. DMS rescue quality
if "rescue_dms_mean" in deep.columns:
    valid = deep["rescue_dms_mean"].dropna()
    valid = valid[valid != 0]
    if len(valid):
        axes[1].hist(valid, bins=40, color="#10B981", alpha=0.7, edgecolor="white")
        axes[1].axvline(0, color="red", linestyle="--", alpha=0.7)
        axes[1].set_title("DMS Rescue Z-score")
        axes[1].set_xlabel("Z-score")
        axes[1].set_ylabel("Count")

# 3. Scores by target
if "target_label" in deep.columns and "score" in deep.columns:
    # Get single-hotspot targets only for readability
    singles = deep[~deep["target_label"].str.contains(r"\\+", na=False)].copy()
    if len(singles):
        target_order = singles.groupby("target_label")["score"].median().sort_values(ascending=False).index
        sns.boxplot(data=singles, x="target_label", y="score", order=target_order,
                    ax=axes[2], palette="Blues_d")
        axes[2].set_title("Score by Hotspot")
        axes[2].set_xlabel("")
        axes[2].tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.show()"""
)

# ── 5c: Mutation Position Heatmap ──
md("### 5c. Mutation Position Heatmap")

code(
    """\
from p53cad.results.visualization import build_candidate_position_heatmap

# Prepare candidate dicts for heatmap
cand_dicts = []
for _, row in top.iterrows():
    muts = json.loads(str(row.get("mutations_json", "[]")))
    cand_dicts.append({
        "candidate_id": row.get("presentation_rank", 0),
        "profile": row.get("profile", "Unknown"),
        "score": float(row.get("score", 0)),
        "mutations": muts,
    })

matrix_df, freq_df = build_candidate_position_heatmap(cand_dicts)

if not matrix_df.empty:
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={"height_ratios": [3, 1]})

    # Heatmap
    pos_cols = [c for c in matrix_df.columns if c.startswith("pos_")]
    heatmap_data = matrix_df[pos_cols].values
    pos_labels = [c.replace("pos_", "") for c in pos_cols]

    axes[0].imshow(heatmap_data, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    axes[0].set_xticks(range(len(pos_labels)))
    axes[0].set_xticklabels(pos_labels, rotation=90, fontsize=8)
    axes[0].set_ylabel("Candidate")
    axes[0].set_title("Mutation Position Heatmap (Top-30 Shortlist)")

    # Frequency bar chart
    if not freq_df.empty:
        axes[1].bar(range(len(freq_df)), freq_df["frequency"], color="#2563EB", alpha=0.7)
        axes[1].set_xticks(range(len(freq_df)))
        axes[1].set_xticklabels(freq_df["position"].astype(str), rotation=90, fontsize=8)
        axes[1].set_ylabel("Frequency")
        axes[1].set_xlabel("Position")
        axes[1].set_title("Mutation Frequency by Position")

    plt.tight_layout()
    plt.show()
else:
    print("No mutation position data available for heatmap")"""
)

# ── 5d: Clinical Impact ──
md("### 5d. Clinical Impact")

code(
    """\
from p53cad.analysis.clinical_impact import TCGA_P53_MUTATIONS, CANCER_INCIDENCE, P53_MUTATION_RATE
from p53cad.results.schema import BIG8_HOTSPOTS

# Show TCGA hotspot frequencies
print("TCGA p53 Hotspot Mutation Frequencies")
print("=" * 55)
for mut in BIG8_HOTSPOTS:
    info = TCGA_P53_MUTATIONS.get(mut, {})
    freq = info.get("frequency", 0)
    cancers = info.get("cancer_types", [])
    print(f"  {mut:<8}  {freq:>4.1f}%  ->  {', '.join(cancers)}")

# Estimate affected patients
print(f"\\nEstimated Annual Patient Impact (US):")
print(f"  {'Cancer':<15} {'Incidence':>10} {'p53 mut rate':>12} {'p53 mutant':>12}")
print(f"  {'-'*15} {'-'*10} {'-'*12} {'-'*12}")
total_p53 = 0
for cancer, incidence in sorted(CANCER_INCIDENCE.items(), key=lambda x: -x[1]):
    rate = P53_MUTATION_RATE.get(cancer, 0)
    p53_cases = int(incidence * rate / 100)
    total_p53 += p53_cases
    print(f"  {cancer:<15} {incidence:>10,} {rate:>11}% {p53_cases:>12,}")
print(f"  {'TOTAL':<15} {'':<10} {'':<12} {total_p53:>12,}")"""
)

# ── 5e: Clinical Impact per Candidate ──
md("### 5e. Clinical Impact per Shortlisted Candidate")

code(
    """\
# Load clinical impact data if available
clinical_path = os.path.join(run_dir, "clinical.parquet")
if os.path.exists(clinical_path):
    clinical_df = pd.read_parquet(clinical_path)
    print(f"Clinical impact records: {len(clinical_df)}")
    clinical_df.head(10)
else:
    print("Clinical impact data not generated in this run.")
    print("(Set with_clinical=True and budget >= 'medium' for clinical analysis)")"""
)

# ── Section 6: Export Results ──
md("## 6. Export Results")

code(
    """\
# Export shortlist to CSV
csv_path = os.path.join(run_dir, "top30.csv")
if os.path.exists(csv_path):
    print(f"Shortlist CSV already saved: {csv_path}")
else:
    shortlist_df.to_csv(csv_path, index=False)
    print(f"Shortlist saved to: {csv_path}")

# Summary
summary_path = os.path.join(run_dir, "summary.md")
if os.path.exists(summary_path):
    with open(summary_path) as f:
        print(f.read())
else:
    print("\\nNo summary.md generated. Key results:")
    print(f"  Run dir: {run_dir}")
    print(f"  Candidates: {result['n_candidates']}")
    print(f"  Shortlist: {result['n_shortlist']}")"""
)

code(
    """\
# List all output artifacts
print(f"\\nAll artifacts in {run_dir}:")
for f in sorted(os.listdir(run_dir)):
    fpath = os.path.join(run_dir, f)
    if os.path.isfile(fpath):
        size = os.path.getsize(fpath)
        unit = "KB" if size > 1024 else "B"
        size_str = f"{size/1024:.1f} KB" if size > 1024 else f"{size} B"
        print(f"  {f:<35} {size_str}")
    elif os.path.isdir(fpath):
        n_files = len(os.listdir(fpath))
        print(f"  {f + '/':<35} ({n_files} files)")"""
)

code(
    """\
# Download results (Colab only)
IN_COLAB = 'google.colab' in str(get_ipython()) if 'get_ipython' in dir() else False

if IN_COLAB:
    from google.colab import files
    import shutil

    # Create a zip of all results
    results_zip = shutil.make_archive('p53_campaign_results', 'zip', run_dir)
    print(f"\\nDownloading results archive: {results_zip}")
    files.download(results_zip)
else:
    print(f"\\nResults saved to: {run_dir}")"""
)

md(
    """\
---

**Done!** The campaign has completed. Key outputs:
- `candidates.parquet` — All evaluated candidates
- `top30.parquet` / `top30.csv` — Shortlisted rescue mutations
- `clinical.parquet` — Patient impact estimates
- `trajectories.parquet` — Optimization trajectories
- `summary.md` — Human-readable report

On Colab, results are automatically zipped and downloaded."""
)

# ── Assemble the notebook ────────────────────────────────────────────────

notebook = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {
            "name": "python",
            "version": "3.10.0",
            "mimetype": "text/x-python",
            "file_extension": ".py",
        },
        "colab": {
            "provenance": [],
            "gpuType": "T4",
        },
        "accelerator": "GPU",
    },
    "cells": cells,
}

out_path = PROJECT_ROOT / "run_simulation.ipynb"
out_path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"Wrote {len(cells)} cells to {out_path}")
print(f"  Markdown cells: {sum(1 for c in cells if c['cell_type'] == 'markdown')}")
print(f"  Code cells:     {sum(1 for c in cells if c['cell_type'] == 'code')}")
