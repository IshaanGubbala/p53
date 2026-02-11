"""
Generate publication-quality figures for p53-proteoMgCAD campaign results.
Style: Nature/Science journal format with clean typography.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib import patheffects
from pathlib import Path

# ── Style ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.2,
    "patch.linewidth": 0.5,
})

# ── Color palette (Nature-inspired) ─────────────────────────────────────
C_ORACLE  = "#2166AC"   # Blue
C_DMS     = "#B2182B"   # Red
C_FUNC    = "#4DAF4A"   # Green
C_NONFUNC = "#984EA3"   # Purple
C_GENE    = "#2166AC"   # Blue
C_MRNA    = "#E08214"   # Orange
C_PROTEIN = "#1B7837"   # Green
C_BG      = "#F7F7F7"
C_GRID    = "#E0E0E0"

DELIVERY_COLORS = {
    "gene_therapy": C_GENE,
    "mrna_therapy": C_MRNA,
    "protein_therapy": C_PROTEIN,
}
DELIVERY_LABELS = {
    "gene_therapy": "Gene Therapy",
    "mrna_therapy": "mRNA Therapy",
    "protein_therapy": "Protein Therapy",
}

# ── Data ─────────────────────────────────────────────────────────────────
CAMPAIGN_DIR = Path("/Users/ishaangubbala/Documents/p53/data/campaigns/campaign_20260211_130341")
OUT_DIR = CAMPAIGN_DIR / "figures"
OUT_DIR.mkdir(exist_ok=True)

top = pd.read_csv(CAMPAIGN_DIR / "top30.csv")
cand = pd.read_parquet(CAMPAIGN_DIR / "candidates.parquet")
deep = cand[cand["pass_name"] == "deep"].copy()


def parse_rescue(row):
    """Extract rescue mutations (mutations - targets)."""
    muts = json.loads(str(row.get("mutations_json", "[]")))
    targets = json.loads(str(row.get("targets_json", "[]")))
    return [m for m in muts if m not in targets]


# =====================================================================
# FIGURE 1: Top 15 Candidates — Oracle Score + DMS Z-Score (dual bars)
# =====================================================================
def fig1_top_candidates():
    fig, ax1 = plt.subplots(figsize=(7.2, 4.5))

    df = top.head(15).copy()
    df["rescue"] = df.apply(parse_rescue, axis=1).apply(lambda x: "+".join(x))
    df["label"] = df["target_label"] + "\n" + df["rescue"]

    y = np.arange(len(df))
    bar_h = 0.35

    # Oracle score bars (normalized to fit visual range)
    oracle_vals = df["score"].values
    bars1 = ax1.barh(y + bar_h/2, oracle_vals, bar_h, color=C_ORACLE, alpha=0.85,
                     edgecolor="white", linewidth=0.3, label="Oracle Score", zorder=3)

    ax1.set_xlabel("Oracle Score", color=C_ORACLE, fontweight="bold")
    ax1.set_xlim(1.55, 1.75)
    ax1.set_yticks(y)
    ax1.set_yticklabels(df["label"], fontsize=7.5)
    ax1.invert_yaxis()
    ax1.tick_params(axis="x", colors=C_ORACLE)

    # DMS Z-score on secondary axis
    ax2 = ax1.twiny()
    dms_vals = df["rescue_dms_mean"].values
    colors = [C_FUNC if v < 0 else C_DMS for v in dms_vals]
    bars2 = ax2.barh(y - bar_h/2, dms_vals, bar_h, color=colors, alpha=0.75,
                     edgecolor="white", linewidth=0.3, zorder=3)

    ax2.set_xlabel("DMS Z-Score (negative = functional)", color=C_DMS, fontweight="bold")
    ax2.axvline(0, color="#666666", linewidth=0.6, linestyle="--", zorder=2)
    ax2.set_xlim(-2.5, 2.0)
    ax2.tick_params(axis="x", colors=C_DMS)

    # Delivery method indicators
    for i, (_, row) in enumerate(df.iterrows()):
        dm = row["delivery_method"]
        color = DELIVERY_COLORS.get(dm, "#999999")
        ax1.plot(1.555, y[i], "s", color=color, markersize=5, zorder=4, clip_on=False)

    # Legend
    oracle_patch = mpatches.Patch(color=C_ORACLE, alpha=0.85, label="Oracle Score")
    func_patch = mpatches.Patch(color=C_FUNC, alpha=0.75, label="DMS Z < 0 (Functional)")
    nonfunc_patch = mpatches.Patch(color=C_DMS, alpha=0.75, label="DMS Z > 0 (LoF)")
    gene_patch = mpatches.Patch(color=C_GENE, label="Gene Therapy")
    mrna_patch = mpatches.Patch(color=C_MRNA, label="mRNA Therapy")
    prot_patch = mpatches.Patch(color=C_PROTEIN, label="Protein Therapy")

    ax1.legend(handles=[oracle_patch, func_patch, nonfunc_patch, gene_patch, mrna_patch, prot_patch],
               loc="lower right", frameon=True, fancybox=False, edgecolor="#CCCCCC",
               fontsize=7, ncol=2)

    ax1.set_title("Top 15 Rescue Candidates: Oracle Score & DMS Fitness", fontweight="bold", pad=25)
    ax1.set_facecolor(C_BG)
    ax1.grid(axis="x", color=C_GRID, linewidth=0.5, zorder=0)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_top_candidates.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig1_top_candidates.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 1: Top 15 candidates")


# =====================================================================
# FIGURE 2: DMS Z-Score Distribution — Shortlist vs All Deep
# =====================================================================
def fig2_dms_distribution():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [2, 1]})

    # Panel A: Histogram comparison
    deep_dms = deep["rescue_dms_mean"].dropna()
    deep_dms = deep_dms[deep_dms != 0]
    short_dms = top["rescue_dms_mean"].dropna()
    short_dms = short_dms[short_dms != 0]

    bins = np.linspace(-3, 3, 30)
    ax1.hist(deep_dms, bins=bins, density=True, alpha=0.4, color="#AAAAAA",
             edgecolor="white", linewidth=0.3, label=f"All deep (n={len(deep_dms)})", zorder=2)
    ax1.hist(short_dms, bins=bins, density=True, alpha=0.7, color=C_FUNC,
             edgecolor="white", linewidth=0.3, label=f"Shortlisted (n={len(short_dms)})", zorder=3)

    ax1.axvline(0, color="#666666", linewidth=0.8, linestyle="--", zorder=4)
    ax1.axvline(short_dms.mean(), color=C_DMS, linewidth=1.2, linestyle="-", zorder=4)
    ax1.text(short_dms.mean() - 0.05, ax1.get_ylim()[1] * 0.9,
             f"Shortlist\nmean={short_dms.mean():.2f}",
             ha="right", va="top", fontsize=7, color=C_DMS, fontweight="bold")

    ax1.set_xlabel("Rescue DMS Z-Score")
    ax1.set_ylabel("Density")
    ax1.set_title("A", fontweight="bold", loc="left", fontsize=12)
    ax1.legend(frameon=True, fancybox=False, edgecolor="#CCCCCC")
    ax1.set_facecolor(C_BG)
    ax1.grid(axis="y", color=C_GRID, linewidth=0.5, zorder=0)

    # Panel B: Functional rescue fraction
    categories = ["All Deep", "Shortlisted"]
    functional_fracs = [
        (deep_dms < 0).sum() / len(deep_dms) * 100,
        (short_dms < 0).sum() / len(short_dms) * 100,
    ]
    nonfunc_fracs = [100 - f for f in functional_fracs]

    bars_func = ax2.bar(categories, functional_fracs, color=C_FUNC, alpha=0.8,
                        edgecolor="white", linewidth=0.5, label="Functional (Z<0)", zorder=3)
    bars_nonfunc = ax2.bar(categories, nonfunc_fracs, bottom=functional_fracs,
                           color="#DDDDDD", edgecolor="white", linewidth=0.5,
                           label="Non-functional", zorder=3)

    for bar, frac in zip(bars_func, functional_fracs):
        ax2.text(bar.get_x() + bar.get_width()/2, frac/2,
                 f"{frac:.0f}%", ha="center", va="center",
                 fontweight="bold", fontsize=10, color="white",
                 path_effects=[patheffects.withStroke(linewidth=2, foreground="black")])

    ax2.set_ylabel("Percentage")
    ax2.set_ylim(0, 105)
    ax2.set_title("B", fontweight="bold", loc="left", fontsize=12)
    ax2.legend(frameon=True, fancybox=False, edgecolor="#CCCCCC", fontsize=7)
    ax2.set_facecolor(C_BG)
    ax2.grid(axis="y", color=C_GRID, linewidth=0.5, zorder=0)

    fig.suptitle("DMS-Aware Selection Enriches for Functional Rescue Mutations", fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig2_dms_distribution.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig2_dms_distribution.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 2: DMS distribution")


# =====================================================================
# FIGURE 3: Oracle Score vs DMS Z-Score Scatter (colored by delivery)
# =====================================================================
def fig3_oracle_vs_dms():
    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    # Plot all deep candidates as background
    deep_valid = deep[deep["rescue_dms_mean"].notna() & (deep["rescue_dms_mean"] != 0)].copy()
    ax.scatter(deep_valid["score"], deep_valid["rescue_dms_mean"],
               s=8, alpha=0.15, color="#BBBBBB", edgecolors="none", zorder=2, label="_nolegend_")

    # Shortlisted candidates with delivery color
    for dm, color in DELIVERY_COLORS.items():
        subset = top[top["delivery_method"] == dm]
        ax.scatter(subset["score"], subset["rescue_dms_mean"],
                   s=50, alpha=0.85, color=color, edgecolors="white", linewidth=0.5,
                   zorder=4, label=DELIVERY_LABELS[dm])

    # Quadrant shading
    ax.axhline(0, color="#666666", linewidth=0.6, linestyle="--", zorder=1)
    ax.axvspan(1.65, 1.75, ymin=0, ymax=0.5, alpha=0.05, color=C_FUNC, zorder=0)

    # Label top 5
    for _, row in top.head(5).iterrows():
        rescue = parse_rescue(row)
        label = "+".join(rescue)
        ax.annotate(label, (row["score"], row["rescue_dms_mean"]),
                    xytext=(8, 5), textcoords="offset points",
                    fontsize=6.5, fontweight="bold", color="#333333",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))

    ax.set_xlabel("Oracle Score (predicted function)", fontweight="bold")
    ax.set_ylabel("Rescue DMS Z-Score (negative = functional)", fontweight="bold")
    ax.set_title("Predicted Function vs. Experimental Fitness of Rescue Mutations", fontweight="bold")
    ax.legend(frameon=True, fancybox=False, edgecolor="#CCCCCC", markerscale=1.2)
    ax.set_facecolor(C_BG)
    ax.grid(color=C_GRID, linewidth=0.5, zorder=0)

    # Annotation for ideal quadrant
    ax.text(1.71, -1.8, "Ideal: High oracle,\nfunctional rescue",
            fontsize=7, ha="center", color=C_FUNC, fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_FUNC, alpha=0.8))

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_oracle_vs_dms.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig3_oracle_vs_dms.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 3: Oracle vs DMS scatter")


# =====================================================================
# FIGURE 4: Target Retention Heatmap
# =====================================================================
def fig4_retention_heatmap():
    # Build retention matrix
    hotspots = ["R175H", "G245S", "R248Q", "R249S", "R273H", "R273C", "R282W", "R248W"]

    retention_data = {}
    for _, row in deep.iterrows():
        tl = row["target_label"]
        targets = json.loads(str(row.get("targets_json", "[]")))
        muts = json.loads(str(row.get("mutations_json", "[]")))
        if tl not in retention_data:
            retention_data[tl] = [0, 0]
        retention_data[tl][1] += 1
        if set(targets).issubset(set(muts)):
            retention_data[tl][0] += 1

    # Build matrix for pairs
    n = len(hotspots)
    matrix = np.full((n, n), np.nan)
    counts = np.full((n, n), "", dtype=object)

    for tl, (ret, total) in retention_data.items():
        parts = tl.split("+")
        if len(parts) == 1:
            i = hotspots.index(parts[0]) if parts[0] in hotspots else -1
            if i >= 0:
                matrix[i, i] = ret / total * 100
                counts[i, i] = f"{ret}/{total}"
        elif len(parts) == 2:
            i = hotspots.index(parts[0]) if parts[0] in hotspots else -1
            j = hotspots.index(parts[1]) if parts[1] in hotspots else -1
            if i >= 0 and j >= 0:
                matrix[i, j] = ret / total * 100
                matrix[j, i] = ret / total * 100
                counts[i, j] = f"{ret}/{total}"
                counts[j, i] = f"{ret}/{total}"

    fig, ax = plt.subplots(figsize=(6, 5))

    # Custom colormap: red -> yellow -> green
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("retention",
        ["#D73027", "#FC8D59", "#FEE08B", "#D9EF8B", "#91CF60", "#1A9850"])

    im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=100, aspect="equal")

    # Annotate
    for i in range(n):
        for j in range(n):
            if not np.isnan(matrix[i, j]):
                val = matrix[i, j]
                txt = f"{val:.0f}%\n{counts[i,j]}"
                color = "white" if val < 30 else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.5,
                        fontweight="bold", color=color)

    ax.set_xticks(range(n))
    ax.set_xticklabels(hotspots, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(hotspots)
    ax.set_title("Target Mutation Retention Rate (%)\nacross Cancer Hotspot Combinations", fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="Retention %")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_retention_heatmap.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig4_retention_heatmap.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 4: Retention heatmap")


# =====================================================================
# FIGURE 5: Mutation Landscape — which positions are used for rescue
# =====================================================================
def fig5_mutation_landscape():
    fig, ax = plt.subplots(figsize=(7.2, 3.5))

    # Count rescue mutation positions across shortlist
    position_counts = {}
    position_dms = {}

    for _, row in top.iterrows():
        rescue = parse_rescue(row)
        dms = row.get("rescue_dms_mean", 0)
        for mut in rescue:
            # Extract position number
            pos = int("".join(c for c in mut if c.isdigit()))
            position_counts[pos] = position_counts.get(pos, 0) + 1
            if pos not in position_dms:
                position_dms[pos] = []
            position_dms[pos].append(dms)

    positions = sorted(position_counts.keys())
    counts = [position_counts[p] for p in positions]
    mean_dms = [np.mean(position_dms[p]) for p in positions]

    # Color by mean DMS
    colors = [C_FUNC if d < 0 else C_DMS for d in mean_dms]

    bars = ax.bar(range(len(positions)), counts, color=colors, alpha=0.8,
                  edgecolor="white", linewidth=0.5, zorder=3)

    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([str(p) for p in positions], rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("Residue Position")
    ax.set_ylabel("Frequency in Shortlist")
    ax.set_title("Rescue Mutation Landscape: Residue Positions Used Across Top 30 Candidates", fontweight="bold")
    ax.set_facecolor(C_BG)
    ax.grid(axis="y", color=C_GRID, linewidth=0.5, zorder=0)

    # Annotate most frequent
    for i, (pos, count) in enumerate(zip(positions, counts)):
        if count >= 3:
            ax.text(i, count + 0.15, f"pos {pos}", ha="center", va="bottom",
                    fontsize=6.5, fontweight="bold", rotation=0)

    # Legend
    func_p = mpatches.Patch(color=C_FUNC, alpha=0.8, label="Functional (mean DMS Z < 0)")
    nonfunc_p = mpatches.Patch(color=C_DMS, alpha=0.8, label="Non-functional (mean DMS Z > 0)")
    ax.legend(handles=[func_p, nonfunc_p], frameon=True, fancybox=False, edgecolor="#CCCCCC")

    # Domain annotations
    ax.axvspan(-0.5, -0.5, alpha=0)  # dummy
    # Mark DNA-binding domain range
    dbd_positions = [i for i, p in enumerate(positions) if 94 <= p <= 292]
    tet_positions = [i for i, p in enumerate(positions) if 323 <= p <= 356]

    if dbd_positions:
        ax.axvspan(min(dbd_positions)-0.5, max(dbd_positions)+0.5, alpha=0.08,
                   color="#457B9D", zorder=0, label="DNA-binding domain")
    if tet_positions:
        ax.axvspan(min(tet_positions)-0.5, max(tet_positions)+0.5, alpha=0.08,
                   color="#E63946", zorder=0, label="Tetramerization domain")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_mutation_landscape.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig5_mutation_landscape.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 5: Mutation landscape")


# =====================================================================
# FIGURE 6: Campaign Summary Panel (4-panel composite)
# =====================================================================
def fig6_summary_panel():
    fig = plt.figure(figsize=(7.2, 6.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    # Panel A: Delivery method distribution
    ax_a = fig.add_subplot(gs[0, 0])
    delivery_counts = top["delivery_method"].value_counts()
    wedges, texts, autotexts = ax_a.pie(
        delivery_counts.values,
        labels=[DELIVERY_LABELS.get(d, d) for d in delivery_counts.index],
        colors=[DELIVERY_COLORS.get(d, "#999") for d in delivery_counts.index],
        autopct=lambda p: f"{p:.0f}%\n({int(p*len(top)/100)})",
        startangle=90,
        pctdistance=0.65,
        wedgeprops=dict(edgecolor="white", linewidth=1.5),
    )
    for t in autotexts:
        t.set_fontsize(8)
        t.set_fontweight("bold")
    ax_a.set_title("A  Delivery Methods", fontweight="bold", loc="left", fontsize=10)

    # Panel B: Score distribution violin
    ax_b = fig.add_subplot(gs[0, 1])
    shortlist_scores = top["score"].values
    deep_scores = deep["score"].values
    parts = ax_b.violinplot([deep_scores, shortlist_scores], positions=[1, 2],
                            showmeans=True, showmedians=True, widths=0.6)
    for pc in parts["bodies"]:
        pc.set_facecolor(C_ORACLE)
        pc.set_alpha(0.4)
    parts["cmeans"].set_color(C_DMS)
    parts["cmedians"].set_color("#333333")
    ax_b.set_xticks([1, 2])
    ax_b.set_xticklabels(["All Deep\n(n=528)", "Shortlisted\n(n=30)"])
    ax_b.set_ylabel("Oracle Score")
    ax_b.set_title("B  Score Distribution", fontweight="bold", loc="left", fontsize=10)
    ax_b.set_facecolor(C_BG)
    ax_b.grid(axis="y", color=C_GRID, linewidth=0.5, zorder=0)

    # Panel C: Target combinations in shortlist
    ax_c = fig.add_subplot(gs[1, 0])
    target_counts = top["target_label"].value_counts().head(10)
    colors_c = [C_FUNC if i < 5 else "#AAAAAA" for i in range(len(target_counts))]
    bars_c = ax_c.barh(range(len(target_counts)), target_counts.values,
                       color=colors_c, alpha=0.8, edgecolor="white", linewidth=0.5, zorder=3)
    ax_c.set_yticks(range(len(target_counts)))
    ax_c.set_yticklabels(target_counts.index, fontsize=7.5)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Candidates in Shortlist")
    ax_c.set_title("C  Target Hotspot Coverage", fontweight="bold", loc="left", fontsize=10)
    ax_c.set_facecolor(C_BG)
    ax_c.grid(axis="x", color=C_GRID, linewidth=0.5, zorder=0)

    # Panel D: Key metrics summary
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis("off")

    valid_dms = top["rescue_dms_mean"].dropna()
    valid_dms = valid_dms[valid_dms != 0]

    metrics = [
        ("Campaign Scenarios", "108 screen + 44 deep"),
        ("Total Candidates", f"{len(cand):,}"),
        ("Deep Candidates", f"{len(deep):,}"),
        ("Shortlist Size", f"{len(top)}/30"),
        ("Best Oracle Score", f"{top['score'].max():.3f}"),
        ("Best DMS Z-Score", f"{valid_dms.min():.2f}"),
        ("Functional Rescues", f"{(valid_dms < 0).sum()}/{len(valid_dms)} ({(valid_dms < 0).mean()*100:.0f}%)"),
        ("Mean Shortlist DMS Z", f"{valid_dms.mean():.3f}"),
        ("Target Retention", f"206/528 (39.0%)"),
        ("Unique Targets", f"{top['target_label'].nunique()}"),
    ]

    for i, (label, value) in enumerate(metrics):
        y = 0.95 - i * 0.095
        ax_d.text(0.02, y, label, fontsize=8, fontweight="bold", va="top",
                  transform=ax_d.transAxes, color="#333333")
        ax_d.text(0.98, y, value, fontsize=8, va="top", ha="right",
                  transform=ax_d.transAxes, color=C_ORACLE, fontweight="bold")

    ax_d.set_title("D  Campaign Metrics", fontweight="bold", loc="left", fontsize=10)
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    # Add border
    for spine in ax_d.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("#CCCCCC")
        spine.set_linewidth(0.8)
    ax_d.set_facecolor(C_BG)

    fig.suptitle("p53-proteoMgCAD: DMS-Aware Campaign Summary",
                 fontweight="bold", fontsize=12, y=1.01)

    fig.savefig(OUT_DIR / "fig6_summary_panel.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig6_summary_panel.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 6: Summary panel")


# =====================================================================
# FIGURE 7: Evidence Score Waterfall
# =====================================================================
def fig7_evidence_waterfall():
    fig, ax = plt.subplots(figsize=(7.2, 3.5))

    df = top.copy()
    # Compute evidence score: oracle_norm * dms_component * clinical_norm
    # Already in selection_score or we approximate
    df["rescue"] = df.apply(parse_rescue, axis=1).apply(lambda x: "+".join(x[:2]))
    df = df.sort_values("score", ascending=False).head(20)

    x = np.arange(len(df))

    # Stacked components: oracle contribution, DMS contribution, clinical contribution
    oracle_norm = (df["score"] - 1.5) / 0.25  # normalize to ~[0,1]
    oracle_norm = oracle_norm.clip(0, 1)

    dms_component = (-df["rescue_dms_mean"].fillna(0)).clip(-1, 2) / 2  # normalize
    dms_component = dms_component.clip(0, 1)

    clinical_norm = df["clinical_score"].replace(-np.inf, 0) / 100
    clinical_norm = clinical_norm.clip(0, 1)

    # Composite
    composite = (oracle_norm * 0.4 + dms_component * 0.35 + clinical_norm * 0.25).values

    colors = [DELIVERY_COLORS.get(dm, "#999") for dm in df["delivery_method"]]

    bars = ax.bar(x, composite, color=colors, alpha=0.85, edgecolor="white",
                  linewidth=0.5, zorder=3)

    ax.set_xticks(x)
    labels = [f"{row['target_label']}\n{row['rescue']}" for _, row in df.iterrows()]
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Composite Evidence Score")
    ax.set_title("Evidence-Weighted Candidate Ranking (Top 20)", fontweight="bold")
    ax.set_facecolor(C_BG)
    ax.grid(axis="y", color=C_GRID, linewidth=0.5, zorder=0)

    # Delivery legend
    handles = [mpatches.Patch(color=c, label=DELIVERY_LABELS[d])
               for d, c in DELIVERY_COLORS.items()]
    ax.legend(handles=handles, frameon=True, fancybox=False, edgecolor="#CCCCCC",
              fontsize=7, loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig7_evidence_waterfall.png", facecolor="white")
    fig.savefig(OUT_DIR / "fig7_evidence_waterfall.pdf", facecolor="white")
    plt.close(fig)
    print("  [+] Figure 7: Evidence waterfall")


# =====================================================================
# RUN ALL
# =====================================================================
if __name__ == "__main__":
    print(f"\nGenerating publication figures -> {OUT_DIR}/\n")
    fig1_top_candidates()
    fig2_dms_distribution()
    fig3_oracle_vs_dms()
    fig4_retention_heatmap()
    fig5_mutation_landscape()
    fig6_summary_panel()
    fig7_evidence_waterfall()
    print(f"\nDone! 7 figures (PNG + PDF) in {OUT_DIR}/")
