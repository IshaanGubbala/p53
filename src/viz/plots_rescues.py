from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb


_MUT_POS_RE = re.compile(r"[A-Za-z\\*](\\d+)[A-Za-z\\*]")


def _gradient(ax, top: str, bottom: str, alpha: float = 1.0) -> None:
    top_rgb = to_rgb(top)
    bottom_rgb = to_rgb(bottom)
    cmap = LinearSegmentedColormap.from_list("bg", [bottom_rgb, top_rgb])
    gradient = np.linspace(0, 1, 256).reshape(256, 1)
    ax.imshow(
        gradient,
        aspect="auto",
        cmap=cmap,
        extent=[0, 1, 0, 1],
        transform=ax.transAxes,
        zorder=0,
        alpha=alpha,
    )


def _style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(color="#ffffff", alpha=0.35, linewidth=0.6)
    ax.tick_params(axis="both", length=0)


def _to_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    df = df.copy()
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def plot_rescue_bubble(
    df: pd.DataFrame,
    out_path: Path,
    title: str | None = None,
    annotate_n: int = 6,
) -> None:
    if not {"ddg_gain", "risk", "rescue_mutations"}.issubset(df.columns):
        raise ValueError("Rescue dataframe missing required columns for bubble plot")

    data = _to_numeric(df, ["ddg_gain", "risk"]).dropna(subset=["ddg_gain", "risk"]).copy()
    if data.empty:
        raise ValueError("No rescue data available for bubble plot")
    data = data.reset_index(drop=True)

    score = (-data["ddg_gain"]).clip(lower=0.0)
    score_min = float(score.min())
    score_max = float(score.max())
    score_norm = (score - score_min) / max(score_max - score_min, 1e-6)
    size = 40 + 240 * (score_norm**0.6)

    risk_values = data["risk"].to_numpy()
    risk_span = float(risk_values.max() - risk_values.min()) if len(risk_values) else 0.0
    jitter_scale = max(risk_span * 0.08, 0.002)
    rng = np.random.default_rng(42)
    risk_plot = np.clip(
        risk_values + rng.normal(0.0, jitter_scale, size=len(risk_values)),
        risk_values.min() - jitter_scale,
        risk_values.max() + jitter_scale,
    )

    cmap = LinearSegmentedColormap.from_list("heat", ["#2a9d8f", "#e9c46a", "#e76f51"])

    with plt.rc_context(
        {
            "font.family": ["DejaVu Serif", "STIXGeneral"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
        }
    ):
        fig, ax = plt.subplots(figsize=(8.4, 6.2))
        _gradient(ax, "#fef6e4", "#edf6ff", alpha=1.0)
        _style_axes(ax)

        scatter = ax.scatter(
            risk_plot,
            data["ddg_gain"],
            s=size,
            c=score,
            cmap=cmap,
            alpha=0.65,
            edgecolor="white",
            linewidth=0.5,
            zorder=2,
        )

        if "is_pareto" in data.columns:
            pareto = data[data["is_pareto"] == True]
            if not pareto.empty:
                ax.scatter(
                    pareto["risk"],
                    pareto["ddg_gain"],
                    s=180,
                    marker="*",
                    color="#ffb703",
                    edgecolor="#ffffff",
                    linewidth=0.8,
                    zorder=3,
                    label="Pareto",
                )

        ax.axhline(0.0, color="#1d3557", linewidth=1.0, alpha=0.5)
        ax.set_xlabel("Risk score (lower is safer)")
        ax.set_ylabel("DDG gain vs seed (negative is better)")
        if title:
            ax.set_title(title)
        ax.margins(x=0.08, y=0.1)

        top = data.nsmallest(max(annotate_n, 1), "ddg_gain")
        for idx, row in top.iterrows():
            ax.annotate(
                row["rescue_mutations"],
                (risk_plot[idx], row["ddg_gain"]),
                xytext=(8, 6),
                textcoords="offset points",
                fontsize=8,
                color="#1d3557",
                weight="bold",
            )

        cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label("Stability gain magnitude")
        cbar.outline.set_visible(False)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220)
        plt.close(fig)


def plot_risk_breakdown(
    df: pd.DataFrame,
    out_path: Path,
    top_n: int = 12,
) -> None:
    if "risk_components" not in df.columns:
        raise ValueError("Rescue dataframe missing risk_components column")

    data = _to_numeric(df, ["ddg_gain", "risk"]).copy()
    data = data.dropna(subset=["ddg_gain"])
    if data.empty:
        raise ValueError("No rescue data available for risk breakdown plot")

    data = data.nsmallest(top_n, "ddg_gain").copy()

    components = {"functional": [], "conservation": [], "burial": []}
    for raw in data["risk_components"]:
        try:
            parsed = json.loads(raw)
        except Exception:
            parsed = {}
        components["functional"].append(float(parsed.get("functional", 0.0)))
        components["conservation"].append(float(parsed.get("conservation", 0.0)))
        components["burial"].append(float(parsed.get("burial", 0.0)))

    labels = data["rescue_mutations"].tolist()
    colors = {
        "functional": "#264653",
        "conservation": "#e76f51",
        "burial": "#e9c46a",
    }

    with plt.rc_context(
        {
            "font.family": ["DejaVu Serif", "STIXGeneral"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
        }
    ):
        fig, ax = plt.subplots(figsize=(8.6, 6.2))
        _gradient(ax, "#fdf0f5", "#f0f7f9", alpha=1.0)
        _style_axes(ax)

        left = np.zeros(len(labels))
        for key in ["functional", "conservation", "burial"]:
            values = np.array(components[key])
            ax.barh(
                labels,
                values,
                left=left,
                color=colors[key],
                alpha=0.85,
                edgecolor="white",
                linewidth=0.5,
                label=key.capitalize(),
                zorder=2,
            )
            left += values

        for idx, (total, gain) in enumerate(zip(left, data["ddg_gain"].tolist())):
            ax.text(
                total + 0.02,
                idx,
                f"{gain:.2f}",
                va="center",
                fontsize=8,
                color="#1d3557",
                weight="bold",
            )

        ax.set_xlabel("Risk component score (lower is safer)")
        ax.set_title("Risk component mix for top rescue sets")
        ax.legend(frameon=False, loc="lower right")
        ax.invert_yaxis()
        total_max = float(left.max()) if len(left) else 0.0
        pad = max(total_max * 0.12, 0.15)
        ax.set_xlim(0, total_max + pad)
        ax.margins(x=0.02)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220)
        plt.close(fig)


def plot_position_frequency(
    df: pd.DataFrame,
    out_path: Path,
    top_n: int | None = None,
) -> None:
    if "rescue_mutations" not in df.columns:
        raise ValueError("Rescue dataframe missing rescue_mutations column")

    data = _to_numeric(df, ["ddg_gain"]).copy()
    if top_n:
        data = data.nsmallest(top_n, "ddg_gain")

    counts: dict[int, int] = {}
    for muts in data["rescue_mutations"]:
        if pd.isna(muts):
            continue
        for token in str(muts).split(","):
            match = _MUT_POS_RE.search(token.strip())
            if not match:
                continue
            pos = int(match.group(1))
            counts[pos] = counts.get(pos, 0) + 1

    if not counts:
        return

    positions = sorted(counts.keys())
    values = np.array([counts[pos] for pos in positions])
    pos_min = float(min(positions))
    pos_max = float(max(positions))
    pos_span = pos_max - pos_min
    x_pad = max(pos_span * 0.05, 1.0)
    y_max = float(values.max())
    y_pad = max(y_max * 0.12, 1.0)

    cmap = LinearSegmentedColormap.from_list("freq", ["#2a9d8f", "#f4a261", "#e76f51"])
    colors = cmap(values / max(values.max(), 1))

    with plt.rc_context(
        {
            "font.family": ["DejaVu Serif", "STIXGeneral"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
        }
    ):
        fig, ax = plt.subplots(figsize=(9.0, 4.5))
        _gradient(ax, "#f6fff8", "#f7f2ff", alpha=1.0)
        _style_axes(ax)

        ax.vlines(positions, 0, values, color=colors, alpha=0.9, linewidth=2.2, zorder=2)
        ax.scatter(positions, values, s=50 + 30 * (values / max(values.max(), 1)), color=colors, zorder=3)

        ax.set_xlabel("Residue position")
        ax.set_ylabel("Frequency in rescue sets")
        ax.set_title("Mutation position hotspots across top rescues")
        ax.set_xlim(pos_min - x_pad, pos_max + x_pad)
        ax.set_ylim(0, y_max + y_pad)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220)
        plt.close(fig)


def plot_ddg_gain_by_target(
    target_map: dict[str, pd.DataFrame],
    out_path: Path,
) -> None:
    if not target_map:
        raise ValueError("No target data provided for comparison plot")

    labels = []
    series = []
    for target, df in target_map.items():
        if "ddg_gain" not in df.columns:
            continue
        values = pd.to_numeric(df["ddg_gain"], errors="coerce").dropna()
        if values.empty:
            continue
        labels.append(target)
        series.append(values.to_numpy())

    if not series:
        raise ValueError("No ddg_gain values available for comparison plot")

    with plt.rc_context(
        {
            "font.family": ["DejaVu Serif", "STIXGeneral"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
        }
    ):
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        _gradient(ax, "#eef7ff", "#fff4e6", alpha=1.0)
        _style_axes(ax)

        parts = ax.violinplot(series, showmeans=True, showextrema=False, widths=0.7)
        palette = ["#2a9d8f", "#e9c46a", "#e76f51", "#264653"]
        for idx, body in enumerate(parts["bodies"]):
            body.set_facecolor(palette[idx % len(palette)])
            body.set_alpha(0.8)
            body.set_edgecolor("#ffffff")
            body.set_linewidth(0.6)

        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        ax.set_ylabel("DDG gain vs seed")
        ax.set_title("Rescue stability gains by target")
        ax.axhline(0.0, color="#1d3557", linewidth=1.0, alpha=0.5)
        ax.margins(y=0.12)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
