from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.logging import get_logger
from src.eval.variant_separation import bootstrap_ci, compute_auc, load_labelled_scores
from src.viz.plots_pareto import plot_pareto_front
from src.viz.plots_variants import plot_ddg_by_label


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    processed_dir = Path(data_cfg.get("processed", "data/processed"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
        "processed": processed_dir,
        "reports": reports_dir,
    }


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))
    report_cfg = configs.get("report", {})

    output_dir = Path(report_cfg.get("output_dir", paths["reports"]))
    tables_dir = Path(report_cfg.get("tables_dir", output_dir / "tables"))
    figures_dir = Path(report_cfg.get("figures_dir", output_dir / "figures"))
    fig_format = str(report_cfg.get("fig_format", "png"))
    top_rescues = int(report_cfg.get("top_rescues", 20))

    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    scores_path = paths["processed"] / "variant_scores.parquet"
    labels_dir = paths["processed"] / "labels"
    if scores_path.exists() and labels_dir.exists():
        try:
            labelled = load_labelled_scores(scores_path, labels_dir)
        except RuntimeError as exc:
            logger.warning("Skipping variant separation: %s", exc)
            labelled = None

        if labelled is not None:
            auc = compute_auc(labelled["label"], labelled["ddg"])
            lo, hi = bootstrap_ci(labelled["label"], labelled["ddg"], n=2000)

            summary = {
                "auc": auc,
                "auc_ci_lower": lo,
                "auc_ci_upper": hi,
                "n_variants": int(len(labelled)),
            }
            (tables_dir / "variant_separation.json").write_text(
                json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
            )
            pd.DataFrame([summary]).to_csv(tables_dir / "variant_separation.csv", index=False)
            plot_ddg_by_label(labelled, figures_dir / f"variant_ddg_by_label.{fig_format}")
            logger.info("Variant separation report written")
    else:
        logger.warning("Variant scores or labels not found; skipping variant report")

    rescues_root = paths["processed"] / "rescues"
    if rescues_root.exists():
        for target_dir in sorted(rescues_root.iterdir()):
            if not target_dir.is_dir():
                continue
            candidates_path = target_dir / "candidates.parquet"
            pareto_path = target_dir / "pareto.parquet"
            if not candidates_path.exists():
                continue

            candidates = pd.read_parquet(candidates_path)
            target = target_dir.name
            plot_pareto_front(
                candidates,
                figures_dir / f"pareto_{target}.{fig_format}",
            )

            # Export Pareto front first, then fill with best non-Pareto
            # Sort by: Pareto status (desc), then n_rescue (asc), then ddg_gain (asc), then risk (asc)
            ranked = candidates.sort_values(
                ["is_pareto", "n_rescue", "ddg_gain", "risk"],
                ascending=[False, True, True, True]
            )
            ranked.head(top_rescues).to_csv(tables_dir / f"rescues_{target}.csv", index=False)

        logger.info("Rescue report outputs written")
    else:
        logger.warning("No rescue outputs found; skipping rescue report")

    return 0
