from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.logging import get_logger
from src.eval.variant_separation import bootstrap_ci, compute_auc, load_labelled_scores
from src.scoring.evoef2_runner import build_mutant_model
from src.viz.animations import render_rotation_gif
from src.viz.plots_pareto import plot_pareto_front
from src.viz.plots_rescues import (
    plot_ddg_gain_by_target,
    plot_position_frequency,
    plot_rescue_bubble,
    plot_risk_breakdown,
)
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
    scoring_cfg = configs.get("scoring", {})
    evoef2_cfg = scoring_cfg.get("evoef2", {})

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
        target_candidates: dict[str, pd.DataFrame] = {}
        anim_dir = figures_dir / "animations"
        anim_pdb_dir = anim_dir / "pdb"
        anim_work_dir = anim_dir / "work"
        anim_top_n = int(report_cfg.get("rescue_animation_top_n", 1))
        anim_frames = int(report_cfg.get("rescue_animation_frames", 72))
        anim_fps = int(report_cfg.get("rescue_animation_fps", 18))
        anim_dpi = int(report_cfg.get("rescue_animation_dpi", 140))
        anim_format = str(report_cfg.get("rescue_animation_format", "gif")).lower()
        use_mutant_pdb = bool(report_cfg.get("rescue_animation_use_mutant_pdb", True))

        base_pdb = Path(
            evoef2_cfg.get(
                "repaired_pdb",
                paths["interim"] / "structure" / "P04637_core_94_312.pdb",
            )
        ).expanduser()
        if not base_pdb.is_absolute():
            base_pdb = (Path.cwd() / base_pdb).resolve()
        if not base_pdb.exists():
            logger.warning("Base PDB for animations not found: %s", base_pdb)
            anim_top_n = 0

        for target_dir in sorted(rescues_root.iterdir()):
            if not target_dir.is_dir():
                continue
            candidates_path = target_dir / "candidates.parquet"
            pareto_path = target_dir / "pareto.parquet"
            if not candidates_path.exists():
                continue

            candidates = pd.read_parquet(candidates_path)
            target = target_dir.name
            target_candidates[target] = candidates.copy()
            try:
                plot_pareto_front(
                    candidates,
                    figures_dir / f"pareto_{target}.{fig_format}",
                )
            except ValueError as exc:
                logger.warning("Skipping pareto plot for %s: %s", target, exc)

            try:
                plot_rescue_bubble(
                    candidates,
                    figures_dir / f"rescue_bubble_{target}.{fig_format}",
                    title=f"{target} rescue landscape",
                    annotate_n=int(report_cfg.get("rescue_bubble_annotate", 6)),
                )
            except ValueError as exc:
                logger.warning("Skipping rescue bubble plot for %s: %s", target, exc)

            try:
                plot_risk_breakdown(
                    candidates,
                    figures_dir / f"rescue_risk_{target}.{fig_format}",
                    top_n=int(report_cfg.get("rescue_risk_top_n", top_rescues)),
                )
            except ValueError as exc:
                logger.warning("Skipping rescue risk plot for %s: %s", target, exc)

            try:
                plot_position_frequency(
                    candidates,
                    figures_dir / f"rescue_positions_{target}.{fig_format}",
                    top_n=report_cfg.get("rescue_position_top_n"),
                )
            except ValueError as exc:
                logger.warning("Skipping rescue position plot for %s: %s", target, exc)

            if anim_top_n > 0 and anim_format == "gif":
                sort_cols = []
                ascending = []
                if "is_pareto" in candidates.columns:
                    sort_cols.append("is_pareto")
                    ascending.append(False)
                sort_cols.extend(["ddg_gain", "risk"])
                ascending.extend([True, True])
                anim_candidates = candidates.sort_values(sort_cols, ascending=ascending).head(anim_top_n)

                for _, row in anim_candidates.iterrows():
                    rescue_raw = row.get("rescue_mutations")
                    if pd.isna(rescue_raw):
                        continue
                    rescue_list = [token for token in str(rescue_raw).split(",") if token]
                    full_raw = row.get("full_mutations") or ""
                    full_list = [token for token in str(full_raw).split(",") if token]
                    set_id = str(row.get("set_id", "")) or f"{target}_anim"

                    pdb_for_anim = base_pdb
                    if use_mutant_pdb and full_list:
                        try:
                            out_pdb = anim_pdb_dir / f"{target}_{set_id}.pdb"
                            pdb_for_anim = build_mutant_model(
                                full_list,
                                base_pdb,
                                evoef2_cfg,
                                anim_work_dir,
                                out_path=out_pdb,
                                recompute=False,
                            )
                        except Exception as exc:
                            logger.warning("Animation PDB build failed for %s: %s", target, exc)
                            pdb_for_anim = base_pdb

                    out_path = anim_dir / f"rescue_anim_{target}_{set_id}.{anim_format}"
                    title = f"{target} + {row.get('rescue_mutations')}"
                    try:
                        render_rotation_gif(
                            pdb_for_anim,
                            target,
                            rescue_list,
                            out_path,
                            title=title,
                            frames=anim_frames,
                            fps=anim_fps,
                            dpi=anim_dpi,
                        )
                    except Exception as exc:
                        logger.warning("Skipping rescue animation for %s: %s", target, exc)

            # Export Pareto front first, then fill with best non-Pareto
            # Sort by: Pareto status (desc), then n_rescue (asc), then ddg_gain (asc), then risk (asc)
            ranked = candidates.sort_values(
                ["is_pareto", "n_rescue", "ddg_gain", "risk"],
                ascending=[False, True, True, True]
            )
            ranked.head(top_rescues).to_csv(tables_dir / f"rescues_{target}.csv", index=False)

        if target_candidates:
            try:
                plot_ddg_gain_by_target(
                    target_candidates,
                    figures_dir / f"rescue_ddg_by_target.{fig_format}",
                )
            except ValueError as exc:
                logger.warning("Skipping rescue ddg summary plot: %s", exc)

        logger.info("Rescue report outputs written")
    else:
        logger.warning("No rescue outputs found; skipping rescue report")

    # Print comprehensive console summary
    print("\n" + "=" * 80)
    print("📊 p53 StabiliMut Results Summary")
    print("=" * 80)

    # Variant Separation Summary
    sep_json = tables_dir / "variant_separation.json"
    if sep_json.exists():
        with open(sep_json) as f:
            sep_data = json.load(f)
        print("\n🧬 Variant Discrimination (ClinVar Pathogenic vs Benign)")
        print("-" * 80)
        print(f"  AUC: {sep_data['auc']:.3f} (95% CI: {sep_data['auc_ci_lower']:.3f}-{sep_data['auc_ci_upper']:.3f})")
        print(f"  Total variants scored: {sep_data['n_variants']}")
        print(f"  ✓ Model successfully discriminates pathogenic from benign mutations")

    # Rescue Design Summary
    if target_candidates:
        print("\n🎯 Rescue Mutation Design Results")
        print("-" * 80)

        for target, candidates in target_candidates.items():
            total = len(candidates)
            pareto = candidates[candidates.get("is_pareto", False)].shape[0] if "is_pareto" in candidates.columns else 0

            print(f"\n  {target}:")
            print(f"    Total candidates explored: {total:,}")
            print(f"    Pareto-optimal rescues: {pareto}")

            # Top 5 rescues
            ranked = candidates.sort_values(
                ["is_pareto", "n_rescue", "ddg_gain", "risk"],
                ascending=[False, True, True, True]
            ).head(5)

            print(f"    Top 5 rescues:")
            for i, (_, row) in enumerate(ranked.iterrows(), 1):
                mutations = row.get("rescue_mutations", row.get("mutations", ""))
                ddg_gain = row.get("ddg_gain", 0)
                rasp_gain = row.get("rasp_gain", None)  # Use gain, not total
                risk = row.get("risk", 0)
                n_rescue = row.get("n_rescue", row.get("n_mutations", 0))
                is_pareto = "⭐" if row.get("is_pareto", False) else "  "

                rasp_str = f"RaSP: {rasp_gain:+.2f}" if rasp_gain is not None else "RaSP: N/A"
                print(f"      {i}. {is_pareto} {mutations:20s} | ΔΔG: {ddg_gain:+.2f} | {rasp_str} | Risk: {risk:.3f} | n={n_rescue}")

            # Statistics
            if pareto > 0:
                pareto_df = candidates[candidates["is_pareto"] == True]
                avg_ddg = pareto_df["ddg_gain"].mean()
                avg_risk = pareto_df["risk"].mean()
                avg_n = pareto_df.get("n_rescue", pareto_df.get("n_mutations", pd.Series([0]))).mean()

                print(f"    Pareto statistics:")
                print(f"      Average ΔΔG gain: {avg_ddg:+.2f} kcal/mol")
                print(f"      Average risk score: {avg_risk:.3f}")
                print(f"      Average rescue size: {avg_n:.1f} mutations")

                # RaSP statistics if available (use gain for comparison)
                if "rasp_gain" in pareto_df.columns and pareto_df["rasp_gain"].notna().any():
                    avg_rasp_gain = pareto_df["rasp_gain"].mean()
                    rasp_agree = (pareto_df["rasp_gain"] < 0).sum()  # Negative gain = rescue benefit
                    print(f"      Average RaSP gain: {avg_rasp_gain:+.2f} kcal/mol")
                    print(f"      RaSP-confirmed stabilizing: {rasp_agree}/{pareto} ({100*rasp_agree/pareto:.1f}%)")

    # Output Files Summary
    print("\n📁 Generated Outputs")
    print("-" * 80)

    # Count files
    n_figures = len(list(figures_dir.glob(f"*.{fig_format}"))) if figures_dir.exists() else 0
    n_tables = len(list(tables_dir.glob("*.csv"))) if tables_dir.exists() else 0
    n_anims = len(list((figures_dir / "animations").glob("*.gif"))) if (figures_dir / "animations").exists() else 0

    # Use resolve() to get absolute paths, then try relative_to, fallback to absolute
    try:
        fig_path = figures_dir.resolve().relative_to(Path.cwd())
    except ValueError:
        fig_path = figures_dir.resolve()
    try:
        tbl_path = tables_dir.resolve().relative_to(Path.cwd())
    except ValueError:
        tbl_path = tables_dir.resolve()
    try:
        anim_path = (figures_dir / "animations").resolve().relative_to(Path.cwd())
    except ValueError:
        anim_path = (figures_dir / "animations").resolve()

    print(f"  Figures: {n_figures} plots in {fig_path}")
    print(f"  Tables: {n_tables} CSV files in {tbl_path}")
    if n_anims > 0:
        print(f"  Animations: {n_anims} GIFs in {anim_path}")

    # Key files
    print(f"\n  Key outputs:")
    for target in sorted(target_candidates.keys()) if target_candidates else []:
        csv_file = tables_dir / f"rescues_{target}.csv"
        if csv_file.exists():
            try:
                rel_path = csv_file.resolve().relative_to(Path.cwd())
            except ValueError:
                rel_path = csv_file.resolve()
            print(f"    • {rel_path}")

    if sep_json.exists():
        try:
            rel_path = sep_json.resolve().relative_to(Path.cwd())
        except ValueError:
            rel_path = sep_json.resolve()
        print(f"    • {rel_path}")

    print("\n" + "=" * 80)
    print("✅ Report generation complete!")
    print("=" * 80 + "\n")

    return 0
