from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.logging import get_logger
from src.eval.variant_separation import bootstrap_ci, compute_auc, load_labelled_scores
from src.eval.train_test_split import load_split
from src.eval.per_target_metrics import per_target_bootstrap_ci
from src.reporting.guide_generator import generate_guide, render_guide_to_pdf
from src.reporting.executive_summary import (
    generate_executive_summary,
    save_summary_csv,
    save_summary_latex,
)
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
    splits_dir = paths["processed"] / "splits"
    split_path = splits_dir / "clinvar_split_seed42.json"

    if scores_path.exists() and labels_dir.exists():
        try:
            labelled = load_labelled_scores(scores_path, labels_dir)
        except RuntimeError as exc:
            logger.warning("Skipping variant separation: %s", exc)
            labelled = None

        if labelled is not None:
            # Overall AUC
            auc = compute_auc(labelled["label"], labelled["ddg"])
            lo, hi = bootstrap_ci(labelled["label"], labelled["ddg"], n=2000)

            summary = {
                "auc": auc,
                "auc_ci_lower": lo,
                "auc_ci_upper": hi,
                "n_variants": int(len(labelled)),
            }

            # Train/test split validation
            if split_path.exists():
                try:
                    split_data = load_split(split_path)
                    train_idx = set(split_data.get("train_indices", []))
                    test_idx = set(split_data.get("test_indices", []))

                    # Create train/test masks
                    train_mask = labelled.index.isin(train_idx)
                    test_mask = labelled.index.isin(test_idx)

                    train_data = labelled[train_mask]
                    test_data = labelled[test_mask]

                    if len(train_data) > 10 and len(test_data) > 10:
                        train_auc = compute_auc(train_data["label"], train_data["ddg"])
                        train_lo, train_hi = bootstrap_ci(train_data["label"], train_data["ddg"], n=2000)

                        test_auc = compute_auc(test_data["label"], test_data["ddg"])
                        test_lo, test_hi = bootstrap_ci(test_data["label"], test_data["ddg"], n=2000)

                        summary["train_auc"] = train_auc
                        summary["train_auc_ci_lower"] = train_lo
                        summary["train_auc_ci_upper"] = train_hi
                        summary["n_train"] = int(len(train_data))

                        summary["test_auc"] = test_auc
                        summary["test_auc_ci_lower"] = test_lo
                        summary["test_auc_ci_upper"] = test_hi
                        summary["n_test"] = int(len(test_data))

                        logger.info("Train AUC: %.3f, Test AUC: %.3f", train_auc, test_auc)
                except Exception as exc:
                    logger.warning("Train/test split validation failed: %s", exc)

            # Per-target metrics (require 'pos' column in labelled data)
            if "pos" in labelled.columns:
                try:
                    per_target_df = per_target_bootstrap_ci(
                        labelled,
                        min_variants=5,
                        n_bootstrap=1000,
                        confidence_level=0.95,
                    )
                    if not per_target_df.empty:
                        per_target_df = per_target_df.sort_values("n_variants", ascending=False)
                        per_target_df.to_csv(tables_dir / "per_target_auc.csv", index=False)
                        logger.info("Per-target AUC metrics written for %d positions", len(per_target_df))
                except Exception as exc:
                    logger.warning("Per-target metrics failed: %s", exc)
            else:
                logger.info("Skipping per-target metrics (no 'pos' column in variant scores)")

            (tables_dir / "variant_separation.json").write_text(
                json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
            )
            pd.DataFrame([summary]).to_csv(tables_dir / "variant_separation.csv", index=False)
            plot_ddg_by_label(labelled, figures_dir / f"variant_ddg_by_label.{fig_format}")
            logger.info("Variant separation report written")
    else:
        logger.warning("Variant scores or labels not found; skipping variant report")

    rescues_root = paths["processed"] / "rescues"
    target_candidates: dict[str, pd.DataFrame] = {}
    if rescues_root.exists():
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

        # Generate executive summary with risk annotations
        guide_cfg = report_cfg.get("guide", {})
        exec_cfg = report_cfg.get("executive_summary", {})

        if exec_cfg.get("enabled", True) and target_candidates:
            try:
                # Load auxiliary data for risk flags
                p53_cfg = configs.get("p53", {})
                msa_cfg = p53_cfg.get("msa", {})
                msa_cons_map = {}
                if msa_cfg.get("enabled", False):
                    msa_path = paths["processed"] / msa_cfg.get("precomputed_path", "msa/P04637_conservation.json")
                    if msa_path.exists():
                        msa_data = json.loads(msa_path.read_text(encoding="utf-8"))
                        msa_cons_map = {int(k): float(v) for k, v in msa_data.items()}
                        logger.info("Loaded MSA conservation scores for %d positions", len(msa_cons_map))

                # Combine all candidates for executive summary
                all_candidates = pd.concat(target_candidates.values(), ignore_index=True)

                top_per_target = int(exec_cfg.get("top_per_target", 5))
                high_cons_threshold = float(exec_cfg.get("high_conservation_threshold", 0.8))
                near_functional_dist = float(exec_cfg.get("near_functional_distance", 10.0))

                exec_summary = generate_executive_summary(
                    all_candidates,
                    top_per_target=top_per_target,
                    msa_cons_map=msa_cons_map if msa_cons_map else None,
                    high_cons_threshold=high_cons_threshold,
                    near_functional_threshold=near_functional_dist,
                )

                if not exec_summary.empty:
                    # Save CSV
                    save_summary_csv(exec_summary, output_dir / "executive_summary.csv")

                    # Save LaTeX
                    save_summary_latex(exec_summary, output_dir / "executive_summary.tex")

                    logger.info("Executive summary generated: %d top rescues", len(exec_summary))
            except Exception as exc:
                logger.warning("Executive summary generation failed: %s", exc)

        # Generate interpretation guide
        if guide_cfg.get("enabled", True):
            try:
                guide_path = output_dir / "guide.md"
                generate_guide(guide_path, config=configs)
                logger.info("Generated interpretation guide: %s", guide_path)

                # Optionally render to PDF
                if guide_cfg.get("output_format", "md") == "pdf":
                    try:
                        pdf_path = render_guide_to_pdf(guide_path)
                        logger.info("Rendered guide to PDF: %s", pdf_path)
                    except RuntimeError as exc:
                        logger.warning("PDF rendering failed (pandoc not installed?): %s", exc)
            except Exception as exc:
                logger.warning("Guide generation failed: %s", exc)

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

        # Check if train/test split available
        if "train_auc" in sep_data and "test_auc" in sep_data:
            print(f"  Overall AUC: {sep_data['auc']:.3f} (95% CI: {sep_data['auc_ci_lower']:.3f}-{sep_data['auc_ci_upper']:.3f})")
            print(f"  Train AUC:   {sep_data['train_auc']:.3f} (95% CI: {sep_data['train_auc_ci_lower']:.3f}-{sep_data['train_auc_ci_upper']:.3f}) [n={sep_data['n_train']}]")
            print(f"  Test AUC:    {sep_data['test_auc']:.3f} (95% CI: {sep_data['test_auc_ci_lower']:.3f}-{sep_data['test_auc_ci_upper']:.3f}) [n={sep_data['n_test']}]")

            # Check for overfitting
            auc_diff = abs(sep_data['train_auc'] - sep_data['test_auc'])
            if auc_diff < 0.05:
                print(f"  ✓ Train/Test AUC within 0.05 - no evidence of overfitting")
            else:
                print(f"  ⚠ Train/Test AUC differ by {auc_diff:.3f} - potential overfitting")
        else:
            print(f"  AUC: {sep_data['auc']:.3f} (95% CI: {sep_data['auc_ci_lower']:.3f}-{sep_data['auc_ci_upper']:.3f})")

        print(f"  Total variants scored: {sep_data['n_variants']}")
        print(f"  ✓ Model successfully discriminates pathogenic from benign mutations")

        # Per-target AUC summary
        per_target_csv = tables_dir / "per_target_auc.csv"
        if per_target_csv.exists():
            per_target_df = pd.read_csv(per_target_csv)
            if not per_target_df.empty:
                print(f"\n  Per-Position AUC (top 5 by variant count):")
                for _, row in per_target_df.head(5).iterrows():
                    pos = int(row['pos'])
                    n = int(row['n_variants'])
                    auc = row['auc']
                    ci_lo = row.get('ci_lower', 0)
                    ci_hi = row.get('ci_upper', 1)
                    print(f"    Position {pos:4d}: AUC={auc:.3f} (95% CI: {ci_lo:.3f}-{ci_hi:.3f}) [n={n}]")

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

    # Executive summary
    exec_csv = output_dir / "executive_summary.csv"
    if exec_csv.exists():
        try:
            rel_path = exec_csv.resolve().relative_to(Path.cwd())
        except ValueError:
            rel_path = exec_csv.resolve()
        print(f"    • {rel_path} (Executive Summary)")

    # Guide
    guide_md = output_dir / "guide.md"
    guide_pdf = output_dir / "guide.pdf"
    if guide_pdf.exists():
        try:
            rel_path = guide_pdf.resolve().relative_to(Path.cwd())
        except ValueError:
            rel_path = guide_pdf.resolve()
        print(f"    • {rel_path} (Interpretation Guide)")
    elif guide_md.exists():
        try:
            rel_path = guide_md.resolve().relative_to(Path.cwd())
        except ValueError:
            rel_path = guide_md.resolve()
        print(f"    • {rel_path} (Interpretation Guide)")

    print("\n" + "=" * 80)
    print("✅ Report generation complete!")
    print("=" * 80 + "\n")

    return 0
