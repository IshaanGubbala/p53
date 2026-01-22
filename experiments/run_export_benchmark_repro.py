"""Export benchmark reproducibility artifacts.

Saves the exact variant lists, scores, and metadata needed to reproduce
the variant separation benchmark (AUC = 0.844).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.logging import get_logger
from src.eval.variant_separation import load_labelled_scores


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
    logs_dir = output_dir / "logs"

    tables_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    scores_path = paths["processed"] / "variant_scores.parquet"
    labels_dir = paths["processed"] / "labels"

    # Load raw label sets (before merging with scores)
    benign = pd.read_parquet(labels_dir / "benign.parquet")
    pathogenic = pd.read_parquet(labels_dir / "pathogenic.parquet")

    # Create labeled input dataset
    benign_input = benign.copy()
    pathogenic_input = pathogenic.copy()
    benign_input["label"] = "benign"
    pathogenic_input["label"] = "pathogenic"

    input_variants = pd.concat([benign_input, pathogenic_input], ignore_index=True)
    input_variants = input_variants.sort_values(["pos", "ref", "alt"])

    # Save input variants (before scoring)
    input_path = tables_dir / "variant_benchmark_input.csv"
    input_variants.to_csv(input_path, index=False)
    logger.info("Saved benchmark input variants to %s", input_path)

    # Load scored and labeled variants (after merging)
    labelled = load_labelled_scores(scores_path, labels_dir)
    labelled["label_name"] = labelled["label"].map({0: "benign", 1: "pathogenic"})

    # Save scored variants
    scored_path = tables_dir / "variant_benchmark_scored.csv"
    scored_cols = ["pos", "ref", "alt", "mutation", "ddg", "label", "label_name"]
    labelled[scored_cols].to_csv(scored_path, index=False)
    logger.info("Saved benchmark scored variants to %s", scored_path)

    # Record metadata
    scores = pd.read_parquet(scores_path)

    # Count unique variants in labels (before conflict removal)
    benign_positions = set(zip(benign["pos"], benign["ref"], benign["alt"]))
    pathogenic_positions = set(zip(pathogenic["pos"], pathogenic["ref"], pathogenic["alt"]))
    overlapping = benign_positions & pathogenic_positions

    metadata = {
        "dataset_info": {
            "total_variants_scored": int(len(scores)),
            "unique_mutations_scored": int(scores["mutation"].nunique()),
            "benign_labels_raw": int(len(benign)),
            "pathogenic_labels_raw": int(len(pathogenic)),
            "conflicting_variants_removed": int(len(overlapping)),
            "final_benchmark_size": int(len(labelled)),
            "final_benign": int((labelled["label"] == 0).sum()),
            "final_pathogenic": int((labelled["label"] == 1).sum()),
        },
        "filters": {
            "gene": "TP53",
            "consequence": "missense_variant",
            "label_conflict_resolution": "remove_conflicting",
            "duplicate_resolution": "average_ddg",
        },
        "random_seeds": {
            "bootstrap_ci": 1337,
            "n_bootstrap_samples": 2000,
        },
        "files": {
            "input": str(input_path.relative_to(output_dir)),
            "scored": str(scored_path.relative_to(output_dir)),
            "scores_parquet": str(scores_path),
            "labels_benign": str(labels_dir / "benign.parquet"),
            "labels_pathogenic": str(labels_dir / "pathogenic.parquet"),
        },
    }

    metadata_path = logs_dir / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    logger.info("Saved benchmark metadata to %s", metadata_path)

    # Print summary
    print("\n=== Benchmark Reproducibility Summary ===")
    print(f"Input variants: {len(input_variants)} ({len(benign)} benign + {len(pathogenic)} pathogenic)")
    print(f"Conflicting labels removed: {len(overlapping)}")
    print(f"Final benchmark size: {len(labelled)} ({(labelled['label'] == 0).sum()} benign + {(labelled['label'] == 1).sum()} pathogenic)")
    print(f"\nArtifacts saved:")
    print(f"  - {input_path}")
    print(f"  - {scored_path}")
    print(f"  - {metadata_path}")

    return 0
