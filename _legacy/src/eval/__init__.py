"""Evaluation utilities for variant scoring and rescue designs."""

from .variant_separation import (
    load_labelled_scores,
    compute_auc,
    bootstrap_ci,
)

from .per_target_metrics import (
    compute_per_target_auc,
    per_target_bootstrap_ci,
    filter_hotspot_positions,
    compare_train_test_per_target,
    identify_challenging_positions,
)

from .train_test_split import (
    stratified_train_test_split,
    create_split_indices,
    save_split,
    load_split,
    apply_split,
)

from .druggability_analysis import (
    calculate_pocket_score,
    analyze_rescue_druggability,
    summarize_druggability,
)

__all__ = [
    # Variant separation
    "load_labelled_scores",
    "compute_auc",
    "bootstrap_ci",
    # Per-target metrics
    "compute_per_target_auc",
    "per_target_bootstrap_ci",
    "filter_hotspot_positions",
    "compare_train_test_per_target",
    "identify_challenging_positions",
    # Train/test splitting
    "stratified_train_test_split",
    "create_split_indices",
    "save_split",
    "load_split",
    "apply_split",
    # Druggability analysis
    "calculate_pocket_score",
    "analyze_rescue_druggability",
    "summarize_druggability",
]
