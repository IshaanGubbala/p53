"""Result artifact schema, storage, and campaign helpers."""

from p53cad.results.schema import (
    BIG8_HOTSPOTS,
    DEFAULT_DELIVERY_METHODS,
    SCHEMA_VERSION,
    build_run_id,
    build_scenario_matrix,
    select_presentation_shortlist,
)
from p53cad.results.store import CampaignStore
from p53cad.results.presentation import hydrate_display_candidates
from p53cad.results.visualization import (
    build_candidate_position_heatmap,
    build_loss_mesh,
    build_multi_metric_frame,
)

__all__ = [
    "BIG8_HOTSPOTS",
    "DEFAULT_DELIVERY_METHODS",
    "SCHEMA_VERSION",
    "build_run_id",
    "build_scenario_matrix",
    "select_presentation_shortlist",
    "CampaignStore",
    "hydrate_display_candidates",
    "build_candidate_position_heatmap",
    "build_multi_metric_frame",
    "build_loss_mesh",
]
