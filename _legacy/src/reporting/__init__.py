"""Reporting utilities for generating analysis reports and summaries."""

from .guide_generator import generate_guide, render_guide_to_pdf
from .risk_annotator import compute_risk_flags, format_risk_flags
from .executive_summary import (
    generate_executive_summary,
    save_summary_csv,
    save_summary_latex,
)

__all__ = [
    "generate_guide",
    "render_guide_to_pdf",
    "compute_risk_flags",
    "format_risk_flags",
    "generate_executive_summary",
    "save_summary_csv",
    "save_summary_latex",
]
