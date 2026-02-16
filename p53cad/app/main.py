from p53cad.core.runtime import bootstrap_runtime, get_runtime_capabilities, log_runtime_capabilities
from p53cad.core.logging import setup_logging

bootstrap_runtime()
setup_logging()
_RUNTIME_CAPS = log_runtime_capabilities("p53cad.app.runtime")

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from pathlib import Path
import torch
import torch.nn.functional as F
import requests
import time
import json
import logging
import io
import zipfile
from typing import Any, Dict, Optional

import importlib
from p53cad.engine.latent import ManifoldEmbedder, ManifoldWalker
from p53cad.engine.oracle import FunctionalOracle
from p53cad.engine.explain import SaliencyMap
from p53cad.analysis.grassmann import GrassmannMetric
from p53cad.data import dms as dms_module
from p53cad.viz.pymol import PyMolGenerator
from p53cad.results.store import CampaignStore
from p53cad.results.presentation import hydrate_display_candidates
from p53cad.results.visualization import (
    build_candidate_position_heatmap,
    build_multi_metric_frame,
    build_loss_mesh,
)

# DMS core helpers (required)
P53_WT = dms_module.P53_WT
apply_mutation = dms_module.apply_mutation
get_dms_data = dms_module.get_dms_data
parse_single_mutation = dms_module.parse_single_mutation

# DMS optional helpers (guarded)
physics_based_score_fn = getattr(dms_module, "physics_based_score", None)
get_coverage_report_fn = getattr(dms_module, "get_coverage_report", None)

logger = logging.getLogger(__name__)
if not callable(physics_based_score_fn):
    available_exports = sorted(name for name in dir(dms_module) if not name.startswith("_"))
    logger.warning(
        "Optional DMS helper 'physics_based_score' unavailable in module %s. "
        "Validation will run in DMS-only fail-closed mode. "
        "Available exports: %s",
        getattr(dms_module, "__file__", "unknown"),
        ", ".join(available_exports),
    )

# Optional enhancement modules loaded independently so one failure
# does not disable every analysis section.
MODULE_IMPORT_ERRORS = {}

DrugGeneratorEngine = None
P53_BINDING_POCKETS = {}
ExplainabilityEngine = None
P53_FUNCTIONAL_SITES = {}
ExperimentalPipeline = None
ClinicalImpactEngine = None
PatientStratifier = None
MultiTargetPlatform = None
TUMOR_SUPPRESSORS = {}
MDValidationEngine = None
ExplainabilityDependencyError = RuntimeError

try:
    from p53cad.engine.drug_generator import DrugGeneratorEngine, P53_BINDING_POCKETS
except Exception as e:
    MODULE_IMPORT_ERRORS["drug"] = str(e)

try:
    from p53cad.engine.explainability import (
        ExplainabilityDependencyError,
        ExplainabilityEngine,
        P53_FUNCTIONAL_SITES,
    )
except Exception as e:
    MODULE_IMPORT_ERRORS["explainability"] = str(e)

try:
    from p53cad.analysis.clinical_impact import ClinicalImpactEngine, PatientStratifier
except Exception as e:
    MODULE_IMPORT_ERRORS["clinical"] = str(e)

try:
    from p53cad.engine.experimental import ExperimentalPipeline
except Exception as e:
    MODULE_IMPORT_ERRORS["experimental"] = str(e)

try:
    from p53cad.engine.multi_target import MultiTargetPlatform, TUMOR_SUPPRESSORS
except Exception as e:
    MODULE_IMPORT_ERRORS["multi_target"] = str(e)

try:
    from p53cad.engine.md_validation import MDValidationEngine
except Exception as e:
    MODULE_IMPORT_ERRORS["md_validation"] = str(e)

DRUG_MODULE_AVAILABLE = DrugGeneratorEngine is not None and bool(P53_BINDING_POCKETS)
EXPLAINABILITY_MODULE_AVAILABLE = ExplainabilityEngine is not None
CLINICAL_MODULE_AVAILABLE = ClinicalImpactEngine is not None
ENHANCED_MODULES_AVAILABLE = any(
    [
        DRUG_MODULE_AVAILABLE,
        EXPLAINABILITY_MODULE_AVAILABLE,
        CLINICAL_MODULE_AVAILABLE,
        ExperimentalPipeline is not None,
        MultiTargetPlatform is not None,
        MDValidationEngine is not None,
    ]
)

if MODULE_IMPORT_ERRORS:
    logging.getLogger(__name__).warning("Optional module import issues: %s", MODULE_IMPORT_ERRORS)


@st.cache_data(show_spinner=False)
def get_runtime_capabilities_snapshot() -> Dict[str, Any]:
    """Return lightweight runtime capability probe for dashboard diagnostics."""
    return get_runtime_capabilities()


def get_explainability_backend_status() -> Dict[str, Any]:
    """
    Probe explainability attention availability once per session.
    Fail-closed: if attention tensors cannot be produced, mark unavailable.
    """
    cached = st.session_state.get("_attention_backend_status")
    if isinstance(cached, dict):
        return cached

    status = {"ready": False, "reason": "Explainability module not loaded."}
    if not EXPLAINABILITY_MODULE_AVAILABLE:
        st.session_state["_attention_backend_status"] = status
        return status

    embedder_obj = globals().get("embedder")
    oracle_obj = globals().get("oracle")
    if embedder_obj is None or oracle_obj is None:
        status["reason"] = "Oracle/embedder model stack is unavailable."
        st.session_state["_attention_backend_status"] = status
        return status

    try:
        engine = ExplainabilityEngine(
            model=embedder_obj.model,
            tokenizer=embedder_obj.tokenizer,
            oracle=oracle_obj,
            embedder=embedder_obj,
        )
        probe_seq = P53_WT[:80]
        engine.attention_extractor.extract_attention(probe_seq)
        status = {"ready": True, "reason": "Attention tensors available."}
    except ExplainabilityDependencyError as err:
        status = {"ready": False, "reason": str(err)}
    except Exception as err:  # pragma: no cover - defensive runtime path
        status = {"ready": False, "reason": f"Runtime probe failed: {err}"}

    st.session_state["_attention_backend_status"] = status
    return status


def get_drug_mode_capability(
    pocket_key: str,
    method: str,
    allow_wt_receptor_fallback: bool = False,
) -> Dict[str, Any]:
    """
    Probe per-mode drug readiness for strict UI gating.
    """
    if not DRUG_MODULE_AVAILABLE:
        return {"ready": False, "reason": MODULE_IMPORT_ERRORS.get("drug", "Drug module unavailable.")}
    try:
        engine = DrugGeneratorEngine()
        return engine.get_mode_capabilities(
            pocket_key=pocket_key,
            method=method,
            allow_wt_receptor_fallback=allow_wt_receptor_fallback,
        )
    except Exception as err:
        return {"ready": False, "reason": str(err)}


# === LIVE SIMULATION FUNCTIONS ===
def predict_structure_esmfold(sequence: str) -> str:
    """Call ESMFold API to predict structure. Returns PDB string."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=120)
        if response.status_code == 200:
            return response.text
        else:
            return None
    except Exception as e:
        return None


def get_real_structure_checkpoint(sequence: str) -> Optional[str]:
    """
    Get a real structure checkpoint from ESMFold API with per-session cache.
    """
    if not sequence:
        return None
    cache = st.session_state.setdefault("_esmfold_checkpoint_cache", {})
    if sequence in cache:
        return cache[sequence]

    pdb = predict_structure_esmfold(sequence)
    if pdb and "ATOM" in pdb:
        cache[sequence] = pdb
        return pdb
    return None

def get_optimization_stage(step: float, total_steps: int) -> dict:
    """Return stage metadata for the live optimization timeline."""
    progress = float(step) / max(float(total_steps), 1.0)
    if progress < 0.25:
        return {
            "name": "Stage 1: Exploration",
            "description": "Sampling broad structural alternatives.",
            "color": "#0EA5E9",
        }
    if progress < 0.55:
        return {
            "name": "Stage 2: Growth",
            "description": "Amplifying beneficial mutation motifs.",
            "color": "#2563EB",
        }
    if progress < 0.82:
        return {
            "name": "Stage 3: Refinement",
            "description": "Pruning unstable and low-signal changes.",
            "color": "#10B981",
        }
    return {
        "name": "Stage 4: Convergence",
        "description": "Locking in the final rescue architecture.",
        "color": "#F59E0B",
    }


def compute_camera_focus_positions(
    current_positions: list,
    previous_positions: Optional[list] = None,
    max_points: int = 8,
) -> list:
    """Prioritize newly changing residues for live camera focus."""
    current = sorted({int(p) for p in current_positions if p})
    if not current:
        return []

    previous = sorted({int(p) for p in (previous_positions or []) if p})
    current_set = set(current)
    previous_set = set(previous)

    # Focus first on residues that changed since the last viewer refresh.
    changed = sorted((current_set - previous_set) | (previous_set - current_set))
    if changed:
        return changed[:max_points]
    return current[:max_points]

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="p53-proteoMgCAD",
    page_icon="https://em-content.zobj.net/source/apple/391/dna_1f9ec.png",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ============================================================================
# CLEAN, MODERN DESIGN - Professional Scientific Application
# ============================================================================
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Manrope:wght@400;500;600;700;800&family=Sora:wght@500;600;700;800&display=swap');

    :root {
        --bg-base: #EEF4FA;
        --bg-tint-1: #D8E7FF;
        --bg-tint-2: #D2F5EA;
        --bg-tint-3: #DDF1FF;
        --surface: #F9FCFF;
        --surface-muted: #EEF5FF;
        --surface-strong: #E5F0FF;
        --border-soft: #CEDDF2;
        --border-strong: #9EB8E0;
        --text-main: #0F213F;
        --text-muted: #4A5D7C;
        --accent: #2563EB;
        --accent-strong: #1E40AF;
        --accent-soft: #DBEAFE;
        --accent-2: #0891B2;
        --accent-3: #0F766E;
        --accent-dark: #0A1736;
    }

    html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"] {
        background:
            radial-gradient(1000px 620px at 12% -10%, var(--bg-tint-1), transparent 65%),
            radial-gradient(860px 620px at 90% -20%, var(--bg-tint-2), transparent 63%),
            radial-gradient(640px 460px at 48% 2%, var(--bg-tint-3), transparent 62%),
            linear-gradient(180deg, #F5F9FF 0%, var(--bg-base) 100%) !important;
        font-family: 'Manrope', -apple-system, BlinkMacSystemFont, sans-serif !important;
        color: var(--text-main) !important;
    }

    [data-testid="stAppViewContainer"] {
        position: relative;
        overflow-x: hidden;
    }

    [data-testid="stAppViewContainer"]::before {
        content: "";
        position: fixed;
        right: -180px;
        top: 14vh;
        width: 420px;
        height: 420px;
        border-radius: 999px;
        background: radial-gradient(circle, rgba(56, 189, 248, 0.22), rgba(56, 189, 248, 0.0) 68%);
        pointer-events: none;
        filter: blur(1px);
        animation: floatOrb 9s ease-in-out infinite;
        z-index: 0;
    }

    [data-testid="stAppViewContainer"]::after {
        content: "";
        position: fixed;
        left: -180px;
        bottom: 5vh;
        width: 360px;
        height: 360px;
        border-radius: 999px;
        background: radial-gradient(circle, rgba(37, 99, 235, 0.18), rgba(37, 99, 235, 0.0) 70%);
        pointer-events: none;
        animation: floatOrb 11s ease-in-out infinite reverse;
        z-index: 0;
    }

    @keyframes floatOrb {
        0% { transform: translateY(0px) translateX(0px); }
        50% { transform: translateY(-16px) translateX(10px); }
        100% { transform: translateY(0px) translateX(0px); }
    }

    @keyframes fadeLift {
        from { opacity: 0; transform: translateY(8px); }
        to { opacity: 1; transform: translateY(0); }
    }

    .main .block-container {
        max-width: 1320px !important;
        padding: 2.05rem 2.7rem 3.25rem 2.7rem !important;
        position: relative;
        z-index: 1;
        animation: fadeLift 0.45s ease-out both;
    }

    #MainMenu, footer, header, .stDeployButton {display: none !important; visibility: hidden !important;}

    html, body {
        font-size: 18px !important;
    }

    h1, h2, h3, h4, h5, h6, .stMarkdown h1, .stMarkdown h2, .stMarkdown h3 {
        font-family: 'Sora', sans-serif !important;
        color: var(--text-main) !important;
        font-weight: 700 !important;
        letter-spacing: -0.02em !important;
    }

    h1, .stMarkdown h1 { font-size: 2.74rem !important; }
    h2, .stMarkdown h2 { font-size: 1.96rem !important; }
    h3, .stMarkdown h3 { font-size: 1.42rem !important; }
    h4, .stMarkdown h4 { font-size: 1.16rem !important; }

    p, span, label, .stMarkdown p, .stMarkdown li {
        font-family: 'Manrope', sans-serif !important;
        color: var(--text-muted) !important;
        line-height: 1.68 !important;
        font-size: 1.08rem !important;
    }

    .hero-container {
        text-align: center;
        padding: 2.6rem 1.6rem 1.85rem 1.6rem;
        margin-bottom: 1.9rem;
        background: linear-gradient(130deg, #FFFFFF 0%, #ECF4FF 52%, #D9EEFF 100%);
        border: 1px solid var(--border-soft);
        border-radius: 20px;
        box-shadow: 0 20px 42px rgba(17, 47, 106, 0.16), 0 6px 14px rgba(7, 23, 63, 0.08);
        position: relative;
        overflow: hidden;
    }

    .hero-container::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        height: 6px;
        width: 100%;
        background: linear-gradient(90deg, var(--accent-dark), var(--accent), #38BDF8, var(--accent-3));
    }

    .hero-container::after {
        content: "";
        position: absolute;
        right: -50px;
        top: -80px;
        width: 250px;
        height: 250px;
        background: radial-gradient(circle at 40% 45%, rgba(34, 211, 238, 0.42), transparent 68%);
        filter: blur(2px);
    }

    .hero-badge {
        display: inline-block;
        background: linear-gradient(140deg, var(--accent-dark) 0%, var(--accent) 100%);
        color: #E7F1FF !important;
        font-size: 0.79rem !important;
        font-weight: 700 !important;
        padding: 0.46rem 1.1rem;
        border-radius: 999px;
        margin-bottom: 0.75rem;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        box-shadow: 0 11px 24px rgba(10, 23, 54, 0.3);
    }

    .hero-title {
        font-family: 'Sora', sans-serif !important;
        font-size: 3.08rem !important;
        font-weight: 800 !important;
        color: var(--accent-dark) !important;
        letter-spacing: -0.034em;
        margin: 0 0 0.5rem 0 !important;
        line-height: 1.15;
    }

    .hero-subtitle {
        font-size: 1.14rem !important;
        color: #324A70 !important;
        max-width: 780px;
        margin: 0 auto 0.5rem auto;
        line-height: 1.48;
        font-weight: 600 !important;
    }

    .hero-description {
        font-size: 1.08rem !important;
        color: #4F678D !important;
        max-width: 620px;
        margin: 0 auto;
    }

    .studio-banner {
        padding: 1.2rem 1.35rem;
        border: 1px solid var(--border-soft);
        border-radius: 16px;
        background: linear-gradient(135deg, #F6FBFF 0%, #EBF4FF 60%, #DEF7F2 100%);
        box-shadow: 0 14px 24px rgba(37, 99, 235, 0.1);
        margin-bottom: 1rem;
        position: relative;
        overflow: hidden;
    }

    .studio-banner::before {
        content: "";
        position: absolute;
        left: 0;
        top: 0;
        height: 100%;
        width: 5px;
        background: linear-gradient(180deg, var(--accent-dark), var(--accent), var(--accent-2));
    }

    .studio-kicker {
        margin: 0 0 0.25rem 0 !important;
        text-transform: uppercase;
        letter-spacing: 0.08em;
        font-size: 0.72rem !important;
        font-weight: 800 !important;
        color: #315586 !important;
    }

    .studio-subtitle {
        margin: 0 !important;
        color: #406089 !important;
        font-size: 0.96rem !important;
    }

    .lovable-card {
        background: linear-gradient(180deg, #FFFFFF 0%, var(--surface) 100%);
        border: 1px solid var(--border-soft);
        border-radius: 16px;
        padding: 1.28rem;
        margin-bottom: 1rem;
        box-shadow: 0 12px 28px rgba(28, 66, 133, 0.11), 0 2px 8px rgba(15, 23, 42, 0.06);
        position: relative;
        overflow: hidden;
        transition: transform 0.22s ease, box-shadow 0.22s ease;
    }

    .lovable-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 16px 32px rgba(28, 66, 133, 0.15), 0 5px 12px rgba(15, 23, 42, 0.08);
    }

    .lovable-card::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        height: 4px;
        width: 100%;
        background: linear-gradient(90deg, var(--accent-dark), var(--accent), #38BDF8);
    }

    .lovable-card h2, .lovable-card h3, .lovable-card h4 {
        color: var(--text-main) !important;
        margin: 0 0 0.45rem 0 !important;
        font-weight: 700 !important;
    }

    .lovable-card p {
        color: var(--text-muted) !important;
        margin: 0 !important;
    }

    .lovable-card-accent {
        background: linear-gradient(135deg, #EAF3FF 0%, #DDF7FF 55%, #D9F8EE 100%);
        border: 1px solid #8EB9EB;
        border-left: 5px solid var(--accent-dark);
        border-radius: 14px;
        padding: 1.2rem;
        margin-bottom: 1rem;
        box-shadow: 0 12px 24px rgba(30, 85, 166, 0.13);
    }

    .lovable-card-accent h3 {
        color: #133966 !important;
        margin: 0 0 0.28rem 0 !important;
        font-size: 1.02rem !important;
        font-weight: 700 !important;
    }

    .lovable-card-accent p {
        color: #2F5A84 !important;
        font-size: 0.87rem !important;
    }

    .data-highlight {
        background: linear-gradient(120deg, #E1F1FF 0%, #D4FAF3 100%);
        border: 1px solid #56A4E8;
        border-radius: 10px;
        padding: 0.74rem 1rem;
        font-weight: 700;
        color: #124A78;
    }

    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #F8FCFF 0%, #EFF6FF 100%) !important;
        border-right: 1px solid var(--border-soft) !important;
        backdrop-filter: blur(8px);
    }

    [data-testid="stSidebar"] [data-testid="stMarkdown"] {
        color: #2F4A70 !important;
    }

    .section-header {
        font-size: 0.72rem !important;
        font-weight: 800 !important;
        color: #6A80A4 !important;
        text-transform: uppercase !important;
        letter-spacing: 0.1em !important;
        margin-bottom: 0.7rem !important;
    }

    .feature-pill {
        display: inline-block;
        background: linear-gradient(120deg, #E7F1FF 0%, #DCF9F5 100%);
        color: #2E4B72;
        font-size: 0.71rem;
        font-weight: 700;
        padding: 0.2rem 0.62rem;
        border-radius: 999px;
        margin-right: 0.4rem;
        border: 1px solid #C1D8F4;
    }

    .stTabs [data-baseweb="tab-list"] {
        background: linear-gradient(135deg, #EAF3FF 0%, #E2EEFF 100%) !important;
        border-radius: 14px !important;
        padding: 6px !important;
        gap: 6px !important;
        border: 1px solid #D4E4FB !important;
    }

    .stTabs [data-baseweb="tab"] {
        background: transparent !important;
        border-radius: 10px !important;
        color: #3D567A !important;
        font-weight: 700 !important;
        font-size: 1.08rem !important;
        padding: 0.86rem 1.32rem !important;
        border: none !important;
        transition: all 0.15s ease !important;
    }

    .stTabs [data-baseweb="tab"]:hover {
        color: #1F3F70 !important;
        background: rgba(255, 255, 255, 0.45) !important;
    }

    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, var(--accent-dark) 0%, var(--accent) 100%) !important;
        color: #ECF4FF !important;
        font-weight: 800 !important;
        box-shadow: 0 8px 18px rgba(30, 64, 175, 0.3) !important;
    }

    .stTabs [data-baseweb="tab-highlight"] {
        display: none !important;
    }

    .stRadio > div {
        flex-direction: row !important;
        flex-wrap: wrap !important;
        gap: 0.55rem !important;
    }

    .stRadio > div > label {
        background: var(--surface) !important;
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 10px !important;
        padding: 0.62rem 1.15rem !important;
        margin: 0 !important;
        font-size: 1.04rem !important;
        color: #34547D !important;
        cursor: pointer !important;
        transition: all 0.16s ease !important;
    }

    .stRadio > div > label:hover {
        border-color: var(--accent) !important;
        background: #F0F6FF !important;
    }

    .stRadio > div > label[data-checked="true"],
    .stRadio > div > label:has(input:checked) {
        background: linear-gradient(135deg, #E6F0FF 0%, #DBF4FF 100%) !important;
        border-color: var(--accent) !important;
        color: #1C4686 !important;
        font-weight: 700 !important;
    }

    .stRadio [role="radiogroup"] > label > div:first-child {
        display: none !important;
    }

    .stButton > button {
        background: linear-gradient(135deg, var(--accent-dark) 0%, var(--accent) 100%) !important;
        color: #ECF4FF !important;
        border: none !important;
        border-radius: 12px !important;
        padding: 0.88rem 1.62rem !important;
        font-weight: 700 !important;
        font-size: 1.05rem !important;
        transition: transform 0.18s ease, box-shadow 0.18s ease !important;
        box-shadow: 0 10px 18px rgba(30, 64, 175, 0.28) !important;
    }

    /* Force high-contrast text/icons inside blue buttons */
    .stButton > button,
    .stButton > button span,
    .stButton > button p,
    .stButton > button div,
    .stButton > button svg {
        color: #F8FBFF !important;
        fill: #F8FBFF !important;
        stroke: #F8FBFF !important;
    }

    .stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0 16px 24px rgba(30, 64, 175, 0.33) !important;
    }

    .stButton > button[kind="primary"] {
        background: linear-gradient(135deg, var(--accent-dark) 0%, var(--accent) 100%) !important;
    }

    .stButton > button[kind="primary"]:hover {
        background: linear-gradient(135deg, #112F77 0%, var(--accent-strong) 100%) !important;
    }

    .stDownloadButton > button {
        background: linear-gradient(180deg, #F9FCFF 0%, #F1F7FF 100%) !important;
        color: #27476F !important;
        border: 1.4px solid var(--border-soft) !important;
        font-size: 1rem !important;
        font-weight: 700 !important;
        padding: 0.62rem 1.2rem !important;
        border-radius: 11px !important;
    }

    .stDownloadButton > button:hover {
        background: #EEF5FF !important;
        border-color: #9BB8E1 !important;
    }

    .stTextInput > div > div > input,
    .stNumberInput input {
        background: #FAFCFF !important;
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 10px !important;
        padding: 0.82rem 1.05rem !important;
        font-size: 1.05rem !important;
        color: var(--text-main) !important;
    }

    .stTextInput > div > div > input:focus,
    .stNumberInput input:focus {
        border-color: var(--accent) !important;
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.18) !important;
    }

    .stSelectbox > div > div,
    .stMultiSelect > div > div {
        background: #FAFCFF !important;
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 10px !important;
        font-size: 1.03rem !important;
    }

    .stMultiSelect [data-baseweb="tag"] {
        background: linear-gradient(120deg, #E3F0FF 0%, #D9F7F1 100%) !important;
        border: 1px solid #B8D2F2 !important;
        border-radius: 8px !important;
        color: #1F4D82 !important;
    }

    [data-testid="stMetric"] {
        background: linear-gradient(180deg, #FFFFFF 0%, #F4F9FF 100%) !important;
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 14px !important;
        padding: 1.22rem 1.24rem !important;
        box-shadow: 0 8px 16px rgba(23, 57, 116, 0.1) !important;
        position: relative;
        overflow: hidden;
    }

    [data-testid="stMetric"]::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 4px;
        background: linear-gradient(90deg, var(--accent), #38BDF8);
    }

    [data-testid="stMetricLabel"] {
        font-size: 0.81rem !important;
        font-weight: 800 !important;
        color: #5D779D !important;
        text-transform: uppercase !important;
        letter-spacing: 0.08em !important;
    }

    [data-testid="stMetricValue"] {
        font-size: 1.98rem !important;
        font-weight: 800 !important;
        color: #0E2A50 !important;
    }

    [data-testid="stMetricDelta"] {
        font-weight: 700 !important;
        font-size: 0.98rem !important;
    }

    [data-testid="stMetricDelta"] svg {
        stroke-width: 2.7px !important;
    }

    [data-testid="stAlert"] > div {
        border-radius: 12px !important;
        padding: 1rem 1.2rem !important;
        font-size: 1.04rem !important;
        border: 1px solid #C9D9EF !important;
    }

    .streamlit-expanderHeader {
        background: linear-gradient(135deg, #F3F8FF 0%, #EAF2FF 100%) !important;
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 12px !important;
        font-weight: 700 !important;
        font-size: 1.07rem !important;
        padding: 0.9rem 1rem !important;
    }

    .stDataFrame {
        border: 1.4px solid var(--border-soft) !important;
        border-radius: 12px !important;
        overflow: hidden !important;
        font-size: 1rem !important;
        box-shadow: 0 8px 16px rgba(19, 53, 103, 0.08) !important;
    }

    hr, .stDivider {
        border: none !important;
        border-top: 1.5px solid #D7E5F8 !important;
        margin: 1.8rem 0 !important;
    }

    .stSlider [data-baseweb="slider"] [role="slider"] {
        background: linear-gradient(135deg, var(--accent-dark), var(--accent)) !important;
        width: 20px !important;
        height: 20px !important;
        box-shadow: 0 0 0 5px rgba(37, 99, 235, 0.14) !important;
    }

    .stSlider [data-testid="stTickBarMin"],
    .stSlider [data-testid="stTickBarMax"] {
        font-size: 0.85rem !important;
        color: #607B9E !important;
    }

    .js-plotly-plot {
        border-radius: 14px !important;
        border: 1px solid var(--border-soft) !important;
        box-shadow: 0 12px 22px rgba(22, 57, 116, 0.09) !important;
        background: linear-gradient(180deg, #FFFFFF 0%, #F7FBFF 100%) !important;
    }

    .js-plotly-plot .plotly {
        border-radius: 14px !important;
    }

    .design-card {
        background: linear-gradient(180deg, #FFFFFF 0%, var(--surface) 100%);
        border: 1.4px solid var(--border-soft);
        border-radius: 15px;
        padding: 1.2rem;
        margin-bottom: 1rem;
        box-shadow: 0 12px 22px rgba(21, 76, 143, 0.1), 0 2px 7px rgba(14, 25, 42, 0.06);
        position: relative;
        overflow: hidden;
        transition: transform 0.18s ease;
    }

    .design-card:hover {
        transform: translateY(-2px);
    }

    .design-card::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        height: 4px;
        width: 100%;
        background: linear-gradient(90deg, var(--accent-dark), var(--accent), #38BDF8);
    }

    .design-card b {
        color: #102947;
        font-weight: 700;
        font-size: 1.03rem;
    }

    .constraint-card {
        background: linear-gradient(180deg, #FFFFFF 0%, #F1F9FF 100%);
        border: 1.4px solid var(--border-soft);
        border-radius: 15px;
        padding: 1.22rem;
        box-shadow: 0 12px 22px rgba(36, 101, 187, 0.11), 0 2px 7px rgba(15, 23, 42, 0.05);
        position: relative;
        overflow: hidden;
    }

    .constraint-card::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        height: 4px;
        width: 100%;
        background: linear-gradient(90deg, var(--accent-dark), var(--accent-2), #22D3EE);
    }

    .accent-box {
        background: linear-gradient(135deg, #EAF3FF 0%, #DAF6FF 60%, #D5F7EE 100%);
        border: 1px solid #A5C4EA;
        border-left: 5px solid var(--accent-dark);
        border-radius: 14px;
        padding: 1rem 1.2rem;
        box-shadow: 0 11px 20px rgba(37, 99, 235, 0.13);
    }

    .candidate-card {
        background: linear-gradient(135deg, #F7FBFF 0%, #ECF5FF 100%);
        border: 1px solid var(--border-soft);
        border-left: 5px solid var(--accent);
        border-radius: 12px;
        padding: 0.9rem 1rem;
        margin-bottom: 0.75rem;
        box-shadow: 0 10px 18px rgba(28, 76, 146, 0.1);
    }

    .candidate-card h4 {
        margin: 0 0 0.35rem 0;
        color: #113156;
        font-size: 1rem;
        font-family: 'Sora', sans-serif !important;
    }

    .candidate-card p {
        margin: 0;
        color: #2E557F !important;
        font-size: 0.88rem !important;
    }

    .gen-design-header {
        font-family: 'Sora', sans-serif !important;
        font-size: 1.86rem;
        font-weight: 800;
        color: #0F2A4D;
        margin: 0 0 0.25rem 0;
        letter-spacing: -0.02em;
    }

    .stSpinner > div {
        border-color: var(--accent) !important;
    }

    .stCheckbox label {
        font-size: 0.87rem !important;
        color: #3A587F !important;
    }

    @media (max-width: 900px) {
        .main .block-container {
            padding: 1rem 1rem 2rem 1rem !important;
        }

        .hero-container {
            padding: 1.65rem 1rem 1.2rem 1rem;
            border-radius: 16px;
        }

        .hero-title {
            font-size: 2.32rem !important;
        }

        .hero-subtitle {
            font-size: 1.11rem !important;
        }

        .stTabs [data-baseweb="tab"] {
            padding: 0.7rem 0.95rem !important;
            font-size: 0.98rem !important;
        }
    }
</style>
""", unsafe_allow_html=True)

# --- PLOTLY VISUAL SYSTEM ---
P53CAD_PLOTLY_TEMPLATE = go.layout.Template(
    layout=go.Layout(
        paper_bgcolor="rgba(248, 252, 255, 0.92)",
        plot_bgcolor="rgba(241, 247, 255, 0.88)",
        colorway=[
            "#1D4ED8", "#0284C7", "#06B6D4", "#14B8A6", "#0F766E", "#F59E0B", "#EF4444"
        ],
        font=dict(family="Manrope, sans-serif", size=14, color="#18385F"),
        title=dict(
            x=0.02, xanchor="left",
            font=dict(family="Sora, sans-serif", size=19, color="#0F2A4D")
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            bgcolor="rgba(255,255,255,0.72)",
            bordercolor="rgba(158, 184, 224, 0.55)",
            borderwidth=1,
            font=dict(size=12, color="#284C79")
        ),
        hoverlabel=dict(
            bgcolor="rgba(10, 23, 54, 0.93)",
            font=dict(color="#ECF4FF", size=12)
        ),
        margin=dict(l=24, r=20, t=52, b=42),
        xaxis=dict(
            showgrid=True,
            gridcolor="rgba(148, 180, 221, 0.35)",
            zeroline=False,
            linecolor="rgba(122, 152, 191, 0.42)",
            tickcolor="rgba(122, 152, 191, 0.42)",
            title_font=dict(color="#1F4878", size=14),
            tickfont=dict(color="#3A5D89", size=12)
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor="rgba(148, 180, 221, 0.35)",
            zeroline=False,
            linecolor="rgba(122, 152, 191, 0.42)",
            tickcolor="rgba(122, 152, 191, 0.42)",
            title_font=dict(color="#1F4878", size=14),
            tickfont=dict(color="#3A5D89", size=12)
        ),
    )
)
pio.templates["p53cad"] = P53CAD_PLOTLY_TEMPLATE
px.defaults.template = "p53cad"
px.defaults.color_discrete_sequence = ["#1D4ED8", "#0284C7", "#06B6D4", "#14B8A6", "#0F766E", "#F59E0B", "#EF4444"]


def style_plotly_figure(fig: go.Figure, chart_kind: str = "default") -> go.Figure:
    """Apply a shared chart style so every panel looks visually consistent."""
    if fig is None:
        return fig

    # Force light surfaces even if an external/default template is dark in the runtime.
    fig.update_layout(
        template="p53cad",
        transition=dict(duration=250, easing="cubic-in-out"),
        paper_bgcolor="rgba(248, 252, 255, 0.97)",
        plot_bgcolor="rgba(241, 247, 255, 0.93)",
    )

    # Consistent larger chart canvases for readability.
    if fig.layout.height is None:
        default_heights = {
            "bar": 500,
            "scatter": 520,
            "pareto": 520,
            "heatmap": 430,
            "trajectory": 420,
            "multiline": 430,
            "mini": 340,
            "3d": 520,
            "default": 460,
        }
        fig.update_layout(height=default_heights.get(chart_kind, 460))

    if chart_kind in {"trajectory", "multiline", "mini"}:
        fig.update_layout(hovermode="x unified")
        fig.update_traces(marker=dict(line=dict(width=1, color="rgba(255,255,255,0.85)")), selector=dict(type="scatter"))
        fig.update_traces(line=dict(width=2.8), selector=dict(type="scatter"))

    if chart_kind in {"scatter", "pareto"}:
        fig.update_layout(hovermode="closest")
        fig.update_traces(marker=dict(line=dict(width=1.2, color="rgba(255,255,255,0.9)")), selector=dict(type="scatter"))

    if chart_kind in {"bar"}:
        fig.update_layout(
            bargap=0.24,
            margin=dict(l=24, r=20, t=64, b=106),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.23,
                xanchor="left",
                x=0,
                bgcolor="rgba(255,255,255,0.84)",
            ),
        )
        fig.update_traces(
            marker=dict(line=dict(width=0.8, color="rgba(255,255,255,0.82)")),
            selector=dict(type="bar")
        )

    if chart_kind in {"scatter", "pareto"}:
        fig.update_layout(
            margin=dict(l=24, r=20, t=64, b=96),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.2,
                xanchor="left",
                x=0,
                bgcolor="rgba(255,255,255,0.84)",
            ),
        )

    if chart_kind in {"heatmap"}:
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
        fig.update_layout(margin=dict(l=10, r=10, t=44, b=16))
        fig.update_traces(
            xgap=1,
            ygap=1,
            selector=dict(type="heatmap")
        )

    if chart_kind in {"3d"}:
        fig.update_traces(
            marker=dict(line=dict(width=0.8, color="rgba(255,255,255,0.85)")),
            selector=dict(type="scatter3d")
        )
        fig.update_layout(
            scene=dict(
                bgcolor="rgba(241,247,255,0.70)",
                xaxis=dict(
                    showbackground=True,
                    backgroundcolor="rgba(223,237,255,0.55)",
                    gridcolor="rgba(146,177,218,0.28)",
                    zeroline=False,
                    color="#2E5888"
                ),
                yaxis=dict(
                    showbackground=True,
                    backgroundcolor="rgba(221,247,241,0.48)",
                    gridcolor="rgba(146,177,218,0.28)",
                    zeroline=False,
                    color="#2E5888"
                ),
                zaxis=dict(
                    showbackground=True,
                    backgroundcolor="rgba(230,240,255,0.60)",
                    gridcolor="rgba(146,177,218,0.28)",
                    zeroline=False,
                    color="#2E5888"
                ),
            )
        )

    # Dynamic margin guard: keep legends and long labels from overlapping plot/title.
    has_visible_legend = any(getattr(trace, "showlegend", True) is not False for trace in fig.data)
    if has_visible_legend:
        current_margin = fig.layout.margin or go.layout.Margin()
        min_top = 64
        min_bottom = 86 if chart_kind in {"bar", "scatter", "pareto", "trajectory", "multiline", "mini"} else 32
        fig.update_layout(
            margin=dict(
                l=max(int(getattr(current_margin, "l", 24) or 24), 20),
                r=max(int(getattr(current_margin, "r", 20) or 20), 16),
                t=max(int(getattr(current_margin, "t", min_top) or min_top), min_top),
                b=max(int(getattr(current_margin, "b", min_bottom) or min_bottom), min_bottom),
            )
        )

    return fig


def _normalize_plotly_kwargs(kwargs: dict) -> dict:
    if "use_container_width" in kwargs:
        use_container = bool(kwargs.pop("use_container_width"))
        kwargs.setdefault("width", "stretch" if use_container else "content")
    return kwargs


def render_plotly(fig: go.Figure, chart_kind: str = "default", **kwargs):
    return st.plotly_chart(style_plotly_figure(fig, chart_kind=chart_kind), theme=None, **_normalize_plotly_kwargs(kwargs))


def render_plotly_in(container, fig: go.Figure, chart_kind: str = "default", **kwargs):
    return container.plotly_chart(style_plotly_figure(fig, chart_kind=chart_kind), theme=None, **_normalize_plotly_kwargs(kwargs))


def safe_physics_based_score(mutation: str) -> Optional[Dict[str, Any]]:
    """Safely execute optional DMS physics fallback scoring."""
    if not callable(physics_based_score_fn):
        return None
    try:
        result = physics_based_score_fn(mutation)
    except Exception as err:
        logger.warning("physics_based_score failed for %s: %s", mutation, err)
        return None
    if not isinstance(result, dict):
        logger.warning("physics_based_score returned non-dict for %s", mutation)
        return None
    return result


def parse_mutation_tokens(raw_text: str) -> list:
    """Parse comma-separated mutation strings into validated point mutations."""
    if not raw_text:
        return []

    mutations = []
    for token in str(raw_text).replace(";", ",").split(","):
        candidate = token.strip().upper()
        parsed = parse_single_mutation(candidate)
        if parsed is None:
            continue
        wt_aa, pos, mut_aa = parsed
        mutations.append(f"{wt_aa}{pos}{mut_aa}")
    return mutations


def build_sequence_from_mutations(target_mutation: str, rescue_mutations: list) -> str:
    """Apply target + rescue mutations sequentially and return final sequence."""
    sequence = apply_mutation(P53_WT, target_mutation) if target_mutation else None
    if sequence is None:
        sequence = P53_WT

    for mut in rescue_mutations:
        updated = apply_mutation(sequence, mut)
        if updated is not None:
            sequence = updated

    return sequence


@st.cache_data(show_spinner=False)
def get_dms_lookup_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load single-mutation DMS data and return lookup tables by mutation and position."""
    try:
        dms_df = get_dms_data()
    except Exception:
        return pd.DataFrame(columns=["mutation", "score"]), pd.DataFrame(columns=["pos", "pos_mean_score", "pos_count"])

    if dms_df is None or dms_df.empty or "mutation" not in dms_df.columns or "score" not in dms_df.columns:
        return pd.DataFrame(columns=["mutation", "score"]), pd.DataFrame(columns=["pos", "pos_mean_score", "pos_count"])

    dms_work = dms_df.copy()
    if "n_mutations" in dms_work.columns:
        dms_work = dms_work[dms_work["n_mutations"] == 1].copy()

    dms_work["mutation"] = dms_work["mutation"].astype(str).str.strip().str.upper()
    dms_work["score"] = pd.to_numeric(dms_work["score"], errors="coerce")
    dms_work = dms_work.dropna(subset=["mutation", "score"])

    by_mutation = dms_work.groupby("mutation", as_index=False)["score"].mean()

    if "pos" in dms_work.columns:
        by_position = dms_work.groupby("pos", as_index=False).agg(
            pos_mean_score=("score", "mean"),
            pos_count=("score", "count"),
        )
    else:
        by_position = pd.DataFrame(columns=["pos", "pos_mean_score", "pos_count"])

    return by_mutation, by_position


@st.cache_data(show_spinner=False)
def get_score_calibration_profile() -> Dict[str, float]:
    """Build a robust score calibration profile from observed single-mutation DMS labels."""
    fallback = {
        "clip_low": -3.0,
        "clip_high": 4.0,
        "center": 0.0,
        "scale": 1.0,
        "n_scores": 0.0,
    }
    try:
        dms_df = get_dms_data()
    except Exception:
        return fallback
    if dms_df is None or dms_df.empty or "score" not in dms_df.columns:
        return fallback

    dms_work = dms_df.copy()
    if "n_mutations" in dms_work.columns:
        dms_work = dms_work[dms_work["n_mutations"] == 1].copy()

    scores = pd.to_numeric(dms_work["score"], errors="coerce").dropna().to_numpy(dtype=float)
    if scores.size < 20:
        return fallback

    q1, q5, q95, q99 = np.percentile(scores, [1, 5, 95, 99])
    center = float(np.median(scores))
    scale = float(max((q95 - q5) / 2.0, 1e-3))
    clip_low = float(min(q1, q5))
    clip_high = float(max(q99, q95))
    if clip_high <= clip_low:
        clip_low, clip_high = -3.0, 4.0

    return {
        "clip_low": clip_low,
        "clip_high": clip_high,
        "center": center,
        "scale": scale,
        "n_scores": float(scores.size),
    }


def calibrate_oracle_score(raw_score: float, profile: Dict[str, float]) -> float:
    """Map raw oracle output to observed DMS-like range with smooth saturation + clipping."""
    center = float(profile.get("center", 0.0))
    scale = float(max(profile.get("scale", 1.0), 1e-6))
    clip_low = float(profile.get("clip_low", -3.0))
    clip_high = float(profile.get("clip_high", 4.0))
    squashed = center + scale * np.tanh((float(raw_score) - center) / scale)
    return float(np.clip(squashed, clip_low, clip_high))


def estimate_oracle_uncertainty(pooled: torch.Tensor, mc_samples: int = 8) -> float:
    """Estimate predictive uncertainty via MC-dropout standard deviation."""
    if oracle is None or pooled is None:
        return 0.0
    n_samples = max(int(mc_samples), 2)
    was_training = oracle.model.training
    preds = []
    try:
        oracle.model.train()
        with torch.no_grad():
            for _ in range(n_samples):
                preds.append(float(oracle.model(pooled).squeeze(-1).mean().item()))
    except Exception:
        return 0.0
    finally:
        if not was_training:
            oracle.model.eval()
    return float(np.std(np.asarray(preds, dtype=float)))


def build_trust_adjusted_score(
    raw_score: float,
    pooled: Optional[torch.Tensor],
    pooled_ref: Optional[torch.Tensor],
    calibration_profile: Dict[str, float],
    uncertainty_weight: float = 0.8,
    ood_rank_weight: float = 1.25,
    ood_radius: float = 1.75,
    mc_samples: int = 8,
) -> Dict[str, float]:
    """Compute calibrated and trust-adjusted score used for ranking/filtering."""
    calibrated = calibrate_oracle_score(raw_score, calibration_profile)
    uncertainty = estimate_oracle_uncertainty(pooled, mc_samples=mc_samples)
    ood_distance = 0.0
    if pooled is not None and pooled_ref is not None:
        with torch.no_grad():
            ood_distance = float(torch.norm(pooled - pooled_ref, p=2, dim=-1).mean().item())
    ood_excess = max(0.0, ood_distance - float(ood_radius))
    adjusted = calibrated - float(uncertainty_weight) * uncertainty - float(ood_rank_weight) * ood_excess
    adjusted = float(
        np.clip(
            adjusted,
            float(calibration_profile.get("clip_low", -3.0)),
            float(calibration_profile.get("clip_high", 4.0)),
        )
    )
    return {
        "score_raw": float(raw_score),
        "score_calibrated": float(calibrated),
        "score_adjusted": float(adjusted),
        "uncertainty": float(uncertainty),
        "ood_distance": float(ood_distance),
        "ood_excess": float(ood_excess),
    }


def get_active_analysis_candidate(default_target_mutation: str) -> dict | None:
    """Build a normalized candidate payload from session state for analysis panels."""
    selected = st.session_state.get("selected_candidate")
    if isinstance(selected, dict) and selected.get("sequence"):
        clean_muts = [
            m.strip().upper()
            for m in selected.get("mutations", [])
            if parse_single_mutation(str(m).strip().upper()) is not None
        ]
        return {
            "target_mutation": st.session_state.get("target_mut_saved", default_target_mutation),
            "sequence": str(selected.get("sequence", P53_WT)),
            "mutations": clean_muts,
            "score": float(selected.get("score", 0.0)),
            "stability": float(selected.get("stability", 0.0)),
            "binding": float(selected.get("binding", 0.0)),
            "identity": float(selected.get("identity", 0.0)),
            "n_mutations": int(selected.get("n_mutations", len(clean_muts))),
            "source": "selected_candidate",
            "mut_summary_raw": ", ".join(clean_muts),
        }

    results = st.session_state.get("results")
    if isinstance(results, pd.DataFrame) and not results.empty and "Score" in results.columns:
        best = results.loc[results["Score"].idxmax()]
        parsed_muts = parse_mutation_tokens(best.get("MutSummary", ""))
        sequence = str(best.get("Sequence", ""))
        if not sequence:
            sequence = build_sequence_from_mutations(default_target_mutation, parsed_muts)
        return {
            "target_mutation": st.session_state.get("target_mut_saved", default_target_mutation),
            "sequence": sequence,
            "mutations": parsed_muts,
            "score": float(best.get("Score", 0.0)),
            "stability": float(best.get("Stability", 0.0)),
            "binding": float(best.get("DNARecruitment", 0.0)),
            "identity": float(best.get("Identity", 0.0)),
            "n_mutations": int(best.get("MutationCount", len(parsed_muts))),
            "source": "results_table",
            "mut_summary_raw": str(best.get("MutSummary", "")),
        }

    return None


def sync_analysis_inputs_from_candidate(active_candidate: Optional[Dict[str, Any]]) -> None:
    """Sync analysis text inputs when a different candidate is selected."""
    if not active_candidate:
        return

    target = str(active_candidate.get("target_mutation", "")).strip().upper()
    rescue_list = [str(m).strip().upper() for m in active_candidate.get("mutations", []) if str(m).strip()]
    rescue_text = ", ".join(rescue_list[:6])
    candidate_key = f"{target}|{','.join(rescue_list)}"

    if st.session_state.get("_analysis_last_candidate_key") == candidate_key:
        return

    st.session_state["_analysis_last_candidate_key"] = candidate_key
    st.session_state["exp_cancer_main"] = target or "R175H"
    st.session_state["exp_rescue_main"] = rescue_text or "N239Y, T284R"
    st.session_state["ci_cancer_main"] = target or "R175H"
    st.session_state["ci_rescue_main"] = rescue_text or "N239Y, T284R"


@st.cache_data(show_spinner=False)
def load_campaign_bundle_from_disk(run_id: str) -> Dict[str, Any]:
    """Load a campaign artifact bundle from local storage."""
    store = CampaignStore()
    bundle = store.load_run_bundle(run_id)
    bundle["run_id"] = str(run_id)
    return bundle


@st.cache_data(show_spinner=False)
def load_campaign_bundle_from_zip(raw_bytes: bytes) -> Dict[str, Any]:
    """Load a campaign artifact bundle from uploaded zip bytes."""
    default_bundle: Dict[str, Any] = {
        "manifest": {},
        "scenarios": pd.DataFrame(),
        "candidates": pd.DataFrame(),
        "trajectories": pd.DataFrame(),
        "clinical": pd.DataFrame(),
        "top30": pd.DataFrame(),
        "run_id": "uploaded_bundle",
    }
    if not raw_bytes:
        return default_bundle

    with zipfile.ZipFile(io.BytesIO(raw_bytes)) as zf:
        names = zf.namelist()

        def _member_for(filename: str) -> Optional[str]:
            for name in names:
                if name.endswith("/"):
                    continue
                if Path(name).name == filename:
                    return name
            return None

        manifest_member = _member_for("manifest.json")
        if manifest_member:
            try:
                default_bundle["manifest"] = json.loads(zf.read(manifest_member).decode("utf-8"))
            except Exception:
                default_bundle["manifest"] = {}

        if default_bundle["manifest"].get("run_id"):
            default_bundle["run_id"] = str(default_bundle["manifest"]["run_id"])
        elif manifest_member:
            parent = Path(manifest_member).parent
            default_bundle["run_id"] = parent.name or "uploaded_bundle"

        for key, filename in [
            ("scenarios", "scenarios.parquet"),
            ("candidates", "candidates.parquet"),
            ("trajectories", "trajectories.parquet"),
            ("clinical", "clinical.parquet"),
            ("top30", "top30.parquet"),
        ]:
            member = _member_for(filename)
            if not member:
                continue
            try:
                default_bundle[key] = pd.read_parquet(io.BytesIO(zf.read(member)))
            except Exception:
                default_bundle[key] = pd.DataFrame()

    return default_bundle


def apply_campaign_bundle_to_session(bundle: Dict[str, Any], source_label: str) -> bool:
    """Hydrate Streamlit candidate state from campaign artifact tables."""
    candidates_df = bundle.get("candidates", pd.DataFrame())
    trajectories_df = bundle.get("trajectories", pd.DataFrame())
    top30_df = bundle.get("top30", pd.DataFrame())
    display_candidates = hydrate_display_candidates(
        candidates_df=candidates_df,
        trajectories_df=trajectories_df,
        top_df=top30_df,
        max_candidates=30,
    )
    if not display_candidates:
        return False

    st.session_state["gd_candidates"] = display_candidates
    st.session_state["gd_source_mode"] = "campaign_artifacts"
    st.session_state["gd_artifact_source"] = str(source_label)
    st.session_state["gd_artifact_manifest"] = bundle.get("manifest", {})
    st.session_state["gd_artifact_run_id"] = str(bundle.get("run_id", source_label))
    st.session_state["gd_artifact_scenarios"] = int(len(bundle.get("scenarios", pd.DataFrame())))

    # Keep downstream analysis synced to a concrete selected candidate.
    st.session_state["selected_candidate"] = display_candidates[0]
    first_target = str(display_candidates[0].get("target_label", "")).strip()
    if first_target:
        st.session_state["target_mut_saved"] = first_target.split("+")[0]
    return True


def compute_target_baseline_metrics(target_mutation: str) -> Dict[str, float]:
    """Compute baseline oracle metrics for the target-only mutant."""
    if oracle is None:
        return {
            "score": 0.0,
            "score_raw": 0.0,
            "score_calibrated": 0.0,
            "stability": 0.0,
            "binding": 0.0,
            "uncertainty": 0.0,
            "ood_distance": 0.0,
        }

    seq = apply_mutation(P53_WT, target_mutation) if target_mutation else None
    if seq is None:
        seq = P53_WT

    calibration_profile = get_score_calibration_profile()
    aa_ids = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    with torch.no_grad():
        emb = embedder.get_embeddings(seq).detach()
        z, logits, probs = embedder.latent_forward_ascent(emb)
        pooled = z.mean(dim=1)
        if pooled.shape[-1] != oracle.input_dim:
            pooled = pooled[:, :oracle.input_dim]
        raw_score = float(oracle.model(pooled).item())
        logits_aa = logits[:, :, aa_ids]
        stability = float(F.log_softmax(logits_aa, dim=-1).max(dim=-1).values.mean().item())
        binding = float(embedder.get_dna_contact_prob(z, logits, probs=probs).item())

    score_bundle = build_trust_adjusted_score(
        raw_score=raw_score,
        pooled=pooled,
        pooled_ref=pooled,
        calibration_profile=calibration_profile,
        uncertainty_weight=0.8,
        ood_rank_weight=1.25,
        ood_radius=1.75,
        mc_samples=8,
    )
    return {
        "score": float(score_bundle["score_adjusted"]),
        "score_raw": float(score_bundle["score_raw"]),
        "score_calibrated": float(score_bundle["score_calibrated"]),
        "stability": stability,
        "binding": binding,
        "uncertainty": float(score_bundle["uncertainty"]),
        "ood_distance": float(score_bundle["ood_distance"]),
    }


def filter_candidate_set(
    candidates: list[dict],
    target_mutation: str,
    baseline_score: float,
    min_score_gain: float = 0.05,
    min_rescue_mutations: int = 1,
) -> tuple[list[dict], Dict[str, Any]]:
    """Drop baseline-like candidates that do not meaningfully improve over target-only baseline."""
    target_mut = str(target_mutation).strip().upper()
    kept = []
    removed = 0
    rejected_reasons: Dict[str, int] = {}

    for cand in candidates:
        muts = [str(m).strip().upper() for m in cand.get("mutations", []) if str(m).strip()]
        rescue_muts = [m for m in muts if m != target_mut]
        score_primary = float(cand.get("score_adjusted", cand.get("score", 0.0)))
        score_gain = score_primary - float(baseline_score)
        has_rescue_count = len(rescue_muts) >= int(min_rescue_mutations)
        has_min_gain = score_gain >= float(min_score_gain)
        is_baseline_like = (not has_rescue_count) or (not has_min_gain)
        if not has_rescue_count:
            rejection_reason = "insufficient_rescue_mutations"
        elif not has_min_gain:
            rejection_reason = "below_min_gain"
        else:
            rejection_reason = "passed_quality_filter"

        cand["score_for_selection"] = float(score_primary)
        cand["rescue_mutations"] = rescue_muts
        cand["score_gain_vs_target"] = score_gain
        cand["is_baseline_like"] = is_baseline_like
        cand["selection_reason"] = rejection_reason

        if is_baseline_like:
            removed += 1
            rejected_reasons[rejection_reason] = rejected_reasons.get(rejection_reason, 0) + 1
        else:
            kept.append(cand)

    fallback_used = False
    if not kept and candidates:
        fallback_used = True
        kept = sorted(candidates, key=lambda c: float(c.get("score", 0.0)), reverse=True)[:1]
        kept[0]["selection_reason"] = "fallback_top_score_no_quality_pass"

    summary = {
        "total": len(candidates),
        "kept": len(kept),
        "removed": removed,
        "fallback_used": fallback_used,
        "min_score_gain": float(min_score_gain),
        "baseline_score": float(baseline_score),
        "rejection_reasons": rejected_reasons,
    }
    return kept, summary


def mutation_set_for_candidate(candidate: Dict[str, Any], target_mutation: str) -> set[str]:
    target = str(target_mutation).strip().upper()
    muts = [str(m).strip().upper() for m in candidate.get("mutations", []) if str(m).strip()]
    return {m for m in muts if m != target}


def jaccard_overlap(set_a: set[str], set_b: set[str]) -> float:
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 0.0
    return float(len(set_a & set_b) / len(union))


def candidate_rank_value(candidate: Dict[str, Any]) -> float:
    if "selection_score" in candidate:
        return float(candidate["selection_score"])
    if "score_gain_vs_target" in candidate:
        return float(candidate["score_gain_vs_target"])
    if "score_adjusted" in candidate:
        return float(candidate["score_adjusted"])
    return float(candidate.get("score", 0.0))


def select_diverse_candidates(
    candidates: list[dict],
    target_mutation: str,
    desired_count: int,
    novelty_weight: float = 0.35,
    overlap_penalty: float = 0.20,
) -> tuple[list[dict], Dict[str, Any]]:
    """Greedy re-ranker using score gain + novelty and explicit mutation overlap penalty."""
    if not candidates:
        return [], {
            "pool_size": 0,
            "selected": 0,
            "novelty_weight": float(novelty_weight),
            "overlap_penalty": float(overlap_penalty),
            "fallback_used": False,
        }

    target = str(target_mutation).strip().upper()
    enriched: list[Dict[str, Any]] = []
    for cand in candidates:
        mut_set = mutation_set_for_candidate(cand, target)
        score_gain = float(cand.get("score_gain_vs_target", cand.get("score_adjusted", cand.get("score", 0.0))))
        enriched.append(
            {
                "candidate": cand,
                "mut_set": mut_set,
                "profile": str(cand.get("profile", "Unknown")),
                "score_gain": score_gain,
            }
        )

    # Stage 1: keep the best restart per profile.
    best_by_profile: Dict[str, Dict[str, Any]] = {}
    for row in enriched:
        profile = row["profile"]
        prev = best_by_profile.get(profile)
        if prev is None or row["score_gain"] > prev["score_gain"]:
            best_by_profile[profile] = row

    selected_rows: list[Dict[str, Any]] = []
    seed_rows = sorted(best_by_profile.values(), key=lambda r: r["score_gain"], reverse=True)
    for row in seed_rows:
        if len(selected_rows) >= desired_count:
            break
        selected_rows.append(row)

    selected_ids = {id(r["candidate"]) for r in selected_rows}
    remaining = [r for r in enriched if id(r["candidate"]) not in selected_ids]

    # Stage 2: greedy utility = score_gain + novelty bonus - overlap penalty.
    while len(selected_rows) < desired_count and remaining:
        best_idx = -1
        best_utility = -float("inf")
        for idx, row in enumerate(remaining):
            if selected_rows:
                overlaps = [jaccard_overlap(row["mut_set"], s["mut_set"]) for s in selected_rows]
                max_overlap = float(max(overlaps))
                novelty = float(1.0 - np.mean(overlaps))
            else:
                max_overlap = 0.0
                novelty = 1.0
            utility = row["score_gain"] + (novelty_weight * novelty) - (overlap_penalty * max_overlap)
            if utility > best_utility:
                best_utility = utility
                best_idx = idx
        if best_idx < 0:
            break
        chosen = remaining.pop(best_idx)
        chosen["selection_utility"] = float(best_utility)
        selected_rows.append(chosen)

    if not selected_rows:
        top = max(enriched, key=lambda r: r["score_gain"])
        selected_rows = [top]

    selected: list[dict] = []
    for rank, row in enumerate(selected_rows[:desired_count], start=1):
        cand = row["candidate"]
        if selected:
            prior_sets = [mutation_set_for_candidate(c, target) for c in selected]
            overlaps = [jaccard_overlap(row["mut_set"], s) for s in prior_sets]
            max_overlap = float(max(overlaps))
            novelty = float(1.0 - np.mean(overlaps))
        else:
            max_overlap = 0.0
            novelty = 1.0
        selection_score = row["score_gain"] + (novelty_weight * novelty) - (overlap_penalty * max_overlap)
        cand["selection_rank"] = rank
        cand["selection_score"] = float(selection_score)
        cand["diversity_novelty"] = float(novelty)
        cand["max_mutation_overlap"] = float(max_overlap)
        cand["selection_reason"] = (
            f"diversity_rank_{rank}_utility={selection_score:.3f}_"
            f"gain={row['score_gain']:.3f}_novelty={novelty:.3f}_overlap={max_overlap:.3f}"
        )
        selected.append(cand)

    summary = {
        "pool_size": len(candidates),
        "selected": len(selected),
        "novelty_weight": float(novelty_weight),
        "overlap_penalty": float(overlap_penalty),
        "fallback_used": len(selected_rows) == 1 and len(candidates) > 1,
    }
    return selected, summary


# Initialize models with Streamlit resource cache to avoid repeated heavy reloads.
@st.cache_resource(show_spinner=False)
def load_models_v8():
    from p53cad.engine.latent import ManifoldEmbedder
    from p53cad.engine.oracle import FunctionalOracle
    from p53cad.engine.explain import SaliencyMap
    from p53cad.analysis.grassmann import GrassmannMetric
    from p53cad.viz.pymol import PyMolGenerator
    
    embedder = ManifoldEmbedder()
    oracle_path = Path("data/models/functional_oracle.pt")
    oracle = FunctionalOracle(model_path=oracle_path) if oracle_path.exists() else None
    explainer = SaliencyMap(oracle, embedder) if oracle else None
    grassmann = GrassmannMetric(embedder)
    viz = PyMolGenerator()
    return embedder, oracle, explainer, grassmann, viz


def clear_model_cache() -> None:
    load_models_v8.clear()

embedder, oracle, explainer, grassmann, viz = load_models_v8()

# --- HERO HEADER ---
st.markdown("""
<div class="hero-container">
    <div class="hero-badge">ISEF 2026</div>
    <h1 class="hero-title">p53-proteoMgCAD</h1>
    <p class="hero-subtitle">Mutative Generative Computer-Assisted Design of Second-Site Rescues for p53</p>
    <p class="hero-description">Constraint-based protein engineering inspired by mechanical topology optimization</p>
</div>
""", unsafe_allow_html=True)

# --- CONFIG ---
# Big 8 Hotspots (~28% of all p53 mutations)
HOTSPOT_FLEET = ["R175H", "R248Q", "R248W", "R273H", "R273C", "G245S", "R249S", "R282W"]

# Extended mutation database (Top 50 most frequent in cancer)
P53_MUTATIONS = {
    "Structural (Zinc Region)": [
        "R175H", "R175G", "R175L",  # Zinc coordination
        "C176F", "C176Y",           # Zinc ligand
        "H179R", "H179Y",           # Zinc ligand
        "C242F", "C242S",           # Zinc ligand
        "G245S", "G245D", "G245C",  # β-sandwich core
    ],
    "DNA Contact": [
        "R248Q", "R248W", "R248L",  # Major groove contact
        "R273H", "R273C", "R273L",  # DNA backbone contact
        "R280K", "R280T",           # Minor groove
        "R282W", "R282Q",           # DNA positioning
    ],
    "Loop Regions": [
        "R249S", "R249M",           # L3 loop
        "G266E", "G266R",           # L3 loop
        "R213Q", "R213X",           # L2 loop
        "Y220C", "Y220S",           # β-sandwich/loop junction
    ],
    "β-Sandwich Core": [
        "V157F", "V157G",           # Hydrophobic core
        "I195T", "I195F",           # Core packing
        "Y163C", "Y163H",           # Core stability
        "I255F", "I255T",           # Core packing
        "F270L", "F270S",           # Aromatic core
    ],
    "Other Frequent": [
        "P151S", "P152L",           # Proline turns
        "M237I", "M237V",           # Core
        "S241F", "S241C",           # Near zinc
        "N239D", "N239S",           # Loop
        "H168R", "H168Y",           # Interface
        "E286K", "E286G",           # Surface
    ]
}

# Flatten for easy access
ALL_MUTATIONS = [m for muts in P53_MUTATIONS.values() for m in muts]

# --- SIDEBAR CONTROLS ---
with st.sidebar:
    st.markdown("""
    <div style="padding: 1.5rem 0; border-bottom: 1px solid #E5E7EB; margin-bottom: 1.5rem;">
        <h2 style="font-size: 1.1rem; font-weight: 600; color: #1A1A1A; margin: 0;">Design Parameters</h2>
        <p style="font-size: 0.8rem; color: #6B7280; margin: 0.25rem 0 0 0;">Configure your rescue design</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<p class="section-header">Target Selection</p>', unsafe_allow_html=True)

    study_mode = st.radio("Laboratory Scope", [
        "Individual Target",
        "Fleet Study (Big 8)",
        "Extended Fleet (Top 50)",
        "Universal Design"
    ], label_visibility="collapsed")

    if study_mode == "Individual Target":
        mut_category = st.selectbox("Mutation Category",
            ["Custom Input"] + list(P53_MUTATIONS.keys()),
            help="Select a category or enter custom mutation",
            label_visibility="collapsed")

        if mut_category == "Custom Input":
            target_mut = st.text_input("Target Mutation", value="R175H",
                help="Enter any p53 mutation (e.g., R175H, Y220C, V157F)",
                placeholder="e.g., R175H")
        else:
            target_mut = st.selectbox("Select Mutation", P53_MUTATIONS[mut_category],
                help=f"Common {mut_category.lower()} mutations")

        if target_mut:
            pos = int(''.join(filter(str.isdigit, target_mut)))
            domain = "DNA-binding domain" if 94 <= pos <= 292 else "N-terminal" if pos < 94 else "C-terminal"
            st.markdown(f'<span class="feature-pill">Position {pos}</span><span class="feature-pill">{domain}</span>', unsafe_allow_html=True)

    elif study_mode == "Fleet Study (Big 8)":
        st.info(f"Targeting {len(HOTSPOT_FLEET)} hotspot mutations")
        target_mut = "R175H"

    elif study_mode == "Extended Fleet (Top 50)":
        selected_categories = st.multiselect("Categories",
            list(P53_MUTATIONS.keys()), default=["Structural (Zinc Region)", "DNA Contact"])
        fleet_muts = [m for cat in selected_categories for m in P53_MUTATIONS.get(cat, [])]
        target_mut = fleet_muts[0] if fleet_muts else "R175H"

    else:
        target_mut = "UNIVERSAL"

    st.markdown('<p class="section-header" style="margin-top: 1.5rem;">Strategy</p>', unsafe_allow_html=True)

    strategy = st.selectbox("Rescue Strategy",
                          ["Gradient Steering (Adaptive)", "Linear Manifold Interpolation"],
                          help="Adaptive uses AI optimization; Linear is geometric baseline.",
                          label_visibility="collapsed")

    st.markdown('<p class="section-header" style="margin-top: 1.5rem;">Constraints</p>', unsafe_allow_html=True)

    lock_res = st.text_input("Locked Residues", value="248, 273",
                            help="Critical sites protected from mutation.",
                            placeholder="e.g., 248, 273")

    similarity_weight = st.slider("Identity Preservation", 0.0, 50.0, 35.0,
                                help="Higher = fewer mutations. Target >90% for therapeutic viability.")

    stability_bias = st.slider("Stability Bias", 0.0, 1.0, 0.2,
                             help="Prioritize structural stability.")

    st.markdown('<p class="section-header" style="margin-top: 1.5rem;">Sampling</p>', unsafe_allow_html=True)

    rescue_steps = st.slider("Resolution", 20, 500, 300, label_visibility="collapsed")
    deep_manifold = st.checkbox("Deep Manifold Sampling", value=False,
                                help="Reveals true fitness landscape topology.")

    st.markdown("""
    <div style="margin-top: 2rem; padding-top: 1rem; border-top: 1px solid #E5E7EB;">
        <p style="font-size: 0.75rem; color: #9CA3AF; margin: 0;">GPU Accelerated</p>
    </div>
    """, unsafe_allow_html=True)

# --- DATA PROCESSING ---
def run_search(target_mut_override=None):
    current_mut = target_mut_override if target_mut_override else target_mut
    walker = ManifoldWalker(embedder)
    
    # Handle Universal Mode initialization
    init_mut = HOTSPOT_FLEET[0] if current_mut == "UNIVERSAL" else current_mut
    cancer_seq = apply_mutation(P53_WT, init_mut)
    if cancer_seq is None:
        cancer_seq = P53_WT # Fallback to WT if mutation parsing fails
    if not cancer_seq: return None
    
    # Standard Amino Acid Vocabulary IDs for ESM-2
    AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    
    if "Gradient" in strategy:
        # EMBEDDING-SPACE OPTIMIZATION Core
        emb = embedder.get_embeddings(cancer_seq).detach().requires_grad_(True)
        emb_wt = embedder.get_embeddings(P53_WT).detach()
        
        optimizer = torch.optim.Adam([emb], lr=0.04) # Optimal from grid search
        locked_indices = [int(i.strip())-1 for i in lock_res.split(",") if i.strip()]
        
        history = []
        
        # RECORD BASELINE (Step 0)
        with torch.no_grad():
            z_init, logits_init, probs_init = embedder.latent_forward_ascent(emb)
            score_init = oracle.model(z_init.mean(dim=1)).item()
            history.append({
                "Step": 0,
                "Score": score_init,
                "Stability": F.log_softmax(logits_init[:, :, AA_IDS], dim=-1).max(dim=-1).values.mean().item(),
                "DNARecruitment": embedder.get_dna_contact_prob(z_init, logits_init, probs=probs_init).item(),
                "SurfaceCharge": embedder.get_surface_charge_density(logits_init, probs=probs_init).item(),
                "HydroPacking": embedder.get_hydrophobic_packing(logits_init, probs=probs_init).item(),
                "GrassmannNovelty": 0.0,
                "Identity": 100.0 * (1.0 - (1.0 / len(P53_WT))),
                "LatentIdentity": 100.0,
                "LatentExcitation": 0.0,
                "LX": z_init.mean(dim=1)[0,0].item() * 5.0,
                "LY": z_init.mean(dim=1)[0,1].item() * 5.0,
                "MutationCount": 1,
                "MutSummary": current_mut,
                "Sequence": cancer_seq,
                "Phase": "Start"
            })

        for step_idx in range(1, rescue_steps + 1):
            optimizer.zero_grad()
            
            # Unpack the optimized 3-value signature: (h, logits, probabilities)
            z, logits, probs = embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            
            # Validation: Force dimension consistency for Oracle
            if pooled.shape[-1] != oracle.input_dim:
                 pooled = pooled[:, :oracle.input_dim] 

            score = oracle.model(pooled)
            
            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()
            
            dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
            hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)
            surface_charge = embedder.get_surface_charge_density(logits, probs=probs)
            
            # LOSS FUNCTION: Balanced for curriculum learning
            # Functional terms boosted to compete with identity during late-phase
            if current_mut == "UNIVERSAL":
                loss = -score * 5.0  # Boosted for multi-target
            else:
                loss = -score * 4.0  # Boosted from 2.0

            if stability_bias > 0:
                loss -= (stability_bias * 8.0) * stability  # Slightly boosted

            loss -= 2.5 * dna_force         # Mechanistic DNA recruitment (boosted)
            loss -= 3.0 * hydro_packing     # Structural core integrity (boosted)
            loss -= 2.5 * surface_charge    # Interface electrostatic affinity (boosted)
                
            if similarity_weight > 0:
                # === DIRECT SEQUENCE IDENTITY CONSTRAINT ===
                # Use probability-weighted mutation penalty (differentiable!)
                # probs shape: (L, V_aa) where V_aa = 20 amino acids

                # Get WT token indices in AA-restricted space
                wt_aa_indices = []
                for aa in P53_WT:
                    aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
                    if aa_id in AA_IDS:
                        wt_aa_indices.append(AA_IDS.index(aa_id))
                    else:
                        wt_aa_indices.append(0)  # fallback
                wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb.device)

                # Probability of NOT being the WT amino acid at each position
                # This is differentiable unlike argmax
                probs_aa = F.softmax(logits_aa[0], dim=-1)  # (L, 20)
                wt_probs = probs_aa[torch.arange(len(P53_WT), device=emb.device), wt_aa_tensor]
                mutation_prob = 1.0 - wt_probs  # Probability of mutation at each position

                # Expected number of mutations (differentiable!)
                expected_mutations = mutation_prob.sum()

                # Also get hard count for barriers
                with torch.no_grad():
                    decoded_ids = torch.argmax(probs_aa, dim=-1)
                    n_mutations = (decoded_ids != wt_aa_tensor).sum().item()
                    seq_identity = 100.0 * (1.0 - n_mutations / len(P53_WT))

                # DIFFERENTIABLE PENALTY: Penalize expected mutations
                max_mutations = int(len(P53_WT) * 0.10)  # ~40 for 90% identity
                loss += similarity_weight * 3.0 * F.relu(expected_mutations - max_mutations)

                # DIVERSITY BONUS: Encourage spreading mutations, not concentrating
                # High entropy = mutations spread across many positions (good)
                # Low entropy = mutations concentrated in few positions (bad - stuck in rut)
                if expected_mutations > 1.0:
                    # Normalize mutation_prob to get distribution
                    mut_dist = mutation_prob / (mutation_prob.sum() + 1e-8)
                    # Entropy of mutation distribution (higher = more diverse)
                    entropy = -torch.sum(mut_dist * torch.log(mut_dist + 1e-8))
                    # Max entropy for uniform distribution over ~40 positions
                    max_entropy = torch.log(torch.tensor(40.0, device=emb.device))
                    # Penalize low diversity (reward high entropy)
                    diversity_penalty = (max_entropy - entropy) * 0.5
                    loss += diversity_penalty

                # HARD BARRIERS (non-differentiable but prevents catastrophe)
                if seq_identity < 85.0:
                    loss += 300.0 * (85.0 - seq_identity)

                if seq_identity < 80.0:
                    loss += 1000.0 * (80.0 - seq_identity)

                # Light embedding regularization for smoothness
                dist_l1 = torch.norm(emb - emb_wt, p=1) / emb.numel()
                loss += (similarity_weight * 15.0) * dist_l1
                
                # AA composition constraint (stronger)
                with torch.no_grad():
                    _, logits_wt, probs_wt = embedder.latent_forward_ascent(emb_wt)

                aa_dist_current = probs.mean(dim=0)
                aa_dist_wt = probs_wt.mean(dim=0)

                kl_div = F.kl_div(
                    torch.log(aa_dist_current + 1e-10),
                    aa_dist_wt,
                    reduction='batchmean'
                )
                loss += 15.0 * kl_div  # Stronger AA composition constraint

                # STABILITY FLOOR: Penalize if PLL drops too much
                # The stability crash in your graph shows PLL dropping from ~0 to ~-0.4
                # Healthy proteins should have stability > -0.2
                if stability.item() < -0.2:
                    loss += 50.0 * (-0.2 - stability)
                
            if locked_indices:
                loss += 500.0 * F.mse_loss(emb[:, locked_indices, :], emb_wt[:, locked_indices, :])
            
            loss.backward()
            optimizer.step()

            # EXPLORATION NOISE: Add small perturbation to escape local minima
            # Temperature schedule: high early (exploration), low late (exploitation)
            progress = step_idx / rescue_steps
            if progress < 0.5:  # Only add noise in first half
                noise_scale = 0.02 * (1.0 - 2 * progress)  # 0.02 → 0.0
                with torch.no_grad():
                    noise = torch.randn_like(emb) * noise_scale
                    emb.data += noise

            if step_idx % 10 == 0 or step_idx == rescue_steps:
                with torch.no_grad():
                    # High-fidelity AA-restricted decoding
                    top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                    top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                    tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                    seq = "".join(tokens)
                    
                    if len(seq) != len(P53_WT):
                        seq = seq[:len(P53_WT)].ljust(len(P53_WT), "X")
                        
                    muts = [f"{P53_WT[j]}{j+1}{seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != seq[j]]
                    g_dist = grassmann.grassmann_distance(seq, P53_WT)
                    
                    # Smooth Latent Indicators
                    z_norm = torch.norm(emb - emb_wt).item()
                    cos_sim = F.cosine_similarity(emb.view(1, -1), emb_wt.view(1, -1)).item()
                    
                    # Track curriculum phase for visualization
                    curr_progress = step_idx / rescue_steps
                    phase = "Exploring" if curr_progress < 0.4 else "Constraining"

                    history.append({
                        "Step": step_idx,
                        "Score": score.item(),
                        "Stability": stability.item(),
                        "DNARecruitment": dna_force.item(),
                        "Phase": phase,
                        "SurfaceCharge": surface_charge.item(),
                        "HydroPacking": hydro_packing.item(),
                        "GrassmannNovelty": g_dist,
                        "Identity": 100.0 * (1.0 - len(muts)/len(P53_WT)),
                        "LatentIdentity": cos_sim * 100.0,
                        "LatentExcitation": z_norm,
                        "LX": pooled[0, 0].item() * 5.0, # Real Latent Coordinate (Dim 1)
                        "LY": pooled[0, 1].item() * 5.0, # Real Latent Coordinate (Dim 2)
                        "MutationCount": len(muts),
                        "MutSummary": ", ".join(muts[:3]) + ("..." if len(muts)>3 else ""),
                        "Sequence": seq
                    })
        df_res = pd.DataFrame(history)
        df_res.to_csv("results_raw.csv", index=False)
        return df_res
    else:
        # Interpolation Logic
        path = walker.interpolate(cancer_seq, P53_WT, steps=20)
        data = []
        for i, seq in enumerate(path):
            z = embedder.encode(seq)
            muts = [f"{P53_WT[j]}{j+1}{seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != seq[j]]
            g_dist = grassmann.grassmann_distance(seq, P53_WT)
            data.append({
                "Step": i, "Score": oracle.predict(z), "Stability": walker.stability_score(seq),
                "GrassmannNovelty": g_dist,
                "Identity": 100.0 * (1.0 - len(muts)/len(P53_WT)), "MutationCount": len(muts),
                "LX": z.mean(dim=1)[0,0].item() * 5.0 if hasattr(z, 'mean') else 0.0,
                "LY": z.mean(dim=1)[0,1].item() * 5.0 if hasattr(z, 'mean') else 0.0,
                "MutSummary": "WT" if not muts else ", ".join(muts[:2]), "Sequence": seq
            })
        return pd.DataFrame(data)

# === GENERATIVE DESIGN ENGINE (LIVE VERSION) ===
def run_generative_design_live(constraints: dict, n_candidates: int = 6,
                                n_steps: Optional[int] = None,
                                progress_callback=None, structure_callback=None):
    """
    Batched generative design engine with:
    - trust score cadence (reduced MC-dropout calls),
    - per-batch early stopping,
    - optional successive halving pre-screen,
    - trial batching for faster throughput.
    """
    target_mut = constraints.get("target_mutation", "R175H")
    cancer_seq = apply_mutation(P53_WT, target_mut)
    if cancer_seq is None:
        cancer_seq = P53_WT

    AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    min_identity = float(constraints.get("min_identity", 90.0))
    min_stability = float(constraints.get("min_stability", -0.3))
    min_binding = float(constraints.get("min_binding", 5.0))
    locked_positions = constraints.get("locked_positions", [248, 273])
    delivery_method = constraints.get("delivery_method", "gene_therapy")
    exploration_diversity = float(constraints.get("diversity", 0.5))

    calibration_profile = get_score_calibration_profile()
    clip_low = float(calibration_profile.get("clip_low", -3.0))
    clip_high = float(calibration_profile.get("clip_high", 4.0))
    uncertainty_weight = float(constraints.get("uncertainty_penalty_weight", 0.8))
    ood_rank_weight = float(constraints.get("ood_rank_penalty_weight", 1.25))
    ood_radius = float(constraints.get("ood_trust_radius", 1.75))
    ood_loss_weight = float(constraints.get("ood_loss_weight", 12.0))
    mc_dropout_samples = int(max(2, constraints.get("mc_dropout_samples", 8)))
    trust_eval_stride = int(max(5, constraints.get("trust_eval_stride", 20)))
    early_stop_patience = int(max(1, constraints.get("early_stop_patience", 8)))
    early_stop_min_delta = float(max(0.0, constraints.get("early_stop_min_delta", 1e-3)))
    early_stop_warmup_steps = int(max(20, constraints.get("early_stop_warmup_steps", 40)))
    batch_trial_size = int(max(1, constraints.get("batch_trial_size", 4)))
    optimizer_lr = float(max(1e-3, constraints.get("optimizer_lr", 0.03)))

    enable_successive_halving = bool(constraints.get("enable_successive_halving", True))
    halving_keep_ratio = float(np.clip(constraints.get("successive_halving_keep_ratio", 0.5), 0.25, 1.0))
    if n_steps is None:
        n_steps = int(constraints.get("optimization_steps", 100))
    n_steps = int(np.clip(n_steps, 20, 500))
    halving_warmup_steps = int(
        np.clip(
            constraints.get("successive_halving_warmup_steps", max(20, n_steps // 4)),
            20,
            n_steps,
        )
    )
    capture_checkpoints = bool(constraints.get("capture_structure_checkpoints", False))
    checkpoint_stride = int(max(10, constraints.get("structure_checkpoint_stride", 40)))

    if delivery_method == "protein_therapy":
        min_identity = max(min_identity, 95.0)
    elif delivery_method == "mrna_therapy":
        min_identity = max(min_identity, 92.0)

    emb_target_ref = embedder.get_embeddings(cancer_seq).detach()
    emb_wt = embedder.get_embeddings(P53_WT).detach()
    with torch.no_grad():
        z_target_ref, _, _ = embedder.latent_forward_ascent(emb_target_ref)
        pooled_target_ref = z_target_ref.mean(dim=1)
        if pooled_target_ref.shape[-1] != oracle.input_dim:
            pooled_target_ref = pooled_target_ref[:, :oracle.input_dim]

    wt_aa_indices = []
    for aa in P53_WT:
        aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
        if aa_id in AA_IDS:
            wt_aa_indices.append(AA_IDS.index(aa_id))
        else:
            wt_aa_indices.append(0)
    wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb_target_ref.device)
    seq_len = len(P53_WT)
    max_mutations = int(seq_len * (100 - min_identity) / 100)
    locked_indices = [int(p) - 1 for p in locked_positions if p]

    weight_profiles = [
        {"function": 4.0, "stability": 8.0, "binding": 2.5, "name": "Balanced", "color": "#2563EB"},
        {"function": 2.0, "stability": 15.0, "binding": 2.0, "name": "Stability-First", "color": "#0EA5E9"},
        {"function": 3.0, "stability": 5.0, "binding": 8.0, "name": "Binding-Optimized", "color": "#10B981"},
        {"function": 8.0, "stability": 4.0, "binding": 3.0, "name": "Function-Maximized", "color": "#14B8A6"},
        {"function": 5.0, "stability": 10.0, "binding": 5.0, "name": "Conservative", "color": "#1D4ED8"},
        {"function": 6.0, "stability": 6.0, "binding": 6.0, "name": "Experimental", "color": "#0891B2"},
    ]

    profile_attempt_counts = {p["name"]: 0 for p in weight_profiles}
    trial_specs = []
    for trial_idx in range(int(max(n_candidates, 1))):
        profile = weight_profiles[trial_idx % len(weight_profiles)]
        profile_attempt_counts[profile["name"]] += 1
        trial_specs.append(
            {
                "trial_idx": trial_idx + 1,
                "orig_idx": trial_idx,
                "seed": 42 + trial_idx * 17,
                "restart": profile_attempt_counts[profile["name"]],
                "profile": profile,
            }
        )

    def _batch_dna_force(z_batch: torch.Tensor, probs_full: torch.Tensor) -> torch.Tensor:
        hotspots = [119, 174, 240, 247, 272, 279]
        hotspots = [i for i in hotspots if i < z_batch.shape[1]]
        if not hotspots:
            return torch.zeros(z_batch.shape[0], device=z_batch.device)
        pos_charge_ids = [10, 15, 21]
        latent_force = z_batch[:, hotspots, :].norm(dim=-1).mean(dim=-1)
        charge_prob = probs_full[:, hotspots][:, :, pos_charge_ids].sum(dim=-1).mean(dim=-1)
        return 0.5 * latent_force + 5.0 * charge_prob

    def _batch_hydrophobic_packing(probs_full: torch.Tensor) -> torch.Tensor:
        hydro_ids = [4, 12, 7, 18, 22, 20]
        core_res = [i for i in range(93, 312) if i < probs_full.shape[1]]
        loops = set(range(111, 124)) | set(range(162, 195)) | set(range(235, 251))
        true_core = [i for i in core_res if i not in loops]
        if not true_core:
            return torch.zeros(probs_full.shape[0], device=probs_full.device)
        return probs_full[:, true_core][:, :, hydro_ids].sum(dim=-1).mean(dim=-1)

    def _init_batch_embeddings(batch_specs: list[Dict[str, Any]]) -> torch.Tensor:
        emb = emb_target_ref.repeat(len(batch_specs), 1, 1).clone().detach()
        noise_bank = []
        for spec in batch_specs:
            torch.manual_seed(int(spec["seed"]))
            noise_bank.append(torch.randn_like(emb_target_ref) * 0.05 * exploration_diversity)
        emb += torch.cat(noise_bank, dim=0)
        return emb.requires_grad_(True)

    def _score_with_cached_uncertainty(
        raw_score: float,
        ood_distance: float,
        cached_uncertainty: float,
    ) -> Dict[str, float]:
        calibrated = calibrate_oracle_score(raw_score, calibration_profile)
        adjusted = calibrated - uncertainty_weight * cached_uncertainty - ood_rank_weight * max(
            0.0, ood_distance - ood_radius
        )
        adjusted = float(np.clip(adjusted, clip_low, clip_high))
        return {
            "score_raw": float(raw_score),
            "score_calibrated": float(calibrated),
            "score_adjusted": float(adjusted),
            "uncertainty": float(cached_uncertainty),
            "ood_distance": float(ood_distance),
        }

    def _screen_trial_specs(specs: list[Dict[str, Any]]) -> list[Dict[str, Any]]:
        if not specs:
            return specs
        scored_rows = []
        for batch_start in range(0, len(specs), batch_trial_size):
            batch_specs = specs[batch_start: batch_start + batch_trial_size]
            emb = _init_batch_embeddings(batch_specs)
            batch_size = emb.shape[0]
            emb_wt_batch = emb_wt.expand(batch_size, -1, -1)
            pooled_ref_batch = pooled_target_ref.expand(batch_size, -1)
            func_w = torch.tensor([s["profile"]["function"] for s in batch_specs], device=emb.device)
            stab_w = torch.tensor([s["profile"]["stability"] for s in batch_specs], device=emb.device)
            bind_w = torch.tensor([s["profile"]["binding"] for s in batch_specs], device=emb.device)
            optimizer = torch.optim.Adam([emb], lr=optimizer_lr)

            for _ in range(halving_warmup_steps):
                optimizer.zero_grad()
                z, logits, _ = embedder.latent_forward_ascent(emb)
                pooled = z.mean(dim=1)
                if pooled.shape[-1] != oracle.input_dim:
                    pooled = pooled[:, :oracle.input_dim]
                score = oracle.model(pooled).squeeze(-1)
                probs_full = torch.softmax(logits, dim=-1)
                logits_aa = logits[:, :, AA_IDS]
                log_probs = F.log_softmax(logits_aa, dim=-1)
                stability = log_probs.max(dim=-1).values.mean(dim=-1)
                dna_force = _batch_dna_force(z, probs_full)
                hydro_packing = _batch_hydrophobic_packing(probs_full)
                ood_distance_t = torch.norm(pooled - pooled_ref_batch, p=2, dim=-1)

                probs_aa = torch.softmax(logits_aa, dim=-1)
                wt_probs = probs_aa[:, torch.arange(seq_len, device=emb.device), wt_aa_tensor]
                expected_mutations = (1.0 - wt_probs).sum(dim=-1)
                expected_identity = 100.0 * (1.0 - expected_mutations / float(seq_len))

                loss_vec = -score * func_w
                loss_vec -= stab_w * stability
                loss_vec -= bind_w * dna_force
                loss_vec -= 3.0 * hydro_packing
                loss_vec += ood_loss_weight * F.relu(ood_distance_t - ood_radius)
                loss_vec += 50.0 * F.relu(expected_mutations - max_mutations)
                loss_vec += 500.0 * F.relu((min_identity - 5.0) - expected_identity)
                loss_vec += 100.0 * F.relu(min_stability - stability)
                loss_vec += 80.0 * F.relu(min_binding - dna_force)
                if locked_indices:
                    lock_pen = (emb[:, locked_indices, :] - emb_wt_batch[:, locked_indices, :]).pow(2).mean(dim=(1, 2))
                    loss_vec += 500.0 * lock_pen
                loss_vec += 40.0 * (emb - emb_wt_batch).abs().mean(dim=(1, 2))

                loss = loss_vec.mean()
                loss.backward()
                torch.nn.utils.clip_grad_norm_([emb], max_norm=1.0)
                optimizer.step()

            with torch.no_grad():
                z, logits, _ = embedder.latent_forward_ascent(emb)
                pooled = z.mean(dim=1)
                if pooled.shape[-1] != oracle.input_dim:
                    pooled = pooled[:, :oracle.input_dim]
                raw_scores = oracle.model(pooled).squeeze(-1)
                ood_distances = torch.norm(pooled - pooled_ref_batch, p=2, dim=-1)
                for i, spec in enumerate(batch_specs):
                    calibrated = calibrate_oracle_score(float(raw_scores[i].item()), calibration_profile)
                    screen_score = calibrated - ood_rank_weight * max(0.0, float(ood_distances[i].item()) - ood_radius)
                    scored_rows.append((spec, float(screen_score)))

        keep_count = max(1, int(np.ceil(len(specs) * halving_keep_ratio)))
        sorted_rows = sorted(scored_rows, key=lambda x: x[1], reverse=True)

        # Keep at least one seed per profile before filling by score.
        kept: list[Dict[str, Any]] = []
        seen_profiles = set()
        for spec, sc in sorted_rows:
            pname = spec["profile"]["name"]
            if pname in seen_profiles:
                continue
            row = dict(spec)
            row["screen_score"] = float(sc)
            kept.append(row)
            seen_profiles.add(pname)
            if len(kept) >= keep_count:
                break
        if len(kept) < keep_count:
            kept_ids = {k["trial_idx"] for k in kept}
            for spec, sc in sorted_rows:
                if spec["trial_idx"] in kept_ids:
                    continue
                row = dict(spec)
                row["screen_score"] = float(sc)
                kept.append(row)
                kept_ids.add(spec["trial_idx"])
                if len(kept) >= keep_count:
                    break
        return kept

    effective_specs = list(trial_specs)
    if enable_successive_halving and len(trial_specs) > 4 and halving_warmup_steps < n_steps:
        effective_specs = _screen_trial_specs(trial_specs)
        logger.info(
            "Successive halving retained %d/%d trial plans (warmup steps=%d, keep_ratio=%.2f).",
            len(effective_specs),
            len(trial_specs),
            halving_warmup_steps,
            halving_keep_ratio,
        )

    effective_specs = sorted(effective_specs, key=lambda s: s["orig_idx"])
    total_effective = len(effective_specs)
    all_candidates = []

    for batch_start in range(0, total_effective, batch_trial_size):
        batch_specs = effective_specs[batch_start: batch_start + batch_trial_size]
        emb = _init_batch_embeddings(batch_specs)
        batch_size = emb.shape[0]
        emb_wt_batch = emb_wt.expand(batch_size, -1, -1)
        pooled_ref_batch = pooled_target_ref.expand(batch_size, -1)
        func_w = torch.tensor([s["profile"]["function"] for s in batch_specs], device=emb.device)
        stab_w = torch.tensor([s["profile"]["stability"] for s in batch_specs], device=emb.device)
        bind_w = torch.tensor([s["profile"]["binding"] for s in batch_specs], device=emb.device)
        optimizer = torch.optim.Adam([emb], lr=optimizer_lr)

        trajectories = [[] for _ in range(batch_size)]
        best_valid_states: list[Optional[Dict[str, Any]]] = [None for _ in range(batch_size)]
        best_valid_scores = np.full(batch_size, -float("inf"), dtype=float)
        cached_uncertainty = np.zeros(batch_size, dtype=float)
        checks_without_improve = 0

        for step_idx in range(1, n_steps + 1):
            optimizer.zero_grad()

            z, logits, _ = embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            if pooled.shape[-1] != oracle.input_dim:
                pooled = pooled[:, :oracle.input_dim]
            raw_score_t = oracle.model(pooled).squeeze(-1)
            probs_full = torch.softmax(logits, dim=-1)
            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability_t = log_probs.max(dim=-1).values.mean(dim=-1)
            dna_force_t = _batch_dna_force(z, probs_full)
            hydro_packing_t = _batch_hydrophobic_packing(probs_full)
            ood_distance_t = torch.norm(pooled - pooled_ref_batch, p=2, dim=-1)

            probs_aa = torch.softmax(logits_aa, dim=-1)
            wt_probs = probs_aa[:, torch.arange(seq_len, device=emb.device), wt_aa_tensor]
            expected_mutations = (1.0 - wt_probs).sum(dim=-1)
            expected_identity = 100.0 * (1.0 - expected_mutations / float(seq_len))

            loss_vec = -raw_score_t * func_w
            loss_vec -= stab_w * stability_t
            loss_vec -= bind_w * dna_force_t
            loss_vec -= 3.0 * hydro_packing_t
            loss_vec += ood_loss_weight * F.relu(ood_distance_t - ood_radius)
            loss_vec += 50.0 * F.relu(expected_mutations - max_mutations)
            loss_vec += 500.0 * F.relu((min_identity - 5.0) - expected_identity)
            loss_vec += 100.0 * F.relu(min_stability - stability_t)
            loss_vec += 80.0 * F.relu(min_binding - dna_force_t)
            if locked_indices:
                lock_pen = (emb[:, locked_indices, :] - emb_wt_batch[:, locked_indices, :]).pow(2).mean(dim=(1, 2))
                loss_vec += 500.0 * lock_pen
            loss_vec += 40.0 * (emb - emb_wt_batch).abs().mean(dim=(1, 2))

            loss = loss_vec.mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_([emb], max_norm=1.0)
            optimizer.step()

            if step_idx % 5 == 0 or step_idx == n_steps:
                improved_any = False
                with torch.no_grad():
                    top_ids_aa_batch = torch.argmax(logits_aa, dim=-1)

                    for i, spec in enumerate(batch_specs):
                        aa_local = [AA_IDS[idx] for idx in top_ids_aa_batch[i].tolist()]
                        tokens = embedder.tokenizer.convert_ids_to_tokens(aa_local)
                        current_seq = "".join(tokens)[:seq_len]
                        muts = [
                            f"{P53_WT[j]}{j+1}{current_seq[j]}"
                            for j in range(seq_len)
                            if P53_WT[j] != current_seq[j]
                        ]
                        mut_positions = [int("".join(filter(str.isdigit, m))) for m in muts if m]
                        seq_identity = 100.0 * (1.0 - len(muts) / float(seq_len))
                        raw_score = float(raw_score_t[i].item())
                        current_ood = float(ood_distance_t[i].item())
                        do_full_trust_eval = (
                            step_idx == n_steps
                            or step_idx == 5
                            or (step_idx % trust_eval_stride == 0)
                        )
                        if do_full_trust_eval:
                            unc = estimate_oracle_uncertainty(pooled[i: i + 1], mc_samples=mc_dropout_samples)
                            cached_uncertainty[i] = float(unc)
                            score_bundle = _score_with_cached_uncertainty(
                                raw_score=raw_score,
                                ood_distance=current_ood,
                                cached_uncertainty=float(cached_uncertainty[i]),
                            )
                        else:
                            score_bundle = _score_with_cached_uncertainty(
                                raw_score=raw_score,
                                ood_distance=current_ood,
                                cached_uncertainty=float(cached_uncertainty[i]),
                            )

                        structure_checkpoint = None
                        first_capture_step = 5
                        should_capture = (
                            capture_checkpoints
                            and (
                                step_idx == n_steps
                                or step_idx == first_capture_step
                                or (step_idx % checkpoint_stride == 0)
                            )
                        )
                        if should_capture:
                            structure_checkpoint = get_real_structure_checkpoint(current_seq)

                        state = {
                            "step": step_idx,
                            "score": float(score_bundle["score_adjusted"]),
                            "score_raw": float(score_bundle["score_raw"]),
                            "score_calibrated": float(score_bundle["score_calibrated"]),
                            "stability": float(stability_t[i].item()),
                            "binding": float(dna_force_t[i].item()),
                            "identity": float(seq_identity),
                            "n_mutations": int(len(muts)),
                            "mutations": muts,
                            "mut_positions": mut_positions,
                            "sequence": current_seq,
                            "uncertainty": float(score_bundle["uncertainty"]),
                            "ood_distance": float(score_bundle["ood_distance"]),
                            "lx": float(pooled[i, 0].item()),
                            "ly": float(pooled[i, 1].item()),
                            "lz": float(pooled[i, 2].item()) if pooled.shape[-1] > 2 else float(score_bundle["score_adjusted"]),
                            "loss_total": float(loss_vec[i].item()),
                            "structure_checkpoint": structure_checkpoint,
                        }
                        trajectories[i].append(state)

                        is_valid = (
                            seq_identity >= min_identity
                            and float(stability_t[i].item()) >= min_stability
                            and float(dna_force_t[i].item()) >= min_binding
                        )
                        if is_valid and state["score"] > (best_valid_scores[i] + early_stop_min_delta):
                            best_valid_scores[i] = float(state["score"])
                            best_valid_states[i] = state.copy()
                            improved_any = True

                        if progress_callback:
                            progress_callback(
                                {
                                    "candidate_idx": int(spec["orig_idx"]),
                                    "candidate_total": int(total_effective),
                                    "profile": spec["profile"],
                                    "step": int(step_idx),
                                    "total_steps": int(n_steps),
                                    "current_state": state,
                                    "trajectory": trajectories[i],
                                }
                            )

                if improved_any:
                    checks_without_improve = 0
                else:
                    checks_without_improve += 1
                if step_idx >= early_stop_warmup_steps and checks_without_improve >= early_stop_patience:
                    break

        for i, spec in enumerate(batch_specs):
            trial_trajectory = trajectories[i]
            final_state = best_valid_states[i] if best_valid_states[i] is not None else (trial_trajectory[-1] if trial_trajectory else None)
            if final_state is None:
                # Defensive fallback in case no checkpoint state was recorded.
                with torch.no_grad():
                    z, logits, _ = embedder.latent_forward_ascent(emb[i: i + 1])
                    pooled = z.mean(dim=1)
                    if pooled.shape[-1] != oracle.input_dim:
                        pooled = pooled[:, :oracle.input_dim]
                    raw_score = float(oracle.model(pooled).item())
                    tokens = embedder.tokenizer.convert_ids_to_tokens(
                        [AA_IDS[idx] for idx in torch.argmax(logits[:, :, AA_IDS], dim=-1)[0].tolist()]
                    )
                    seq = "".join(tokens)[:seq_len]
                    muts = [f"{P53_WT[j]}{j+1}{seq[j]}" for j in range(seq_len) if P53_WT[j] != seq[j]]
                    final_state = _score_with_cached_uncertainty(
                        raw_score=raw_score,
                        ood_distance=float(torch.norm(pooled - pooled_target_ref, p=2, dim=-1).item()),
                        cached_uncertainty=0.0,
                    )
                    final_state.update(
                        {
                            "step": 0,
                            "stability": float(min_stability),
                            "binding": float(min_binding),
                            "identity": 100.0 * (1.0 - len(muts) / float(seq_len)),
                            "n_mutations": len(muts),
                            "mutations": muts,
                            "mut_positions": [int("".join(filter(str.isdigit, m))) for m in muts if m],
                            "sequence": seq,
                            "lx": float(pooled[0, 0].item()),
                            "ly": float(pooled[0, 1].item()),
                            "lz": float(pooled[0, 2].item()) if pooled.shape[-1] > 2 else float(final_state["score_adjusted"]),
                            "loss_total": float("nan"),
                            "structure_checkpoint": None,
                        }
                    )

            candidate = {
                "candidate_id": int(spec["trial_idx"]),
                "profile": spec["profile"]["name"],
                "color": spec["profile"]["color"],
                "restart": int(spec["restart"]),
                "trial_idx": int(spec["trial_idx"]),
                "sequence": final_state["sequence"],
                "score": float(final_state["score"]),
                "score_raw": float(final_state.get("score_raw", final_state["score"])),
                "score_calibrated": float(final_state.get("score_calibrated", final_state["score"])),
                "stability": float(final_state["stability"]),
                "binding": float(final_state["binding"]),
                "identity": float(final_state["identity"]),
                "n_mutations": int(final_state["n_mutations"]),
                "mutations": final_state["mutations"],
                "mut_positions": final_state.get("mut_positions", []),
                "uncertainty": float(final_state.get("uncertainty", 0.0)),
                "ood_distance": float(final_state.get("ood_distance", 0.0)),
                "trajectory": trial_trajectory,
                "meets_constraints": (
                    float(final_state["identity"]) >= min_identity
                    and float(final_state["stability"]) >= min_stability
                    and float(final_state["binding"]) >= min_binding
                ),
            }
            all_candidates.append(candidate)

    all_candidates.sort(key=lambda c: int(c.get("trial_idx", c.get("candidate_id", 0))))
    for idx, cand in enumerate(all_candidates, start=1):
        cand["candidate_id"] = idx
    return all_candidates


# === GENERATIVE DESIGN ENGINE ===
def run_generative_design(constraints: dict, n_candidates: int = 6, n_steps: Optional[int] = None):
    """
    Generative Design Mode: Like mechanical CAD topology optimization.

    User specifies CONSTRAINTS (not solutions):
    - Physics: min stability, binding thresholds
    - Geometry: locked residues, protected regions
    - Material: identity level (how "human-like")
    - Manufacturing: delivery method constraints

    AI generates MULTIPLE diverse solutions exploring the constraint space.
    """
    walker = ManifoldWalker(embedder)
    target_mut = constraints.get('target_mutation', 'R175H')

    cancer_seq = apply_mutation(P53_WT, target_mut)
    if cancer_seq is None:
        cancer_seq = P53_WT

    AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

    # Extract constraints
    min_identity = constraints.get('min_identity', 90.0)
    min_stability = constraints.get('min_stability', -0.3)
    min_binding = constraints.get('min_binding', 5.0)
    locked_positions = constraints.get('locked_positions', [248, 273])
    delivery_method = constraints.get('delivery_method', 'gene_therapy')
    exploration_diversity = constraints.get('diversity', 0.5)
    calibration_profile = get_score_calibration_profile()
    uncertainty_weight = float(constraints.get("uncertainty_penalty_weight", 0.8))
    ood_rank_weight = float(constraints.get("ood_rank_penalty_weight", 1.25))
    ood_radius = float(constraints.get("ood_trust_radius", 1.75))
    ood_loss_weight = float(constraints.get("ood_loss_weight", 12.0))
    mc_dropout_samples = int(max(2, constraints.get("mc_dropout_samples", 8)))
    if n_steps is None:
        n_steps = int(constraints.get("optimization_steps", 150))
    n_steps = int(np.clip(n_steps, 20, 500))

    # Adjust identity based on delivery method
    if delivery_method == 'protein_therapy':
        min_identity = max(min_identity, 95.0)  # Stricter for direct protein
    elif delivery_method == 'mrna_therapy':
        min_identity = max(min_identity, 92.0)

    all_candidates = []
    emb_target_ref = embedder.get_embeddings(cancer_seq).detach()
    emb_wt = embedder.get_embeddings(P53_WT).detach()
    with torch.no_grad():
        z_target_ref, _, _ = embedder.latent_forward_ascent(emb_target_ref)
        pooled_target_ref = z_target_ref.mean(dim=1)
        if pooled_target_ref.shape[-1] != oracle.input_dim:
            pooled_target_ref = pooled_target_ref[:, :oracle.input_dim]

    wt_aa_indices = []
    for aa in P53_WT:
        aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
        if aa_id in AA_IDS:
            wt_aa_indices.append(AA_IDS.index(aa_id))
        else:
            wt_aa_indices.append(0)
    wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb_target_ref.device)

    for candidate_idx in range(n_candidates):
        # DIVERSITY: Each candidate explores different regions via:
        # 1. Different random seeds for exploration noise
        # 2. Different weighting between objectives
        # 3. Different starting perturbations

        torch.manual_seed(42 + candidate_idx * 17)
        np.random.seed(42 + candidate_idx * 17)

        # Vary objective weights for diversity on Pareto frontier
        # Candidate 0: Balanced
        # Candidate 1: Favor stability
        # Candidate 2: Favor binding
        # Candidate 3: Favor function
        # etc.
        weight_profiles = [
            {'function': 4.0, 'stability': 8.0, 'binding': 2.5, 'name': 'Balanced'},
            {'function': 2.0, 'stability': 15.0, 'binding': 2.0, 'name': 'Stability-First'},
            {'function': 3.0, 'stability': 5.0, 'binding': 8.0, 'name': 'Binding-Optimized'},
            {'function': 8.0, 'stability': 4.0, 'binding': 3.0, 'name': 'Function-Maximized'},
            {'function': 5.0, 'stability': 10.0, 'binding': 5.0, 'name': 'Conservative'},
            {'function': 6.0, 'stability': 6.0, 'binding': 6.0, 'name': 'Experimental'},
        ]

        profile = weight_profiles[candidate_idx % len(weight_profiles)]

        emb = emb_target_ref.clone().detach().requires_grad_(True)

        # Add initial perturbation for diversity
        with torch.no_grad():
            perturbation = torch.randn_like(emb) * 0.05 * exploration_diversity
            emb.data += perturbation

        optimizer = torch.optim.Adam([emb], lr=0.04)
        locked_indices = [int(p) - 1 for p in locked_positions if p]

        best_valid_state = None
        best_valid_score = -float('inf')

        for step_idx in range(1, n_steps + 1):
            optimizer.zero_grad()

            z, logits, probs = embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            if pooled.shape[-1] != oracle.input_dim:
                pooled = pooled[:, :oracle.input_dim]

            score = oracle.model(pooled)

            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()

            dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
            hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)
            ood_distance_t = torch.norm(pooled - pooled_target_ref, p=2, dim=-1).mean()

            # LOSS with profile-specific weights
            loss = -score * profile['function']
            loss -= profile['stability'] * stability
            loss -= profile['binding'] * dna_force
            loss -= 3.0 * hydro_packing
            loss += ood_loss_weight * F.relu(ood_distance_t - ood_radius)

            # CONSTRAINT ENFORCEMENT (hard constraints from user)
            probs_aa = F.softmax(logits_aa[0], dim=-1)
            wt_probs = probs_aa[torch.arange(len(P53_WT), device=emb.device), wt_aa_tensor]
            mutation_prob = 1.0 - wt_probs
            expected_mutations = mutation_prob.sum()

            with torch.no_grad():
                decoded_ids = torch.argmax(probs_aa, dim=-1)
                n_mutations = (decoded_ids != wt_aa_tensor).sum().item()
                seq_identity = 100.0 * (1.0 - n_mutations / len(P53_WT))

            # HARD CONSTRAINT: Identity must meet minimum
            max_mutations = int(len(P53_WT) * (100 - min_identity) / 100)
            loss += 50.0 * F.relu(expected_mutations - max_mutations)

            # HARD CONSTRAINT: Identity barrier
            if seq_identity < min_identity - 5:
                loss += 500.0 * (min_identity - 5 - seq_identity)

            # HARD CONSTRAINT: Stability floor
            if stability.item() < min_stability:
                loss += 100.0 * (min_stability - stability)
            if dna_force.item() < min_binding:
                loss += 80.0 * (min_binding - dna_force)

            # LOCKED POSITIONS
            if locked_indices:
                loss += 500.0 * F.mse_loss(emb[:, locked_indices, :], emb_wt[:, locked_indices, :])

            # Regularization
            dist_l1 = torch.norm(emb - emb_wt, p=1) / emb.numel()
            loss += 40.0 * dist_l1

            loss.backward()
            optimizer.step()

            # Track best valid state
            with torch.no_grad():
                if (
                    seq_identity >= min_identity
                    and stability.item() >= min_stability
                    and dna_force.item() >= min_binding
                ):
                    score_bundle = build_trust_adjusted_score(
                        raw_score=float(score.item()),
                        pooled=pooled,
                        pooled_ref=pooled_target_ref,
                        calibration_profile=calibration_profile,
                        uncertainty_weight=uncertainty_weight,
                        ood_rank_weight=ood_rank_weight,
                        ood_radius=ood_radius,
                        mc_samples=mc_dropout_samples,
                    )
                    if score_bundle["score_adjusted"] > best_valid_score:
                        best_valid_score = score_bundle["score_adjusted"]
                        top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                        top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                        tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                        best_valid_state = {
                            'sequence': "".join(tokens)[:len(P53_WT)],
                            'score': score_bundle["score_adjusted"],
                            'score_raw': score_bundle["score_raw"],
                            'score_calibrated': score_bundle["score_calibrated"],
                            'stability': stability.item(),
                            'binding': dna_force.item(),
                            'identity': seq_identity,
                            'n_mutations': n_mutations,
                            'uncertainty': score_bundle["uncertainty"],
                            'ood_distance': score_bundle["ood_distance"],
                        }

        # Final decode
        with torch.no_grad():
            z, logits, probs = embedder.latent_forward_ascent(emb)
            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()
            pooled = z.mean(dim=1)
            if pooled.shape[-1] != oracle.input_dim:
                pooled = pooled[:, :oracle.input_dim]
            score = oracle.model(pooled)
            dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)

            probs_aa = F.softmax(logits_aa[0], dim=-1)
            decoded_ids = torch.argmax(probs_aa, dim=-1)
            n_mutations = (decoded_ids != wt_aa_tensor).sum().item()
            seq_identity = 100.0 * (1.0 - n_mutations / len(P53_WT))

            top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
            top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
            tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
            final_seq = "".join(tokens)[:len(P53_WT)]

            muts = [f"{P53_WT[j]}{j+1}{final_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != final_seq[j]]
            final_score_bundle = build_trust_adjusted_score(
                raw_score=float(score.item()),
                pooled=pooled,
                pooled_ref=pooled_target_ref,
                calibration_profile=calibration_profile,
                uncertainty_weight=uncertainty_weight,
                ood_rank_weight=ood_rank_weight,
                ood_radius=ood_radius,
                mc_samples=mc_dropout_samples,
            )

        # Use best valid state if final doesn't meet constraints
        if best_valid_state and (
            seq_identity < min_identity
            or stability.item() < min_stability
            or dna_force.item() < min_binding
        ):
            candidate = {
                'candidate_id': candidate_idx + 1,
                'profile': profile['name'],
                'sequence': best_valid_state['sequence'],
                'score': best_valid_state['score'],
                'score_raw': best_valid_state.get('score_raw', best_valid_state['score']),
                'score_calibrated': best_valid_state.get('score_calibrated', best_valid_state['score']),
                'stability': best_valid_state['stability'],
                'binding': best_valid_state['binding'],
                'identity': best_valid_state['identity'],
                'n_mutations': best_valid_state['n_mutations'],
                'uncertainty': float(best_valid_state.get('uncertainty', 0.0)),
                'ood_distance': float(best_valid_state.get('ood_distance', 0.0)),
                'mutations': [f"{P53_WT[j]}{j+1}{best_valid_state['sequence'][j]}"
                             for j in range(len(P53_WT))
                             if j < len(best_valid_state['sequence']) and P53_WT[j] != best_valid_state['sequence'][j]],
                'meets_constraints': True
            }
        else:
            candidate = {
                'candidate_id': candidate_idx + 1,
                'profile': profile['name'],
                'sequence': final_seq,
                'score': final_score_bundle['score_adjusted'],
                'score_raw': final_score_bundle['score_raw'],
                'score_calibrated': final_score_bundle['score_calibrated'],
                'stability': stability.item(),
                'binding': dna_force.item(),
                'identity': seq_identity,
                'n_mutations': n_mutations,
                'uncertainty': float(final_score_bundle.get('uncertainty', 0.0)),
                'ood_distance': float(final_score_bundle.get('ood_distance', 0.0)),
                'mutations': muts,
                'meets_constraints': (
                    seq_identity >= min_identity
                    and stability.item() >= min_stability
                    and dna_force.item() >= min_binding
                )
            }

        all_candidates.append(candidate)

    return all_candidates

# --- KNOWN EXPERIMENTAL RESCUES (Literature-validated) ---
KNOWN_RESCUES = {
    "R175H": {
        "published": ["N239Y", "N268D", "H178Y"],
        "source": "## Cell Reports, 2020",
        "notes": "Second-site suppressors identified via yeast screening"
    },
    "Y220C": {
        "published": ["PhiKan083", "PK7088"],  # Small molecules, not mutations
        "source": "Nature, 2008 (Fersht lab)",
        "notes": "Small molecule stabilizers bind Y220C cavity"
    },
    "G245S": {
        "published": ["R249M", "T123A"],
        "source": "PNAS, 2019",
        "notes": "Computational + experimental validation"
    },
    "R248Q": {
        "published": ["H115N", "T123P"],
        "source": "Oncogene, 2017",
        "notes": "Restores DNA binding via allosteric mechanism"
    },
    "R273H": {
        "published": ["T284R", "S240R"],
        "source": "JBC, 2018",
        "notes": "Charge compensation rescues"
    },
    "R249S": {
        "published": ["H168R", "T123A"],
        "source": "Cancer Research, 2016",
        "notes": "Loop stabilization rescues"
    }
}

# --- SESSION STATE INITIALIZATION ---
if 'results' not in st.session_state:
    st.session_state['results'] = None
if 'target_mut_saved' not in st.session_state:
    st.session_state['target_mut_saved'] = None

# --- TABS (2 focused views) ---
tab1, tab2 = st.tabs(["🧬 Design Studio", "📊 Analysis & Drugs"])

# Old Quick Design tab removed - using Studio as main interface
# === TAB 1: DESIGN STUDIO (Main Design Interface) ===
with tab1:
    # === GENERATIVE DESIGN CAD MODE ===
    st.markdown("""
    <div class="studio-banner">
        <p class="studio-kicker">Generative Protein Engineering</p>
        <p class="gen-design-header">proteoMgCAD Studio</p>
        <p class="studio-subtitle">Define constraints. AI generates optimal second-site rescues. Topology optimization for proteins.</p>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("### 📦 Result Source")
    result_source_mode = st.radio(
        "Choose how to populate candidate results",
        [
            "Local campaign artifacts (auto latest)",
            "Live session",
            "Upload run bundle",
        ],
        horizontal=True,
        key="gd_result_source_mode",
    )

    if result_source_mode == "Local campaign artifacts (auto latest)":
        store = CampaignStore()
        runs_df = store.list_runs()
        if runs_df.empty:
            st.info("No local campaign artifacts found yet in `data/campaigns`. Run `p53cad campaign-run` first.")
        else:
            run_options = runs_df["run_id"].astype(str).tolist()
            default_idx = 0
            previous_run = st.session_state.get("gd_artifact_run_id")
            if previous_run in run_options:
                default_idx = run_options.index(previous_run)
            selected_run = st.selectbox(
                "Campaign run",
                options=run_options,
                index=default_idx,
                key="gd_selected_campaign_run",
            )
            load_now = st.button("Load selected campaign results", key="gd_load_campaign_results_btn")
            if load_now or st.session_state.get("_gd_loaded_campaign_run") != selected_run:
                try:
                    bundle = load_campaign_bundle_from_disk(selected_run)
                    if apply_campaign_bundle_to_session(bundle, source_label=f"local:{selected_run}"):
                        st.session_state["_gd_loaded_campaign_run"] = selected_run
                        st.success(f"Loaded campaign artifacts: `{selected_run}`")
                    else:
                        st.warning(f"Run `{selected_run}` contains no candidates to display.")
                except Exception as err:
                    st.error(f"Failed to load campaign artifacts: {err}")
    elif result_source_mode == "Upload run bundle":
        upload = st.file_uploader(
            "Upload campaign run bundle (.zip)",
            type=["zip"],
            key="gd_upload_bundle",
        )
        if upload is not None and st.button("Load uploaded bundle", key="gd_load_uploaded_bundle_btn"):
            try:
                bundle = load_campaign_bundle_from_zip(upload.getvalue())
                if apply_campaign_bundle_to_session(bundle, source_label=f"upload:{upload.name}"):
                    st.success(f"Loaded bundle `{upload.name}` ({bundle.get('run_id', 'uploaded_bundle')}).")
                else:
                    st.warning("Uploaded bundle did not contain displayable candidates.")
            except Exception as err:
                st.error(f"Failed to parse uploaded bundle: {err}")
    else:
        st.caption("Live session mode uses in-memory candidates generated during this app run.")

    # === CONSTRAINT SPECIFICATION PANEL ===
    st.markdown("##  Define Design Constraints")
    st.markdown("*Specify what your design MUST achieve. The AI explores the solution space to find optimal rescue mutations.*")

    const_col1, const_col2, const_col3 = st.columns(3)

    with const_col1:
        st.markdown('<div class="constraint-card">', unsafe_allow_html=True)
        st.markdown("###  Physics Constraints")
        st.caption("*Like load/stress requirements in mechanical CAD*")

        gd_min_stability = st.slider(
            "Minimum Stability (PLL)",
            min_value=-0.5, max_value=0.0, value=-0.2,
            help="Minimum folding stability. Higher = stricter (more stable required)"
        )

        gd_min_binding = st.slider(
            "Minimum DNA Binding",
            min_value=0.0, max_value=10.0, value=5.0,
            help="Minimum DNA recruitment force. Higher = must bind DNA better"
        )

        gd_min_function = st.slider(
            "Minimum Function Score",
            min_value=-1.0, max_value=0.5, value=-0.2,
            help="Minimum rescue score to be considered valid"
        )
        st.markdown('</div>', unsafe_allow_html=True)

    with const_col2:
        st.markdown('<div class="constraint-card">', unsafe_allow_html=True)
        st.markdown("###  Geometry Constraints")
        st.caption("*Like fixed supports in mechanical CAD*")

        gd_target = st.selectbox(
            "Target Cancer Mutation",
            ALL_MUTATIONS[:20],
            index=0,
            help="The cancer mutation to rescue"
        )

        gd_locked = st.multiselect(
            "Locked Positions (Cannot Mutate)",
            options=list(range(94, 293)),
            default=[248, 273, 175],
            help="Critical positions that must remain unchanged"
        )

        gd_protected_regions = st.multiselect(
            "Protected Regions",
            ["L1 Loop (112-124)", "L2 Loop (163-195)", "L3 Loop (236-251)", "Zinc Site (176-179, 238-242)"],
            default=["Zinc Site (176-179, 238-242)"],
            help="Entire regions to protect from mutation"
        )
        st.markdown('</div>', unsafe_allow_html=True)

    with const_col3:
        st.markdown('<div class="constraint-card">', unsafe_allow_html=True)
        st.markdown("###  Material Constraints")
        st.caption("*Like material selection in mechanical CAD*")

        gd_identity = st.slider(
            "Minimum Sequence Identity",
            min_value=80.0, max_value=99.0, value=92.0,
            help="How 'human-like' the design must be. Higher = fewer mutations allowed"
        )

        gd_delivery = st.selectbox(
            "Delivery Method",
            ["Gene Therapy (AAV)", "mRNA Therapy", "Protein Therapy (Direct)"],
            index=0,
            help="Manufacturing constraint: affects identity requirements"
        )

        gd_diversity = st.slider(
            "Exploration Diversity",
            min_value=0.0, max_value=1.0, value=0.5,
            help="How different should candidates be from each other?"
        )
        st.markdown('</div>', unsafe_allow_html=True)

    # === GENERATION CONTROLS ===
    st.markdown("---")
    gen_col1, gen_col2 = st.columns([1, 3])

    with gen_col1:
        gd_n_candidates = st.selectbox(
            "Number of Candidates",
            [3, 4, 5, 6, 8, 10],
            index=3,
            help="How many diverse solutions to generate"
        )
        gd_opt_steps = st.slider(
            "Optimization Steps",
            min_value=40,
            max_value=400,
            value=120,
            step=10,
            help="Gradient steps per candidate. Higher can improve quality but increases runtime."
        )
        gd_min_gain = st.slider(
            "Min Score Gain vs Target",
            min_value=0.00,
            max_value=1.00,
            value=0.05,
            step=0.01,
            help="Filter out baseline-like candidates unless they exceed this score gain over target-only baseline."
        )
        gd_restarts = st.slider(
            "Restarts Per Profile",
            min_value=1,
            max_value=4,
            value=2,
            step=1,
            help="Independent optimization restarts per objective profile to escape local minima."
        )
        gd_novelty_weight = st.slider(
            "Novelty Weight",
            min_value=0.0,
            max_value=1.0,
            value=0.35,
            step=0.05,
            help="Higher values prioritize mutation-pattern diversity in final reranking."
        )

    with gen_col2:
        st.info(f"""
        **Generation Preview:**
        - Target: **{gd_target}** rescue
        - Identity: ≥{gd_identity:.0f}% (max {int(len(P53_WT) * (100 - gd_identity) / 100)} mutations)
        - Locked: {len(gd_locked)} positions
        - Delivery: {gd_delivery}
        - Candidates: {gd_n_candidates} diverse solutions
        - Steps per candidate: {gd_opt_steps}
        - Min gain filter: +{gd_min_gain:.2f} vs target baseline
        - Restarts/profile: {gd_restarts}
        - Novelty weight: {gd_novelty_weight:.2f}
        - Scoring: calibrated + uncertainty/OOD trust penalty
        """)

    # === VISUALIZATION MODE ===
    gd_live_mode = st.toggle(" Watch Live Optimization", value=True,
                              help="Like watching CAD topology optimization - see proteins being built in real-time")
    gd_structure_checkpoints = st.checkbox(
        "Capture real structure checkpoints (slow, strict mode)",
        value=False,
        help="When enabled, the viewer requests real ESMFold checkpoints during optimization. "
        "If unavailable, structure trajectory is hidden instead of using synthetic morphing.",
    )

    # === GENERATE BUTTON ===
    if st.button(" GENERATE CANDIDATE DESIGNS", type="primary", width="stretch"):
        # Parse constraints
        delivery_map = {
            "Gene Therapy (AAV)": "gene_therapy",
            "mRNA Therapy": "mrna_therapy",
            "Protein Therapy (Direct)": "protein_therapy"
        }

        constraints = {
            'target_mutation': gd_target,
            'min_identity': gd_identity,
            'min_stability': gd_min_stability,
            'min_binding': gd_min_binding,
            'locked_positions': gd_locked,
            'delivery_method': delivery_map[gd_delivery],
            'diversity': gd_diversity,
            'optimization_steps': gd_opt_steps,
            'min_score_gain': gd_min_gain,
            'restarts_per_profile': gd_restarts,
            'novelty_weight': gd_novelty_weight,
            'ood_trust_radius': 1.75,
            'ood_loss_weight': 12.0,
            'ood_rank_penalty_weight': 1.25,
            'uncertainty_penalty_weight': 0.8,
            'mc_dropout_samples': 8,
            'trust_eval_stride': 20,
            'early_stop_patience': 8,
            'early_stop_min_delta': 1e-3,
            'early_stop_warmup_steps': 40,
            'batch_trial_size': 4,
            'optimizer_lr': 0.03,
            'enable_successive_halving': True,
            'successive_halving_warmup_steps': int(np.clip(max(20, gd_opt_steps // 4), 20, gd_opt_steps)),
            'successive_halving_keep_ratio': 0.5,
            'capture_structure_checkpoints': bool(gd_structure_checkpoints),
            'structure_checkpoint_stride': 40,
        }

        if gd_live_mode:
            # === LIVE VISUALIZATION MODE ===
            st.markdown("---")
            st.markdown("##  Live Optimization Viewer")
            st.markdown("*Watch the AI build rescue proteins in real-time*")

            # Two-column layout: protein viewer on left, metrics on right
            viewer_col, metrics_col = st.columns([3, 2])

            with viewer_col:
                st.markdown("### 🧬 Live Protein Structure")
                structure_placeholder = st.empty()

            with metrics_col:
                status_placeholder = st.empty()
                metrics_placeholder = st.empty()
                mutations_placeholder = st.empty()
                progress_placeholder = st.empty()

            st.markdown("---")
            trajectory_placeholder = st.empty()

            # Show initial protein structure immediately
            try:
                with open("data/raw/p53_wt.pdb", "r") as f:
                    initial_pdb = f.read()

                initial_html = f"""
                <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                <div id="init_viewer" style="width:100%; height:400px; border-radius:12px; border:2px solid #6366F1; background:linear-gradient(135deg,#F8FCFF 0%,#EAF3FF 100%);"></div>
                <script>
                    let viewer = $3Dmol.createViewer('init_viewer', {{backgroundColor: '0xF8FCFF'}});
                    let pdb = `{initial_pdb.replace(chr(10), chr(92) + 'n').replace('`', '')}`;
                    viewer.addModel(pdb, 'pdb');
                    viewer.setStyle({{}}, {{cartoon: {{color: '0x64748B', opacity: 0.8}}}});
                    viewer.addStyle({{resi: [112,113,114,115,116,117,118,119,120,121,122,123,124]}},
                                   {{cartoon: {{color: '0x0EA5E9', opacity: 0.95}}}});
                    viewer.addStyle({{resi: [236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251]}},
                                   {{cartoon: {{color: '0x10B981', opacity: 0.95}}}});
                    viewer.zoomTo();
                    viewer.spin('y', 0.5);
                    viewer.render();
                </script>
                <p style="text-align:center; color:#6366F1; font-weight:600; margin-top:8px;">
                    ⏳ Starting optimization... Watch the structure evolve!
                </p>
                """
                import streamlit.components.v1 as components
                with structure_placeholder.container():
                    components.html(initial_html, height=450)
            except Exception as e:
                structure_placeholder.warning(f"Could not load initial structure: {e}")

            all_candidates = []
            weight_profiles = [
                {'function': 4.0, 'stability': 8.0, 'binding': 2.5, 'name': 'Balanced', 'color': '#2563EB'},
                {'function': 2.0, 'stability': 15.0, 'binding': 2.0, 'name': 'Stability-First', 'color': '#0EA5E9'},
                {'function': 3.0, 'stability': 5.0, 'binding': 8.0, 'name': 'Binding-Optimized', 'color': '#10B981'},
                {'function': 8.0, 'stability': 4.0, 'binding': 3.0, 'name': 'Function-Maximized', 'color': '#14B8A6'},
                {'function': 5.0, 'stability': 10.0, 'binding': 5.0, 'name': 'Conservative', 'color': '#1D4ED8'},
                {'function': 6.0, 'stability': 6.0, 'binding': 6.0, 'name': 'Experimental', 'color': '#0891B2'},
            ]
            profile_attempt_counts = {p["name"]: 0 for p in weight_profiles}
            total_trials = int(max(gd_n_candidates, 1) * max(gd_restarts, 1))
            score_calibration = get_score_calibration_profile()
            uncertainty_weight = float(constraints.get("uncertainty_penalty_weight", 0.8))
            ood_rank_weight = float(constraints.get("ood_rank_penalty_weight", 1.25))
            ood_radius = float(constraints.get("ood_trust_radius", 1.75))
            ood_loss_weight = float(constraints.get("ood_loss_weight", 12.0))
            mc_dropout_samples = int(max(2, constraints.get("mc_dropout_samples", 8)))
            trust_eval_stride = int(max(5, constraints.get("trust_eval_stride", 20)))
            early_stop_patience = int(max(1, constraints.get("early_stop_patience", 8)))
            early_stop_min_delta = float(max(0.0, constraints.get("early_stop_min_delta", 1e-3)))
            early_stop_warmup_steps = int(max(20, constraints.get("early_stop_warmup_steps", 40)))
            capture_checkpoints = bool(constraints.get("capture_structure_checkpoints", False))
            checkpoint_stride = int(max(10, constraints.get("structure_checkpoint_stride", 40)))
            cancer_seq_base = apply_mutation(P53_WT, gd_target)
            if cancer_seq_base is None:
                cancer_seq_base = P53_WT
            AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
            emb_target_ref = embedder.get_embeddings(cancer_seq_base).detach()
            emb_wt = embedder.get_embeddings(P53_WT).detach()
            with torch.no_grad():
                z_target_ref, _, _ = embedder.latent_forward_ascent(emb_target_ref)
                pooled_target_ref = z_target_ref.mean(dim=1)
                if pooled_target_ref.shape[-1] != oracle.input_dim:
                    pooled_target_ref = pooled_target_ref[:, :oracle.input_dim]

            wt_aa_indices = []
            for aa in P53_WT:
                aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
                if aa_id in AA_IDS:
                    wt_aa_indices.append(AA_IDS.index(aa_id))
                else:
                    wt_aa_indices.append(0)
            wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb_target_ref.device)

            trial_order = list(range(total_trials))
            if constraints.get("enable_successive_halving", True) and total_trials > gd_n_candidates:
                screen_steps = int(np.clip(max(20, gd_opt_steps // 4), 20, gd_opt_steps))
                keep_ratio = float(np.clip(constraints.get("successive_halving_keep_ratio", 0.5), 0.25, 1.0))
                keep_count = max(int(max(gd_n_candidates, 1)), int(np.ceil(total_trials * keep_ratio)))
                screen_constraints = dict(constraints)
                screen_constraints["capture_structure_checkpoints"] = False
                screen_constraints["enable_successive_halving"] = False
                screen_constraints["mc_dropout_samples"] = 2
                with st.spinner(
                    f"Speed pre-screen: evaluating {total_trials} trial plans for {screen_steps} warmup steps..."
                ):
                    screened = run_generative_design_live(
                        screen_constraints,
                        n_candidates=total_trials,
                        n_steps=screen_steps,
                    )
                screened_sorted = sorted(screened, key=candidate_rank_value, reverse=True)
                trial_order = sorted(
                    max(0, int(c.get("trial_idx", c.get("candidate_id", 1))) - 1)
                    for c in screened_sorted[:keep_count]
                )
                st.info(
                    f"Successive halving kept {len(trial_order)}/{total_trials} trial plans "
                    f"(warmup {screen_steps} steps, keep ratio {keep_ratio:.2f})."
                )
            active_trials = max(len(trial_order), 1)

            # Run generation with live updates
            for live_idx, trial_idx in enumerate(trial_order):
                torch.manual_seed(42 + trial_idx * 17)
                np.random.seed(42 + trial_idx * 17)
                profile = weight_profiles[trial_idx % len(weight_profiles)]
                profile_attempt_counts[profile["name"]] += 1
                restart_idx = profile_attempt_counts[profile["name"]]

                status_placeholder.markdown(f"""
                ###  Building Trial {live_idx + 1}/{active_trials}
                **Strategy:** {profile['name']} | Restart {restart_idx}/{gd_restarts}
                """)

                emb = emb_target_ref.clone().detach().requires_grad_(True)

                with torch.no_grad():
                    perturbation = torch.randn_like(emb) * 0.05 * gd_diversity
                    emb.data += perturbation

                optimizer = torch.optim.Adam([emb], lr=float(constraints.get("optimizer_lr", 0.03)))
                locked_indices = [int(p) - 1 for p in gd_locked if p]

                n_steps = int(np.clip(gd_opt_steps, 20, 500))
                trajectory = []
                best_valid_state = None
                best_valid_score = -float('inf')
                checks_without_improve = 0
                cached_uncertainty = 0.0
                mut_positions = []
                previous_render_mut_positions = []
                last_render_step = 0
                render_stride = 10

                for step_idx in range(1, n_steps + 1):
                    optimizer.zero_grad()

                    z, logits, probs = embedder.latent_forward_ascent(emb)
                    pooled = z.mean(dim=1)
                    if pooled.shape[-1] != oracle.input_dim:
                        pooled = pooled[:, :oracle.input_dim]

                    score = oracle.model(pooled)

                    logits_aa = logits[:, :, AA_IDS]
                    log_probs = F.log_softmax(logits_aa, dim=-1)
                    stability = log_probs.max(dim=-1).values.mean()
                    dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
                    hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)
                    ood_distance_t = torch.norm(pooled - pooled_target_ref, p=2, dim=-1).mean()

                    loss = -score * profile['function']
                    loss -= profile['stability'] * stability
                    loss -= profile['binding'] * dna_force
                    loss -= 3.0 * hydro_packing
                    loss += ood_loss_weight * F.relu(ood_distance_t - ood_radius)

                    probs_aa = F.softmax(logits_aa[0], dim=-1)
                    wt_probs = probs_aa[torch.arange(len(P53_WT), device=emb.device), wt_aa_tensor]
                    mutation_prob = 1.0 - wt_probs
                    expected_mutations = mutation_prob.sum()

                    with torch.no_grad():
                        decoded_ids = torch.argmax(probs_aa, dim=-1)
                        n_mutations = (decoded_ids != wt_aa_tensor).sum().item()
                        seq_identity = 100.0 * (1.0 - n_mutations / len(P53_WT))

                    max_mutations = int(len(P53_WT) * (100 - gd_identity) / 100)
                    loss += 50.0 * F.relu(expected_mutations - max_mutations)

                    if seq_identity < gd_identity - 5:
                        loss += 500.0 * (gd_identity - 5 - seq_identity)
                    if stability.item() < gd_min_stability:
                        loss += 100.0 * (gd_min_stability - stability)
                    if dna_force.item() < gd_min_binding:
                        loss += 80.0 * (gd_min_binding - dna_force)
                    if locked_indices:
                        loss += 500.0 * F.mse_loss(emb[:, locked_indices, :], emb_wt[:, locked_indices, :])

                    dist_l1 = torch.norm(emb - emb_wt, p=1) / emb.numel()
                    loss += 40.0 * dist_l1

                    loss.backward()
                    optimizer.step()

                    # === LIVE UPDATE every 5 steps ===
                    if step_idx % 5 == 0 or step_idx == n_steps:
                        improved_this_checkpoint = False
                        with torch.no_grad():
                            top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                            top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                            tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                            current_seq = "".join(tokens)[:len(P53_WT)]
                            muts = [f"{P53_WT[j]}{j+1}{current_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != current_seq[j]]
                            mut_positions = [int(''.join(filter(str.isdigit, m))) for m in muts if m]
                            do_full_trust_eval = (
                                step_idx == n_steps
                                or step_idx == 5
                                or (step_idx % trust_eval_stride == 0)
                            )
                            raw_score = float(score.item())
                            current_ood = float(ood_distance_t.item())
                            calibrated_score = calibrate_oracle_score(raw_score, score_calibration)
                            if do_full_trust_eval:
                                cached_uncertainty = estimate_oracle_uncertainty(
                                    pooled,
                                    mc_samples=mc_dropout_samples,
                                )
                            adjusted_score = calibrated_score - uncertainty_weight * cached_uncertainty - ood_rank_weight * max(
                                0.0, current_ood - ood_radius
                            )
                            score_bundle = {
                                "score_raw": raw_score,
                                "score_calibrated": float(calibrated_score),
                                "score_adjusted": float(
                                    np.clip(
                                        adjusted_score,
                                        float(score_calibration.get("clip_low", -3.0)),
                                        float(score_calibration.get("clip_high", 4.0)),
                                    )
                                ),
                                "uncertainty": float(cached_uncertainty),
                                "ood_distance": float(current_ood),
                            }
                            structure_checkpoint = None
                            first_capture_step = 5  # first recorded step in this loop cadence
                            should_capture = (
                                capture_checkpoints
                                and (
                                    step_idx == n_steps
                                    or step_idx == first_capture_step
                                    or (step_idx % checkpoint_stride == 0)
                                )
                            )
                            if should_capture:
                                structure_checkpoint = get_real_structure_checkpoint(current_seq)

                            trajectory.append({
                                'step': step_idx, 'score': score_bundle['score_adjusted'], 'stability': stability.item(),
                                'score_raw': score_bundle['score_raw'],
                                'score_calibrated': score_bundle['score_calibrated'],
                                'binding': dna_force.item(), 'identity': seq_identity, 'n_mutations': n_mutations,
                                'mutations': muts, 'mut_positions': mut_positions, 'sequence': current_seq,
                                'uncertainty': score_bundle['uncertainty'],
                                'ood_distance': score_bundle['ood_distance'],
                                'lx': pooled[0, 0].item(), 'ly': pooled[0, 1].item(),
                                'lz': pooled[0, 2].item() if pooled.shape[-1] > 2 else score.item(),  # 3D latent
                                'loss_total': float(loss.item()),
                                'structure_checkpoint': structure_checkpoint,
                            })

                            if (
                                seq_identity >= gd_identity
                                and stability.item() >= gd_min_stability
                                and dna_force.item() >= gd_min_binding
                            ):
                                if score_bundle['score_adjusted'] > (best_valid_score + early_stop_min_delta):
                                    best_valid_score = score_bundle['score_adjusted']
                                    best_valid_state = trajectory[-1].copy()
                                    improved_this_checkpoint = True

                            stage_meta = get_optimization_stage(step_idx, n_steps)
                            # Update live metrics display
                            with metrics_placeholder.container():
                                st.markdown(
                                    f"""
                                    <div style="background:linear-gradient(135deg,#EAF3FF 0%,#DBF4FF 100%);
                                                border:1px solid #B7D4F7; border-left:5px solid {stage_meta['color']};
                                                border-radius:10px; padding:0.55rem 0.75rem; margin-bottom:0.7rem;">
                                        <p style="margin:0; color:#1E3A5F; font-weight:700; font-size:0.85rem;">{stage_meta['name']}</p>
                                        <p style="margin:0.15rem 0 0 0; color:#3B5C85; font-size:0.78rem;">{stage_meta['description']}</p>
                                    </div>
                                    """,
                                    unsafe_allow_html=True,
                                )
                                st.markdown(f"**Step {step_idx}/{n_steps}**")
                                mc1, mc2 = st.columns(2)
                                mc1.metric("Score", f"{score_bundle['score_adjusted']:.3f}")
                                mc2.metric("Identity", f"{seq_identity:.1f}%")
                                mc3, mc4 = st.columns(2)
                                mc3.metric("Stability", f"{stability.item():.3f}")
                                mc4.metric("Uncertainty", f"{score_bundle['uncertainty']:.3f}")

                            # Update mutations list
                            with mutations_placeholder.container():
                                st.markdown("**Current Mutations:**")
                                for m in muts[:4]:
                                    st.write(f"• {m}")

                            # Update progress
                            overall_progress = (live_idx * n_steps + step_idx) / (active_trials * n_steps)
                            progress_placeholder.progress(
                                overall_progress,
                                text=f"Overall: {overall_progress*100:.0f}% | {stage_meta['name']}",
                            )

                            # Update trajectory plot (only every 10 steps to reduce flicker)
                            if len(trajectory) > 1 and step_idx % render_stride == 0:
                                traj_df = pd.DataFrame(trajectory)
                                fig_traj = go.Figure()
                                fig_traj.add_trace(go.Scatter(
                                    x=traj_df['step'],
                                    y=traj_df['score'],
                                    mode='lines+markers',
                                    name='Rescue Score',
                                    line=dict(color=profile['color'], width=3, shape='spline', smoothing=0.7),
                                    marker=dict(size=5),
                                    fill='tozeroy',
                                    fillcolor='rgba(37, 99, 235, 0.14)'
                                ))
                                fig_traj.add_trace(go.Scatter(
                                    x=traj_df['step'],
                                    y=traj_df['stability'],
                                    mode='lines',
                                    name='Stability',
                                    line=dict(color='#38BDF8', width=2.3, dash='dot')
                                ))
                                fig_traj.add_trace(go.Scatter(
                                    x=traj_df['step'],
                                    y=traj_df['identity'],
                                    mode='lines',
                                    name='Identity (%)',
                                    line=dict(color='#1D4ED8', width=2, dash='dash'),
                                    yaxis='y2'
                                ))
                                fig_traj.add_hline(
                                    y=gd_min_stability,
                                    line_dash='dot',
                                    line_width=1.2,
                                    line_color='#0EA5E9',
                                    annotation_text='Stability floor',
                                    annotation_position='top left'
                                )
                                fig_traj.update_layout(
                                    height=320,
                                    margin=dict(l=24, r=20, t=70, b=98),
                                    title=f"Trial {live_idx+1} Optimization",
                                    xaxis_title="Step",
                                    yaxis_title="Score / Stability",
                                    yaxis2=dict(
                                        title="Identity %",
                                        overlaying='y',
                                        side='right',
                                        range=[75, 101],
                                        showgrid=False
                                    ),
                                    legend=dict(orientation="h", yanchor="top", y=-0.2, x=0)
                                )
                                render_plotly_in(trajectory_placeholder, fig_traj, chart_kind="mini", width="stretch")

                        if improved_this_checkpoint:
                            checks_without_improve = 0
                        else:
                            checks_without_improve += 1
                        if step_idx >= early_stop_warmup_steps and checks_without_improve >= early_stop_patience:
                            with metrics_placeholder.container():
                                st.caption(
                                    f"Early stop triggered after {checks_without_improve} non-improving checkpoints."
                                )
                            break

                # === 3D STRUCTURE - strict mode: real checkpoint structures only ===
                if step_idx % render_stride == 0 or step_idx == n_steps:
                    try:
                        stage_meta = get_optimization_stage(step_idx, n_steps)
                        progress_pct = int((step_idx / n_steps) * 100)
                        is_final = step_idx == n_steps
                        checkpoint_frames = [
                            row.get("structure_checkpoint")
                            for row in trajectory
                            if row.get("structure_checkpoint")
                        ]
                        # De-duplicate while preserving order.
                        dedup_frames: list[str] = []
                        seen_frames = set()
                        for frame in checkpoint_frames:
                            if frame and frame not in seen_frames:
                                seen_frames.add(frame)
                                dedup_frames.append(frame)

                        if not dedup_frames:
                            with structure_placeholder.container():
                                if not capture_checkpoints:
                                    st.info(
                                        "Live structure trajectory is disabled in strict mode. "
                                        "Enable `Capture real structure checkpoints` to render real frames."
                                    )
                                elif step_idx < 5:
                                    st.info("Warming up optimization before first checkpoint capture...")
                                elif step_idx < checkpoint_stride:
                                    st.info(
                                        f"Waiting for first real checkpoint (capturing every {checkpoint_stride} steps)."
                                    )
                                else:
                                    st.warning(
                                        "No real structure checkpoint returned yet. "
                                        "Check ESMFold API/network access and keep checkpoint capture enabled."
                                    )
                            last_render_step = step_idx
                            previous_render_mut_positions = list(mut_positions)
                            continue

                        frame_payload = json.dumps(dedup_frames)
                        mut_payload = json.dumps(mut_positions[:10])
                        focus_positions = compute_camera_focus_positions(
                            mut_positions,
                            previous_positions=previous_render_mut_positions,
                            max_points=10,
                        )
                        focus_payload = json.dumps(focus_positions)
                        grow_radius = 0.35 + (progress_pct / 100.0) * 0.35
                        base_opacity = 0.60 + (progress_pct / 100.0) * 0.25
                        viewer_id = f"morph_view_{trial_idx}_{step_idx}"
                        status_text = (
                            f"Finalized {profile['name']} candidate"
                            if is_final
                            else f"{profile['name']} | Step {step_idx}/{n_steps}"
                        )

                        live_3d_html = f"""
                        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px; gap:8px; flex-wrap:wrap;">
                        <div style="background:linear-gradient(135deg,#EEF5FF 0%,#DFF0FF 100%); color:#1E3A5F; border-left:4px solid {stage_meta['color']};
                                        padding:6px 10px; border-radius:8px; font-size:12px; font-weight:700;">
                                {stage_meta['name']}
                        </div>
                            <div style="color:#475569; font-size:12px; font-weight:600;">
                                {stage_meta['description']}
                            </div>
                        </div>
                        <div id="{viewer_id}" style="width:100%; height:382px; border-radius:12px; border:3px solid {profile["color"]}; background:linear-gradient(135deg,#F8FCFF 0%,#EAF3FF 100%);"></div>
                        <div style="height:10px; background:#D6E5F8; border-radius:6px; margin-top:8px; overflow:hidden;">
                            <div style="height:100%; width:{progress_pct}%; background:linear-gradient(90deg, {profile['color']} 0%, {stage_meta['color']} 100%);"></div>
                        </div>
                        <script>
                            (() => {{
                                const frames = {frame_payload};
                                const mutPositions = {mut_payload};
                                const focusPositions = {focus_payload};
                                const viewer = $3Dmol.createViewer('{viewer_id}', {{backgroundColor: '0xF8FCFF'}});
                                const profileColor = '{profile["color"]}';
                                const growRadius = {grow_radius:.3f};
                                const baseOpacity = {base_opacity:.3f};
                                const hasFocus = focusPositions.length > 0 || mutPositions.length > 0;
                                const focusSelection = focusPositions.length > 0
                                    ? {{resi: focusPositions}}
                                    : (mutPositions.length > 0 ? {{resi: mutPositions}} : {{}});

                                function applyStyles() {{
                                    viewer.setStyle({{}}, {{cartoon: {{color: '0x64748B', opacity: baseOpacity}}}});
                                    viewer.addStyle({{resi: [112,113,114,115,116,117,118,119,120,121,122,123,124]}},
                                                   {{cartoon: {{color: '0x0EA5E9', opacity: 0.95}}}});
                                    viewer.addStyle({{resi: [236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251]}},
                                                   {{cartoon: {{color: '0x10B981', opacity: 0.95}}}});
                                    mutPositions.forEach((p) => {{
                                        viewer.addStyle({{resi: p}}, {{
                                            cartoon: {{color: profileColor, opacity: 1.0}},
                                            stick: {{color: profileColor, radius: 0.28 + growRadius * 0.22}}
                                        }});
                                        viewer.addStyle({{resi: p, atom: 'CA'}}, {{
                                            sphere: {{color: profileColor, radius: growRadius, opacity: 0.95}}
                                        }});
                                    }});
                                }}

                                function focusCamera(durationMs) {{
                                    if (hasFocus) {{
                                        viewer.zoomTo(focusSelection, durationMs);
                                    }} else {{
                                        viewer.zoomTo({{}}, durationMs);
                                    }}
                                }}

                                function renderFrame(frameIndex) {{
                                    viewer.removeAllModels();
                                    viewer.addModel(frames[frameIndex], 'pdb');
                                    applyStyles();
                                    if (frameIndex === 0) {{
                                        focusCamera(320);
                                    }} else if (frameIndex % 2 === 0 && hasFocus) {{
                                        focusCamera(120);
                                    }}
                                    viewer.render();
                                }}

                                let frameIndex = 0;
                                renderFrame(frameIndex);
                                viewer.spin('y', hasFocus ? 0.55 : 1.0);

                                if (frames.length > 1) {{
                                    const interval = setInterval(() => {{
                                        frameIndex = (frameIndex + 1) % frames.length;
                                        renderFrame(frameIndex);
                                    }}, 95);
                                    setTimeout(() => clearInterval(interval), Math.max(2200, frames.length * 100));
                                }}
                            }})();
                        </script>
                        <p style="text-align:center; color:{profile['color']}; font-weight:600; margin-top:10px; font-size:14px;">
                            {status_text} | {len(mut_positions)} active mutation sites | {len(dedup_frames)} real structure checkpoints
                        </p>
                        """
                        with structure_placeholder.container():
                            components.html(live_3d_html, height=474)
                        last_render_step = step_idx
                        previous_render_mut_positions = list(mut_positions)
                    except Exception as e:
                        structure_placeholder.error(f"Could not render protein: {e}")

                # Store final candidate
                final_state = best_valid_state if best_valid_state else trajectory[-1]
                candidate = {
                    'candidate_id': live_idx + 1,
                    'profile': profile['name'],
                    'color': profile['color'],
                    'restart': restart_idx,
                    'trial_idx': trial_idx + 1,
                    'sequence': final_state['sequence'],
                    'score': final_state['score'],
                    'score_raw': final_state.get('score_raw', final_state['score']),
                    'score_calibrated': final_state.get('score_calibrated', final_state['score']),
                    'stability': final_state['stability'],
                    'binding': final_state['binding'],
                    'identity': final_state['identity'],
                    'n_mutations': final_state['n_mutations'],
                    'mutations': final_state['mutations'],
                    'mut_positions': final_state.get('mut_positions', []),
                    'uncertainty': float(final_state.get('uncertainty', 0.0)),
                    'ood_distance': float(final_state.get('ood_distance', 0.0)),
                    'trajectory': trajectory,
                    'meets_constraints': (
                        final_state['identity'] >= gd_identity
                        and final_state['stability'] >= gd_min_stability
                        and final_state['binding'] >= gd_min_binding
                    )
                }
                all_candidates.append(candidate)

            candidates = all_candidates
            if not candidates:
                logger.warning(
                    "Live generation produced zero candidates; attempting non-live fallback. "
                    "target=%s trials=%s steps=%s",
                    gd_target,
                    total_trials,
                    gd_opt_steps,
                )
                st.warning(
                    "Live pass returned zero candidates. Retrying once with robust fallback generation."
                )
                try:
                    candidates = run_generative_design_live(
                        constraints,
                        n_candidates=total_trials,
                        n_steps=gd_opt_steps,
                    )
                except Exception as fallback_exc:
                    logger.exception("Fallback generation failed: %s", fallback_exc)
                    candidates = []

            if not candidates:
                logger.error(
                    "Generation fallback still produced zero candidates; emitting target baseline candidate."
                )
                baseline_for_emit = compute_target_baseline_metrics(gd_target)
                fallback_seq = apply_mutation(P53_WT, gd_target) or P53_WT
                fallback_muts = [
                    f"{P53_WT[j]}{j+1}{fallback_seq[j]}"
                    for j in range(len(P53_WT))
                    if P53_WT[j] != fallback_seq[j]
                ]
                fallback_positions = [
                    int("".join(filter(str.isdigit, m)))
                    for m in fallback_muts
                    if any(ch.isdigit() for ch in m)
                ]
                candidates = [
                    {
                        "candidate_id": 1,
                        "profile": "Baseline",
                        "color": "#64748B",
                        "restart": 1,
                        "trial_idx": 1,
                        "sequence": fallback_seq,
                        "score": float(baseline_for_emit.get("score", 0.0)),
                        "score_raw": float(baseline_for_emit.get("score_raw", 0.0)),
                        "score_calibrated": float(baseline_for_emit.get("score_calibrated", 0.0)),
                        "stability": float(baseline_for_emit.get("stability", 0.0)),
                        "binding": float(baseline_for_emit.get("binding", 0.0)),
                        "identity": 100.0 * (1.0 - len(fallback_muts) / len(P53_WT)),
                        "n_mutations": len(fallback_muts),
                        "mutations": fallback_muts,
                        "mut_positions": fallback_positions,
                        "uncertainty": float(baseline_for_emit.get("uncertainty", 0.0)),
                        "ood_distance": float(baseline_for_emit.get("ood_distance", 0.0)),
                        "trajectory": [],
                        "meets_constraints": False,
                        "selection_reason": "baseline_emit_after_zero_candidate_failure",
                    }
                ]
                st.error(
                    "Both live and fallback generation returned no candidates. "
                    "Showing target baseline only so downstream analysis remains visible."
                )

            status_placeholder.success(f"Generated {len(candidates)} raw trial designs. Applying quality/diversity selection...")

        else:
            # Non-live mode (faster): show compact numeric progress without live structure rendering.
            total_trials = int(max(gd_n_candidates, 1) * max(gd_restarts, 1))
            compact_status = st.empty()
            compact_progress = st.progress(
                0.0,
                text=f"Generating {total_trials} trial designs ({gd_restarts} restarts/profile)...",
            )
            compact_metrics = st.empty()

            def _non_live_progress(update: Dict[str, Any]) -> None:
                cand_idx = int(update.get("candidate_idx", 0)) + 1
                cand_total = int(update.get("candidate_total", total_trials))
                step_idx = int(update.get("step", 0))
                total_steps = int(update.get("total_steps", max(gd_opt_steps, 1)))
                state = update.get("current_state", {}) or {}
                profile = (update.get("profile") or {}).get("name", "Unknown")

                denom = max(cand_total * total_steps, 1)
                overall = ((cand_idx - 1) * total_steps + step_idx) / denom
                compact_progress.progress(
                    float(np.clip(overall, 0.0, 1.0)),
                    text=(
                        f"Candidate {cand_idx}/{cand_total} ({profile}) "
                        f"- step {step_idx}/{total_steps}"
                    ),
                )
                compact_status.markdown(
                    f"**Progress:** candidate {cand_idx}/{cand_total} | profile `{profile}` | step {step_idx}/{total_steps}"
                )
                with compact_metrics.container():
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Score", f"{float(state.get('score', 0.0)):.3f}")
                    m2.metric("Identity", f"{float(state.get('identity', 0.0)):.1f}%")
                    m3.metric("Stability", f"{float(state.get('stability', 0.0)):.3f}")
                    m4.metric("Uncertainty", f"{float(state.get('uncertainty', 0.0)):.3f}")

            candidates = run_generative_design_live(
                constraints,
                n_candidates=total_trials,
                n_steps=gd_opt_steps,
                progress_callback=_non_live_progress,
            )
            compact_progress.progress(1.0, text="Generation complete.")

        baseline_metrics = compute_target_baseline_metrics(gd_target)
        candidates, filter_summary = filter_candidate_set(
            candidates=candidates,
            target_mutation=gd_target,
            baseline_score=baseline_metrics["score"],
            min_score_gain=gd_min_gain,
            min_rescue_mutations=1,
        )
        candidates, diversity_summary = select_diverse_candidates(
            candidates=candidates,
            target_mutation=gd_target,
            desired_count=gd_n_candidates,
            novelty_weight=gd_novelty_weight,
            overlap_penalty=0.20,
        )
        for rank_idx, cand in enumerate(candidates, start=1):
            cand["candidate_id"] = rank_idx

        if filter_summary["removed"] > 0:
            st.info(
                "Quality filter removed "
                f"{filter_summary['removed']} baseline-like candidates "
                f"(min gain +{gd_min_gain:.2f} vs trust-adjusted baseline {baseline_metrics['score']:.3f})."
            )
        if filter_summary["fallback_used"]:
            st.warning(
                "No candidate passed the quality filter; showing the top-scoring design only. "
                "Lower `Min Score Gain vs Target` or increase `Optimization Steps`."
            )
        st.info(
            "Diversity rerank selected "
            f"{len(candidates)} final candidates from {diversity_summary['pool_size']} filtered designs "
            f"(novelty weight {diversity_summary['novelty_weight']:.2f}, overlap penalty {diversity_summary['overlap_penalty']:.2f})."
        )
        st.caption(
            "Scores shown below are trust-adjusted: calibrated to DMS range, then penalized by MC-dropout "
            "uncertainty and latent OOD distance."
        )
        rejection_reason_text = ", ".join(
            f"{k}={v}" for k, v in filter_summary.get("rejection_reasons", {}).items()
        ) or "none"
        debug_df = pd.DataFrame(
            [
                {
                    "Baseline trust score": f"{baseline_metrics['score']:.3f}",
                    "Min gain threshold": f"+{gd_min_gain:.2f}",
                    "Filtered kept/total": f"{filter_summary['kept']}/{filter_summary['total']}",
                    "Fallback used": "yes" if filter_summary["fallback_used"] else "no",
                    "Rejected reasons": rejection_reason_text,
                }
            ]
        )
        st.dataframe(debug_df, hide_index=True, width="stretch")

        # Store in session state
        st.session_state['gd_candidates'] = candidates
        st.session_state['gd_constraints'] = constraints
        st.session_state['gd_filter_summary'] = filter_summary
        st.session_state['gd_target_baseline'] = baseline_metrics
        st.session_state['gd_diversity_summary'] = diversity_summary

    # === CANDIDATE GALLERY ===
    if 'gd_candidates' in st.session_state and st.session_state['gd_candidates']:
        candidates = st.session_state['gd_candidates']
        constraints = st.session_state.get('gd_constraints', {})

        st.markdown("---")
        st.markdown("##  Design Candidates Gallery")
        st.markdown("*AI-generated solutions exploring your constraint space. Like CAD design exploration.*")

        # Summary metrics
        valid_count = sum(1 for c in candidates if c['meets_constraints'])
        avg_score = np.mean([c['score'] for c in candidates])
        avg_identity = np.mean([c['identity'] for c in candidates])

        sum_col1, sum_col2, sum_col3, sum_col4 = st.columns(4)
        sum_col1.metric("Valid Designs", f"{valid_count}/{len(candidates)}")
        sum_col2.metric("Avg Trust Score", f"{avg_score:.3f}")
        sum_col3.metric("Avg Identity", f"{avg_identity:.1f}%")
        sum_col4.metric("Design Space Coverage", f"{len(set(c['profile'] for c in candidates))} profiles")

        if st.session_state.get("gd_source_mode") == "campaign_artifacts":
            run_id = st.session_state.get("gd_artifact_run_id", "unknown")
            source = st.session_state.get("gd_artifact_source", "local")
            scenario_count = st.session_state.get("gd_artifact_scenarios", "n/a")
            st.caption(
                f"Artifact-backed presentation mode: run `{run_id}` from `{source}` "
                f"(scenarios logged: {scenario_count})."
            )

        filter_summary = st.session_state.get("gd_filter_summary")
        target_baseline = st.session_state.get("gd_target_baseline")
        diversity_summary = st.session_state.get("gd_diversity_summary")
        if isinstance(filter_summary, dict) and isinstance(target_baseline, dict):
            st.caption(
                "Quality filter: "
                f"kept {filter_summary.get('kept', len(candidates))}/"
                f"{filter_summary.get('total', len(candidates))}, "
                f"removed {filter_summary.get('removed', 0)} baseline-like designs. "
                f"Target baseline trust score={float(target_baseline.get('score', 0.0)):.3f}, "
                f"min gain=+{float(filter_summary.get('min_score_gain', 0.0)):.2f}."
            )
        if isinstance(diversity_summary, dict):
            st.caption(
                "Diversification: "
                f"pool {diversity_summary.get('pool_size', len(candidates))} → "
                f"selected {diversity_summary.get('selected', len(candidates))}; "
                f"novelty weight={float(diversity_summary.get('novelty_weight', 0.0)):.2f}, "
                f"overlap penalty={float(diversity_summary.get('overlap_penalty', 0.0)):.2f}."
            )

        # Pareto frontier visualization
        st.markdown("###  Pareto Frontier (Trade-offs)")

        pareto_df = pd.DataFrame(candidates)

        fig_pareto = go.Figure()

        # Color by profile
        colors = {'Balanced': '#2563EB', 'Stability-First': '#0EA5E9', 'Binding-Optimized': '#10B981',
                  'Function-Maximized': '#14B8A6', 'Conservative': '#1D4ED8', 'Experimental': '#0891B2'}
        min_identity_gate = constraints.get('min_identity', 90)
        min_function_gate = constraints.get('min_function', -0.2)

        bind_min = float(pareto_df['binding'].min()) if len(pareto_df) else 0.0
        bind_span = float(pareto_df['binding'].max() - bind_min) if len(pareto_df) else 1.0
        bind_span = bind_span if bind_span > 1e-8 else 1.0
        pareto_df['binding_size'] = 12 + 18 * ((pareto_df['binding'] - bind_min) / bind_span)

        for profile in pareto_df['profile'].unique():
            df_profile = pareto_df[pareto_df['profile'] == profile]
            customdata = np.stack(
                [
                    df_profile['binding'].to_numpy(),
                    df_profile['stability'].to_numpy(),
                    df_profile['n_mutations'].to_numpy()
                ],
                axis=-1
            )
            fig_pareto.add_trace(go.Scatter(
                x=df_profile['identity'],
                y=df_profile['score'],
                mode='markers+text',
                name=profile,
                text=df_profile['candidate_id'].astype(str),
                textposition='top center',
                customdata=customdata,
                marker=dict(
                    size=df_profile['binding_size'],
                    color=colors.get(profile, '#FFFFFF'),
                    line=dict(width=1.4, color='rgba(255,255,255,0.9)'),
                    opacity=0.9,
                    symbol='circle'
                ),
                hovertemplate=f"<b>{profile}</b><br>" +
                              "Identity: %{x:.1f}%<br>" +
                              "Trust score: %{y:.3f}<br>" +
                              "DNA Binding: %{customdata[0]:.2f}<br>" +
                              "Stability: %{customdata[1]:.3f}<br>" +
                              "Mutations: %{customdata[2]:.0f}<extra></extra>"
            ))

        y_top = float(pareto_df['score'].max()) + 0.12 if len(pareto_df) else 0.5
        x_right = float(pareto_df['identity'].max()) + 0.7 if len(pareto_df) else 100
        fig_pareto.add_shape(
            type="rect",
            x0=min_identity_gate,
            x1=x_right,
            y0=min_function_gate,
            y1=y_top,
            fillcolor="rgba(16, 185, 129, 0.10)",
            line=dict(width=0),
            layer="below"
        )
        fig_pareto.add_vline(
            x=min_identity_gate,
            line_dash="dash",
            line_width=1.7,
            line_color="#EF4444",
            annotation_text=f"Identity floor {min_identity_gate:.0f}%"
        )
        fig_pareto.add_hline(
            y=min_function_gate,
            line_dash="dash",
            line_width=1.7,
            line_color="#F59E0B",
            annotation_text=f"Function floor {min_function_gate:.2f}"
        )

        best_idx = pareto_df['score'].idxmax()
        best_row = pareto_df.loc[best_idx]
        fig_pareto.add_trace(go.Scatter(
            x=[best_row['identity']],
            y=[best_row['score']],
            mode='markers+text',
            text=["Best"],
            textposition='top right',
            name='Best Candidate',
            marker=dict(size=24, color="#F59E0B", symbol='star', line=dict(width=1.3, color="#FFFFFF")),
            hovertemplate="<b>Best Candidate</b><br>Identity: %{x:.1f}%<br>Trust score: %{y:.3f}<extra></extra>"
        ))

        fig_pareto.update_layout(
            title="Design Space Exploration (Size = DNA Binding, Green = Feasible Region)",
            xaxis_title="Sequence Identity (%)",
            yaxis_title="Trust-Adjusted Rescue Score",
            height=520,
            margin=dict(l=24, r=20, t=70, b=96),
            legend=dict(orientation="h", yanchor="top", y=-0.2, x=0),
            xaxis=dict(range=[max(79, pareto_df['identity'].min() - 1), min(100, x_right)]),
            yaxis=dict(range=[min(min_function_gate - 0.15, pareto_df['score'].min() - 0.08), y_top])
        )
        render_plotly(fig_pareto, chart_kind="pareto", width="stretch")

        # === 3D STRUCTURE GALLERY ===
        st.markdown("### ️ 3D Structure Gallery")
        st.markdown("*Compare candidate structures side-by-side. Each shows unique mutations highlighted in profile color.*")

        # Read PDB once
        try:
            with open("data/raw/p53_wt.pdb", "r") as f:
                pdb_content = f.read()
            pdb_available = True
        except:
            pdb_available = False
            st.warning("PDB file not found. 3D gallery unavailable.")

        if pdb_available:
            # Display 3 candidates per row
            sorted_candidates = sorted(candidates, key=candidate_rank_value, reverse=True)
            n_cols = 3

            for row_start in range(0, min(len(sorted_candidates), 6), n_cols):
                gallery_cols = st.columns(n_cols)

                for i, col in enumerate(gallery_cols):
                    cand_idx = row_start + i
                    if cand_idx < len(sorted_candidates):
                        cand = sorted_candidates[cand_idx]
                        profile_color = cand.get('color', '#2563EB')

                        # Get mutation positions
                        mut_positions = cand.get('mut_positions', [])
                        if not mut_positions and cand.get('mutations'):
                            mut_positions = [int(''.join(filter(str.isdigit, m))) for m in cand['mutations'] if any(c.isdigit() for c in m)]

                        mut_str = "+".join([str(p) for p in mut_positions[:10]]) if mut_positions else "none"

                        # Find unique mutations (mutations this candidate has that others don't)
                        all_other_muts = set()
                        for other in sorted_candidates:
                            if other['candidate_id'] != cand['candidate_id']:
                                other_positions = other.get('mut_positions', [])
                                if not other_positions and other.get('mutations'):
                                    other_positions = [int(''.join(filter(str.isdigit, m))) for m in other['mutations'] if any(c.isdigit() for c in m)]
                                all_other_muts.update(other_positions)

                        unique_muts = [p for p in mut_positions if p not in all_other_muts]
                        unique_str = "+".join([str(p) for p in unique_muts[:5]]) if unique_muts else "none"

                        with col:
                            # Validity badge
                            badge_color = "#10B981" if cand['meets_constraints'] else "#EF4444"
                            badge_text = "" if cand['meets_constraints'] else "️"

                            st.markdown(f"""
                            <div style="text-align:center; padding:5px; background:linear-gradient(135deg, #EEF5FF 0%, #E0EEFF 100%);
                                        border-radius:10px; border:2px solid {profile_color}; margin-bottom:5px;">
                                <span style="color:{profile_color}; font-weight:bold; font-size:1.1em;">
                                    #{cand['candidate_id']} {cand['profile']} {badge_text}
                                </span>
                            </div>
                            """, unsafe_allow_html=True)

                            # 3D viewer for this candidate
                            viewer_id = f"gallery_view_{cand['candidate_id']}"
                            gallery_3d_html = f"""
                            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                            <div id="{viewer_id}" style="width:100%; height:280px; border-radius:8px; overflow:hidden; background:#FAFAFA;"></div>
                            <script>
                                (function() {{
                                    let viewer = $3Dmol.createViewer('{viewer_id}', {{backgroundColor: '0xFAFAFA'}});
                                    let pdb = `{pdb_content.replace(chr(10), chr(92) + 'n').replace('`', '')}`;
                                    viewer.addModel(pdb, 'pdb');

                                    // Base structure - dark slate for good contrast
                                    viewer.setStyle({{}}, {{cartoon: {{color: '0x64748B', opacity: 0.8}}}});

                                    // DNA binding loops - vibrant colors
                                    viewer.addStyle({{resi: [112,113,114,115,116,117,118,119,120,121,122,123,124]}},
                                                   {{cartoon: {{color: '0x0EA5E9', opacity: 0.95}}}});
                                    viewer.addStyle({{resi: [236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251]}},
                                                   {{cartoon: {{color: '0x10B981', opacity: 0.95}}}});

                                    // ALL mutations for this candidate (profile color)
                                    let allMuts = "{mut_str}".split("+");
                                    allMuts.forEach(pos => {{
                                        if (pos && pos !== 'none') {{
                                            viewer.addStyle({{resi: parseInt(pos)}}, {{
                                                stick: {{color: '{profile_color}', radius: 0.4}},
                                                cartoon: {{color: '{profile_color}'}}
                                            }});
                                        }}
                                    }});

                                    // UNIQUE mutations (gold glow effect - what sets this candidate apart)
                                    let uniqueMuts = "{unique_str}".split("+");
                                    uniqueMuts.forEach(pos => {{
                                        if (pos && pos !== 'none') {{
                                            viewer.addStyle({{resi: parseInt(pos)}}, {{
                                                sphere: {{color: '0xFBBF24', radius: 1.0, opacity: 0.5}}
                                            }});
                                        }}
                                    }});

                                    viewer.zoomTo();
                                    viewer.spin('y', 0.5);
                                    viewer.render();

                                    // Click to pause rotation
                                    document.getElementById('{viewer_id}').addEventListener('click', function() {{
                                        viewer.spin(false);
                                    }});
                                }})();
                            </script>
                            """
                            import streamlit.components.v1 as components
                            components.html(gallery_3d_html, height=290)

                            # Metrics under viewer
                            m1, m2, m3 = st.columns(3)
                            m1.metric("Score", f"{cand['score']:.2f}", label_visibility="collapsed")
                            m2.metric("ID%", f"{cand['identity']:.0f}%", label_visibility="collapsed")
                            m3.metric("Unc", f"{float(cand.get('uncertainty', 0.0)):.2f}", label_visibility="collapsed")

                            # Show unique mutations
                            if unique_muts:
                                st.caption(f" **Unique:** {', '.join(cand['mutations'][:3] if len(unique_muts) > 0 else ['None'])}")
                            else:
                                st.caption("*No unique mutations*")

            # === MUTATION COMPARISON HEATMAP ===
            st.markdown("###  Mutation Comparison Heatmap")
            st.markdown("*Rows are candidates, columns are informative residue positions; darker cells indicate mutation presence.*")

            matrix_df, freq_df = build_candidate_position_heatmap(sorted_candidates, max_positions=60)
            if not matrix_df.empty:
                pos_cols = [c for c in matrix_df.columns if c.startswith("pos_")]
                positions = [int(col.split("_", 1)[1]) for col in pos_cols]
                z_matrix = matrix_df[pos_cols].to_numpy(dtype=float)
                heat_hover = np.where(z_matrix > 0.5, "Yes", "No")

                fig_heatmap = go.Figure(
                    data=go.Heatmap(
                        z=z_matrix,
                        x=[f"Pos {p}" for p in positions],
                        y=matrix_df["candidate_label"].tolist(),
                        customdata=heat_hover,
                        colorscale=[[0, "#EDF4FF"], [1, "#1D4ED8"]],
                        zmin=0,
                        zmax=1,
                        showscale=True,
                        colorbar=dict(
                            title="Mutated",
                            tickvals=[0, 1],
                            ticktext=["No", "Yes"],
                            len=0.72,
                        ),
                        hovertemplate="Candidate: %{y}<br>Position: %{x}<br>Mutated: %{customdata}<extra></extra>",
                    )
                )
                fig_heatmap.update_layout(
                    height=max(420, 58 + 36 * len(matrix_df)),
                    margin=dict(l=16, r=16, t=64, b=96),
                    title="Mutation Presence Matrix",
                    xaxis_title="Residue Position",
                    yaxis_title="Candidate",
                )
                fig_heatmap.update_xaxes(tickangle=-45)
                render_plotly(fig_heatmap, chart_kind="heatmap", width="stretch")
                st.caption(
                    "Interpretation: vertical bands indicate convergent hotspot positions repeatedly selected across candidates."
                )

                # Frequency summary highlights consensus positions and mutation prevalence.
                freq_sorted = freq_df.sort_values(["count", "position"], ascending=[False, True]).copy()
                top_freq = freq_sorted.head(25)
                fig_hotspots = go.Figure(
                    go.Bar(
                        x=[f"Pos {int(p)}" for p in top_freq["position"]],
                        y=top_freq["count"],
                        marker=dict(
                            color=top_freq["frequency"],
                            colorscale="Blues",
                            line=dict(width=0.8, color="rgba(255,255,255,0.8)"),
                            colorbar=dict(title="Frequency", tickformat=".0%"),
                        ),
                        text=[f"{v:.0%}" for v in top_freq["frequency"]],
                        textposition="outside",
                        hovertemplate=(
                            "Position: %{x}<br>"
                            "Candidates mutating: %{y}<br>"
                            "Frequency: %{marker.color:.1%}<extra></extra>"
                        ),
                    )
                )
                fig_hotspots.update_layout(
                    height=380,
                    margin=dict(l=24, r=20, t=64, b=96),
                    title="Consensus Mutation Frequency (Top informative positions)",
                    xaxis_title="Residue Position",
                    yaxis_title="Mutating Candidates",
                )
                fig_hotspots.update_xaxes(tickangle=-45)
                render_plotly(fig_hotspots, chart_kind="bar", width="stretch")
                st.caption(
                    "Interpretation: higher bars are stronger consensus loci; prioritize these for mechanism and wet-lab follow-up."
                )

                consensus = top_freq[top_freq["count"] >= 2]["position"].astype(int).tolist()
                if consensus:
                    st.success(
                        "**Consensus positions (mutated by 2+ candidates):** "
                        + ", ".join(str(p) for p in consensus)
                    )

        # Candidate cards
        st.markdown("### Candidate Details")

        # Sort by score (best first)
        sorted_candidates = sorted(candidates, key=candidate_rank_value, reverse=True)

        for i in range(0, len(sorted_candidates), 2):
            card_cols = st.columns(2)

            for j, col in enumerate(card_cols):
                if i + j < len(sorted_candidates):
                    cand = sorted_candidates[i + j]

                    # Color based on validity
                    border_color = "#10B981" if cand['meets_constraints'] else "#EF4444"
                    validity_badge = " VALID" if cand['meets_constraints'] else "️ CONSTRAINT VIOLATION"

                    with col:
                        st.markdown(f"""
                        <div class="candidate-card" style="border-left-color: {border_color};">
                            <h4>Candidate #{cand['candidate_id']} - {cand['profile']}</h4>
                            <p><b>{validity_badge}</b></p>
                        </div>
                        """, unsafe_allow_html=True)

                        # Metrics row
                        m1, m2, m3, m4 = st.columns(4)
                        m1.metric("Score", f"{cand['score']:.3f}")
                        m2.metric("Gain vs target", f"{float(cand.get('score_gain_vs_target', 0.0)):+.3f}")
                        m3.metric("Identity", f"{cand['identity']:.1f}%")
                        m4.metric("Mutations", cand['n_mutations'])

                        # Show mutations
                        muts_display = ", ".join(cand['mutations'][:5])
                        if len(cand['mutations']) > 5:
                            muts_display += f" (+{len(cand['mutations']) - 5} more)"
                        if not muts_display:
                            muts_display = "None (no sequence change from baseline)"
                        st.caption(f"**Mutations:** {muts_display}")

                        # Physics scores (older autoregressive campaigns stored 0.0 before backfill)
                        _ar_zero = cand.get('profile') == 'Autoregressive' and cand['stability'] == 0.0 and cand['binding'] == 0.0
                        stab_str = "N/A" if _ar_zero else f"{cand['stability']:.3f}"
                        bind_str = "N/A" if _ar_zero else f"{cand['binding']:.2f}"
                        st.caption(f"Stability: {stab_str} | Binding: {bind_str}")
                        st.caption(
                            f"Raw score: {float(cand.get('score_raw', cand['score'])):.3f} | "
                            f"Calibrated: {float(cand.get('score_calibrated', cand['score'])):.3f} | "
                            f"Uncertainty: {float(cand.get('uncertainty', 0.0)):.3f} | "
                            f"OOD dist: {float(cand.get('ood_distance', 0.0)):.3f}"
                        )
                        if "selection_score" in cand:
                            def _fmt_optional_metric(raw_value, precision):
                                try:
                                    value = float(raw_value)
                                except Exception:
                                    return "n/a"
                                if not np.isfinite(value):
                                    return "n/a"
                                return f"{value:.{precision}f}"

                            st.caption(
                                f"Selection score: {_fmt_optional_metric(cand.get('selection_score'), 3)} | "
                                f"Novelty: {_fmt_optional_metric(cand.get('diversity_novelty'), 2)} | "
                                f"Max overlap: {_fmt_optional_metric(cand.get('max_mutation_overlap'), 2)}"
                            )
                        st.caption(f"Selection reason: {cand.get('selection_reason', 'n/a')}")

                        # Action buttons
                        btn_col1, btn_col2 = st.columns(2)
                        with btn_col1:
                            if st.button(f" Select #{cand['candidate_id']}", key=f"select_{cand['candidate_id']}"):
                                st.session_state['selected_candidate'] = cand
                                st.session_state['results'] = pd.DataFrame([{
                                    'Step': 0, 'Score': cand['score'], 'Stability': cand['stability'],
                                    'DNARecruitment': cand['binding'], 'Identity': cand['identity'],
                                    'MutationCount': cand['n_mutations'], 'MutSummary': ", ".join(cand['mutations'][:5]),
                                    'Sequence': cand['sequence'], 'SurfaceCharge': 0.5, 'HydroPacking': 0.4,
                                    'GrassmannNovelty': 0.1, 'LatentIdentity': cand['identity'],
                                    'LX': 0, 'LY': 0, 'LatentExcitation': 0, 'Phase': 'Final'
                                }])
                                st.session_state['target_mut_saved'] = constraints.get('target_mutation', 'R175H')
                                st.success(f"Selected Candidate #{cand['candidate_id']}! Go to Validation tab.")

                        with btn_col2:
                            fasta = f">p53_rescue_cand{cand['candidate_id']}\n{cand['sequence']}"
                            st.download_button(
                                f" FASTA",
                                fasta,
                                file_name=f"candidate_{cand['candidate_id']}.fasta",
                                key=f"dl_{cand['candidate_id']}"
                            )

        # === COMPARISON TABLE ===
        st.markdown("###  Comparison Matrix")

        compare_df = pd.DataFrame([{
            'ID': c['candidate_id'],
            'Profile': c['profile'],
            'Score': f"{c['score']:.3f}",
            'Raw Score': f"{float(c.get('score_raw', c['score'])):.3f}",
            'Gain vs Target': f"{float(c.get('score_gain_vs_target', 0.0)):+.3f}",
            'Uncertainty': f"{float(c.get('uncertainty', 0.0)):.3f}",
            'OOD Dist': f"{float(c.get('ood_distance', 0.0)):.3f}",
            'Novelty': f"{float(c.get('diversity_novelty', 0.0)):.2f}",
            'Stability': f"{c['stability']:.3f}",
            'Binding': f"{c['binding']:.2f}",
            'Identity': f"{c['identity']:.1f}%",
            'Mutations': c['n_mutations'],
            'Selection Reason': str(c.get('selection_reason', 'n/a')),
            'Valid': '' if c['meets_constraints'] else ''
        } for c in sorted_candidates])

        st.dataframe(
            compare_df,
            hide_index=True,
            width="stretch",
            column_config={
                "ID": st.column_config.NumberColumn("Candidate"),
                "Profile": st.column_config.TextColumn("Strategy"),
                "Valid": st.column_config.TextColumn("Meets Constraints")
            }
        )

        # Export all candidates
        st.download_button(
            " Export All Candidates (CSV)",
            pd.DataFrame(candidates).to_csv(index=False),
            file_name="generative_design_candidates.csv",
            width="stretch"
        )

        # === BEST PERFORMING PROTEIN RESULTS ===
        st.markdown("---")
        st.markdown("## 🏆 Best Performing Protein")
        st.markdown("*Top-scoring rescue design from this optimization run*")

        # Get the best candidate (already sorted by score)
        best_candidate = sorted_candidates[0]

        # Hero metrics row
        hero_col1, hero_col2, hero_col3, hero_col4 = st.columns(4)
        with hero_col1:
            st.metric("🎯 Trust Score", f"{best_candidate['score']:.3f}",
                     delta="Best" if best_candidate['meets_constraints'] else "Constraint Issue")
        with hero_col2:
            st.metric("🧬 Identity", f"{best_candidate['identity']:.1f}%",
                     delta=f"{int(best_candidate['n_mutations'])} mutations")
        with hero_col3:
            _hero_ar_zero = best_candidate.get('profile') == 'Autoregressive' and best_candidate['stability'] == 0.0 and best_candidate['binding'] == 0.0
            st.metric("⚡ Stability", "N/A" if _hero_ar_zero else f"{best_candidate['stability']:.3f}")
        with hero_col4:
            st.metric("📉 Uncertainty", f"{float(best_candidate.get('uncertainty', 0.0)):.3f}")
        _hero_bind_str = "N/A" if _hero_ar_zero else f"{best_candidate['binding']:.2f}"
        st.caption(
            f"Raw oracle score: {float(best_candidate.get('score_raw', best_candidate['score'])):.3f} | "
            f"Calibrated score: {float(best_candidate.get('score_calibrated', best_candidate['score'])):.3f} | "
            f"OOD distance: {float(best_candidate.get('ood_distance', 0.0)):.3f} | "
            f"DNA Binding: {_hero_bind_str}"
        )

        # Two column layout: Mutations + 3D Structure | Graphs
        result_col1, result_col2 = st.columns([1, 1])

        with result_col1:
            # Mutations applied
            st.markdown("### 🔬 Applied Mutations")
            st.markdown(f"**Strategy:** {best_candidate['profile']}")

            if best_candidate['mutations']:
                # Display mutations in a clean format
                mut_display = ""
                for i, mut in enumerate(best_candidate['mutations']):
                    if i < 10:  # Show first 10
                        mut_display += f"• **{mut}**\n"
                if len(best_candidate['mutations']) > 10:
                    mut_display += f"*...and {len(best_candidate['mutations']) - 10} more*"
                st.markdown(mut_display)
            else:
                st.info("No mutations detected")

            # Mutation summary (positions)
            if best_candidate.get('mut_positions'):
                positions = best_candidate['mut_positions']
                st.caption(f"**Positions modified:** {', '.join(str(p) for p in sorted(positions)[:15])}")

        with result_col2:
            # Optimization trajectory graph
            st.markdown("### 📈 Optimization Trajectory")

            if best_candidate.get('trajectory') and len(best_candidate['trajectory']) > 1:
                traj_df = pd.DataFrame(best_candidate['trajectory'])
                min_identity_target = constraints.get('min_identity', 90)
                min_function_target = constraints.get('min_function', -0.2)
                min_stability_target = constraints.get('min_stability', -0.2)

                fig_result = go.Figure()
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['score'],
                    mode='lines+markers', name='Rescue Score',
                    line=dict(color=best_candidate.get('color', '#2563EB'), width=3),
                    marker=dict(size=6),
                    fill='tozeroy',
                    fillcolor='rgba(37, 99, 235, 0.12)'
                ))
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['stability'],
                    mode='lines', name='Stability',
                    line=dict(color='#38BDF8', width=2, dash='dot')
                ))
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['binding'] / 10,  # Scale for visibility
                    mode='lines', name='DNA Binding (scaled)',
                    line=dict(color='#10B981', width=2, dash='dash')
                ))
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['identity'],
                    mode='lines', name='Identity (%)',
                    line=dict(color='#1D4ED8', width=2, dash='dashdot'),
                    yaxis='y2'
                ))

                fig_result.add_hline(
                    y=min_function_target,
                    line_dash='dot',
                    line_width=1.2,
                    line_color='#F59E0B',
                    annotation_text='Function floor',
                    annotation_position='bottom right'
                )
                fig_result.add_hline(
                    y=min_stability_target,
                    line_dash='dot',
                    line_width=1.2,
                    line_color='#0EA5E9',
                    annotation_text='Stability floor',
                    annotation_position='top right'
                )
                fig_result.add_shape(
                    type='line',
                    x0=float(traj_df['step'].min()),
                    x1=float(traj_df['step'].max()),
                    y0=min_identity_target,
                    y1=min_identity_target,
                    xref='x',
                    yref='y2',
                    line=dict(color='#DC2626', width=1.2, dash='dot')
                )
                fig_result.add_annotation(
                    x=float(traj_df['step'].max()),
                    y=min_identity_target,
                    xref='x',
                    yref='y2',
                    text='Identity floor',
                    showarrow=False,
                    xanchor='right',
                    yanchor='bottom',
                    font=dict(size=10, color='#DC2626'),
                    bgcolor='rgba(255,255,255,0.65)'
                )

                fig_result.update_layout(
                    height=420,
                    margin=dict(l=24, r=20, t=70, b=96),
                    title=f"Candidate #{best_candidate['candidate_id']} Evolution",
                    xaxis_title="Optimization Step",
                    yaxis_title="Score / Stability / Binding (scaled)",
                    yaxis2=dict(
                        title="Identity (%)",
                        overlaying='y',
                        side='right',
                        range=[75, 101],
                        showgrid=False
                    ),
                    legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="left", x=0.0)
                )
                render_plotly(fig_result, chart_kind="trajectory", width="stretch")
                st.caption(
                    "Rescue score should trend upward while stability stays above the floor and identity remains above the identity line."
                )
            else:
                st.info("No trajectory data available")

        # Protein sequence display
        st.markdown("### 🧬 Rescued Protein Sequence")

        seq = best_candidate['sequence']

        # Create a styled sequence view highlighting mutations
        if best_candidate.get('mut_positions'):
            mut_positions = set(best_candidate['mut_positions'])
            # Build colored sequence HTML
            seq_html = '<div style="font-family: monospace; font-size: 12px; line-height: 1.8; background: #F4F8FF; padding: 15px; border-radius: 8px; overflow-x: auto; border:1px solid #D5E5F7;">'
            for i, aa in enumerate(seq):
                pos = i + 1  # 1-indexed
                if pos in mut_positions:
                    seq_html += f'<span style="color: #F43F5E; font-weight: bold; background: rgba(244,63,94,0.18); padding: 2px 1px; border-radius: 3px;">{aa}</span>'
                else:
                    seq_html += f'<span style="color: #3A5C84;">{aa}</span>'
                # Add line break every 60 residues
                if (i + 1) % 60 == 0:
                    seq_html += '<br>'
            seq_html += '</div>'
            st.markdown(seq_html, unsafe_allow_html=True)
        else:
            # Plain sequence display
            st.code(seq, language=None)

        # Legend and download
        leg_col1, leg_col2 = st.columns([2, 1])
        with leg_col1:
            st.caption("🔴 **Red residues** = Mutations applied | ⬜ Gray = Wild-type positions")
            st.caption(f"**Total length:** {len(seq)} residues | **Modified:** {best_candidate['n_mutations']} positions")

        with leg_col2:
            # Download best protein
            fasta_best = f">p53_best_rescue_{best_candidate['profile'].replace(' ', '_')}\n{seq}"
            st.download_button(
                "⬇️ Download Best Protein (FASTA)",
                fasta_best,
                file_name=f"best_rescue_{best_candidate['candidate_id']}.fasta",
                width="stretch"
            )

        # === ADVANCED VISUALIZATIONS ===
        st.markdown("---")
        st.markdown("## 📊 Advanced Analysis")

        if best_candidate.get('trajectory') and len(best_candidate['trajectory']) > 2:
            traj_df = pd.DataFrame(best_candidate['trajectory'])

            viz_col1, viz_col2 = st.columns(2)

            with viz_col1:
                _has_latent = "lx" in traj_df.columns and "ly" in traj_df.columns

                if _has_latent:
                    # === 3D LATENT SPACE PLOT (gradient-based trajectories) ===
                    st.markdown("### 🌐 3D Latent Space Trajectory")
                    st.caption("*Surface = estimated loss landscape from sampled optimization states; line = actual optimization path.*")

                    fig_3d = go.Figure()
                    mesh = build_loss_mesh(traj_df, grid_size=40)
                    if mesh is not None:
                        grid_x, grid_y, grid_z = mesh
                        fig_3d.add_trace(
                            go.Surface(
                                x=grid_x,
                                y=grid_y,
                                z=grid_z,
                                colorscale="Blues",
                                opacity=0.55,
                                name="Estimated loss surface",
                                showscale=False,
                                hovertemplate="Latent X: %{x:.2f}<br>Latent Y: %{y:.2f}<br>Loss: %{z:.3f}<extra></extra>",
                            )
                        )
                    else:
                        st.info(
                            "Loss surface mesh unavailable for this trajectory (need enough points with `loss_total`). "
                            "Showing path-only view."
                        )

                    traj_z = pd.to_numeric(traj_df.get("loss_total", traj_df.get("score")), errors="coerce")
                    if traj_z.isna().all():
                        traj_z = pd.to_numeric(traj_df.get("lz", traj_df.get("score")), errors="coerce")

                    fig_3d.add_trace(
                        go.Scatter3d(
                            x=traj_df["lx"],
                            y=traj_df["ly"],
                            z=traj_z,
                            mode="lines+markers",
                            line=dict(color="#0EA5E9", width=6),
                            marker=dict(
                                size=5,
                                color=traj_df["step"],
                                colorscale=[[0.0, "#1D4ED8"], [0.5, "#0EA5E9"], [1.0, "#10B981"]],
                                colorbar=dict(title="Step", x=1.02, thickness=12),
                                opacity=0.95,
                            ),
                            text=[
                                f"Step {s}<br>Score: {sc:.3f}<br>Loss: {lz:.3f}<br>Identity: {ident:.1f}%"
                                for s, sc, lz, ident in zip(
                                    traj_df["step"],
                                    traj_df["score"],
                                    traj_z,
                                    traj_df["identity"],
                                )
                            ],
                            hoverinfo="text",
                            name="Optimization trajectory",
                        )
                    )
                    fig_3d.add_trace(
                        go.Scatter3d(
                            x=[traj_df["lx"].iloc[0]],
                            y=[traj_df["ly"].iloc[0]],
                            z=[traj_z.iloc[0]],
                            mode="markers",
                            marker=dict(size=13, color="#DC2626", symbol="diamond", line=dict(color="white", width=1)),
                            name="Start (cancer)",
                        )
                    )
                    fig_3d.add_trace(
                        go.Scatter3d(
                            x=[traj_df["lx"].iloc[-1]],
                            y=[traj_df["ly"].iloc[-1]],
                            z=[traj_z.iloc[-1]],
                            mode="markers",
                            marker=dict(size=14, color="#059669", symbol="diamond", line=dict(color="white", width=1)),
                            name="End (selected)",
                        )
                    )
                    fig_3d.update_layout(
                        height=560,
                        margin=dict(l=24, r=20, t=64, b=100),
                        scene=dict(
                            xaxis_title="Latent X",
                            yaxis_title="Latent Y",
                            zaxis_title="Loss (lower is better)",
                            camera=dict(eye=dict(x=1.4, y=1.35, z=1.1)),
                        ),
                        legend=dict(orientation="h", yanchor="top", y=-0.24, xanchor="left", x=0.0),
                    )
                    render_plotly(fig_3d, chart_kind="3d", width="stretch")
                    st.caption("Interpretation: good runs move toward lower-surface basins while preserving identity and constraints.")
                else:
                    # === 2D MUTATION TRAJECTORY (autoregressive trajectories) ===
                    st.markdown("### 📈 Mutation Trajectory")
                    st.caption("*Score improvement as rescue mutations are added one-by-one.*")

                    fig_traj = go.Figure()
                    hover_texts = []
                    for _, row in traj_df.iterrows():
                        mut_label = row.get("mutation_applied", "")
                        step_label = f"Step {int(row['step'])}"
                        if mut_label:
                            step_label += f"<br>Applied: {mut_label}"
                        step_label += f"<br>Score: {row['score']:.3f}"
                        if "n_mutations" in traj_df.columns:
                            step_label += f"<br>Total mutations: {int(row['n_mutations'])}"
                        hover_texts.append(step_label)

                    fig_traj.add_trace(go.Scatter(
                        x=traj_df["step"],
                        y=traj_df["score"],
                        mode="lines+markers",
                        line=dict(color="#0EA5E9", width=3),
                        marker=dict(
                            size=9,
                            color=traj_df["score"],
                            colorscale=[[0.0, "#DC2626"], [0.5, "#F59E0B"], [1.0, "#10B981"]],
                            colorbar=dict(title="Score", thickness=12),
                        ),
                        text=hover_texts,
                        hoverinfo="text",
                        name="Score trajectory",
                    ))
                    # Mark start and end
                    fig_traj.add_trace(go.Scatter(
                        x=[traj_df["step"].iloc[0]],
                        y=[traj_df["score"].iloc[0]],
                        mode="markers",
                        marker=dict(size=14, color="#DC2626", symbol="diamond"),
                        name="Start (cancer baseline)",
                    ))
                    fig_traj.add_trace(go.Scatter(
                        x=[traj_df["step"].iloc[-1]],
                        y=[traj_df["score"].iloc[-1]],
                        mode="markers",
                        marker=dict(size=14, color="#059669", symbol="diamond"),
                        name="End (best rescue)",
                    ))
                    fig_traj.update_layout(
                        height=560,
                        margin=dict(l=24, r=20, t=64, b=100),
                        xaxis_title="Mutation Step",
                        yaxis_title="Oracle Score",
                        legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="left", x=0.0),
                    )
                    render_plotly(fig_traj, chart_kind="trajectory", width="stretch")
                    st.caption("Each step adds one greedy rescue mutation. Higher scores indicate improved p53 function.")

            with viz_col2:
                # === MULTI-METRIC COMPARISON ===
                st.markdown("### 📈 Multi-Metric Evolution")
                normalize_metrics = st.toggle(
                    "Normalize metrics to thresholds",
                    value=True,
                    key=f"metric_norm_toggle_{best_candidate['candidate_id']}",
                    help="Normalized mode maps each metric relative to its acceptance threshold.",
                )
                metric_df = build_multi_metric_frame(
                    traj_df,
                    normalize=normalize_metrics,
                    min_function=float(constraints.get("min_function", -0.2)),
                    min_stability=float(constraints.get("min_stability", -0.2)),
                    min_binding=float(constraints.get("min_binding", 5.0)),
                    min_identity=float(constraints.get("min_identity", 90.0)),
                )

                fig_multi = go.Figure()
                fig_multi.add_trace(
                    go.Scatter(
                        x=metric_df["step"],
                        y=metric_df["score"],
                        mode="lines+markers",
                        name="Rescue Score",
                        line=dict(color="#0EA5E9", width=3),
                        marker=dict(size=4),
                    )
                )
                fig_multi.add_trace(
                    go.Scatter(
                        x=metric_df["step"],
                        y=metric_df["stability"],
                        mode="lines",
                        name="Stability",
                        line=dict(color="#38BDF8", width=2),
                    )
                )
                fig_multi.add_trace(
                    go.Scatter(
                        x=metric_df["step"],
                        y=metric_df["binding"],
                        mode="lines",
                        name="DNA Binding",
                        line=dict(color="#10B981", width=2),
                    )
                )
                fig_multi.add_trace(
                    go.Scatter(
                        x=metric_df["step"],
                        y=metric_df["identity"],
                        mode="lines",
                        name="Identity",
                        line=dict(color="#1D4ED8", width=2, dash="dot"),
                        yaxis="y" if normalize_metrics else "y2",
                    )
                )

                if normalize_metrics:
                    fig_multi.add_hline(y=0.0, line_dash="dot", line_color="rgba(239,68,68,0.45)", line_width=1.1)
                    fig_multi.add_hline(y=0.8, line_dash="dot", line_color="rgba(16,185,129,0.55)", line_width=1.1)
                    fig_multi.update_layout(
                        yaxis_title="Threshold-anchored normalized value [0,1]",
                        yaxis=dict(range=[-0.05, 1.05]),
                    )
                else:
                    fig_multi.add_hline(
                        y=float(constraints.get("min_function", -0.2)),
                        line_dash="dot",
                        line_color="#F59E0B",
                        line_width=1.1,
                        annotation_text="Function floor",
                    )
                    fig_multi.add_hline(
                        y=float(constraints.get("min_stability", -0.2)),
                        line_dash="dot",
                        line_color="#0EA5E9",
                        line_width=1.1,
                        annotation_text="Stability floor",
                    )
                    fig_multi.add_hline(
                        y=float(constraints.get("min_binding", 5.0)),
                        line_dash="dot",
                        line_color="#10B981",
                        line_width=1.1,
                        annotation_text="Binding floor",
                    )
                    fig_multi.update_layout(
                        yaxis_title="Raw score / stability / binding",
                        yaxis2=dict(
                            title="Identity (%)",
                            overlaying="y",
                            side="right",
                            range=[75, 101],
                            showgrid=False,
                        ),
                    )

                fig_multi.update_layout(
                    height=500,
                    margin=dict(l=24, r=20, t=64, b=100),
                    xaxis_title="Optimization Step",
                    legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="left", x=0.0),
                )
                render_plotly(fig_multi, chart_kind="multiline", width="stretch")
                st.caption(
                    "Interpretation: upward score with stable/above-threshold constraints indicates improving rescue quality; "
                    "oscillation early is expected from exploration before convergence."
                )

            # === RESIDUE CONTRIBUTION HEATMAP ===
            st.markdown("### 🧬 Position-wise Mutation Impact")
            st.caption("*Tracks mutation persistence across optimization checkpoints (not generic chunk colors).*")

            def _region_label(pos: int) -> str:
                if pos in set(range(112, 125)) | set(range(164, 175)) | set(range(236, 252)):
                    return "DNA loop"
                if pos in {176, 179, 238, 242}:
                    return "Zinc site"
                if pos in set(range(325, 355)):
                    return "Dimer interface"
                return "Other"

            pos_stats: Dict[int, Dict[str, Any]] = {}
            for _, row in traj_df.iterrows():
                step_val = int(row.get("step", 0) or 0)
                mut_pos = row.get("mut_positions", [])
                if isinstance(mut_pos, str):
                    try:
                        mut_pos = json.loads(mut_pos)
                    except Exception:
                        mut_pos = []
                for p in mut_pos or []:
                    try:
                        pos = int(p)
                    except Exception:
                        continue
                    entry = pos_stats.setdefault(pos, {"count": 0, "first_step": step_val})
                    entry["count"] += 1
                    entry["first_step"] = min(int(entry["first_step"]), step_val)

            # Fallback for sparse trajectories: use final selected mutations.
            if not pos_stats:
                for p in best_candidate.get("mut_positions", []):
                    pos = int(p)
                    pos_stats[pos] = {"count": 1, "first_step": int(traj_df["step"].max() if not traj_df.empty else 0)}

            if pos_stats:
                total_steps = max(int(traj_df["step"].max()), 1) if "step" in traj_df.columns and not traj_df.empty else 1
                impact_df = pd.DataFrame(
                    [
                        {
                            "position": int(pos),
                            "count": int(meta["count"]),
                            "frequency": float(meta["count"] / total_steps),
                            "first_step": int(meta["first_step"]),
                            "region": _region_label(int(pos)),
                        }
                        for pos, meta in sorted(pos_stats.items(), key=lambda kv: kv[0])
                    ]
                )
                impact_df = impact_df.sort_values(["frequency", "position"], ascending=[False, True]).reset_index(drop=True)
                top_impact = impact_df.head(40).sort_values("position")

                fig_heatmap_pos = go.Figure(
                    data=go.Heatmap(
                        z=[top_impact["frequency"].to_numpy(dtype=float)],
                        x=[f"Pos {int(p)}" for p in top_impact["position"]],
                        y=["Mutation persistence"],
                        customdata=np.array(top_impact[["region", "first_step"]]),
                        colorscale="YlOrRd",
                        zmin=0.0,
                        zmax=max(1e-6, float(top_impact["frequency"].max())),
                        colorbar=dict(title="Step frequency", tickformat=".0%"),
                        hovertemplate=(
                            "Position: %{x}<br>"
                            "Persistence: %{z:.1%}<br>"
                            "Region: %{customdata[0]}<br>"
                            "First seen step: %{customdata[1]}<extra></extra>"
                        ),
                    )
                )
                fig_heatmap_pos.update_layout(
                    height=300,
                    margin=dict(l=16, r=16, t=64, b=90),
                    xaxis_title="Residue position",
                    yaxis_title="",
                )
                render_plotly(fig_heatmap_pos, chart_kind="heatmap", width="stretch")

                fig_first_seen = go.Figure(
                    go.Bar(
                        x=[f"Pos {int(p)}" for p in top_impact["position"]],
                        y=top_impact["first_step"],
                        marker=dict(
                            color=top_impact["frequency"],
                            colorscale="Blues",
                            line=dict(width=0.8, color="rgba(255,255,255,0.8)"),
                            colorbar=dict(title="Persistence", tickformat=".0%"),
                        ),
                        hovertemplate=(
                            "Position: %{x}<br>"
                            "First seen step: %{y}<br>"
                            "Persistence: %{marker.color:.1%}<extra></extra>"
                        ),
                    )
                )
                fig_first_seen.update_layout(
                    height=330,
                    margin=dict(l=24, r=20, t=64, b=96),
                    title="Mutation appearance timing (lower is earlier)",
                    xaxis_title="Residue position",
                    yaxis_title="First optimization step observed",
                )
                fig_first_seen.update_xaxes(tickangle=-45)
                render_plotly(fig_first_seen, chart_kind="bar", width="stretch")

                st.caption(
                    "Interpretation: persistent/high-frequency positions are stronger optimization anchors; "
                    "early-appearing positions tend to be core drivers."
                )
                st.dataframe(
                    top_impact.rename(
                        columns={
                            "position": "Position",
                            "count": "Observed checkpoints",
                            "frequency": "Step frequency",
                            "first_step": "First step",
                            "region": "Functional region",
                        }
                    ),
                    hide_index=True,
                    width="stretch",
                )

        else:
            st.info("Run optimization to see advanced visualizations")

    else:
        # Show concept explanation when no candidates
        st.markdown("---")
        st.markdown("##  How Generative Design Works")

        concept_col1, concept_col2 = st.columns(2)

        with concept_col1:
            st.markdown("""
            ### Mechanical CAD (Fusion 360)
            ```
            INPUT:
            - Load forces (where stress applies)
            - Support points (fixed positions)
            - Material (steel, aluminum, etc.)
            - Manufacturing (3D print, CNC, etc.)

            AI GENERATES:
            → Topology-optimized structures
            → Multiple designs on Pareto frontier
            → Non-intuitive, organic shapes
            ```
            """)

        with concept_col2:
            st.markdown("""
            ### p53-proteoMgCAD
            ```
            INPUT:
            - Physics (stability, binding thresholds)
            - Geometry (locked residues, protected regions)
            - Material (identity level = "human-like")
            - Manufacturing (delivery method)

            AI GENERATES:
            → Second-site suppressor mutations
            → Multiple candidates exploring trade-offs
            → Novel, non-obvious rescue designs
            ```
            """)

        st.info("""
        **The Key Insight**: Instead of asking "what mutations should I try?", you specify
        "what must my protein achieve?" and the AI explores the solution space to find
        optimal designs that meet your constraints.

        This is the same paradigm shift that transformed mechanical engineering with
        topology optimization. Now applied to protein engineering.
        """)

# === TAB 2: ANALYSIS & DRUGS ===
with tab2:
    st.markdown("""
    <div class="lovable-card" style="margin-bottom: 1.25rem;">
        <h2 style="margin: 0 0 0.45rem 0; font-size: 1.72rem;">📊 Analysis Dashboard</h2>
        <p style="color: #496487; margin: 0;">Single-page analysis pipeline: validation, mechanism, drug strategy, and clinical impact.</p>
    </div>
    """, unsafe_allow_html=True)

    active_target = st.session_state.get("target_mut_saved", target_mut)
    active_candidate = get_active_analysis_candidate(active_target)
    sync_analysis_inputs_from_candidate(active_candidate)
    runtime_caps = get_runtime_capabilities_snapshot()

    with st.expander("Runtime Capabilities", expanded=False):
        cap1, cap2, cap3, cap4 = st.columns(4)
        cap1.metric("Torch", runtime_caps.get("torch_version", "n/a"))
        cap2.metric("MPS", "yes" if runtime_caps.get("mps_available") else "no")
        cap3.metric("Vina CLI", "yes" if runtime_caps.get("vina_cli_available") else "no")
        md_ready = bool(runtime_caps.get("openmm_installed") and runtime_caps.get("openff_installed"))
        cap4.metric("MD stack", "ready" if md_ready else "missing")

        explain_status = get_explainability_backend_status()
        st.caption(
            "Explainability attention backend: "
            + ("ready" if explain_status.get("ready") else "unavailable")
            + f" ({explain_status.get('reason', 'n/a')})"
        )
        if DRUG_MODULE_AVAILABLE:
            try:
                manifest = DrugGeneratorEngine().get_receptor_manifest(allow_wt_receptor_fallback=False)
                dup_map = manifest.get("duplicate_receptor_paths", {})
                if dup_map:
                    st.warning(
                        "Pocket receptor manifest warning: one receptor file is reused by multiple pockets "
                        f"({dup_map})."
                    )
                else:
                    st.caption("Pocket receptor manifest: unique receptor mapping per pocket.")
            except Exception as err:
                st.caption(f"Drug runtime manifest probe unavailable: {err}")

    if active_candidate:
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Target Mutation", active_candidate["target_mutation"])
        c2.metric("Rescue Mutations", active_candidate["n_mutations"])
        c3.metric("Rescue Score", f"{active_candidate['score']:.3f}")
        c4.metric("Identity", f"{active_candidate['identity']:.1f}%")
        c5.metric("Stability", f"{active_candidate['stability']:.3f}")
    else:
        st.info("Run and select a candidate in the Design Studio to unlock real validation and downstream analysis.")

    # === 1) VALIDATION ===
    st.markdown("---")
    physics_fallback_available = callable(physics_based_score_fn)
    if physics_fallback_available:
        st.markdown("### 1) 🔬 Real Validation (DMS + Physics Fallback)")
        st.caption("Uses Giacomelli 2018 DMS when available, physics-based estimates (BLOSUM62 + AA properties) as fallback. Evidence triage, not wet-lab proof.")
    else:
        st.markdown("### 1) 🔬 Real Validation (DMS only)")
        st.caption("Uses Giacomelli 2018 DMS where available. Evidence triage, not wet-lab proof.")
        st.warning("Running in DMS-only mode; physics fallback unavailable.")

    if active_candidate:
        dms_by_mut, dms_by_pos = get_dms_lookup_tables()
        target_for_eval = str(active_candidate["target_mutation"]).upper()
        rescue_mutations = [m.upper() for m in active_candidate["mutations"]]
        mut_summary_raw = active_candidate.get("mut_summary_raw", "")

        if active_candidate.get("source") == "results_table" and "..." in mut_summary_raw:
            st.warning("Mutation list from results is truncated. Use `Select` on a candidate card for full mutation-level validation.")

        target_row = dms_by_mut[dms_by_mut["mutation"] == target_for_eval] if not dms_by_mut.empty else pd.DataFrame()
        target_dms_score = float(target_row["score"].iloc[0]) if not target_row.empty else np.nan

        dms_rows = []
        for mut in rescue_mutations:
            parsed = parse_single_mutation(mut)
            if parsed is None:
                continue
            _, pos, _ = parsed

            mut_row = dms_by_mut[dms_by_mut["mutation"] == mut] if not dms_by_mut.empty else pd.DataFrame()
            has_dms = not mut_row.empty

            physics_used = False
            if has_dms:
                mut_score = float(mut_row["score"].iloc[0])
                score_method = "DMS"
                confidence = 1.0
            else:
                physics_result = safe_physics_based_score(mut)
                if physics_result is not None and pd.notna(physics_result.get("score")):
                    mut_score = float(physics_result["score"])
                    score_method = "Physics"
                    confidence = float(physics_result.get("confidence", 0.0))
                    physics_used = True
                else:
                    mut_score = np.nan
                    score_method = "Unavailable"
                    confidence = 0.0

            pos_row = dms_by_pos[dms_by_pos["pos"] == pos] if not dms_by_pos.empty else pd.DataFrame()
            pos_mean = float(pos_row["pos_mean_score"].iloc[0]) if not pos_row.empty else np.nan
            pos_count = int(pos_row["pos_count"].iloc[0]) if not pos_row.empty else 0

            # Calculate delta vs target (use physics for both if neither has DMS)
            if has_dms and pd.notna(target_dms_score):
                delta_vs_target = mut_score - target_dms_score
            elif (not has_dms) and pd.isna(target_dms_score) and physics_used:
                target_physics = safe_physics_based_score(target_for_eval)
                if target_physics is not None and pd.notna(target_physics.get("score")):
                    delta_vs_target = mut_score - float(target_physics["score"])
                else:
                    delta_vs_target = np.nan
            else:
                delta_vs_target = np.nan

            # Determine support category
            if has_dms:
                if pd.notna(delta_vs_target) and delta_vs_target > 0:
                    support = "Improves vs target (DMS)"
                elif pd.notna(delta_vs_target):
                    support = "Worse vs target (DMS)"
                else:
                    support = "DMS measured"
            else:
                if score_method == "Physics" and pd.notna(delta_vs_target) and delta_vs_target > 0:
                    support = "Likely improves (physics)"
                elif score_method == "Physics" and pd.notna(delta_vs_target):
                    support = "Likely worse (physics)"
                elif score_method == "Physics":
                    support = "Physics estimate"
                else:
                    support = "No DMS record"

            dms_rows.append({
                "Mutation": mut,
                "Position": pos,
                "DMS Score": mut_score if has_dms else np.nan,
                "Physics Score": mut_score if score_method == "Physics" else np.nan,
                "Combined Score": mut_score,
                "Delta vs target": delta_vs_target,
                "Position mean score": pos_mean,
                "Position record count": pos_count,
                "Support": support,
                "Method": score_method,
                "Confidence": confidence,
            })

        dms_eval_df = pd.DataFrame(dms_rows)

        # Calculate coverage metrics
        dms_measured_count = int((dms_eval_df["Method"] == "DMS").sum()) if not dms_eval_df.empty else 0
        physics_estimated_count = int((dms_eval_df["Method"] == "Physics").sum()) if not dms_eval_df.empty else 0
        unavailable_count = int((dms_eval_df["Method"] == "Unavailable").sum()) if not dms_eval_df.empty else 0
        total_mutations = max(len(dms_eval_df), 1)
        dms_coverage = dms_measured_count / total_mutations if total_mutations else 0.0

        # Mean delta uses Combined Score (DMS where available, physics elsewhere)
        mean_delta = float(dms_eval_df["Delta vs target"].dropna().mean()) if not dms_eval_df.empty and dms_eval_df["Delta vs target"].notna().any() else np.nan

        # Average confidence (DMS=1.0, physics=0.3-0.6)
        avg_confidence = float(dms_eval_df["Confidence"].mean()) if not dms_eval_df.empty else 0.5

        function_norm = float(np.clip((active_candidate["score"] + 1.0) / 1.5, 0, 1))
        identity_norm = float(np.clip((active_candidate["identity"] - 80.0) / 20.0, 0, 1))
        stability_norm = float(np.clip((active_candidate["stability"] + 0.5) / 0.8, 0, 1))
        binding_norm = float(np.clip(active_candidate["binding"] / 10.0, 0, 1))
        physics_norm = (0.35 * function_norm) + (0.25 * identity_norm) + (0.2 * stability_norm) + (0.2 * binding_norm)

        # Evidence score now weights by confidence (DMS-heavy results score higher)
        dms_effect_norm = float(np.clip((mean_delta + 1.5) / 3.0, 0, 1)) if pd.notna(mean_delta) else 0.5
        evidence_score = float(100.0 * ((0.40 * physics_norm) + (0.30 * avg_confidence) + (0.30 * dms_effect_norm)))

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Evidence Score", f"{evidence_score:.0f}/100", delta="weighted by confidence")
        m2_delta = f"{dms_measured_count} DMS / {physics_estimated_count} physics"
        if unavailable_count > 0:
            m2_delta = f"{m2_delta} / {unavailable_count} unavailable"
        m2.metric("DMS Coverage", f"{dms_coverage*100:.0f}%", delta=m2_delta)
        m3.metric("Mean Δ vs target", f"{mean_delta:+.2f}" if pd.notna(mean_delta) else "n/a")
        if pd.notna(target_dms_score):
            target_score_text = f"{target_dms_score:+.2f}"
        elif physics_fallback_available:
            target_score_text = "physics fallback"
        else:
            target_score_text = "not in DMS"
        m4.metric("Target DMS Score", target_score_text)
        confidence_delta = "DMS=100%, Physics=30-60%" if physics_fallback_available else "DMS=100%, unavailable=0%"
        m5.metric("Avg Confidence", f"{avg_confidence*100:.0f}%", delta=confidence_delta)

        if target_for_eval in KNOWN_RESCUES:
            published = KNOWN_RESCUES[target_for_eval]["published"]
            overlap = sorted(set(rescue_mutations) & set([m.upper() for m in published]))
            if overlap:
                st.success(f"Published overlap for {target_for_eval}: {', '.join(overlap)}")
            else:
                st.info(f"Published rescues for {target_for_eval}: {', '.join(published)}")

        vcol1, vcol2 = st.columns([1.9, 1.1])
        with vcol1:
            if not dms_eval_df.empty:
                # Color map for DMS + physics categories
                color_map = {
                    "Improves vs target (DMS)": "#0F766E",      # Green (best - DMS validated)
                    "Worse vs target (DMS)": "#DC2626",         # Red (DMS shows worse)
                    "DMS measured": "#0284C7",                  # Blue (DMS, no target comparison)
                    "Likely improves (physics)": "#86EFAC",     # Light green (physics estimate)
                    "Likely worse (physics)": "#FCA5A5",        # Light red (physics estimate)
                    "Physics estimate": "#94A3B8",              # Gray (physics, no comparison)
                    "No DMS record": "#94A3B8",                 # Gray (no scoring data available)
                }

                # Always show delta vs target when available
                if dms_eval_df["Delta vs target"].notna().any():
                    dms_plot_df = dms_eval_df[dms_eval_df["Delta vs target"].notna()].copy()
                    fig_dms = px.bar(
                        dms_plot_df,
                        x="Mutation",
                        y="Delta vs target",
                        color="Support",
                        text="Delta vs target",
                        color_discrete_map=color_map,
                        title=f"Mutation effect relative to {target_for_eval} (DMS + Physics)",
                    )
                    fig_dms.add_hline(y=0, line_dash="dot", line_color="#64748B", line_width=1.3)
                    fig_dms.update_traces(texttemplate="%{text:.2f}", textposition="outside")
                    fig_dms.update_layout(height=540)
                    dms_count = len(dms_plot_df[dms_plot_df["Method"] == "DMS"])
                    physics_count = len(dms_plot_df[dms_plot_df["Method"] == "Physics"])
                    unavailable_plot_count = len(dms_plot_df[dms_plot_df["Method"] == "Unavailable"])
                    caption_msg = f"Showing {len(dms_plot_df)} mutations: {dms_count} DMS experimental, {physics_count} physics-based estimates"
                    if unavailable_plot_count > 0:
                        caption_msg = f"{caption_msg}, {unavailable_plot_count} unavailable"
                    st.caption(caption_msg)
                elif dms_eval_df["Combined Score"].notna().any():
                    # Show combined scores when no delta available
                    fig_dms = px.bar(
                        dms_eval_df,
                        x="Mutation",
                        y="Combined Score",
                        color="Method",
                        text="Combined Score",
                        color_discrete_map={"DMS": "#0284C7", "Physics": "#94A3B8", "Unavailable": "#94A3B8"},
                        title="Mutation scores (DMS experimental + Physics estimates)",
                    )
                    fig_dms.update_traces(texttemplate="%{text:.2f}", textposition="outside")
                    fig_dms.update_layout(height=540)
                else:
                    support_counts = (
                        dms_eval_df["Support"]
                        .value_counts(dropna=False)
                        .rename_axis("Support")
                        .reset_index(name="Count")
                    )
                    fig_dms = px.bar(
                        support_counts,
                        x="Support",
                        y="Count",
                        color="Support",
                        text="Count",
                        color_discrete_map=color_map,
                        title="DMS evidence availability for rescue mutations",
                    )
                    fig_dms.update_traces(texttemplate="%{text}", textposition="outside")
                    fig_dms.update_layout(height=540, yaxis_title="Number of rescue mutations")
                render_plotly(fig_dms, chart_kind="bar", width="stretch")
                st.caption(
                    "Interpretation: bars above 0 indicate rescue mutations predicted/measured to outperform the target mutation; "
                    "bar color indicates whether evidence is direct DMS measurement or fallback mode."
                )
            else:
                st.warning("No valid point-mutation rescue list was available for DMS validation.")

        with vcol2:
            gauge = go.Figure(
                go.Indicator(
                    mode="gauge+number",
                    value=evidence_score,
                    number={"suffix": "/100"},
                    title={"text": "Validation Evidence"},
                    gauge={
                        "axis": {"range": [0, 100]},
                        "bar": {"color": "#1E40AF"},
                        "steps": [
                            {"range": [0, 40], "color": "#FEE2E2"},
                            {"range": [40, 70], "color": "#FEF3C7"},
                            {"range": [70, 100], "color": "#DCFCE7"},
                        ],
                    },
                )
            )
            gauge.update_layout(height=440, margin=dict(l=20, r=20, t=58, b=12))
            render_plotly(gauge, chart_kind="default", width="stretch")
            st.caption("Evidence score combines model quality, measurement confidence, and rescue improvement signal (0-100).")

            physics_breakdown = pd.DataFrame(
                [
                    {"Component": "Function", "Weight": function_norm},
                    {"Component": "Identity", "Weight": identity_norm},
                    {"Component": "Stability", "Weight": stability_norm},
                    {"Component": "Binding", "Weight": binding_norm},
                ]
            )
            fig_phy = px.bar(
                physics_breakdown,
                x="Component",
                y="Weight",
                text="Weight",
                color="Weight",
                color_continuous_scale=["#D6ECFF", "#7CC7FA", "#0284C7", "#0F766E"],
                title="Physics subscore components",
            )
            fig_phy.update_traces(texttemplate="%{text:.2f}", textposition="outside")
            fig_phy.update_layout(coloraxis_showscale=False, yaxis_range=[0, 1.05], height=440)
            render_plotly(fig_phy, chart_kind="bar", width="stretch")
            st.caption("Each component is normalized to [0,1]; higher bars contribute more strongly to the physics subscore.")

        if not dms_eval_df.empty:
            show_df = dms_eval_df.copy()
            # Format score columns
            show_df["DMS Score"] = show_df["DMS Score"].map(lambda x: f"{x:+.2f}" if pd.notna(x) else "-")
            show_df["Physics Score"] = show_df["Physics Score"].map(lambda x: f"{x:+.2f}" if pd.notna(x) else "-")
            show_df["Combined Score"] = show_df["Combined Score"].map(lambda x: f"{x:+.2f}" if pd.notna(x) else "n/a")
            show_df["Delta vs target"] = show_df["Delta vs target"].map(lambda x: f"{x:+.2f}" if pd.notna(x) else "n/a")
            show_df["Position mean score"] = show_df["Position mean score"].map(lambda x: f"{x:+.2f}" if pd.notna(x) else "-")
            show_df["Confidence"] = show_df["Confidence"].map(lambda x: f"{x*100:.0f}%")
            # Select columns to display
            display_cols = ["Mutation", "Position", "DMS Score", "Physics Score", "Delta vs target", "Method", "Confidence", "Support"]
            st.dataframe(show_df[display_cols], hide_index=True, width="stretch")

            # Add data source info
            with st.expander("📊 About DMS Coverage & Physics Fallback"):
                mode_note = (
                    "**Mode:** DMS + physics fallback enabled.\n"
                    if physics_fallback_available
                    else "**Mode:** DMS-only fail-closed (physics fallback unavailable).\n"
                )
                st.markdown(mode_note + """
**Current DMS Dataset:** Giacomelli 2018 cell line data (~367 entries, 88 positions)

**Physics-Based Fallback** uses amino acid properties when DMS data unavailable:
- **BLOSUM62**: Evolutionary substitution scores
- **Hydropathy**: Kyte-Doolittle scale (burial effects)
- **Volume**: Steric clash potential
- **Charge**: Electrostatic effects

**For full saturation mutagenesis data** (~8000+ variants):
- [MaveDB](https://www.mavedb.org/) - Search "TP53"
- [Nature Genetics 2024](https://www.nature.com/articles/s41588-024-02039-4) - Deep CRISPR study
                """)
    else:
        st.info("Select a generated candidate to run real DMS + physics-based validation.")

    # === 2) EXPLAINABILITY ===
    st.markdown("---")
    st.markdown("### 2) 🧠 Explainability")
    st.caption("Mechanistic breakdown with attribution, counterfactuals, and energy terms.")

    if not EXPLAINABILITY_MODULE_AVAILABLE:
        st.warning("Explainability module is unavailable in this environment.")
        explain_reason = MODULE_IMPORT_ERRORS.get("explainability")
        if explain_reason:
            with st.expander("Explainability import details"):
                st.code(explain_reason, language=None)
    else:
        exp_default_target = active_candidate["target_mutation"] if active_candidate else "R175H"
        exp_default_rescues = ", ".join(active_candidate["mutations"][:6]) if active_candidate and active_candidate["mutations"] else "N239Y, T284R"
        st.caption(
            "Strict mode uses real model attention/occlusion/counterfactual terms only (no synthetic attribution fallback)."
        )
        col_exp1, col_exp2 = st.columns([1, 2])

        with col_exp1:
            if "exp_cancer_main" not in st.session_state:
                st.session_state["exp_cancer_main"] = exp_default_target
            if "exp_rescue_main" not in st.session_state:
                st.session_state["exp_rescue_main"] = exp_default_rescues

            exp_cancer_mut = st.text_input("Cancer mutation", key="exp_cancer_main").strip().upper()
            exp_rescue_muts = st.text_input(
                "Rescue mutations (comma-separated)",
                key="exp_rescue_main",
            )
            explain_status = get_explainability_backend_status()
            explain_ready = bool(explain_status.get("ready"))
            explain_block_reason = explain_status.get("reason", "Explainability dependencies are unavailable.")
            if not explain_ready:
                st.warning(f"Explainability strict mode blocked: {explain_block_reason}")

            if st.button(
                "Run Explainability Analysis",
                key="btn_explain_main",
                width="stretch",
                disabled=not explain_ready,
            ):
                if parse_single_mutation(exp_cancer_mut) is None:
                    st.error("Invalid cancer mutation format. Use format like `R175H`.")
                else:
                    rescue_list = parse_mutation_tokens(exp_rescue_muts)
                    if not rescue_list:
                        st.warning("No valid point-mutation rescue list was parsed.")
                    with st.spinner("Running explainability analysis..."):
                        mutant_seq = apply_mutation(P53_WT, exp_cancer_mut)
                        if mutant_seq is None:
                            mutant_seq = P53_WT
                        rescue_seq = build_sequence_from_mutations(exp_cancer_mut, rescue_list)

                        try:
                            engine = ExplainabilityEngine(
                                model=embedder.model,
                                tokenizer=embedder.tokenizer,
                                oracle=oracle,
                                embedder=embedder,
                            )
                            explanation = engine.explain_rescue(
                                P53_WT,
                                mutant_seq,
                                rescue_seq,
                                exp_cancer_mut,
                                rescue_list,
                            )
                            st.session_state["explanation"] = explanation
                            st.success("Explainability run complete.")
                        except ExplainabilityDependencyError as err:
                            st.session_state.pop("explanation", None)
                            st.error(
                                "Explainability dependencies are unavailable for strict mode "
                                f"({err})."
                            )
                            st.info(
                                "Troubleshooting: run once on CPU and retry. "
                                "If this persists, your transformer backend is not exposing attention tensors."
                            )
                        except Exception as err:
                            st.session_state.pop("explanation", None)
                            st.error(
                                "Explainability runtime error. "
                                "Try `cpu` runtime or rerun after model reload.\n"
                                f"Details: {err}"
                            )

        with col_exp2:
            if "explanation" in st.session_state:
                exp = st.session_state["explanation"]
                mech = exp["mechanism"]
                explain_status = get_explainability_backend_status()
                st.caption(
                    "Explainability backend: "
                    + ("strict-attention" if explain_status.get("ready") else "unavailable")
                    + f" ({explain_status.get('reason', 'n/a')})"
                )

                st.info(exp["summary"])
                em1, em2, em3, em4 = st.columns(4)
                em1.metric("Primary", mech["primary"].replace("_", " ").title())
                em2.metric("Mechanism confidence", f"{mech['confidence']:.0%}")
                em3.metric("Stability contrib", f"{mech['contributions']['stability']:+.2f}")
                em4.metric("Binding contrib", f"{mech['contributions']['binding']:+.2f}")

                evidence_lines = mech.get("evidence", [])[:4]
                if evidence_lines:
                    st.markdown("**Top evidence**")
                    for ev in evidence_lines:
                        st.markdown(f"- {ev}")

                energy = exp["energy"]
                energy_df = pd.DataFrame(
                    [
                        {"Component": "Electrostatic", "Contribution": energy["electrostatic"]},
                        {"Component": "Van der Waals", "Contribution": energy["van_der_waals"]},
                        {"Component": "H-Bonds", "Contribution": energy["hydrogen_bonds"]},
                        {"Component": "Solvation", "Contribution": energy["solvation"]},
                        {"Component": "Backbone", "Contribution": energy["backbone_strain"]},
                        {"Component": "Packing", "Contribution": energy["sidechain_packing"]},
                    ]
                )
                energy_df["Effect"] = np.where(energy_df["Contribution"] < 0, "Stabilizing", "Destabilizing")
                fig_energy = px.bar(
                    energy_df.sort_values("Contribution"),
                    x="Component",
                    y="Contribution",
                    color="Effect",
                    text="Contribution",
                    color_discrete_map={"Stabilizing": "#0891B2", "Destabilizing": "#DC2626"},
                    title="Energy decomposition",
                )
                fig_energy.update_traces(texttemplate="%{text:.2f}", textposition="outside")
                fig_energy.add_hline(y=0, line_dash="dot", line_color="#64748B", line_width=1.3)
                fig_energy.update_layout(height=520)
                render_plotly(fig_energy, chart_kind="bar", width="stretch")
                st.caption(
                    "Negative contributions are stabilizing, positive contributions are destabilizing; the total ddG summarizes net rescue energetics."
                )
                st.metric("Total ddG", f"{energy['total_ddg']:.2f} kcal/mol", delta=energy["rescue_quality"])

                attr_records = exp.get("top_attributions", [])
                if attr_records:
                    attr_df = pd.DataFrame(attr_records)
                    attr_df["Label"] = attr_df.apply(
                        lambda r: f"{r.get('residue', 'X')}{int(r.get('position', 0)) + 1}",
                        axis=1,
                    )
                    attr_df = attr_df.sort_values("composite_importance", ascending=False).head(14)
                    attr_df = attr_df.sort_values("composite_importance", ascending=True)
                    fig_attr = px.bar(
                        attr_df,
                        x="composite_importance",
                        y="Label",
                        color="functional_role",
                        orientation="h",
                        title="Top residue attributions",
                        color_discrete_sequence=["#1D4ED8", "#0284C7", "#06B6D4", "#0F766E", "#64748B"],
                    )
                    fig_attr.update_layout(height=540)
                    render_plotly(fig_attr, chart_kind="bar", width="stretch")
                    st.caption(
                        "Residue attribution combines attention and occlusion terms; larger bars indicate stronger influence on the rescue prediction."
                    )

                cf_rows = []
                for source_mut, entries in exp.get("counterfactuals", {}).items():
                    for item in entries[:3]:
                        cf_rows.append(
                            {
                                "Source": source_mut,
                                "Alternative": item["alternative"],
                                "Delta": item["delta_score"],
                            }
                        )
                if cf_rows:
                    cf_df = pd.DataFrame(cf_rows)
                    fig_cf = px.bar(
                        cf_df,
                        x="Alternative",
                        y="Delta",
                        color="Source",
                        text="Delta",
                        title="Counterfactual alternatives",
                    )
                    fig_cf.update_traces(texttemplate="%{text:.2f}", textposition="outside")
                    fig_cf.add_hline(y=0, line_dash="dot", line_color="#64748B", line_width=1.2)
                    fig_cf.update_layout(height=500)
                    render_plotly(fig_cf, chart_kind="bar", width="stretch")
                    st.caption(
                        "Counterfactual delta compares each suggested substitute against the selected rescue mutation at the same site."
                    )
                else:
                    st.info("No counterfactual alternatives were produced for these rescue mutations.")

                if mech.get("risk_factors"):
                    st.markdown("**Risk factors**")
                    for risk in mech["risk_factors"]:
                        st.warning(risk)
            else:
                st.info("Run explainability to generate mechanism and counterfactual charts.")

    # === 3) DRUG GENERATION ===
    st.markdown("---")
    st.markdown("### 3) 💊 Drug Candidate Generator")
    st.caption(
        "Template/denovo modes use model-estimated affinity. Docking modes use AutoDock Vina when optional "
        "dependencies and receptor PDBQT files are available."
    )

    if DRUG_MODULE_AVAILABLE:
        if active_candidate:
            st.markdown(
                f"""
                <div class="accent-box" style="margin-bottom:0.9rem;">
                    <b>Connected protein rescue:</b> {active_candidate['target_mutation']} with {active_candidate['n_mutations']} rescue mutations.
                </div>
                """,
                unsafe_allow_html=True,
            )

        col_dg1, col_dg2 = st.columns([1, 2])
        with col_dg1:
            pocket_names = list(P53_BINDING_POCKETS.keys())
            selected_pocket = st.selectbox("Target pocket", pocket_names, key="drug_pocket")
            pocket_info = P53_BINDING_POCKETS[selected_pocket]
            st.info(f"**{pocket_info.name}**\n\n{pocket_info.description}")
            st.caption(f"Druggability: {pocket_info.druggability_score:.0%}")

            n_candidates = st.slider("Number of candidates", 5, 30, 15, key="n_drug_cand")
            gen_method = st.selectbox(
                "Generation method",
                ["template", "docking", "docking_md"],
                key="drug_method",
            )
            allow_wt_receptor_fallback = st.checkbox(
                "Allow WT receptor fallback (less realistic)",
                value=False,
                key="allow_wt_receptor_fallback",
            )
            mode_cap = get_drug_mode_capability(
                selected_pocket,
                gen_method,
                allow_wt_receptor_fallback=allow_wt_receptor_fallback,
            )
            if not mode_cap.get("ready", False):
                st.warning(f"{gen_method} blocked: {mode_cap.get('reason', 'runtime unavailable')}")
            if mode_cap.get("receptor_path"):
                st.caption(
                    f"Receptor: {mode_cap.get('receptor_path')} ({mode_cap.get('receptor_source', 'unknown')})"
                )

            if st.button(
                "Generate Drug Candidates",
                key="btn_gen_drugs",
                width="stretch",
                disabled=not mode_cap.get("ready", False),
            ):
                with st.spinner("Generating drug candidates..."):
                    engine = DrugGeneratorEngine(
                        allow_wt_receptor_fallback=allow_wt_receptor_fallback,
                    )
                    protein_context = None
                    if active_candidate:
                        protein_context = {
                            "mutations": active_candidate.get("mutations", []),
                            "score": float(active_candidate.get("score", 0.0)),
                            "target_mutation": active_candidate.get("target_mutation", ""),
                        }
                    try:
                        generated = engine.generate_for_pocket(
                            selected_pocket,
                            n_candidates,
                            gen_method,
                            protein_context=protein_context,
                            allow_wt_receptor_fallback=allow_wt_receptor_fallback,
                        )
                    except Exception as e:
                        st.session_state.pop("drug_candidates", None)
                        st.error(f"Drug generation failed for method `{gen_method}`: {e}")
                    else:
                        st.session_state["drug_candidates"] = generated
                        if active_candidate:
                            st.session_state["drug_target_mutation"] = active_candidate["target_mutation"]
                            st.session_state["drug_rescue_protein"] = {
                                "profile": active_candidate["source"],
                                "mutations": active_candidate["mutations"],
                                "score": active_candidate["score"],
                            }
                        st.success(f"Generated {len(generated)} candidates.")

        with col_dg2:
            if "drug_candidates" in st.session_state and st.session_state["drug_candidates"]:
                candidates = st.session_state["drug_candidates"]
                candidates_sorted = sorted(candidates, key=lambda c: c.binding_affinity)
                props_df = pd.DataFrame(
                    [
                        {
                            "Name": c.name,
                            "MW": c.molecular_weight,
                            "LogP": c.logp,
                            "Drug-likeness": c.drug_likeness,
                            "BindingAffinity": c.binding_affinity,
                            "AffinityMagnitude": abs(c.binding_affinity),
                            "Lipinski": "Pass" if c.passes_lipinski() else "Fail",
                        }
                        for c in candidates_sorted
                    ]
                )
                props_df["Rank"] = np.arange(1, len(props_df) + 1)

                context_applied = any(bool(c.metadata.get("protein_context_applied", False)) for c in candidates)
                run_method = candidates_sorted[0].generation_method if candidates_sorted else "n/a"
                docking_run = run_method in {"docking", "docking_md"}

                s1, s2, s3, s4 = st.columns(4)
                s1.metric("Candidates", len(candidates))
                s2.metric("Best affinity", f"{props_df['BindingAffinity'].min():.2f} kcal/mol")
                s3.metric("Mean drug-likeness", f"{props_df['Drug-likeness'].mean():.2f}")
                s4.metric("Lipinski pass", f"{(props_df['Lipinski'] == 'Pass').mean()*100:.0f}%")
                st.caption(f"Scoring method: `{run_method}`")
                if context_applied:
                    lead = candidates_sorted[0]
                    st.caption(
                        "Protein-conditioned ranking active: lead affinity includes rescue-context shift "
                        f"({lead.metadata.get('protein_context_affinity_shift', 0.0):+.2f} kcal/mol), "
                        f"alignment={lead.metadata.get('protein_context_pocket_alignment', 0.0):.0%}."
                    )
                else:
                    if docking_run:
                        st.caption("Protein-conditioned ranking unavailable; using raw docking scores.")
                    else:
                        st.caption("Protein-conditioned ranking unavailable; using pocket-only heuristic scoring.")
                if docking_run:
                    lead = candidates_sorted[0]
                    st.caption(
                        "Docking backend: "
                        f"{lead.metadata.get('docking_backend', 'unknown')} | "
                        f"raw affinity={lead.metadata.get('docking_raw_affinity', np.nan):.2f} kcal/mol."
                    )
                    if run_method == "docking_md":
                        st.caption(
                            "MD status: "
                            f"{lead.metadata.get('md_status', 'unknown')} - "
                            f"{lead.metadata.get('md_detail', 'no detail')}"
                        )
                        if lead.metadata.get("md_status") != "ran":
                            st.warning("MD refinement did not execute for this run; see status above for missing setup/dependencies.")

                preview = props_df[["Name", "BindingAffinity", "Drug-likeness", "MW", "LogP", "Lipinski"]].copy().head(12)
                preview["BindingAffinity"] = preview["BindingAffinity"].map(lambda x: f"{x:.2f}")
                preview["Drug-likeness"] = preview["Drug-likeness"].map(lambda x: f"{x:.2f}")
                preview["MW"] = preview["MW"].map(lambda x: f"{x:.0f}")
                preview["LogP"] = preview["LogP"].map(lambda x: f"{x:.2f}")
                st.dataframe(preview, hide_index=True, width="stretch")

                props_plot_df = props_df.copy()
                props_plot_df["dup_count"] = props_plot_df.groupby(["MW", "LogP"])["Name"].transform("size")
                props_plot_df["dup_order"] = props_plot_df.groupby(["MW", "LogP"]).cumcount()
                props_plot_df["dup_centered"] = props_plot_df["dup_order"] - (props_plot_df["dup_count"] - 1) / 2.0
                props_plot_df["MW_plot"] = props_plot_df["MW"] + np.where(
                    props_plot_df["dup_count"] > 1,
                    props_plot_df["dup_centered"] * 3.0,
                    0.0,
                )
                props_plot_df["LogP_plot"] = props_plot_df["LogP"] + np.where(
                    props_plot_df["dup_count"] > 1,
                    props_plot_df["dup_centered"] * 0.05,
                    0.0,
                )

                hover_fields = {
                    "MW": ":.0f",
                    "LogP": ":.2f",
                    "BindingAffinity": ":.2f",
                    "MW_plot": False,
                    "LogP_plot": False,
                    "dup_count": False,
                    "dup_order": False,
                    "dup_centered": False,
                }

                if props_df["Drug-likeness"].nunique() > 1:
                    fig_props = px.scatter(
                        props_plot_df,
                        x="MW_plot",
                        y="LogP_plot",
                        color="Drug-likeness",
                        size="AffinityMagnitude",
                        symbol="Lipinski",
                        hover_name="Name",
                        hover_data=hover_fields,
                        color_continuous_scale=["#D8EEFF", "#6FB8F5", "#0284C7", "#0F766E"],
                        title="Drug property landscape",
                    )
                    fig_props.update_layout(coloraxis_colorbar_title_text="Drug-likeness")
                else:
                    fig_props = px.scatter(
                        props_plot_df,
                        x="MW_plot",
                        y="LogP_plot",
                        size="AffinityMagnitude",
                        symbol="Lipinski",
                        hover_name="Name",
                        hover_data=hover_fields,
                        title="Drug property landscape",
                    )
                    fig_props.update_traces(marker=dict(color="#0284C7"))

                fig_props.add_hline(y=5, line_dash="dash", line_color="#EF4444", annotation_text="LogP=5")
                fig_props.add_vline(x=500, line_dash="dash", line_color="#EF4444", annotation_text="MW=500")
                x_min = max(0.0, float(props_plot_df["MW_plot"].min()) - 20.0)
                x_max = float(props_plot_df["MW_plot"].max()) + 20.0
                y_min = float(props_plot_df["LogP_plot"].min()) - 0.5
                y_max = max(5.4, float(props_plot_df["LogP_plot"].max()) + 0.5)
                fig_props.update_xaxes(title="Molecular weight (Da)", range=[x_min, x_max])
                fig_props.update_yaxes(title="LogP", range=[y_min, y_max])
                fig_props.update_layout(height=580)
                render_plotly(fig_props, chart_kind="scatter", width="stretch")
                st.caption(
                    "Each point is a candidate; size tracks affinity magnitude, and dashed thresholds mark Lipinski-style MW/LogP cutoffs."
                )
                if (props_plot_df["dup_count"] > 1).any():
                    st.caption("Overlapping MW/LogP points are slightly jittered so each candidate remains visible.")

                rank_n = min(10, len(props_df))
                rank_df = props_df.sort_values("BindingAffinity", ascending=True).head(rank_n).copy()
                rank_title = (
                    f"Top {rank_n} candidates by docking score"
                    if docking_run
                    else f"Top {rank_n} candidates by predicted affinity"
                )
                fig_rank = px.bar(
                    rank_df,
                    x="Name",
                    y="BindingAffinity",
                    color="Drug-likeness",
                    text="BindingAffinity",
                    color_continuous_scale=["#DAEEFF", "#70C2F8", "#0284C7", "#0F766E"],
                    title=rank_title,
                )
                fig_rank.update_traces(texttemplate="%{text:.2f}", textposition="outside")
                fig_rank.update_xaxes(categoryorder="array", categoryarray=rank_df["Name"].tolist(), tickangle=-20)
                fig_rank.update_layout(height=540, yaxis_title="Binding affinity (kcal/mol)")
                render_plotly(fig_rank, chart_kind="bar", width="stretch")
                if docking_run:
                    st.caption(
                        "More negative bars indicate stronger AutoDock Vina docking scores (or context-adjusted scores if rescue context is active)."
                    )
                else:
                    st.caption(
                        "More negative bars indicate stronger predicted pocket binding under the current heuristic model."
                    )

                top = candidates_sorted[0]
                drug_col1, drug_col2 = st.columns([2, 1])
                with drug_col1:
                    drug_3d_html = f"""
                    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                    <div id="drug_viewer" style="width:100%; height:390px; border-radius:14px; border:2px solid #10B981; background:linear-gradient(135deg,#F8FCFF 0%,#EAF3FF 100%);"></div>
                    <script>
                        (async function() {{
                            let viewer = $3Dmol.createViewer('drug_viewer', {{backgroundColor: '0xF8FCFF'}});
                            let smiles = "{top.smiles}";
                            try {{
                                let response = await fetch(
                                    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' +
                                    encodeURIComponent(smiles) + '/SDF?record_type=3d'
                                );
                                if (response.ok) {{
                                    viewer.addModel(await response.text(), 'sdf');
                                }} else {{
                                    viewer.addModel(smiles, 'smiles');
                                }}
                            }} catch(e) {{
                                viewer.addModel(smiles, 'smiles');
                            }}
                            viewer.setStyle({{}}, {{stick: {{radius: 0.15, colorscheme: 'Jmol'}}, sphere: {{radius: 0.4, colorscheme: 'Jmol'}}}});
                            viewer.zoomTo();
                            viewer.spin('y', 0.7);
                            viewer.render();
                        }})();
                    </script>
                    """
                    components.html(drug_3d_html, height=420)

                with drug_col2:
                    st.metric("Lead candidate", top.name)
                    st.metric("Binding affinity", f"{top.binding_affinity:.2f} kcal/mol")
                    st.metric("Drug-likeness", f"{top.drug_likeness:.2f}")
                    st.metric("Synthetic access", f"{top.synthetic_accessibility:.1f}/10")
                    st.metric("Molecular weight", f"{top.molecular_weight:.0f} Da")
                    st.metric("LogP", f"{top.logp:.1f}")
                    st.code(top.smiles, language=None)
                    st.caption(
                        "3D view is fetched from PubChem when available; otherwise rendered directly from SMILES."
                    )
            else:
                st.info("Generate drug candidates to view redesigned property and ranking charts.")
    else:
        st.warning("Drug generation module is unavailable in this environment.")
        drug_reason = MODULE_IMPORT_ERRORS.get("drug")
        if not drug_reason and not P53_BINDING_POCKETS:
            drug_reason = "No configured binding pockets were loaded."
        if drug_reason:
            with st.expander("Drug module import details"):
                st.code(drug_reason, language=None)

    # === 4) CLINICAL IMPACT ===
    st.markdown("---")
    st.markdown("### 4) 🏥 Clinical Impact")
    st.caption("Clinical outputs are scenario-based estimates from the current model assumptions.")

    if CLINICAL_MODULE_AVAILABLE:
        ci_default_target = active_candidate["target_mutation"] if active_candidate else "R175H"
        ci_default_rescues = ", ".join(active_candidate["mutations"][:6]) if active_candidate and active_candidate["mutations"] else "N239Y, T284R"
        col_ci1, col_ci2 = st.columns([1, 2])

        with col_ci1:
            if "ci_cancer_main" not in st.session_state:
                st.session_state["ci_cancer_main"] = ci_default_target
            if "ci_rescue_main" not in st.session_state:
                st.session_state["ci_rescue_main"] = ci_default_rescues

            ci_cancer = st.text_input("Cancer mutation", key="ci_cancer_main").strip().upper()
            ci_rescue_text = st.text_input("Rescue mutations (comma-separated)", key="ci_rescue_main")

            if st.button("Assess Clinical Impact", key="btn_clinical_main", width="stretch"):
                if parse_single_mutation(ci_cancer) is None:
                    st.error("Invalid cancer mutation format. Use format like `R175H`.")
                else:
                    rescue_list = parse_mutation_tokens(ci_rescue_text)
                    rescue_seq = build_sequence_from_mutations(ci_cancer, rescue_list)
                    with st.spinner("Assessing clinical impact..."):
                        engine = ClinicalImpactEngine()
                        report = engine.generate_report(
                            name=f"p53_{ci_cancer}_rescue",
                            wt_sequence=P53_WT,
                            rescue_sequence=rescue_seq,
                            cancer_mutation=ci_cancer,
                            rescue_mutations=rescue_list,
                        )
                        st.session_state["clinical_report"] = report
                        st.success("Clinical impact assessment complete.")

        with col_ci2:
            if "clinical_report" in st.session_state:
                report = st.session_state["clinical_report"]
                pop = report.patient_population

                cv1, cv2, cv3 = st.columns(3)
                cv1.metric("Clinical score", f"{report.overall_clinical_score:.0f}/100")
                cv2.metric("Viability", report.clinical_viability.upper())
                cv3.metric("US annual patients", f"{pop.total_patients_per_year:,}")

                cp1, cp2 = st.columns(2)
                cp1.metric("Global estimate", f"{pop.global_estimate:,}")
                if hasattr(pop, "mutation_frequency"):
                    cp2.metric("Mutation frequency", f"{float(pop.mutation_frequency):.2f}%")
                else:
                    cp2.metric("Mutation frequency", "n/a")

                delivery_df = pd.DataFrame(
                    [
                        {
                            "Method": d.method,
                            "Feasibility": float(d.feasibility_score),
                            "Cost": d.estimated_cost_per_dose,
                        }
                        for d in report.delivery_options[:5]
                    ]
                )
                show_delivery = delivery_df.copy()
                show_delivery["Feasibility"] = show_delivery["Feasibility"].map(lambda x: f"{x:.0%}")
                st.dataframe(show_delivery, hide_index=True, width="stretch")

                fig_delivery = px.bar(
                    delivery_df.sort_values("Feasibility", ascending=False),
                    x="Method",
                    y="Feasibility",
                    text="Feasibility",
                    color="Feasibility",
                    color_continuous_scale=["#D6ECFF", "#7CC7FA", "#0284C7", "#0F766E"],
                    title="Delivery feasibility",
                )
                fig_delivery.update_traces(texttemplate="%{text:.0%}", textposition="outside")
                fig_delivery.update_layout(coloraxis_showscale=False, yaxis_range=[0, 1.05], height=520)
                render_plotly(fig_delivery, chart_kind="bar", width="stretch")
                st.caption("Higher feasibility scores indicate fewer delivery/manufacturing/regulatory barriers under current assumptions.")

                # Transparent score construction so each dashboard number is auditable.
                if pop.total_patients_per_year >= 10000:
                    population_points = 25.0
                elif pop.total_patients_per_year >= 5000:
                    population_points = 20.0
                elif pop.total_patients_per_year >= 1000:
                    population_points = 15.0
                elif pop.total_patients_per_year >= 100:
                    population_points = 10.0
                else:
                    population_points = 5.0

                ti = report.therapeutic_index
                if ti.sequence_identity >= 98:
                    identity_points = 15.0
                elif ti.sequence_identity >= 95:
                    identity_points = 10.0
                elif ti.sequence_identity >= 92:
                    identity_points = 5.0
                else:
                    identity_points = 0.0

                window_points = 10.0 if ti.therapeutic_window == "wide" else (5.0 if ti.therapeutic_window == "moderate" else 0.0)
                best_delivery = report.delivery_options[0] if report.delivery_options else None
                delivery_points = float(best_delivery.feasibility_score * 25.0) if best_delivery else 0.0
                immuno_points = float(25.0 - report.immunogenicity.risk_score * 25.0)

                breakdown_df = pd.DataFrame(
                    [
                        {
                            "Component": "Patient population",
                            "Value": f"{pop.total_patients_per_year:,} US patients/year",
                            "Points": f"{population_points:.1f}/25",
                            "Computation": "Tiered rule from annual US patient count",
                        },
                        {
                            "Component": "Sequence identity",
                            "Value": f"{ti.sequence_identity:.1f}%",
                            "Points": f"{identity_points:.1f}/15",
                            "Computation": ">=98%=15, >=95%=10, >=92%=5, else 0",
                        },
                        {
                            "Component": "Therapeutic window",
                            "Value": ti.therapeutic_window,
                            "Points": f"{window_points:.1f}/10",
                            "Computation": "wide=10, moderate=5, narrow=0",
                        },
                        {
                            "Component": "Best delivery option",
                            "Value": best_delivery.method if best_delivery else "n/a",
                            "Points": f"{delivery_points:.1f}/25",
                            "Computation": "best feasibility score x 25",
                        },
                        {
                            "Component": "Immunogenicity",
                            "Value": f"risk={report.immunogenicity.risk_score:.2f}",
                            "Points": f"{immuno_points:.1f}/25",
                            "Computation": "25 - (risk_score x 25)",
                        },
                    ]
                )
                with st.expander("How each clinical number is computed"):
                    st.markdown(
                        f"""
                        - `US annual patients` = sum over cancer types of:
                          `incidence x p53 mutation rate x mutation frequency` for `{report.cancer_mutation}`.
                        - `Global estimate` = `US annual patients x 3` (rough multiplier used by this model).
                        - `Clinical score` = sum of component points shown below (clipped to 0-100).
                        - `Viability` thresholds: high >=70, moderate >=50, low >=30, else not viable.
                        """
                    )
                    st.dataframe(breakdown_df, hide_index=True, width="stretch")

                st.success(report.recommended_development_path)
                export_data = {
                    "name": report.rescue_name,
                    "clinical_score": report.overall_clinical_score,
                    "viability": report.clinical_viability,
                    "recommendation": report.recommended_development_path,
                }
                st.download_button(
                    "Download Clinical Report (JSON)",
                    data=json.dumps(export_data, indent=2),
                    file_name=f"{report.rescue_name}_clinical.json",
                    mime="application/json",
                    width="stretch",
                )
            else:
                st.info("Run clinical assessment to populate feasibility and population charts.")
    else:
        st.warning("Clinical impact module is unavailable in this environment.")
        clinical_reason = MODULE_IMPORT_ERRORS.get("clinical")
        if clinical_reason:
            with st.expander("Clinical module import details"):
                st.code(clinical_reason, language=None)
# --- FOOTER ---
st.markdown("""
<div class="footer">
    <p style="margin-bottom: 0.5rem; font-weight: 600; color: #1A1A1A;">
        p53-proteoMgCAD
    </p>
    <p style="color: #9CA3AF; font-size: 0.85rem; margin-bottom: 1rem;">
        Mutative Generative Computer-Assisted Design of Second-Site Rescues
    </p>
    <p style="font-size: 0.75rem; color: #D1D5DB;">
        ISEF 2026 &middot; Advancing Generative Protein Design for Precision Oncology
    </p>
</div>
""", unsafe_allow_html=True)
