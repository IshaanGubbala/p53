import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import torch
import torch.nn.functional as F
import requests
import time
import json

import importlib
from p53cad.engine.latent import ManifoldEmbedder, ManifoldWalker
from p53cad.engine.oracle import FunctionalOracle
from p53cad.engine.explain import SaliencyMap
from p53cad.analysis.grassmann import GrassmannMetric
from p53cad.data.dms import P53_WT, apply_mutation
from p53cad.viz.pymol import PyMolGenerator

# New Enhancement Modules
try:
    from p53cad.engine.drug_generator import DrugGeneratorEngine, P53_BINDING_POCKETS
    from p53cad.engine.explainability import ExplainabilityEngine, P53_FUNCTIONAL_SITES
    from p53cad.engine.experimental import ExperimentalPipeline
    from p53cad.analysis.clinical_impact import ClinicalImpactEngine, PatientStratifier
    from p53cad.engine.multi_target import MultiTargetPlatform, TUMOR_SUPPRESSORS
    from p53cad.engine.md_validation import MDValidationEngine
    ENHANCED_MODULES_AVAILABLE = True
except ImportError as e:
    ENHANCED_MODULES_AVAILABLE = False
    print(f"Enhanced modules not fully loaded: {e}")

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

def morph_pdb_for_mutations(pdb_string: str, mutation_positions: list, step: int, total_steps: int) -> str:
    """
    Morph a PDB structure to show conformational changes around mutation sites.
    This creates a visual effect of the protein structure changing during optimization.
    """
    if not pdb_string or not mutation_positions:
        return pdb_string

    lines = pdb_string.split('\n')
    new_lines = []

    # Calculate morph intensity based on optimization progress
    progress = step / max(total_steps, 1)
    wave_phase = step * 0.2  # Create wave-like motion

    for line in lines:
        if line.startswith('ATOM'):
            try:
                # Parse residue number
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Calculate displacement based on proximity to mutations
                displacement = 0.0
                for mut_pos in mutation_positions:
                    distance_to_mut = abs(res_num - mut_pos)
                    if distance_to_mut < 15:  # Affect residues within 15 positions
                        # Stronger effect closer to mutation
                        effect = (15 - distance_to_mut) / 15.0
                        # Wave-like displacement
                        wave = np.sin(wave_phase + distance_to_mut * 0.3)
                        displacement += effect * wave * 1.5 * progress

                # Apply displacement (primarily in Y direction for visual effect)
                new_x = x + displacement * 0.3 * np.cos(wave_phase)
                new_y = y + displacement * np.sin(wave_phase * 0.7)
                new_z = z + displacement * 0.2 * np.sin(wave_phase * 1.3)

                new_line = f"{line[:30]}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}{line[54:]}"
                new_lines.append(new_line)
            except:
                new_lines.append(line)
        else:
            new_lines.append(line)

    return '\n'.join(new_lines)

def generate_motion_frames(pdb_string: str, n_frames: int = 20) -> list:
    """Generate pseudo-motion frames by perturbing coordinates (simplified NMA-like motion)."""
    import re
    frames = [pdb_string]

    # Parse ATOM lines
    lines = pdb_string.split('\n')
    atom_lines = [l for l in lines if l.startswith('ATOM')]
    other_lines = [l for l in lines if not l.startswith('ATOM')]

    for frame_idx in range(1, n_frames):
        # Sinusoidal breathing motion
        scale = 0.3 * np.sin(2 * np.pi * frame_idx / n_frames)
        new_atoms = []

        for line in atom_lines:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Add harmonic motion (simplified normal mode)
                # Breathing mode - expand/contract from center
                cx, cy, cz = 0, 0, 0  # Approximate center
                dx, dy, dz = x - cx, y - cy, z - cz
                dist = np.sqrt(dx**2 + dy**2 + dz**2) + 0.1

                new_x = x + scale * dx / dist * 2
                new_y = y + scale * dy / dist * 2
                new_z = z + scale * dz / dist * 2

                new_line = f"{line[:30]}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}{line[54:]}"
                new_atoms.append(new_line)
            except:
                new_atoms.append(line)

        frame_pdb = '\n'.join(other_lines[:2] + new_atoms + other_lines[2:])
        frames.append(frame_pdb)

    return frames

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
    /* Import Google Font */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    :root {
        --bg-base: #ECF2FA;
        --bg-tint-1: #D6E7FF;
        --bg-tint-2: #D1FAE5;
        --bg-tint-3: #E0F2FE;
        --surface: #F8FBFF;
        --surface-muted: #EEF4FF;
        --surface-strong: #E6F0FF;
        --border-soft: #D5E1F3;
        --accent: #2563EB;
        --accent-strong: #1D4ED8;
        --accent-soft: #DBEAFE;
        --accent-2: #0EA5E9;
        --accent-2-soft: #CFFAFE;
        --accent-dark: #0B1B3B;
    }

    /* ===== GLOBAL RESET & BASE ===== */
    html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"] {
        background: radial-gradient(1200px 700px at 8% -10%, var(--bg-tint-1), transparent 62%),
                    radial-gradient(900px 600px at 92% -20%, var(--bg-tint-2), transparent 58%),
                    radial-gradient(700px 520px at 50% 0%, var(--bg-tint-3), transparent 60%),
                    var(--bg-base) !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
    }

    .main .block-container {
        max-width: 1200px !important;
        padding: 2rem 2.5rem 3rem 2.5rem !important;
    }

    /* Hide Streamlit chrome */
    #MainMenu, footer, header, .stDeployButton {display: none !important; visibility: hidden !important;}

    /* ===== TYPOGRAPHY (Larger for readability) ===== */
    html, body {
        font-size: 17px !important;
    }

    h1, h2, h3, h4, h5, h6, .stMarkdown h1, .stMarkdown h2, .stMarkdown h3 {
        font-family: 'Inter', sans-serif !important;
        color: #0F172A !important;
        font-weight: 600 !important;
        letter-spacing: -0.025em !important;
    }

    h1, .stMarkdown h1 { font-size: 2.5rem !important; font-weight: 700 !important; }
    h2, .stMarkdown h2 { font-size: 1.85rem !important; }
    h3, .stMarkdown h3 { font-size: 1.4rem !important; }
    h4, .stMarkdown h4 { font-size: 1.2rem !important; }

    p, span, label, .stMarkdown p, .stMarkdown li {
        font-family: 'Inter', sans-serif !important;
        color: #374151 !important;
        line-height: 1.7 !important;
        font-size: 1.05rem !important;
    }

    /* ===== HERO HEADER ===== */
    .hero-container {
        text-align: center;
        padding: 2.5rem 1.5rem 1.75rem 1.5rem;
        margin-bottom: 1.75rem;
        background: linear-gradient(135deg, #FFFFFF 0%, #F1F7FF 55%, #E0F2FE 100%);
        border: 1px solid var(--border-soft);
        border-radius: 18px;
        box-shadow: 0 16px 30px rgba(15,23,42,0.10), 0 4px 10px rgba(30,64,175,0.12);
        position: relative;
        overflow: hidden;
    }

    .hero-container::before {
        content: "";
        position: absolute;
        top: 0;
        left: 0;
        height: 5px;
        width: 100%;
        background: linear-gradient(90deg, var(--accent-dark), var(--accent), #38BDF8);
    }

    .hero-container::after {
        content: "";
        position: absolute;
        right: -40px;
        top: -60px;
        width: 220px;
        height: 220px;
        background: radial-gradient(circle at 40% 40%, rgba(56, 189, 248, 0.35), transparent 65%);
        filter: blur(2px);
    }

    .hero-badge {
        display: inline-block;
        background: linear-gradient(135deg, var(--accent-dark) 0%, var(--accent) 100%);
        color: #E0F2FE;
        font-size: 0.85rem;
        font-weight: 600;
        padding: 0.4rem 1rem;
        border-radius: 50px;
        margin-bottom: 0.75rem;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        box-shadow: 0 10px 20px rgba(11,27,59,0.25);
    }

    .hero-title {
        font-family: 'Inter', sans-serif !important;
        font-size: 2.75rem !important;
        font-weight: 700 !important;
        color: var(--accent-dark) !important;
        letter-spacing: -0.03em;
        margin: 0 0 0.5rem 0 !important;
        line-height: 1.2;
    }

    .hero-subtitle {
        font-size: 1.2rem;
        color: #4B5563;
        max-width: 600px;
        margin: 0 auto 0.5rem auto;
        line-height: 1.5;
        font-weight: 400;
    }

    .hero-description {
        font-size: 1rem;
        color: #6B7280;
        max-width: 500px;
        margin: 0 auto;
    }

    /* ===== CARD COMPONENTS ===== */
    .lovable-card {
        background: linear-gradient(180deg, #FFFFFF 0%, var(--surface) 100%);
        border: 1px solid var(--border-soft);
        border-radius: 14px;
        padding: 1.25rem;
        margin-bottom: 1rem;
        box-shadow: 0 10px 24px rgba(30,64,175,0.08), 0 2px 6px rgba(15,23,42,0.06);
        position: relative;
        overflow: hidden;
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
        color: #0F172A !important;
        margin: 0 0 0.5rem 0 !important;
        font-weight: 600 !important;
    }

    .lovable-card p {
        color: #4B5563 !important;
        margin: 0 !important;
    }

    .lovable-card-accent {
        background: linear-gradient(135deg, var(--accent-soft) 0%, #E0F2FE 55%, #D1FAE5 100%);
        border: 1px solid #93C5FD;
        border-left: 4px solid var(--accent-dark);
        border-radius: 14px;
        padding: 1.25rem;
        margin-bottom: 1rem;
        box-shadow: 0 12px 24px rgba(37,99,235,0.12);
    }

    .lovable-card-accent h3 {
        color: #1E3A5F !important;
        margin: 0 0 0.25rem 0 !important;
        font-size: 1rem !important;
        font-weight: 600 !important;
    }

    .lovable-card-accent p {
        color: #475569 !important;
        font-size: 0.85rem !important;
    }

    /* Important data highlight */
    .data-highlight {
        background: #E0F2FE;
        border: 1px solid #38BDF8;
        border-radius: 8px;
        padding: 0.75rem 1rem;
        font-weight: 600;
        color: #0C4A6E;
    }

    /* ===== SIDEBAR ===== */
    [data-testid="stSidebar"] {
        background: var(--surface) !important;
        border-right: 1px solid var(--border-soft) !important;
    }

    [data-testid="stSidebar"] [data-testid="stMarkdown"] {
        color: #374151 !important;
    }

    .section-header {
        font-size: 0.7rem !important;
        font-weight: 600 !important;
        color: #9CA3AF !important;
        text-transform: uppercase !important;
        letter-spacing: 0.08em !important;
        margin-bottom: 0.75rem !important;
    }

    .feature-pill {
        display: inline-block;
        background: var(--surface-strong);
        color: #334155;
        font-size: 0.7rem;
        font-weight: 500;
        padding: 0.2rem 0.6rem;
        border-radius: 50px;
        margin-right: 0.4rem;
    }

    /* ===== TABS (Larger) ===== */
    .stTabs [data-baseweb="tab-list"] {
        background: var(--surface-muted) !important;
        border-radius: 12px !important;
        padding: 6px !important;
        gap: 6px !important;
        border: none !important;
    }

    .stTabs [data-baseweb="tab"] {
        background: transparent !important;
        border-radius: 10px !important;
        color: #4B5563 !important;
        font-weight: 500 !important;
        font-size: 1rem !important;
        padding: 0.7rem 1.25rem !important;
        border: none !important;
    }

    .stTabs [data-baseweb="tab"]:hover {
        color: #1F2937 !important;
    }

    .stTabs [aria-selected="true"] {
        background: var(--surface) !important;
        color: #0F172A !important;
        font-weight: 600 !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.08) !important;
    }

    /* Tab highlight line removal */
    .stTabs [data-baseweb="tab-highlight"] {
        display: none !important;
    }

    /* ===== RADIO BUTTONS (Larger) ===== */
    .stRadio > div {
        flex-direction: row !important;
        flex-wrap: wrap !important;
        gap: 0.6rem !important;
    }

    .stRadio > div > label {
        background: var(--surface) !important;
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 10px !important;
        padding: 0.65rem 1.2rem !important;
        margin: 0 !important;
        font-size: 1rem !important;
        color: #374151 !important;
        cursor: pointer !important;
        transition: all 0.15s ease !important;
    }

    .stRadio > div > label:hover {
        border-color: var(--accent) !important;
        background: var(--surface-strong) !important;
    }

    .stRadio > div > label[data-checked="true"],
    .stRadio > div > label:has(input:checked) {
        background: var(--accent-soft) !important;
        border-color: var(--accent) !important;
        color: #1E3A8A !important;
        font-weight: 500 !important;
    }

    /* Hide radio circles */
    .stRadio [role="radiogroup"] > label > div:first-child {
        display: none !important;
    }

    /* ===== BUTTONS (Larger) ===== */
    .stButton > button {
        background: #0F172A !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 10px !important;
        padding: 0.75rem 1.5rem !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        transition: all 0.15s ease !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1) !important;
    }

    .stButton > button:hover {
        background: #1E293B !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 4px 8px rgba(0,0,0,0.15) !important;
    }

    .stButton > button[kind="primary"] {
        background: var(--accent) !important;
    }

    .stButton > button[kind="primary"]:hover {
        background: var(--accent-strong) !important;
    }

    .stDownloadButton > button {
        background: var(--surface) !important;
        color: #1F2937 !important;
        border: 1.5px solid var(--border-soft) !important;
        font-size: 0.95rem !important;
        padding: 0.65rem 1.25rem !important;
    }

    .stDownloadButton > button:hover {
        background: var(--surface-muted) !important;
        border-color: #9CA3AF !important;
    }

    /* ===== INPUT FIELDS (Larger) ===== */
    .stTextInput > div > div > input {
        background: var(--surface) !important;
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 10px !important;
        padding: 0.75rem 1rem !important;
        font-size: 1rem !important;
        color: #0F172A !important;
    }

    .stTextInput > div > div > input:focus {
        border-color: var(--accent) !important;
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.18) !important;
    }

    .stSelectbox > div > div {
        background: var(--surface) !important;
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 10px !important;
        font-size: 1rem !important;
    }

    /* ===== METRICS (Larger) ===== */
    [data-testid="stMetric"] {
        background: var(--surface) !important;
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 12px !important;
        padding: 1.25rem !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.06) !important;
    }

    [data-testid="stMetricLabel"] {
        font-size: 0.85rem !important;
        font-weight: 600 !important;
        color: #6B7280 !important;
        text-transform: uppercase !important;
        letter-spacing: 0.06em !important;
    }

    [data-testid="stMetricValue"] {
        font-size: 1.85rem !important;
        font-weight: 700 !important;
        color: #0F172A !important;
    }

    /* Metric delta (positive/negative indicator) */
    [data-testid="stMetricDelta"] {
        font-weight: 600 !important;
        font-size: 0.95rem !important;
    }

    [data-testid="stMetricDelta"] svg {
        stroke-width: 3px !important;
    }

    /* ===== ALERTS (Larger) ===== */
    [data-testid="stAlert"] > div {
        border-radius: 12px !important;
        padding: 1.25rem !important;
        font-size: 1rem !important;
    }

    /* ===== EXPANDERS (Larger) ===== */
    .streamlit-expanderHeader {
        background: var(--surface-muted) !important;
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 12px !important;
        font-weight: 600 !important;
        font-size: 1.05rem !important;
        padding: 1rem !important;
    }

    /* ===== DATAFRAMES (Larger) ===== */
    .stDataFrame {
        border: 1.5px solid var(--border-soft) !important;
        border-radius: 12px !important;
        overflow: hidden !important;
        font-size: 0.95rem !important;
    }

    /* ===== DIVIDERS ===== */
    hr, .stDivider {
        border: none !important;
        border-top: 1.5px solid #E5E7EB !important;
        margin: 2rem 0 !important;
    }

    /* ===== SLIDERS (Larger) ===== */
    .stSlider [data-baseweb="slider"] [role="slider"] {
        background: var(--accent) !important;
        width: 20px !important;
        height: 20px !important;
    }

    .stSlider [data-testid="stTickBarMin"],
    .stSlider [data-testid="stTickBarMax"] {
        font-size: 0.9rem !important;
        color: #6B7280 !important;
    }

    /* ===== PLOTLY CHARTS ===== */
    .js-plotly-plot .plotly {
        border-radius: 12px !important;
    }

    /* ===== CUSTOM UTILITY CLASSES (Larger) ===== */
    .design-card {
        background: linear-gradient(180deg, #FFFFFF 0%, var(--surface) 100%);
        border: 1.5px solid var(--border-soft);
        border-radius: 14px;
        padding: 1.25rem;
        margin-bottom: 1rem;
        box-shadow: 0 10px 20px rgba(14,116,144,0.08), 0 2px 6px rgba(15,23,42,0.06);
        position: relative;
        overflow: hidden;
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
        color: #0F172A;
        font-weight: 600;
        font-size: 1.1rem;
    }

    .constraint-card {
        background: linear-gradient(180deg, #FFFFFF 0%, #F0F9FF 100%);
        border: 1.5px solid var(--border-soft);
        border-radius: 14px;
        padding: 1.25rem;
        box-shadow: 0 10px 20px rgba(59,130,246,0.10), 0 2px 6px rgba(15,23,42,0.06);
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
        background: linear-gradient(135deg, var(--accent-soft) 0%, #E0F2FE 60%);
        border: 1px solid #93C5FD;
        border-left: 4px solid var(--accent-dark);
        border-radius: 14px;
        padding: 1rem 1.25rem;
        box-shadow: 0 10px 20px rgba(37,99,235,0.12);
    }

    .gen-design-header {
        font-size: 1.75rem;
        font-weight: 700;
        color: #0F172A;
        margin-bottom: 0.5rem;
    }

    /* ===== SPINNER ===== */
    .stSpinner > div {
        border-color: var(--accent) !important;
    }

    /* ===== CHECKBOX ===== */
    .stCheckbox label {
        font-size: 0.875rem !important;
        color: #374151 !important;
    }

    /* ===== MULTISELECT ===== */
    .stMultiSelect > div > div {
        background: var(--surface) !important;
        border: 1px solid #E5E7EB !important;
        border-radius: 8px !important;
    }

    .stMultiSelect [data-baseweb="tag"] {
        background: var(--accent-soft) !important;
        border-radius: 6px !important;
    }
</style>
""", unsafe_allow_html=True)

# Initialize models (NO CACHE - force fresh reload)
def load_models_v8():
    # Force reload of core engine to bypass stale Streamlit modules
    import sys
    if 'p53cad.engine.latent' in sys.modules:
        del sys.modules['p53cad.engine.latent']
    if 'p53cad.engine.oracle' in sys.modules:
        del sys.modules['p53cad.engine.oracle']
    if 'p53cad.engine.explain' in sys.modules:
        del sys.modules['p53cad.engine.explain']
    if 'p53cad.analysis.grassmann' in sys.modules:
        del sys.modules['p53cad.analysis.grassmann']
    
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
            if pooled.shape[-1] != 320:
                 pooled = pooled[:, :320] 

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
                                progress_callback=None, structure_callback=None):
    """
    Generative Design Mode with LIVE VISUALIZATION.
    Like watching mechanical CAD topology optimization in real-time.

    Yields intermediate states for real-time UI updates.
    """
    target_mut = constraints.get('target_mutation', 'R175H')

    cancer_seq = apply_mutation(P53_WT, target_mut)
    if cancer_seq is None:
        cancer_seq = P53_WT

    AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

    min_identity = constraints.get('min_identity', 90.0)
    min_stability = constraints.get('min_stability', -0.3)
    min_binding = constraints.get('min_binding', 5.0)
    locked_positions = constraints.get('locked_positions', [248, 273])
    delivery_method = constraints.get('delivery_method', 'gene_therapy')
    exploration_diversity = constraints.get('diversity', 0.5)

    if delivery_method == 'protein_therapy':
        min_identity = max(min_identity, 95.0)
    elif delivery_method == 'mrna_therapy':
        min_identity = max(min_identity, 92.0)

    all_candidates = []

    weight_profiles = [
        {'function': 4.0, 'stability': 8.0, 'binding': 2.5, 'name': 'Balanced', 'color': '#00D4FF'},
        {'function': 2.0, 'stability': 15.0, 'binding': 2.0, 'name': 'Stability-First', 'color': '#FFD700'},
        {'function': 3.0, 'stability': 5.0, 'binding': 8.0, 'name': 'Binding-Optimized', 'color': '#FF6B6B'},
        {'function': 8.0, 'stability': 4.0, 'binding': 3.0, 'name': 'Function-Maximized', 'color': '#00FF88'},
        {'function': 5.0, 'stability': 10.0, 'binding': 5.0, 'name': 'Conservative', 'color': '#9D00FF'},
        {'function': 6.0, 'stability': 6.0, 'binding': 6.0, 'name': 'Experimental', 'color': '#FF9500'},
    ]

    for candidate_idx in range(n_candidates):
        torch.manual_seed(42 + candidate_idx * 17)
        np.random.seed(42 + candidate_idx * 17)

        profile = weight_profiles[candidate_idx % len(weight_profiles)]

        emb = embedder.get_embeddings(cancer_seq).detach().requires_grad_(True)
        emb_wt = embedder.get_embeddings(P53_WT).detach()

        with torch.no_grad():
            perturbation = torch.randn_like(emb) * 0.05 * exploration_diversity
            emb.data += perturbation

        optimizer = torch.optim.Adam([emb], lr=0.04)
        locked_indices = [int(p) - 1 for p in locked_positions if p]

        n_steps = 100  # Faster for live viz
        trajectory = []  # Store optimization trajectory
        best_valid_state = None
        best_valid_score = -float('inf')

        # Pre-compute WT AA indices
        wt_aa_indices = []
        for aa in P53_WT:
            aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
            if aa_id in AA_IDS:
                wt_aa_indices.append(AA_IDS.index(aa_id))
            else:
                wt_aa_indices.append(0)
        wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb.device)

        for step_idx in range(1, n_steps + 1):
            optimizer.zero_grad()

            z, logits, probs = embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            if pooled.shape[-1] != 320:
                pooled = pooled[:, :320]

            score = oracle.model(pooled)

            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()
            dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
            hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)

            # Loss
            loss = -score * profile['function']
            loss -= profile['stability'] * stability
            loss -= profile['binding'] * dna_force
            loss -= 3.0 * hydro_packing

            probs_aa = F.softmax(logits_aa[0], dim=-1)
            wt_probs = probs_aa[torch.arange(len(P53_WT), device=emb.device), wt_aa_tensor]
            mutation_prob = 1.0 - wt_probs
            expected_mutations = mutation_prob.sum()

            with torch.no_grad():
                decoded_ids = torch.argmax(probs_aa, dim=-1)
                n_mutations = (decoded_ids != wt_aa_tensor).sum().item()
                seq_identity = 100.0 * (1.0 - n_mutations / len(P53_WT))

            max_mutations = int(len(P53_WT) * (100 - min_identity) / 100)
            loss += 50.0 * F.relu(expected_mutations - max_mutations)

            if seq_identity < min_identity - 5:
                loss += 500.0 * (min_identity - 5 - seq_identity)
            if stability.item() < min_stability:
                loss += 100.0 * (min_stability - stability)
            if locked_indices:
                loss += 500.0 * F.mse_loss(emb[:, locked_indices, :], emb_wt[:, locked_indices, :])

            dist_l1 = torch.norm(emb - emb_wt, p=1) / emb.numel()
            loss += 40.0 * dist_l1

            loss.backward()
            optimizer.step()

            # Record trajectory every 5 steps for visualization
            if step_idx % 5 == 0 or step_idx == n_steps:
                with torch.no_grad():
                    top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                    top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                    tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                    current_seq = "".join(tokens)[:len(P53_WT)]
                    muts = [f"{P53_WT[j]}{j+1}{current_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != current_seq[j]]

                    # Get mutation positions for visualization
                    mut_positions = [int(''.join(filter(str.isdigit, m))) for m in muts if m]

                    trajectory.append({
                        'step': step_idx,
                        'score': score.item(),
                        'stability': stability.item(),
                        'binding': dna_force.item(),
                        'identity': seq_identity,
                        'n_mutations': n_mutations,
                        'mutations': muts[:5],
                        'mut_positions': mut_positions,
                        'sequence': current_seq,
                        'lx': pooled[0, 0].item(),
                        'ly': pooled[0, 1].item()
                    })

                    # Track best valid
                    if seq_identity >= min_identity and stability.item() >= min_stability:
                        if score.item() > best_valid_score:
                            best_valid_score = score.item()
                            best_valid_state = trajectory[-1].copy()

                    # YIELD for live visualization
                    if progress_callback:
                        progress_callback({
                            'candidate_idx': candidate_idx,
                            'candidate_total': n_candidates,
                            'profile': profile,
                            'step': step_idx,
                            'total_steps': n_steps,
                            'current_state': trajectory[-1],
                            'trajectory': trajectory
                        })

        # Final candidate
        final_state = best_valid_state if best_valid_state else trajectory[-1]
        candidate = {
            'candidate_id': candidate_idx + 1,
            'profile': profile['name'],
            'color': profile['color'],
            'sequence': final_state['sequence'],
            'score': final_state['score'],
            'stability': final_state['stability'],
            'binding': final_state['binding'],
            'identity': final_state['identity'],
            'n_mutations': final_state['n_mutations'],
            'mutations': final_state['mutations'],
            'mut_positions': final_state.get('mut_positions', []),
            'trajectory': trajectory,
            'meets_constraints': final_state['identity'] >= min_identity and final_state['stability'] >= min_stability
        }
        all_candidates.append(candidate)

    return all_candidates


# === GENERATIVE DESIGN ENGINE ===
def run_generative_design(constraints: dict, n_candidates: int = 6):
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

    # Adjust identity based on delivery method
    if delivery_method == 'protein_therapy':
        min_identity = max(min_identity, 95.0)  # Stricter for direct protein
    elif delivery_method == 'mrna_therapy':
        min_identity = max(min_identity, 92.0)

    all_candidates = []

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

        emb = embedder.get_embeddings(cancer_seq).detach().requires_grad_(True)
        emb_wt = embedder.get_embeddings(P53_WT).detach()

        # Add initial perturbation for diversity
        with torch.no_grad():
            perturbation = torch.randn_like(emb) * 0.05 * exploration_diversity
            emb.data += perturbation

        optimizer = torch.optim.Adam([emb], lr=0.04)
        locked_indices = [int(p) - 1 for p in locked_positions if p]

        # Shorter optimization for multiple candidates
        n_steps = 150
        best_valid_state = None
        best_valid_score = -float('inf')

        for step_idx in range(1, n_steps + 1):
            optimizer.zero_grad()

            z, logits, probs = embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            if pooled.shape[-1] != 320:
                pooled = pooled[:, :320]

            score = oracle.model(pooled)

            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()

            dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
            hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)

            # LOSS with profile-specific weights
            loss = -score * profile['function']
            loss -= profile['stability'] * stability
            loss -= profile['binding'] * dna_force
            loss -= 3.0 * hydro_packing

            # CONSTRAINT ENFORCEMENT (hard constraints from user)
            wt_aa_indices = []
            for aa in P53_WT:
                aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
                if aa_id in AA_IDS:
                    wt_aa_indices.append(AA_IDS.index(aa_id))
                else:
                    wt_aa_indices.append(0)
            wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb.device)

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
                if seq_identity >= min_identity and stability.item() >= min_stability:
                    if score.item() > best_valid_score:
                        best_valid_score = score.item()
                        top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                        top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                        tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                        best_valid_state = {
                            'sequence': "".join(tokens)[:len(P53_WT)],
                            'score': score.item(),
                            'stability': stability.item(),
                            'binding': dna_force.item(),
                            'identity': seq_identity,
                            'n_mutations': n_mutations
                        }

        # Final decode
        with torch.no_grad():
            z, logits, probs = embedder.latent_forward_ascent(emb)
            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability = log_probs.max(dim=-1).values.mean()
            score = oracle.model(z.mean(dim=1))
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

        # Use best valid state if final doesn't meet constraints
        if best_valid_state and (seq_identity < min_identity or stability.item() < min_stability):
            candidate = {
                'candidate_id': candidate_idx + 1,
                'profile': profile['name'],
                'sequence': best_valid_state['sequence'],
                'score': best_valid_state['score'],
                'stability': best_valid_state['stability'],
                'binding': best_valid_state['binding'],
                'identity': best_valid_state['identity'],
                'n_mutations': best_valid_state['n_mutations'],
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
                'score': score.item(),
                'stability': stability.item(),
                'binding': dna_force.item(),
                'identity': seq_identity,
                'n_mutations': n_mutations,
                'mutations': muts,
                'meets_constraints': seq_identity >= min_identity and stability.item() >= min_stability
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
    st.markdown('<p class="gen-design-header"> proteoMgCAD Studio</p>', unsafe_allow_html=True)
    st.markdown("*Define constraints. AI generates optimal second-site rescues. Topology optimization for proteins.*")
    st.markdown("---")

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

    with gen_col2:
        st.info(f"""
        **Generation Preview:**
        - Target: **{gd_target}** rescue
        - Identity: ≥{gd_identity:.0f}% (max {int(len(P53_WT) * (100 - gd_identity) / 100)} mutations)
        - Locked: {len(gd_locked)} positions
        - Delivery: {gd_delivery}
        - Candidates: {gd_n_candidates} diverse solutions
        """)

    # === VISUALIZATION MODE ===
    gd_live_mode = st.toggle(" Watch Live Optimization", value=True,
                              help="Like watching CAD topology optimization - see proteins being built in real-time")

    # === GENERATE BUTTON ===
    if st.button(" GENERATE CANDIDATE DESIGNS", type="primary", use_container_width=True):
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
            'diversity': gd_diversity
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
                <div id="init_viewer" style="width:100%; height:400px; border-radius:12px; border:2px solid #6366F1; background:#0F172A;"></div>
                <script>
                    let viewer = $3Dmol.createViewer('init_viewer', {{backgroundColor: '0x0F172A'}});
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

            # Run generation with live updates
            for candidate_idx in range(gd_n_candidates):
                torch.manual_seed(42 + candidate_idx * 17)
                np.random.seed(42 + candidate_idx * 17)

                weight_profiles = [
                    {'function': 4.0, 'stability': 8.0, 'binding': 2.5, 'name': 'Balanced', 'color': '#00D4FF'},
                    {'function': 2.0, 'stability': 15.0, 'binding': 2.0, 'name': 'Stability-First', 'color': '#FFD700'},
                    {'function': 3.0, 'stability': 5.0, 'binding': 8.0, 'name': 'Binding-Optimized', 'color': '#FF6B6B'},
                    {'function': 8.0, 'stability': 4.0, 'binding': 3.0, 'name': 'Function-Maximized', 'color': '#00FF88'},
                    {'function': 5.0, 'stability': 10.0, 'binding': 5.0, 'name': 'Conservative', 'color': '#9D00FF'},
                    {'function': 6.0, 'stability': 6.0, 'binding': 6.0, 'name': 'Experimental', 'color': '#FF9500'},
                ]
                profile = weight_profiles[candidate_idx % len(weight_profiles)]

                status_placeholder.markdown(f"""
                ###  Building Candidate {candidate_idx + 1}/{gd_n_candidates}
                **Strategy:** {profile['name']}
                """)

                cancer_seq = apply_mutation(P53_WT, gd_target)
                if cancer_seq is None:
                    cancer_seq = P53_WT

                AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

                emb = embedder.get_embeddings(cancer_seq).detach().requires_grad_(True)
                emb_wt = embedder.get_embeddings(P53_WT).detach()

                with torch.no_grad():
                    perturbation = torch.randn_like(emb) * 0.05 * gd_diversity
                    emb.data += perturbation

                optimizer = torch.optim.Adam([emb], lr=0.04)
                locked_indices = [int(p) - 1 for p in gd_locked if p]

                wt_aa_indices = []
                for aa in P53_WT:
                    aa_id = embedder.tokenizer.convert_tokens_to_ids(aa)
                    if aa_id in AA_IDS:
                        wt_aa_indices.append(AA_IDS.index(aa_id))
                    else:
                        wt_aa_indices.append(0)
                wt_aa_tensor = torch.tensor(wt_aa_indices, device=emb.device)

                n_steps = 80  # Faster for live viz
                trajectory = []
                best_valid_state = None
                best_valid_score = -float('inf')

                for step_idx in range(1, n_steps + 1):
                    optimizer.zero_grad()

                    z, logits, probs = embedder.latent_forward_ascent(emb)
                    pooled = z.mean(dim=1)
                    if pooled.shape[-1] != 320:
                        pooled = pooled[:, :320]

                    score = oracle.model(pooled)

                    logits_aa = logits[:, :, AA_IDS]
                    log_probs = F.log_softmax(logits_aa, dim=-1)
                    stability = log_probs.max(dim=-1).values.mean()
                    dna_force = embedder.get_dna_contact_prob(z, logits, probs=probs)
                    hydro_packing = embedder.get_hydrophobic_packing(logits, probs=probs)

                    loss = -score * profile['function']
                    loss -= profile['stability'] * stability
                    loss -= profile['binding'] * dna_force
                    loss -= 3.0 * hydro_packing

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
                    if locked_indices:
                        loss += 500.0 * F.mse_loss(emb[:, locked_indices, :], emb_wt[:, locked_indices, :])

                    dist_l1 = torch.norm(emb - emb_wt, p=1) / emb.numel()
                    loss += 40.0 * dist_l1

                    loss.backward()
                    optimizer.step()

                    # === LIVE UPDATE every 5 steps ===
                    if step_idx % 5 == 0 or step_idx == n_steps:
                        with torch.no_grad():
                            top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
                            top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
                            tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
                            current_seq = "".join(tokens)[:len(P53_WT)]
                            muts = [f"{P53_WT[j]}{j+1}{current_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != current_seq[j]]
                            mut_positions = [int(''.join(filter(str.isdigit, m))) for m in muts if m]

                            trajectory.append({
                                'step': step_idx, 'score': score.item(), 'stability': stability.item(),
                                'binding': dna_force.item(), 'identity': seq_identity, 'n_mutations': n_mutations,
                                'mutations': muts[:5], 'mut_positions': mut_positions, 'sequence': current_seq,
                                'lx': pooled[0, 0].item(), 'ly': pooled[0, 1].item(),
                                'lz': pooled[0, 2].item() if pooled.shape[-1] > 2 else score.item()  # 3D latent
                            })

                            if seq_identity >= gd_identity and stability.item() >= gd_min_stability:
                                if score.item() > best_valid_score:
                                    best_valid_score = score.item()
                                    best_valid_state = trajectory[-1].copy()

                            # Update live metrics display
                            with metrics_placeholder.container():
                                st.markdown(f"**Step {step_idx}/{n_steps}**")
                                mc1, mc2 = st.columns(2)
                                mc1.metric("Score", f"{score.item():.3f}")
                                mc2.metric("Identity", f"{seq_identity:.1f}%")
                                mc3, mc4 = st.columns(2)
                                mc3.metric("Stability", f"{stability.item():.3f}")
                                mc4.metric("Mutations", n_mutations)

                            # Update mutations list
                            with mutations_placeholder.container():
                                st.markdown("**Current Mutations:**")
                                for m in muts[:4]:
                                    st.write(f"• {m}")

                            # Update progress
                            overall_progress = (candidate_idx * n_steps + step_idx) / (gd_n_candidates * n_steps)
                            progress_placeholder.progress(overall_progress, text=f"Overall: {overall_progress*100:.0f}%")

                            # Update trajectory plot (only every 10 steps to reduce flicker)
                            if len(trajectory) > 1 and step_idx % 10 == 0:
                                traj_df = pd.DataFrame(trajectory)
                                fig_traj = go.Figure()
                                fig_traj.add_trace(go.Scatter(x=traj_df['step'], y=traj_df['score'],
                                                              mode='lines+markers', name='Score',
                                                              line=dict(color=profile['color'], width=3)))
                                fig_traj.add_trace(go.Scatter(x=traj_df['step'], y=traj_df['stability'],
                                                              mode='lines', name='Stability',
                                                              line=dict(color='#FFD700', dash='dot')))
                                fig_traj.update_layout(template='plotly_white', height=200,
                                                       margin=dict(l=0, r=0, t=30, b=0),
                                                       title=f"Candidate {candidate_idx+1} Optimization",
                                                       showlegend=True)
                                trajectory_placeholder.plotly_chart(fig_traj, use_container_width=True)

                # === 3D STRUCTURE - Update every 20 steps for LIVE MORPHING ===
                if step_idx % 20 == 0 or step_idx == n_steps:
                    if mut_positions:
                        mut_str = "+".join([str(p) for p in mut_positions[:10]])
                    else:
                        mut_str = "none"

                    try:
                        with open("data/raw/p53_wt.pdb", "r") as f:
                            pdb_content = f.read()

                        # MORPH the structure based on mutations - coordinates actually change!
                        morphed_pdb = morph_pdb_for_mutations(
                            pdb_content,
                            mut_positions,
                            step_idx,
                            n_steps
                        )

                        progress_pct = int((step_idx / n_steps) * 100)
                        is_final = step_idx == n_steps

                        # Render structure with morphing
                        live_3d_html = f"""
                        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                        <div id="morph_view_{candidate_idx}_{step_idx}" style="width:100%; height:380px; border-radius:12px; border:3px solid {profile["color"]}; background:#0F172A;"></div>
                        <div style="height:8px; background:linear-gradient(90deg, {profile["color"]} {progress_pct}%, #374151 {progress_pct}%); border-radius:4px; margin-top:8px;"></div>
                        <script>
                            let viewer = $3Dmol.createViewer('morph_view_{candidate_idx}_{step_idx}', {{backgroundColor: '0x0F172A'}});
                            let pdb = `{morphed_pdb.replace(chr(10), chr(92) + 'n').replace('`', '')}`;
                            viewer.addModel(pdb, 'pdb');
                            viewer.setStyle({{}}, {{cartoon: {{color: '0x64748B', opacity: 0.7}}}});
                            let muts = "{mut_str}".split("+");
                            muts.forEach(pos => {{
                                if (pos && pos !== 'none') {{
                                    let p = parseInt(pos);
                                    viewer.addStyle({{resi: p}}, {{
                                        cartoon: {{color: '{profile["color"]}', opacity: 1.0}},
                                        stick: {{color: '{profile["color"]}', radius: 0.3}}
                                    }});
                                    viewer.addStyle({{resi: p, atom: 'CA'}}, {{sphere: {{color: '{profile["color"]}', radius: 0.6}}}});
                                }}
                            }});
                            viewer.addStyle({{resi: [112,113,114,115,116,117,118,119,120,121,122,123,124]}}, {{cartoon: {{color: '0x0EA5E9', opacity: 0.9}}}});
                            viewer.addStyle({{resi: [236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251]}}, {{cartoon: {{color: '0x10B981', opacity: 0.9}}}});
                            viewer.zoomTo();
                            viewer.spin('y', 1.5);
                            viewer.render();
                        </script>
                        <p style="text-align:center; color:{profile['color']}; font-weight:600; margin-top:10px; font-size:15px;">
                            {"✅ " + profile['name'] + " Complete" if is_final else "🔄 " + profile['name'] + " Step " + str(step_idx) + "/" + str(n_steps)} | {len(mut_positions)} mutations
                        </p>
                        """
                        with structure_placeholder.container():
                            components.html(live_3d_html, height=450)
                    except Exception as e:
                        structure_placeholder.error(f"Could not render protein: {e}")

                # Store final candidate
                final_state = best_valid_state if best_valid_state else trajectory[-1]
                candidate = {
                    'candidate_id': candidate_idx + 1,
                    'profile': profile['name'],
                    'color': profile['color'],
                    'sequence': final_state['sequence'],
                    'score': final_state['score'],
                    'stability': final_state['stability'],
                    'binding': final_state['binding'],
                    'identity': final_state['identity'],
                    'n_mutations': final_state['n_mutations'],
                    'mutations': final_state['mutations'],
                    'mut_positions': final_state.get('mut_positions', []),
                    'trajectory': trajectory,
                    'meets_constraints': final_state['identity'] >= gd_identity and final_state['stability'] >= gd_min_stability
                }
                all_candidates.append(candidate)

            candidates = all_candidates
            status_placeholder.success(f" Generated {len(candidates)} candidate designs!")

        else:
            # Non-live mode (faster)
            with st.spinner(f" Generating {gd_n_candidates} candidate designs..."):
                candidates = run_generative_design_live(constraints, n_candidates=gd_n_candidates)

        # Store in session state
        st.session_state['gd_candidates'] = candidates
        st.session_state['gd_constraints'] = constraints

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
        sum_col2.metric("Avg Function Score", f"{avg_score:.3f}")
        sum_col3.metric("Avg Identity", f"{avg_identity:.1f}%")
        sum_col4.metric("Design Space Coverage", f"{len(set(c['profile'] for c in candidates))} profiles")

        # Pareto frontier visualization
        st.markdown("###  Pareto Frontier (Trade-offs)")

        pareto_df = pd.DataFrame(candidates)

        fig_pareto = go.Figure()

        # Color by profile
        colors = {'Balanced': '#00D4FF', 'Stability-First': '#FFD700', 'Binding-Optimized': '#FF6B6B',
                  'Function-Maximized': '#00FF88', 'Conservative': '#9D00FF', 'Experimental': '#FF9500'}

        for profile in pareto_df['profile'].unique():
            df_profile = pareto_df[pareto_df['profile'] == profile]
            fig_pareto.add_trace(go.Scatter(
                x=df_profile['identity'],
                y=df_profile['score'],
                mode='markers+text',
                name=profile,
                text=df_profile['candidate_id'].astype(str),
                textposition='top center',
                marker=dict(
                    size=df_profile['binding'] * 3 + 10,
                    color=colors.get(profile, '#FFFFFF'),
                    line=dict(width=2, color='white'),
                    opacity=0.8
                ),
                hovertemplate=f"<b>{profile}</b><br>" +
                              "Identity: %{x:.1f}%<br>" +
                              "Score: %{y:.3f}<br>" +
                              "Binding: %{marker.size:.1f}<extra></extra>"
            ))

        # Add constraint boundaries
        fig_pareto.add_vline(x=constraints.get('min_identity', 90), line_dash="dash", line_color="red",
                            annotation_text=f"Min Identity: {constraints.get('min_identity', 90)}%")
        fig_pareto.add_hline(y=constraints.get('min_function', -0.2), line_dash="dash", line_color="orange",
                            annotation_text="Min Function")

        fig_pareto.update_layout(
            template="plotly_white",
            title="Design Space Exploration (Size = DNA Binding)",
            xaxis_title="Sequence Identity (%)",
            yaxis_title="Rescue Score",
            height=400,
            legend=dict(orientation="h", yanchor="bottom", y=1.02)
        )
        st.plotly_chart(fig_pareto, use_container_width=True)

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
            sorted_candidates = sorted(candidates, key=lambda x: x['score'], reverse=True)
            n_cols = 3

            for row_start in range(0, min(len(sorted_candidates), 6), n_cols):
                gallery_cols = st.columns(n_cols)

                for i, col in enumerate(gallery_cols):
                    cand_idx = row_start + i
                    if cand_idx < len(sorted_candidates):
                        cand = sorted_candidates[cand_idx]
                        profile_color = cand.get('color', '#00D4FF')

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
                            badge_color = "#00FF88" if cand['meets_constraints'] else "#FF6B6B"
                            badge_text = "" if cand['meets_constraints'] else "️"

                            st.markdown(f"""
                            <div style="text-align:center; padding:5px; background:linear-gradient(135deg, #1e2130 0%, #0e1117 100%);
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
                            m3.metric("Muts", cand['n_mutations'], label_visibility="collapsed")

                            # Show unique mutations
                            if unique_muts:
                                st.caption(f" **Unique:** {', '.join(cand['mutations'][:3] if len(unique_muts) > 0 else ['None'])}")
                            else:
                                st.caption("*No unique mutations*")

            # === MUTATION COMPARISON HEATMAP ===
            st.markdown("###  Mutation Comparison Heatmap")
            st.markdown("*Which positions are mutated by each candidate? Brighter = more candidates mutate this position.*")

            # Build mutation matrix
            all_positions = set()
            for c in sorted_candidates:
                muts = c.get('mut_positions', [])
                if not muts and c.get('mutations'):
                    muts = [int(''.join(filter(str.isdigit, m))) for m in c['mutations'] if any(ch.isdigit() for ch in m)]
                all_positions.update(muts)

            if all_positions:
                positions_sorted = sorted(all_positions)
                matrix = []
                for c in sorted_candidates:
                    row = []
                    c_muts = set(c.get('mut_positions', []))
                    if not c_muts and c.get('mutations'):
                        c_muts = set(int(''.join(filter(str.isdigit, m))) for m in c['mutations'] if any(ch.isdigit() for ch in m))
                    for pos in positions_sorted:
                        row.append(1 if pos in c_muts else 0)
                    matrix.append(row)

                fig_heatmap = go.Figure(data=go.Heatmap(
                    z=matrix,
                    x=[f"Pos {p}" for p in positions_sorted],
                    y=[f"#{c['candidate_id']} {c['profile']}" for c in sorted_candidates],
                    colorscale=[[0, '#0e1117'], [1, '#00D4FF']],
                    showscale=False
                ))
                fig_heatmap.update_layout(
                    template='plotly_white',
                    height=250,
                    margin=dict(l=0, r=0, t=30, b=0),
                    title="Mutation Positions by Candidate",
                    xaxis_title="Position",
                    yaxis_title=""
                )
                st.plotly_chart(fig_heatmap, use_container_width=True)

                # Consensus mutations (positions mutated by multiple candidates)
                position_counts = {}
                for pos in all_positions:
                    count = sum(1 for c in sorted_candidates if pos in set(c.get('mut_positions', [])))
                    position_counts[pos] = count

                consensus = {k: v for k, v in position_counts.items() if v >= 2}
                if consensus:
                    st.success(f"** Consensus positions** (mutated by 2+ candidates): {', '.join(str(p) for p in sorted(consensus.keys()))}")

        # Candidate cards
        st.markdown("### Candidate Details")

        # Sort by score (best first)
        sorted_candidates = sorted(candidates, key=lambda x: x['score'], reverse=True)

        for i in range(0, len(sorted_candidates), 2):
            card_cols = st.columns(2)

            for j, col in enumerate(card_cols):
                if i + j < len(sorted_candidates):
                    cand = sorted_candidates[i + j]

                    # Color based on validity
                    border_color = "#00FF88" if cand['meets_constraints'] else "#FF6B6B"
                    validity_badge = " VALID" if cand['meets_constraints'] else "️ CONSTRAINT VIOLATION"

                    with col:
                        st.markdown(f"""
                        <div class="candidate-card" style="border-left-color: {border_color};">
                            <h4>Candidate #{cand['candidate_id']} - {cand['profile']}</h4>
                            <p><b>{validity_badge}</b></p>
                        </div>
                        """, unsafe_allow_html=True)

                        # Metrics row
                        m1, m2, m3 = st.columns(3)
                        m1.metric("Score", f"{cand['score']:.3f}")
                        m2.metric("Identity", f"{cand['identity']:.1f}%")
                        m3.metric("Mutations", cand['n_mutations'])

                        # Show mutations
                        muts_display = ", ".join(cand['mutations'][:5])
                        if len(cand['mutations']) > 5:
                            muts_display += f" (+{len(cand['mutations']) - 5} more)"
                        st.caption(f"**Mutations:** {muts_display}")

                        # Physics scores
                        st.caption(f"Stability: {cand['stability']:.3f} | Binding: {cand['binding']:.2f}")

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
            'Stability': f"{c['stability']:.3f}",
            'Binding': f"{c['binding']:.2f}",
            'Identity': f"{c['identity']:.1f}%",
            'Mutations': c['n_mutations'],
            'Valid': '' if c['meets_constraints'] else ''
        } for c in sorted_candidates])

        st.dataframe(
            compare_df,
            hide_index=True,
            use_container_width=True,
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
            use_container_width=True
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
            st.metric("🎯 Rescue Score", f"{best_candidate['score']:.3f}",
                     delta="Best" if best_candidate['meets_constraints'] else "Constraint Issue")
        with hero_col2:
            st.metric("🧬 Identity", f"{best_candidate['identity']:.1f}%",
                     delta=f"{int(best_candidate['n_mutations'])} mutations")
        with hero_col3:
            st.metric("⚡ Stability", f"{best_candidate['stability']:.3f}")
        with hero_col4:
            st.metric("🔗 DNA Binding", f"{best_candidate['binding']:.2f}")

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

                fig_result = go.Figure()
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['score'],
                    mode='lines+markers', name='Rescue Score',
                    line=dict(color=best_candidate.get('color', '#00D4FF'), width=3),
                    marker=dict(size=6)
                ))
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['stability'],
                    mode='lines', name='Stability',
                    line=dict(color='#FFD700', width=2, dash='dot')
                ))
                fig_result.add_trace(go.Scatter(
                    x=traj_df['step'], y=traj_df['binding'] / 10,  # Scale for visibility
                    mode='lines', name='DNA Binding (scaled)',
                    line=dict(color='#00FF88', width=2, dash='dash')
                ))

                fig_result.update_layout(
                    template='plotly_white',
                    height=250,
                    margin=dict(l=0, r=0, t=30, b=0),
                    title=f"Candidate #{best_candidate['candidate_id']} Evolution",
                    xaxis_title="Optimization Step",
                    yaxis_title="Score",
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5)
                )
                st.plotly_chart(fig_result, use_container_width=True)
            else:
                st.info("No trajectory data available")

        # Protein sequence display
        st.markdown("### 🧬 Rescued Protein Sequence")

        seq = best_candidate['sequence']

        # Create a styled sequence view highlighting mutations
        if best_candidate.get('mut_positions'):
            mut_positions = set(best_candidate['mut_positions'])
            # Build colored sequence HTML
            seq_html = '<div style="font-family: monospace; font-size: 12px; line-height: 1.8; background: #1e2130; padding: 15px; border-radius: 8px; overflow-x: auto;">'
            for i, aa in enumerate(seq):
                pos = i + 1  # 1-indexed
                if pos in mut_positions:
                    seq_html += f'<span style="color: #FF6B6B; font-weight: bold; background: rgba(255,107,107,0.2); padding: 2px 1px; border-radius: 3px;">{aa}</span>'
                else:
                    seq_html += f'<span style="color: #94A3B8;">{aa}</span>'
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
                use_container_width=True
            )

        # === ADVANCED VISUALIZATIONS ===
        st.markdown("---")
        st.markdown("## 📊 Advanced Analysis")

        if best_candidate.get('trajectory') and len(best_candidate['trajectory']) > 2:
            traj_df = pd.DataFrame(best_candidate['trajectory'])

            viz_col1, viz_col2 = st.columns(2)

            with viz_col1:
                # === 3D LATENT SPACE PLOT ===
                st.markdown("### 🌐 3D Latent Space Trajectory")
                st.caption("*Watch how the protein moves through ESM-2 embedding space during optimization*")

                # Create 3D scatter with trajectory line
                fig_3d = go.Figure()

                # Add trajectory line
                fig_3d.add_trace(go.Scatter3d(
                    x=traj_df['lx'], y=traj_df['ly'],
                    z=traj_df.get('lz', traj_df['score']),
                    mode='lines',
                    line=dict(color='rgba(100,100,100,0.5)', width=2),
                    name='Path',
                    showlegend=False
                ))

                # Add points colored by step (time)
                fig_3d.add_trace(go.Scatter3d(
                    x=traj_df['lx'], y=traj_df['ly'],
                    z=traj_df.get('lz', traj_df['score']),
                    mode='markers',
                    marker=dict(
                        size=6,
                        color=traj_df['step'],
                        colorscale='Viridis',
                        colorbar=dict(title='Step', x=1.02),
                        opacity=0.9
                    ),
                    text=[f"Step {s}<br>Score: {sc:.3f}" for s, sc in zip(traj_df['step'], traj_df['score'])],
                    hoverinfo='text',
                    name='Optimization Path'
                ))

                # Mark start (red) and end (green)
                fig_3d.add_trace(go.Scatter3d(
                    x=[traj_df['lx'].iloc[0]], y=[traj_df['ly'].iloc[0]],
                    z=[traj_df.get('lz', traj_df['score']).iloc[0]],
                    mode='markers', marker=dict(size=12, color='#FF6B6B', symbol='diamond'),
                    name='Start (Cancer)'
                ))
                fig_3d.add_trace(go.Scatter3d(
                    x=[traj_df['lx'].iloc[-1]], y=[traj_df['ly'].iloc[-1]],
                    z=[traj_df.get('lz', traj_df['score']).iloc[-1]],
                    mode='markers', marker=dict(size=12, color='#00FF88', symbol='diamond'),
                    name='End (Rescued)'
                ))

                fig_3d.update_layout(
                    template='plotly_white',
                    height=400,
                    margin=dict(l=0, r=0, t=30, b=0),
                    scene=dict(
                        xaxis_title='Latent X',
                        yaxis_title='Latent Y',
                        zaxis_title='Latent Z',
                        camera=dict(eye=dict(x=1.5, y=1.5, z=1.2))
                    ),
                    legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5)
                )
                st.plotly_chart(fig_3d, use_container_width=True)

            with viz_col2:
                # === MULTI-METRIC COMPARISON ===
                st.markdown("### 📈 Multi-Metric Evolution")
                st.caption("*All optimization objectives normalized to [0,1] scale*")

                # Normalize metrics to 0-1 for comparison
                def normalize(series):
                    min_val, max_val = series.min(), series.max()
                    if max_val - min_val < 1e-6:
                        return series * 0 + 0.5
                    return (series - min_val) / (max_val - min_val)

                fig_multi = go.Figure()

                # Score (main objective)
                fig_multi.add_trace(go.Scatter(
                    x=traj_df['step'], y=normalize(traj_df['score']),
                    mode='lines+markers', name='Rescue Score',
                    line=dict(color='#0EA5E9', width=3),
                    marker=dict(size=4)
                ))

                # Stability
                fig_multi.add_trace(go.Scatter(
                    x=traj_df['step'], y=normalize(traj_df['stability']),
                    mode='lines', name='Stability',
                    line=dict(color='#FFD700', width=2)
                ))

                # DNA Binding
                fig_multi.add_trace(go.Scatter(
                    x=traj_df['step'], y=normalize(traj_df['binding']),
                    mode='lines', name='DNA Binding',
                    line=dict(color='#10B981', width=2)
                ))

                # Identity (inverse - lower mutations = higher identity)
                fig_multi.add_trace(go.Scatter(
                    x=traj_df['step'], y=normalize(traj_df['identity']),
                    mode='lines', name='Identity',
                    line=dict(color='#8B5CF6', width=2, dash='dot')
                ))

                fig_multi.update_layout(
                    template='plotly_white',
                    height=400,
                    margin=dict(l=0, r=0, t=30, b=0),
                    xaxis_title='Optimization Step',
                    yaxis_title='Normalized Value (0-1)',
                    yaxis=dict(range=[-0.05, 1.05]),
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5)
                )
                st.plotly_chart(fig_multi, use_container_width=True)

            # === RESIDUE CONTRIBUTION HEATMAP ===
            st.markdown("### 🧬 Position-wise Mutation Impact")
            st.caption("*Heatmap showing which residue positions were modified and their functional regions*")

            # Create position data
            seq_len = len(best_candidate['sequence'])
            mut_positions = set(best_candidate.get('mut_positions', []))

            # Define functional regions
            DNA_BINDING_L1 = set(range(112, 125))  # L1 loop
            DNA_BINDING_L2 = set(range(164, 175))  # L2 loop
            DNA_BINDING_L3 = set(range(236, 252))  # L3 loop
            ZINC_BINDING = {176, 179, 238, 242}    # Zinc coordination
            DIMER_INTERFACE = set(range(325, 355)) # Tetramerization

            # Build heatmap data (chunked by 50 residues for visibility)
            chunk_size = 50
            n_chunks = (seq_len + chunk_size - 1) // chunk_size

            heatmap_data = []
            annotations = []

            for chunk_idx in range(min(n_chunks, 8)):  # Max 8 rows
                start = chunk_idx * chunk_size + 1
                end = min((chunk_idx + 1) * chunk_size, seq_len)
                row = []
                for pos in range(start, end + 1):
                    if pos in mut_positions:
                        value = 3  # Mutation
                    elif pos in DNA_BINDING_L1 or pos in DNA_BINDING_L2 or pos in DNA_BINDING_L3:
                        value = 2  # DNA binding region
                    elif pos in ZINC_BINDING:
                        value = 2.5  # Zinc binding
                    elif pos in DIMER_INTERFACE:
                        value = 1  # Dimer interface
                    else:
                        value = 0  # Normal
                    row.append(value)
                # Pad if needed
                while len(row) < chunk_size:
                    row.append(-1)  # Empty
                heatmap_data.append(row)

            fig_heatmap_pos = go.Figure(data=go.Heatmap(
                z=heatmap_data,
                x=[str(i) for i in range(1, chunk_size + 1)],
                y=[f"{i*chunk_size+1}-{min((i+1)*chunk_size, seq_len)}" for i in range(len(heatmap_data))],
                colorscale=[
                    [0.0, '#1E293B'],      # Normal (dark)
                    [0.25, '#3B82F6'],     # Dimer interface (blue)
                    [0.5, '#10B981'],      # DNA binding (green)
                    [0.75, '#F59E0B'],     # Zinc (orange)
                    [1.0, '#EF4444']       # Mutation (red)
                ],
                showscale=False,
                hovertemplate='Position: %{x}<br>Region: %{y}<br>Value: %{z}<extra></extra>'
            ))

            fig_heatmap_pos.update_layout(
                template='plotly_white',
                height=250,
                margin=dict(l=0, r=0, t=30, b=0),
                xaxis_title='Position (within chunk)',
                yaxis_title='Residue Range'
            )
            st.plotly_chart(fig_heatmap_pos, use_container_width=True)

            # Legend for heatmap
            st.markdown("""
            <div style="display: flex; gap: 20px; justify-content: center; flex-wrap: wrap; margin-top: -10px;">
                <span style="color: #EF4444;">🔴 Mutation</span>
                <span style="color: #10B981;">🟢 DNA Binding Loop</span>
                <span style="color: #F59E0B;">🟠 Zinc Site</span>
                <span style="color: #3B82F6;">🔵 Dimer Interface</span>
                <span style="color: #64748B;">⬛ Other</span>
            </div>
            """, unsafe_allow_html=True)

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
    # === UNIFIED ANALYSIS TAB ===
    st.markdown("""
    <div class="lovable-card" style="margin-bottom: 1.5rem;">
        <h2 style="margin: 0 0 0.5rem 0; font-size: 1.5rem;">📊 Analysis Dashboard</h2>
        <p style="color: #6B7280; margin: 0;">Validate designs, explore mechanisms, and assess clinical potential</p>
    </div>
    """, unsafe_allow_html=True)

    # Sub-section selector
    analysis_section = st.radio(
        "Select Analysis Module",
        ["Validation", "Explainability", "Drug Generator", "Clinical Impact"],
        horizontal=True,
        label_visibility="collapsed"
    )

    st.markdown("<div style='height: 1rem;'></div>", unsafe_allow_html=True)

    # === VALIDATION SECTION ===
    if analysis_section == "Validation":
        st.markdown("### 🔬 Validation Dashboard")
        st.markdown("*Cross-reference AI designs against experimental data and physics*")

        results = st.session_state.get('results', None)
        target = st.session_state.get('target_mut_saved', target_mut)

        if results is not None and len(results) > 0:
            best = results.loc[results['Score'].idxmax()]
            rescued_seq = best['Sequence']
            mut_summary = best['MutSummary']
            our_muts = [m.strip() for m in mut_summary.split(",") if m.strip()]

            val_col1, val_col2 = st.columns([2, 1])

            with val_col1:
                st.markdown("**Literature Cross-Reference**")
                if target in KNOWN_RESCUES:
                    known = KNOWN_RESCUES[target]
                    st.success(f"Published rescues for {target}: {', '.join(known['published'])}")
                    overlap = set(our_muts) & set(known['published'])
                    if overlap:
                        st.success(f"✓ Your design includes validated rescue: {overlap}")
                else:
                    st.warning(f"No published data for {target}. Your design is exploratory.")

                st.markdown("**Physics Scoring**")
                physics_df = pd.DataFrame([
                    {"Metric": "Folding ΔΔG", "Value": f"{best.get('Stability', 0) * -2:.2f} kcal/mol"},
                    {"Metric": "DNA Binding", "Value": f"{best.get('DNARecruitment', 0):.2f}"},
                    {"Metric": "Hydro Packing", "Value": f"{best.get('HydroPacking', 0):.2f}"},
                ])
                st.dataframe(physics_df, hide_index=True, use_container_width=True)

            with val_col2:
                identity_score = min(best['Identity'] / 95.0, 1.0) * 25
                function_score = max(0, (best['Score'] + 0.5) * 25)
                stability_score_val = max(0, (best['Stability'] + 0.3) * 25)
                literature_score = 25 if target in KNOWN_RESCUES else 10
                total_confidence = min(identity_score + function_score + stability_score_val + literature_score, 100)

                st.metric("Confidence", f"{total_confidence:.0f}/100",
                         delta="Validated" if total_confidence > 70 else "Needs Testing")

                st.progress(identity_score / 25, text=f"Identity: {identity_score:.0f}/25")
                st.progress(function_score / 25, text=f"Function: {function_score:.0f}/25")
                st.progress(stability_score_val / 25, text=f"Stability: {stability_score_val:.0f}/25")

            # Summary
            if total_confidence >= 70:
                st.success(f"**High confidence design** - ready for experimental validation")
            elif total_confidence >= 50:
                st.warning(f"**Moderate confidence** - consider optimizing further")
            else:
                st.error(f"**Low confidence** - needs improvement before testing")

        else:
            st.info("💡 Run a design in the **Design** or **Studio** tab first, then return here to validate.")

    # === EXPLAINABILITY SECTION ===
    elif analysis_section == "Explainability":
        if not ENHANCED_MODULES_AVAILABLE:
            st.warning("Enhanced modules not available. Please ensure all dependencies are installed.")
        else:
            st.markdown("""
            <div class="lovable-card-accent">
                <h3 style="margin: 0 0 0.25rem 0; font-size: 1.1rem;">Rescue Mechanism Explainability</h3>
                <p style="color: #6B7280; margin: 0; font-size: 0.85rem;">Understand WHY your rescue mutations work through attention attribution, counterfactuals, and energy decomposition</p>
            </div>
            """, unsafe_allow_html=True)

            col_exp1, col_exp2 = st.columns([1, 2])

            with col_exp1:
                st.markdown("**Configure Analysis**")

                exp_cancer_mut = st.text_input("Cancer Mutation", value="R175H", key="exp_cancer")
                # Default to well-studied rescue mutations from literature
                exp_rescue_muts = st.text_input("Rescue Mutations (comma-separated)", value="N239Y, T284R", key="exp_rescue")

                if st.button("Analyze Rescue Mechanism", key="btn_explain"):
                    rescue_list = [m.strip() for m in exp_rescue_muts.split(",")]

                    with st.spinner("Analyzing rescue mechanism..."):
                        # Apply mutations to get sequences
                        wt_seq = P53_WT

                        # Apply cancer mutation
                        cancer_pos = int(exp_cancer_mut[1:-1]) - 1
                        cancer_aa = exp_cancer_mut[-1]
                        mutant_seq = wt_seq[:cancer_pos] + cancer_aa + wt_seq[cancer_pos+1:]

                        # Apply rescue mutations
                        rescue_seq = mutant_seq
                        for mut in rescue_list:
                            pos = int(mut[1:-1]) - 1
                            new_aa = mut[-1]
                            rescue_seq = rescue_seq[:pos] + new_aa + rescue_seq[pos+1:]

                        # Generate explanation
                        engine = ExplainabilityEngine()
                        explanation = engine.explain_rescue(
                            wt_seq, mutant_seq, rescue_seq,
                            exp_cancer_mut, rescue_list
                        )

                        st.session_state['explanation'] = explanation
                        st.success("Analysis complete!")

            with col_exp2:
                if 'explanation' in st.session_state:
                    exp = st.session_state['explanation']

                    # Summary
                    st.markdown("### Summary")
                    st.info(exp['summary'])

                    # Mechanism
                    st.markdown("### Rescue Mechanism")
                    mech = exp['mechanism']

                    col_m1, col_m2, col_m3, col_m4 = st.columns(4)
                    col_m1.metric("Primary", mech['primary'].replace('_', ' ').title())
                    col_m2.metric("Confidence", f"{mech['confidence']:.0%}")
                    col_m3.metric("Stability", f"{mech['contributions']['stability']:+.2f}")
                    col_m4.metric("Binding", f"{mech['contributions']['binding']:+.2f}")

                    # Evidence
                    st.markdown("**Evidence:**")
                    for ev in mech['evidence'][:4]:
                        st.markdown(f"- {ev}")

                    # Energy decomposition
                    st.markdown("### Energy Decomposition")
                    energy = exp['energy']

                    energy_df = pd.DataFrame({
                        'Component': ['Electrostatic', 'Van der Waals', 'H-Bonds', 'Solvation', 'Backbone', 'Packing'],
                        'Contribution (kcal/mol)': [
                            energy['electrostatic'],
                            energy['van_der_waals'],
                            energy['hydrogen_bonds'],
                            energy['solvation'],
                            energy['backbone_strain'],
                            energy['sidechain_packing']
                        ]
                    })

                    fig_energy = px.bar(energy_df, x='Component', y='Contribution (kcal/mol)',
                                       color='Contribution (kcal/mol)',
                                       color_continuous_scale='RdBu_r',
                                       template='plotly_white')
                    fig_energy.add_hline(y=0, line_dash='dash', line_color='gray')
                    st.plotly_chart(fig_energy, use_container_width=True)

                    st.metric("Total ddG", f"{energy['total_ddg']:.2f} kcal/mol",
                             delta=energy['rescue_quality'])

                    # Risk factors
                    if mech['risk_factors']:
                        st.markdown("### Risk Factors")
                        for risk in mech['risk_factors']:
                            st.warning(risk)

    # === DRUG GENERATOR SECTION ===
    elif analysis_section == "Drug Generator":
            st.subheader("De Novo Drug Generation")
            st.caption("Generate small molecule drug candidates to stabilize p53 mutants")

            # Show connection to protein rescue if available
            if 'gd_candidates' in st.session_state and st.session_state['gd_candidates']:
                best_protein = sorted(st.session_state['gd_candidates'], key=lambda x: x['score'], reverse=True)[0]
                target = st.session_state.get('gd_constraints', {}).get('target_mutation', 'Unknown')

                st.markdown("""
                <div style="background: linear-gradient(135deg, #10B981 0%, #0EA5E9 100%); padding: 12px 16px;
                            border-radius: 10px; margin-bottom: 1rem;">
                    <p style="color: white; margin: 0; font-size: 0.9rem; font-weight: 600;">
                        🔗 Connected to Best Rescue: Candidate #{} ({}) | Target: {}
                    </p>
                    <p style="color: rgba(255,255,255,0.85); margin: 4px 0 0 0; font-size: 0.8rem;">
                        Drug generation will target binding pockets that complement your rescue mutations
                    </p>
                </div>
                """.format(best_protein['candidate_id'], best_protein['profile'], target), unsafe_allow_html=True)

                # Store for later use
                st.session_state['drug_target_mutation'] = target
                st.session_state['drug_rescue_protein'] = best_protein
            else:
                st.info("💡 **Tip:** Run a protein rescue design in the Studio tab first to enable integrated drug-protein combination therapy design.")

            col_dg1, col_dg2 = st.columns([1, 2])

            with col_dg1:
                st.markdown("**Target Binding Pocket**")

                pocket_names = list(P53_BINDING_POCKETS.keys())
                selected_pocket = st.selectbox("Select Pocket", pocket_names, key="drug_pocket")

                pocket_info = P53_BINDING_POCKETS[selected_pocket]
                st.info(f"**{pocket_info.name}**\n\n{pocket_info.description}")
                st.caption(f"Druggability: {pocket_info.druggability_score:.0%}")

                n_candidates = st.slider("Number of Candidates", 5, 30, 15, key="n_drug_cand")
                gen_method = st.selectbox("Generation Method", ["template", "denovo"], key="drug_method")

                if st.button("Generate Drug Candidates", key="btn_gen_drugs"):
                    with st.spinner("Generating drug candidates..."):
                        engine = DrugGeneratorEngine()
                        candidates = engine.generate_for_pocket(selected_pocket, n_candidates, gen_method)
                        st.session_state['drug_candidates'] = candidates
                        st.session_state['drug_pocket'] = selected_pocket  # Store for combo therapy display
                        st.success(f"Generated {len(candidates)} candidates!")

            with col_dg2:
                if 'drug_candidates' in st.session_state:
                    candidates = st.session_state['drug_candidates']

                    # Summary table
                    st.markdown("### Generated Candidates")

                    drug_data = []
                    for c in candidates[:10]:
                        drug_data.append({
                            'Name': c.name,
                            'Affinity (kcal/mol)': f"{c.binding_affinity:.2f}",
                            'Drug-likeness': f"{c.drug_likeness:.2f}",
                            'MW': f"{c.molecular_weight:.0f}",
                            'LogP': f"{c.logp:.1f}",
                            'Lipinski': "Pass" if c.passes_lipinski() else "Fail"
                        })

                    st.dataframe(pd.DataFrame(drug_data), use_container_width=True)

                    # Property distribution
                    st.markdown("### Property Distribution")

                    props_df = pd.DataFrame([{
                        'MW': c.molecular_weight,
                        'LogP': c.logp,
                        'Drug-likeness': c.drug_likeness,
                        'Affinity': -c.binding_affinity
                    } for c in candidates])

                    fig_props = px.scatter(props_df, x='MW', y='LogP',
                                          color='Drug-likeness', size='Affinity',
                                          template='plotly_white',
                                          title='Property Space')

                    # Add Lipinski rule of 5 boundaries
                    fig_props.add_hline(y=5, line_dash='dash', line_color='red', annotation_text='LogP=5')
                    fig_props.add_vline(x=500, line_dash='dash', line_color='red', annotation_text='MW=500')

                    st.plotly_chart(fig_props, use_container_width=True)

                    # Top candidate details with 3D visualization
                    st.markdown("### 🧪 Top Drug Candidate - 3D Structure")
                    top = candidates[0]

                    drug_col1, drug_col2 = st.columns([2, 1])

                    with drug_col1:
                        # 3D Molecule viewer using 3Dmol.js with SMILES
                        drug_3d_html = f"""
                        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                        <script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
                        <style>
                            #drug_viewer {{
                                width: 100%;
                                height: 350px;
                                border-radius: 12px;
                                border: 2px solid #10B981;
                                background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%);
                            }}
                            .drug-label {{
                                text-align: center;
                                color: #10B981;
                                font-weight: 600;
                                margin-top: 8px;
                                font-size: 14px;
                            }}
                        </style>
                        <div id="drug_viewer"></div>
                        <p class="drug-label">💊 {top.name} | Affinity: {top.binding_affinity:.2f} kcal/mol</p>
                        <script>
                            (async function() {{
                                let viewer = $3Dmol.createViewer('drug_viewer', {{
                                    backgroundColor: '0x0F172A'
                                }});

                                // Use PubChem to get 3D structure from SMILES
                                let smiles = "{top.smiles}";
                                try {{
                                    // Fetch 3D conformer from PubChem
                                    let response = await fetch(
                                        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' +
                                        encodeURIComponent(smiles) + '/SDF?record_type=3d'
                                    );
                                    if (response.ok) {{
                                        let sdf = await response.text();
                                        viewer.addModel(sdf, 'sdf');
                                    }} else {{
                                        // Fallback: 2D from SMILES
                                        viewer.addModel(smiles, 'smiles');
                                    }}
                                }} catch(e) {{
                                    viewer.addModel(smiles, 'smiles');
                                }}

                                // Style: stick with element colors
                                viewer.setStyle({{}}, {{
                                    stick: {{radius: 0.15, colorscheme: 'Jmol'}},
                                    sphere: {{radius: 0.4, colorscheme: 'Jmol'}}
                                }});

                                viewer.zoomTo();
                                viewer.spin('y', 0.8);
                                viewer.render();
                            }})();
                        </script>
                        """
                        import streamlit.components.v1 as components
                        components.html(drug_3d_html, height=420)

                    with drug_col2:
                        st.markdown("**Drug Properties**")
                        st.metric("Binding Affinity", f"{top.binding_affinity:.2f} kcal/mol")
                        st.metric("Mechanism", top.mechanism.replace('_', ' ').title())
                        st.metric("Drug-likeness", f"{top.drug_likeness:.2f}")
                        st.metric("Synthetic Access.", f"{top.synthetic_accessibility:.1f}/10")
                        st.metric("Mol. Weight", f"{top.molecular_weight:.0f} Da")
                        st.metric("LogP", f"{top.logp:.1f}")

                        lipinski = "✅ Pass" if top.passes_lipinski() else "❌ Fail"
                        st.markdown(f"**Lipinski Rule of 5:** {lipinski}")

                        st.code(f"{top.smiles}", language=None)
                        st.caption("SMILES notation")

                    # === COMBINATION THERAPY SUMMARY ===
                    if 'drug_rescue_protein' in st.session_state:
                        st.markdown("---")
                        st.markdown("### 🧬💊 Combination Therapy Summary")
                        st.markdown("*Protein rescue + Small molecule stabilization*")

                        protein = st.session_state['drug_rescue_protein']
                        target = st.session_state.get('drug_target_mutation', 'Unknown')

                        combo_col1, combo_col2 = st.columns(2)

                        with combo_col1:
                            st.markdown("""
                            <div style="background: #1E293B; padding: 15px; border-radius: 10px;
                                        border-left: 4px solid #0EA5E9;">
                                <p style="color: #0EA5E9; font-weight: 600; margin: 0 0 8px 0;">🧬 Rescue Protein</p>
                                <p style="color: #E2E8F0; font-size: 0.85rem; margin: 4px 0;">
                                    <b>Target:</b> {}<br>
                                    <b>Strategy:</b> {}<br>
                                    <b>Mutations:</b> {}<br>
                                    <b>Score:</b> {:.3f}
                                </p>
                            </div>
                            """.format(
                                target,
                                protein['profile'],
                                ', '.join(protein['mutations'][:5]) + ('...' if len(protein['mutations']) > 5 else ''),
                                protein['score']
                            ), unsafe_allow_html=True)

                        with combo_col2:
                            st.markdown("""
                            <div style="background: #1E293B; padding: 15px; border-radius: 10px;
                                        border-left: 4px solid #10B981;">
                                <p style="color: #10B981; font-weight: 600; margin: 0 0 8px 0;">💊 Stabilizer Drug</p>
                                <p style="color: #E2E8F0; font-size: 0.85rem; margin: 4px 0;">
                                    <b>Name:</b> {}<br>
                                    <b>Pocket:</b> {}<br>
                                    <b>Affinity:</b> {:.2f} kcal/mol<br>
                                    <b>Mechanism:</b> {}
                                </p>
                            </div>
                            """.format(
                                top.name,
                                st.session_state.get('drug_pocket', 'Y220C Cavity'),
                                top.binding_affinity,
                                top.mechanism.replace('_', ' ').title()
                            ), unsafe_allow_html=True)

                        # Synergy explanation
                        st.info("""
                        **Therapeutic Strategy:**
                        1. **Rescue mutations** restore DNA binding and transcriptional activity
                        2. **Small molecule** provides additional thermodynamic stabilization
                        3. Combined approach increases therapeutic window and reduces required dose

                        This combination therapy approach mirrors successful strategies used in CFTR
                        rescue for cystic fibrosis (correctors + potentiators).
                        """)

    # === CLINICAL IMPACT SECTION ===
    elif analysis_section == "Clinical Impact":
        if not ENHANCED_MODULES_AVAILABLE:
            st.warning("Enhanced modules not available. Please ensure all dependencies are installed.")
        else:
            st.markdown("### 🏥 Clinical Impact Assessment")
            st.markdown("*Quantify the therapeutic potential and clinical viability of your rescue designs*")

            col_ci1, col_ci2 = st.columns([1, 2])

            with col_ci1:
                st.markdown("**Assessment Configuration**")

                ci_cancer = st.text_input("Cancer Mutation", value="R175H", key="ci_cancer")
                ci_rescue = st.text_input("Rescue Mutations (comma-separated)", value="N239Y, T284R", key="ci_rescue")

                if st.button("Assess Clinical Impact", key="btn_clinical"):
                    rescue_list = [m.strip() for m in ci_rescue.split(",")]

                    with st.spinner("Assessing clinical impact..."):
                        rescue_seq = P53_WT

                        # Apply cancer mutation
                        cancer_pos = int(ci_cancer[1:-1]) - 1
                        cancer_aa = ci_cancer[-1]
                        rescue_seq = rescue_seq[:cancer_pos] + cancer_aa + rescue_seq[cancer_pos+1:]

                        # Apply rescue mutations
                        for mut in rescue_list:
                            pos = int(mut[1:-1]) - 1
                            new_aa = mut[-1]
                            rescue_seq = rescue_seq[:pos] + new_aa + rescue_seq[pos+1:]

                        engine = ClinicalImpactEngine()
                        report = engine.generate_report(
                            name=f"p53_{ci_cancer}_rescue",
                            wt_sequence=P53_WT,
                            rescue_sequence=rescue_seq,
                            cancer_mutation=ci_cancer,
                            rescue_mutations=rescue_list
                        )

                        st.session_state['clinical_report'] = report
                        st.success("Assessment complete!")

            with col_ci2:
                if 'clinical_report' in st.session_state:
                    report = st.session_state['clinical_report']

                    st.markdown("**Clinical Viability**")
                    col_v1, col_v2 = st.columns(2)
                    col_v1.metric("Score", f"{report.overall_clinical_score:.0f}/100")
                    col_v2.metric("Viability", report.clinical_viability.upper())

                    st.markdown("**Patient Population**")
                    pop = report.patient_population
                    col_p1, col_p2 = st.columns(2)
                    col_p1.metric("US Annual", f"{pop.total_patients_per_year:,}")
                    col_p2.metric("Global", f"{pop.global_estimate:,}")

                    st.markdown("**Delivery Feasibility**")
                    delivery_data = []
                    for d in report.delivery_options[:3]:
                        delivery_data.append({
                            'Method': d.method,
                            'Feasibility': f"{d.feasibility_score:.0%}",
                            'Cost': d.estimated_cost_per_dose
                        })
                    st.dataframe(pd.DataFrame(delivery_data), use_container_width=True, hide_index=True)

                    st.markdown("**Recommendation**")
                    st.success(report.recommended_development_path)

                    # Export
                    if st.button("Export Report", key="btn_export_clinical"):
                        export_data = {
                            'name': report.rescue_name,
                            'clinical_score': report.overall_clinical_score,
                            'viability': report.clinical_viability,
                            'recommendation': report.recommended_development_path
                        }
                        st.download_button(
                            "Download JSON",
                            data=json.dumps(export_data, indent=2),
                            file_name=f"{report.rescue_name}_clinical.json",
                            mime="application/json"
                        )

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
