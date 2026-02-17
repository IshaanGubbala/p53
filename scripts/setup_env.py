#!/usr/bin/env python3
"""Cross-platform environment setup for p53CAD.

Run once after cloning the repository:

    python scripts/setup_env.py

Handles:
  1. Symlink fixup  — Windows git clones often turn symlinks into text stubs.
     This script detects broken symlinks and replaces them with real copies.
  2. CUDA PyTorch   — Detects NVIDIA GPUs and reinstalls torch with CUDA
     support if the current install is CPU-only.
  3. Dependency check — Validates that all required (and optional) packages
     are importable and prints a clear status matrix.

Requires only Python stdlib — safe to run before any project deps are installed.
"""

from __future__ import annotations

import importlib.util
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Resolve project root (one level up from scripts/)
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
RECEPTORS_DIR = PROJECT_ROOT / "data" / "raw" / "receptors"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _header(msg: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {msg}")
    print(f"{'=' * 60}")


def _ok(msg: str) -> None:
    print(f"  [OK]   {msg}")


def _warn(msg: str) -> None:
    print(f"  [WARN] {msg}")


def _fail(msg: str) -> None:
    print(f"  [FAIL] {msg}")


def _info(msg: str) -> None:
    print(f"  [INFO] {msg}")


def _has_module(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ModuleNotFoundError, ValueError):
        return False


def _run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kwargs)


# ---------------------------------------------------------------------------
# 1. Fix symlinks (Windows git clone problem)
# ---------------------------------------------------------------------------

def fix_symlinks() -> None:
    _header("Checking symlinks")

    if not RECEPTORS_DIR.exists():
        _info(f"Receptor directory not found ({RECEPTORS_DIR}), skipping.")
        return

    canonical = RECEPTORS_DIR / "p53_wt.pdbqt"
    if not canonical.exists():
        _warn(f"Canonical receptor {canonical} missing — nothing to fix against.")
        return

    fixed = 0
    for p in RECEPTORS_DIR.iterdir():
        if p.name == "p53_wt.pdbqt":
            continue
        if not p.suffix == ".pdbqt":
            continue

        needs_fix = False

        if p.is_symlink():
            # Valid symlink — check it resolves
            if p.resolve().exists():
                _ok(f"{p.name} -> valid symlink")
                continue
            else:
                needs_fix = True
                _warn(f"{p.name} -> broken symlink, replacing with copy")
        elif p.is_file():
            # On Windows, git may create a text file containing the symlink target
            try:
                content = p.read_text(encoding="utf-8").strip()
            except (UnicodeDecodeError, OSError):
                content = ""
            if content in ("p53_wt.pdbqt", "./p53_wt.pdbqt"):
                needs_fix = True
                _warn(f"{p.name} -> text stub (not a real symlink), replacing with copy")
            elif p.stat().st_size == canonical.stat().st_size:
                _ok(f"{p.name} -> already a copy")
                continue
            else:
                _info(f"{p.name} -> unknown state ({p.stat().st_size} bytes), skipping")
                continue

        if needs_fix:
            p.unlink()
            shutil.copy2(canonical, p)
            _ok(f"{p.name} -> replaced with copy of p53_wt.pdbqt")
            fixed += 1

    if fixed:
        _info(f"Fixed {fixed} broken symlink(s)")
    else:
        _ok("All receptor symlinks are valid")


# ---------------------------------------------------------------------------
# 2. Ensure CUDA-enabled PyTorch
# ---------------------------------------------------------------------------

def _detect_nvidia_gpu() -> bool:
    """Check if an NVIDIA GPU is present via nvidia-smi."""
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False
    try:
        result = _run([nvidia_smi, "--query-gpu=name", "--format=csv,noheader"])
        if result.returncode == 0 and result.stdout.strip():
            for line in result.stdout.strip().splitlines():
                _info(f"Detected GPU: {line.strip()}")
            return True
    except OSError:
        pass
    return False


def _torch_has_cuda() -> bool:
    """Check if the installed torch has CUDA support."""
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


def _torch_version() -> str | None:
    try:
        import torch
        return torch.__version__
    except ImportError:
        return None


def ensure_cuda_pytorch() -> None:
    _header("Checking PyTorch CUDA support")

    has_gpu = _detect_nvidia_gpu()
    torch_ver = _torch_version()

    if not has_gpu:
        if platform.system() == "Darwin":
            _ok("macOS detected — MPS acceleration will be used (no CUDA needed)")
        else:
            _info("No NVIDIA GPU detected — CPU-only PyTorch is fine")
        return

    # GPU present — check if torch has CUDA
    if torch_ver is None:
        _warn("PyTorch is not installed yet")
        _info("Install with CUDA support:")
        _info("  pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124")
        return

    if _torch_has_cuda():
        import torch
        _ok(f"PyTorch {torch_ver} with CUDA {torch.version.cuda} — GPU ready")
        return

    # torch installed but no CUDA
    _warn(f"PyTorch {torch_ver} is installed WITHOUT CUDA support")
    _info("Your system has an NVIDIA GPU but torch can't use it.")
    _info("")
    _info("To fix, reinstall PyTorch with CUDA:")
    _info("  pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124")
    _info("")

    answer = input("  Reinstall PyTorch with CUDA now? [y/N] ").strip().lower()
    if answer in ("y", "yes"):
        _info("Reinstalling PyTorch with CUDA 12.4 support...")
        result = _run(
            [
                sys.executable, "-m", "pip", "install",
                "--force-reinstall",
                "torch", "torchvision", "torchaudio",
                "--index-url", "https://download.pytorch.org/whl/cu124",
            ],
        )
        if result.returncode == 0:
            _ok("PyTorch reinstalled with CUDA support")
            _info("Restart this script to verify: python scripts/setup_env.py")
        else:
            _fail("pip install failed:")
            for line in (result.stderr or result.stdout).strip().splitlines()[-5:]:
                print(f"         {line}")
    else:
        _info("Skipped. You can run the pip command above manually later.")


# ---------------------------------------------------------------------------
# 3. Validate dependencies
# ---------------------------------------------------------------------------

# (module_name, label, required)
DEPS = [
    ("torch",              "PyTorch",             True),
    ("transformers",       "HuggingFace Transformers", True),
    ("click",              "Click CLI",           True),
    ("pandas",             "Pandas",              True),
    ("numpy",              "NumPy",               True),
    ("matplotlib",         "Matplotlib",          True),
    ("seaborn",            "Seaborn",             True),
    ("scipy",              "SciPy",               True),
    ("pyarrow",            "PyArrow",             True),
    ("streamlit",          "Streamlit",           True),
    ("plotly",             "Plotly",              True),
    ("rdkit",              "RDKit",               False),
    ("meeko",              "Meeko",               False),
    ("openmm",             "OpenMM",              False),
    ("openff.toolkit",     "OpenFF Toolkit",      False),
    ("openmmforcefields",  "OpenMM Forcefields",  False),
    ("pdbfixer",           "PDBFixer",            False),
    ("mdtraj",             "MDTraj",              False),
    ("peft",               "PEFT (LoRA)",         False),
]


def check_dependencies() -> None:
    _header("Dependency check")

    missing_required = []
    missing_optional = []

    for mod, label, required in DEPS:
        found = _has_module(mod)
        tag = "required" if required else "optional"
        if found:
            _ok(f"{label:<28s} ({tag})")
        elif required:
            _fail(f"{label:<28s} ({tag}) — MISSING")
            missing_required.append(mod)
        else:
            _warn(f"{label:<28s} ({tag}) — not installed")
            missing_optional.append(mod)

    # Check vina CLI separately
    vina_path = shutil.which("vina")
    if vina_path:
        _ok(f"{'Vina CLI':<28s} (optional) — {vina_path}")
    else:
        _warn(f"{'Vina CLI':<28s} (optional) — not on PATH")
        missing_optional.append("vina-cli")

    print()
    if missing_required:
        _fail(f"Missing {len(missing_required)} required package(s): {', '.join(missing_required)}")
        _info("Install with: pip install -e .")
    elif missing_optional:
        _ok(f"All required packages present. {len(missing_optional)} optional package(s) missing.")
        _info("For drug docking:   pip install -e '.[drug]'")
        _info("For MD simulations: conda install -c conda-forge openmm openmmforcefields pdbfixer mdtraj")
    else:
        _ok("All packages present!")


# ---------------------------------------------------------------------------
# 4. Device summary
# ---------------------------------------------------------------------------

def device_summary() -> None:
    _header("Device summary")
    try:
        os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
        import torch

        if torch.cuda.is_available():
            name = torch.cuda.get_device_name(0)
            mem = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
            _ok(f"CUDA: {name} ({mem:.1f} GB) — torch {torch.__version__}+cu{torch.version.cuda}")
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            _ok(f"MPS (Apple Silicon) — torch {torch.__version__}")
        else:
            _info(f"CPU only — torch {torch.__version__}")

        _info(f"select_device() will return: ", )
        # Simulate select_device
        if torch.cuda.is_available():
            print("cuda")
        elif torch.backends.mps.is_available():
            print("mps")
        else:
            print("cpu")

    except ImportError:
        _warn("PyTorch not installed — cannot determine device")
    except Exception as exc:
        _warn(f"Device probe failed: {exc}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"p53CAD Environment Setup — {platform.system()} {platform.machine()}")
    print(f"Python {platform.python_version()} at {sys.executable}")

    fix_symlinks()
    ensure_cuda_pytorch()
    check_dependencies()
    device_summary()

    print(f"\n{'=' * 60}")
    print("  Setup complete. Run 'p53cad doctor' for full diagnostics.")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()
