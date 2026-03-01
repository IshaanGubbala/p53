from __future__ import annotations

import importlib.util
import os
import platform
import shutil
from typing import Any, Dict, Optional

from p53cad.core.logging import get_logger


def _has_module(module_name: str) -> bool:
    try:
        return importlib.util.find_spec(module_name) is not None
    except ModuleNotFoundError:
        return False


def select_device(preference: Optional[str] = None) -> "torch.device":
    """Return the best available torch device.

    Priority: explicit *preference* > CUDA > MPS > CPU.
    Passing ``"auto"`` or ``None`` triggers auto-detection.
    """
    import torch  # pylint: disable=import-outside-toplevel

    if preference and preference != "auto":
        return torch.device(preference)
    force_cpu = os.environ.get("P53CAD_FORCE_CPU", "1" if platform.system() == "Windows" else "0") == "1"
    if force_cpu:
        return torch.device("cpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def bootstrap_runtime(
    seed: Optional[int] = None,
    enable_tf32: bool = False,  # Disabled by default - can cause instability on some GPUs
    enable_cudnn_benchmark: bool = True,
    enable_sdpa: bool = True,
    enable_flash_attention: bool = True,
) -> Dict[str, Any]:
    """
    Apply process-level runtime guards before heavy ML/scientific imports.
    
    Parameters
    ----------
    seed : int, optional
        Random seed for reproducibility.
    enable_tf32 : bool
        Enable TF32 tensor cores on Ampere+ GPUs for faster matrix multiplications.
        Safe: gracefully disables on older GPUs.
    enable_cudnn_benchmark : bool
        Enable cuDNN auto-tuner for optimal convolution algorithms.
        Safe: only affects convolution layers.
    enable_sdpa : bool
        Enable scaled dot-product attention (SDPA) for faster transformer inference.
        Falls back to eager mode if unstable.
    enable_flash_attention : bool
        Enable Flash Attention when available. Requires compatible GPU.
    
    Returns
    -------
    Dict[str, Any]
        Summary of enabled optimizations.
    """
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
    if seed is not None:
        os.environ.setdefault("PYTHONHASHSEED", str(int(seed)))
    
    enabled = {
        "tf32": False,
        "cudnn_benchmark": False,
        "sdpa": False,
        "flash_attention": False,
    }

    if enable_tf32 or enable_cudnn_benchmark or enable_sdpa or enable_flash_attention:
        try:
            import torch
            if torch.cuda.is_available():
                if enable_cudnn_benchmark:
                    torch.backends.cudnn.benchmark = True
                    enabled["cudnn_benchmark"] = True
                if enable_tf32:
                    torch.backends.cuda.matmul.allow_tf32 = True
                    torch.backends.cudnn.allow_tf32 = True
                    enabled["tf32"] = True
                if enable_sdpa:
                    try:
                        torch.backends.cuda.enable_flash_sdp(enable_flash_attention)
                        torch.backends.cuda.enable_math_sdp(not enable_flash_attention)
                        torch.backends.cuda.enable_mem_efficient_sdp(enable_flash_attention)
                        enabled["sdpa"] = True
                        enabled["flash_attention"] = enable_flash_attention and torch.cuda.is_available()
                    except Exception as e:
                        import logging as _logging
                        _logging.getLogger("p53cad.runtime").warning(
                            "SDPA enable failed: %s. Falling back to eager attention.", e
                        )
        except Exception:
            pass

    return enabled


def get_runtime_capabilities() -> Dict[str, Any]:
    """
    Collect lightweight runtime capability probes without forcing heavy imports.
    """
    caps: Dict[str, Any] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "kmp_duplicate_lib_ok": os.getenv("KMP_DUPLICATE_LIB_OK", ""),
        "mps_fallback_enabled": os.getenv("PYTORCH_ENABLE_MPS_FALLBACK", ""),
        "torch_installed": _has_module("torch"),
        "streamlit_installed": _has_module("streamlit"),
        "plotly_installed": _has_module("plotly"),
        "rdkit_installed": _has_module("rdkit"),
        "openmm_installed": _has_module("openmm"),
        "openff_installed": _has_module("openff.toolkit"),
        "openmmforcefields_installed": _has_module("openmmforcefields"),
        "meeko_installed": _has_module("meeko"),
        "vina_py_installed": _has_module("vina"),
        "vina_cli_available": shutil.which("vina") is not None,
        "pdbfixer_installed": _has_module("pdbfixer"),
        "mdtraj_installed": _has_module("mdtraj"),
        "esmfold_local_available": _has_module("transformers") and _has_module("torch"),
        "physics_validation_ready": (
            _has_module("openmm")
            and _has_module("pdbfixer")
            and _has_module("mdtraj")
            and _has_module("transformers")
            and _has_module("torch")
        ),
    }

    if caps["torch_installed"]:
        try:
            import torch  # pylint: disable=import-outside-toplevel

            caps["torch_version"] = torch.__version__
            caps["mps_available"] = bool(torch.backends.mps.is_available())
            caps["cuda_available"] = bool(torch.cuda.is_available())
        except Exception as exc:  # pragma: no cover - defensive path
            caps["torch_probe_error"] = str(exc)
            caps["mps_available"] = False
            caps["cuda_available"] = False
    else:
        caps["mps_available"] = False
        caps["cuda_available"] = False

    if _has_module("transformers"):
        try:
            import transformers  # pylint: disable=import-outside-toplevel

            caps["transformers_version"] = transformers.__version__
        except Exception as exc:  # pragma: no cover - defensive path
            caps["transformers_probe_error"] = str(exc)

    return caps


def log_runtime_capabilities(logger_name: str = "p53cad.runtime") -> Dict[str, Any]:
    """
    Log a single compact capability line for workflow integrity debugging.
    """
    logger = get_logger(logger_name)
    caps = get_runtime_capabilities()
    logger.info(
        "Runtime capabilities | python=%s torch=%s mps=%s cuda=%s transformers=%s rdkit=%s "
        "vina_cli=%s openmm=%s openff=%s openmmforcefields=%s",
        caps.get("python"),
        caps.get("torch_version", "n/a"),
        caps.get("mps_available"),
        caps.get("cuda_available"),
        caps.get("transformers_version", "n/a"),
        caps.get("rdkit_installed"),
        caps.get("vina_cli_available"),
        caps.get("openmm_installed"),
        caps.get("openff_installed"),
        caps.get("openmmforcefields_installed"),
    )
    return caps
