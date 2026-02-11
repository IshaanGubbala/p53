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


def bootstrap_runtime(seed: Optional[int] = None) -> None:
    """
    Apply process-level runtime guards before heavy ML/scientific imports.
    """
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
    if seed is not None:
        os.environ.setdefault("PYTHONHASHSEED", str(int(seed)))


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
