from __future__ import annotations

from typing import Any

from src.core.logging import get_logger


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)

    from experiments.run_make_report import run as run_core_report

    status = int(run_core_report(args, configs) or 0)

    fancy_enabled = bool(configs.get("report", {}).get("fancy_visualizations", True))
    if not fancy_enabled:
        logger.info("Fancy visualizations disabled by config")
        return status

    try:
        from visualization.create_visualizations import run as run_fancy_report
    except Exception as exc:
        logger.warning("Skipping fancy visualizations: %s", exc)
        return status

    fancy_status = int(run_fancy_report(args, configs) or 0)
    return max(status, fancy_status)
