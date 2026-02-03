from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

try:
    import yaml
except ImportError as exc:  # pragma: no cover - used at runtime
    raise RuntimeError("PyYAML is required. Install with `pip install pyyaml`.") from exc


def load_yaml(path: Path) -> dict[str, Any]:
    """Load YAML from disk, returning an empty dict if the file is empty."""
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    return payload or {}


def load_configs(config_dir: Path, names: Iterable[str]) -> dict[str, dict[str, Any]]:
    """Load multiple YAML files by basename (without .yaml)."""
    configs: dict[str, dict[str, Any]] = {}
    for name in names:
        path = config_dir / f"{name}.yaml"
        configs[name] = load_yaml(path)
    return configs


def normalize_paths(paths_cfg: dict[str, Any], project_root: Path) -> dict[str, Any]:
    """Convert relative paths in the paths config to absolute Path objects."""
    data_root_value = paths_cfg.get("data_root")
    data_root = _resolve_path(data_root_value, project_root) if data_root_value else project_root

    def _convert(value: Any, key: str | None = None) -> Any:
        if isinstance(value, dict):
            return {child_key: _convert(child_value, child_key) for child_key, child_value in value.items()}
        if isinstance(value, list):
            return [_convert(item, key) for item in value]
        if isinstance(value, str):
            base = project_root
            if key in {"raw", "interim", "processed", "cache_dir"}:
                base = data_root
            return _resolve_path(value, base)
        return value

    return _convert(paths_cfg)


def _resolve_path(value: str | Path | None, base: Path) -> Path:
    if value is None:
        return base
    path = Path(value)
    if path.is_absolute():
        return path.expanduser().resolve()
    return (base / path).expanduser().resolve()
