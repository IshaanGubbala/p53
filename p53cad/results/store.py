from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import pandas as pd

from p53cad.core.logging import get_logger
from p53cad.results.schema import SCHEMA_VERSION


logger = get_logger(__name__)


@dataclass
class RunPaths:
    run_dir: Path
    manifest_path: Path
    scenarios_path: Path
    candidates_path: Path
    trajectories_path: Path
    clinical_path: Path
    top30_path: Path
    top30_csv_path: Path
    summary_path: Path
    checkpoint_state_path: Path
    checkpoint_scenarios_path: Path


class CampaignStore:
    def __init__(self, base_dir: Path | str = "data/campaigns"):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.index_path = self.base_dir / "index.json"

    def run_paths(self, run_id: str) -> RunPaths:
        run_dir = self.base_dir / str(run_id)
        return RunPaths(
            run_dir=run_dir,
            manifest_path=run_dir / "manifest.json",
            scenarios_path=run_dir / "scenarios.parquet",
            candidates_path=run_dir / "candidates.parquet",
            trajectories_path=run_dir / "trajectories.parquet",
            clinical_path=run_dir / "clinical.parquet",
            top30_path=run_dir / "top30.parquet",
            top30_csv_path=run_dir / "top30.csv",
            summary_path=run_dir / "summary.md",
            checkpoint_state_path=run_dir / "checkpoints" / "state.json",
            checkpoint_scenarios_path=run_dir / "checkpoints" / "scenario_status.jsonl",
        )

    def init_run(
        self,
        run_id: str,
        config: Dict[str, Any],
        runtime_caps: Optional[Dict[str, Any]] = None,
        resume: bool = False,
    ) -> RunPaths:
        paths = self.run_paths(run_id)
        paths.run_dir.mkdir(parents=True, exist_ok=True)
        paths.checkpoint_state_path.parent.mkdir(parents=True, exist_ok=True)

        if paths.manifest_path.exists() and resume:
            return paths

        manifest = {
            "run_id": run_id,
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "updated_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "initialized",
            "config": config,
            "runtime_capabilities": runtime_caps or {},
            "artifacts": {
                "scenarios": paths.scenarios_path.name,
                "candidates": paths.candidates_path.name,
                "trajectories": paths.trajectories_path.name,
                "clinical": paths.clinical_path.name,
                "top30": paths.top30_path.name,
                "top30_csv": paths.top30_csv_path.name,
                "summary": paths.summary_path.name,
            },
        }
        self._write_json(paths.manifest_path, manifest)
        self._update_index(run_id, status="initialized")
        return paths

    def update_manifest(self, run_id: str, patch: Dict[str, Any]) -> None:
        paths = self.run_paths(run_id)
        if not paths.manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found for run {run_id}")
        manifest = self._read_json(paths.manifest_path)
        manifest.update(patch)
        manifest["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        self._write_json(paths.manifest_path, manifest)
        if "status" in patch:
            self._update_index(run_id, status=str(patch["status"]))

    def append_scenario_checkpoint(self, run_id: str, payload: Dict[str, Any]) -> None:
        paths = self.run_paths(run_id)
        paths.checkpoint_scenarios_path.parent.mkdir(parents=True, exist_ok=True)
        with paths.checkpoint_scenarios_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(payload) + "\n")

    def load_scenario_checkpoints(self, run_id: str) -> pd.DataFrame:
        paths = self.run_paths(run_id)
        if not paths.checkpoint_scenarios_path.exists():
            return pd.DataFrame()
        rows = []
        with paths.checkpoint_scenarios_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    logger.warning("Skipping malformed checkpoint row in %s", paths.checkpoint_scenarios_path)
        return pd.DataFrame(rows)

    def write_tables(
        self,
        run_id: str,
        scenarios_df: pd.DataFrame,
        candidates_df: pd.DataFrame,
        trajectories_df: pd.DataFrame,
        clinical_df: pd.DataFrame,
    ) -> None:
        paths = self.run_paths(run_id)
        self._write_parquet(scenarios_df, paths.scenarios_path)
        self._write_parquet(candidates_df, paths.candidates_path)
        self._write_parquet(trajectories_df, paths.trajectories_path)
        self._write_parquet(clinical_df, paths.clinical_path)

    def write_top30(self, run_id: str, top30_df: pd.DataFrame) -> None:
        paths = self.run_paths(run_id)
        self._write_parquet(top30_df, paths.top30_path)
        top30_df.to_csv(paths.top30_csv_path, index=False)

    def write_summary(self, run_id: str, text: str) -> None:
        paths = self.run_paths(run_id)
        paths.summary_path.write_text(text, encoding="utf-8")

    def load_run_bundle(self, run_id: str) -> Dict[str, Any]:
        paths = self.run_paths(run_id)
        manifest = self._read_json(paths.manifest_path)

        def _safe_read(path: Path) -> pd.DataFrame:
            if path.exists():
                return pd.read_parquet(path)
            return pd.DataFrame()

        return {
            "manifest": manifest,
            "scenarios": _safe_read(paths.scenarios_path),
            "candidates": _safe_read(paths.candidates_path),
            "trajectories": _safe_read(paths.trajectories_path),
            "clinical": _safe_read(paths.clinical_path),
            "top30": _safe_read(paths.top30_path),
            "run_dir": paths.run_dir,
        }

    def list_runs(self) -> pd.DataFrame:
        if not self.index_path.exists():
            return pd.DataFrame(columns=["run_id", "status", "updated_at_utc", "path"])
        payload = self._read_json(self.index_path)
        runs = payload.get("runs", []) if isinstance(payload, dict) else []
        if not runs:
            return pd.DataFrame(columns=["run_id", "status", "updated_at_utc", "path"])
        df = pd.DataFrame(runs)
        if "updated_at_utc" in df.columns:
            df = df.sort_values("updated_at_utc", ascending=False)
        return df

    def latest_run_id(self) -> Optional[str]:
        df = self.list_runs()
        if df.empty:
            return None
        return str(df.iloc[0]["run_id"])

    def _write_parquet(self, df: pd.DataFrame, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(path, index=False)

    def _read_json(self, path: Path) -> Dict[str, Any]:
        return json.loads(path.read_text(encoding="utf-8"))

    def _write_json(self, path: Path, payload: Dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def _update_index(self, run_id: str, status: str) -> None:
        now = datetime.now(timezone.utc).isoformat()
        payload = {"runs": []}
        if self.index_path.exists():
            try:
                payload = self._read_json(self.index_path)
            except Exception:
                payload = {"runs": []}

        runs = payload.get("runs", []) if isinstance(payload, dict) else []
        runs = [r for r in runs if str(r.get("run_id")) != str(run_id)]
        runs.append(
            {
                "run_id": str(run_id),
                "status": str(status),
                "updated_at_utc": now,
                "path": str((self.base_dir / str(run_id)).resolve()),
            }
        )
        runs = sorted(runs, key=lambda r: str(r.get("updated_at_utc", "")), reverse=True)
        self._write_json(self.index_path, {"runs": runs})
