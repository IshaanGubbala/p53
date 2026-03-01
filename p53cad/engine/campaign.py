from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

import math

from p53cad.analysis.clinical_impact import ClinicalImpactEngine
from p53cad.core.logging import get_logger
from p53cad.data.dms import P53_WT, apply_mutation, get_dms_data, get_pairwise_epistasis_lookup, parse_single_mutation
from p53cad.engine.latent import ManifoldEmbedder
from p53cad.engine.oracle import AttentionPoolingNet, FunctionalOracle, compute_conditional_rescue_scores
from p53cad.engine.structural_awareness import (
    PrecomputedConstants,
    build_precomputed_constants,
    compute_mutation_confidence_penalty,
    compute_interface_floor_penalty,
    compute_crater_penalty,
    compute_contact_consistency_penalty,
    compute_contact_map_regularization,
    AsyncLogger,
    CUDASafetyGuard,
    get_optimization_config,
    PopulationShare,
    ESMFoldGateConfig,
    check_esmfold_gate,
    DNA_BINDING_INTERFACE_RESIDUES,
)
from p53cad.results.schema import (
    BIG8_HOTSPOTS,
    DEFAULT_DELIVERY_METHODS,
    build_run_id,
    build_scenario_matrix,
    scenarios_to_frame,
    select_presentation_shortlist,
)
from p53cad.results.store import CampaignStore


logger = get_logger(__name__)


@contextmanager
def _prevent_sleep():
    """Prevent idle/system sleep during long campaigns.

    macOS: ``caffeinate -is`` subprocess.
    Windows: ``SetThreadExecutionState`` via ctypes (process-scoped, auto-reverts).
    Linux: no-op (most servers don't sleep).
    """
    proc = None
    _win_reset = False
    system = platform.system()

    if system == "Darwin" and shutil.which("caffeinate"):
        try:
            proc = subprocess.Popen(
                ["caffeinate", "-is"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            logger.info("Sleep prevention enabled (caffeinate PID %d)", proc.pid)
        except OSError:
            logger.debug("Failed to start caffeinate, continuing without sleep prevention")

    elif system == "Windows":
        try:
            import ctypes
            ES_CONTINUOUS = 0x80000000
            ES_SYSTEM_REQUIRED = 0x00000001
            ctypes.windll.kernel32.SetThreadExecutionState(
                ES_CONTINUOUS | ES_SYSTEM_REQUIRED
            )
            _win_reset = True
            logger.info("Sleep prevention enabled (SetThreadExecutionState)")
        except Exception:
            logger.debug("Failed to set Windows execution state, continuing without sleep prevention")

    try:
        yield
    finally:
        if proc is not None:
            proc.terminate()
            proc.wait(timeout=5)
            logger.info("Sleep prevention disabled")
        if _win_reset:
            try:
                import ctypes
                ES_CONTINUOUS = 0x80000000
                ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS)
                logger.info("Sleep prevention disabled (SetThreadExecutionState reset)")
            except Exception:
                pass


AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

WEIGHT_PROFILES = [
    {"function": 4.0, "stability": 8.0, "binding": 2.5, "name": "Balanced", "color": "#2563EB"},
    {"function": 2.0, "stability": 15.0, "binding": 2.0, "name": "Stability-First", "color": "#0EA5E9"},
    {"function": 3.0, "stability": 5.0, "binding": 8.0, "name": "Binding-Optimized", "color": "#10B981"},
    {"function": 8.0, "stability": 4.0, "binding": 3.0, "name": "Function-Maximized", "color": "#14B8A6"},
    {"function": 5.0, "stability": 10.0, "binding": 5.0, "name": "Conservative", "color": "#1D4ED8"},
    {"function": 6.0, "stability": 6.0, "binding": 6.0, "name": "Experimental", "color": "#0891B2"},
]


@dataclass
class ScenarioRuntime:
    scenario_id: str
    target_label: str
    targets: List[str]
    delivery_method: str


class CampaignRunner:
    def __init__(
        self,
        *,
        store: Optional[CampaignStore] = None,
        oracle_path: str | Path = "data/models/functional_oracle.pt",
        device: Optional[str] = None,
    ):
        self.store = store or CampaignStore()
        self.oracle_path = Path(oracle_path)
        self.device = device

        self.embedder: Optional[ManifoldEmbedder] = None
        self.oracle: Optional[FunctionalOracle] = None
        self.clinical: Optional[ClinicalImpactEngine] = None

        self.calibration_profile = self._score_calibration_profile()
        self.uncertainty_weight = 0.8
        self.ood_rank_weight = 1.25
        self.ood_radius = 1.75
        self.ood_loss_weight = 12.0

        self._dms_lookup: Dict[Tuple[int, str], float] = {}
        self._dms_quality_tensor: Optional[torch.Tensor] = None
        self._load_dms_lookup()

        self._pairwise_dms: Dict[tuple, float] = get_pairwise_epistasis_lookup()

        self._wt_contacts: Dict[int, List[Tuple[int, float]]] = {}
        self._ca_coords: Dict[int, np.ndarray] = {}
        self._load_contact_map()

        self._wt_hidden: Optional[torch.Tensor] = None

        self._emb_cache: Dict[str, torch.Tensor] = {}
        self._baseline_cache: Dict[str, Dict[str, float]] = {}

        self._optimization_config: Dict[str, bool] = {}
        self._population_share = PopulationShare(max_size=100)
        self._async_logger = AsyncLogger(interval=20)

        self._esmfold_gate_config = ESMFoldGateConfig(
            min_global_plddt=70.0,
            min_dbd_plddt=65.0,
            min_interface_plddt=60.0,
            max_wt_rmsd=5.0,
            enabled=True,
        )

    def _load_dms_lookup(self) -> None:
        """Load raw DMS Z-scores into a (position, amino_acid) → Z-score dict."""
        dms_path = Path(__file__).parent.parent.parent / "data" / "raw" / "p53_DMS_Giacomelli_2018.csv"
        if not dms_path.exists():
            logger.warning("DMS CSV not found at %s; DMS-aware optimization disabled", dms_path)
            return
        try:
            raw = pd.read_csv(dms_path)
            score_col = "A549_p53WT_Nutlin-3_Z-score"
            if score_col not in raw.columns:
                logger.warning("DMS score column '%s' not found; DMS-aware optimization disabled", score_col)
                return
            for _, row in raw.iterrows():
                pos = int(row.get("Position", 0))
                aa = str(row.get("AA_variant", ""))
                z = row.get(score_col)
                if pd.isna(z) or pos < 1 or pos > len(P53_WT) or len(aa) != 1:
                    continue
                self._dms_lookup[(pos, aa.upper())] = float(z)
            logger.info("Loaded %d DMS entries for evidence-aware optimization", len(self._dms_lookup))
        except Exception as exc:
            logger.warning("Failed to load DMS lookup: %s", exc)

    def _load_contact_map(self) -> None:
        """Parse WT PDB and build contact map for structural penalties."""
        pdb_path = Path(__file__).parent.parent.parent / "data" / "raw" / "p53_wt.pdb"
        if not pdb_path.exists():
            logger.warning("WT PDB not found at %s; contact/epistasis penalties disabled", pdb_path)
            return
        try:
            from p53cad.engine.explainability import EnergyDecomposer
            self._ca_coords = EnergyDecomposer._parse_ca_coordinates(str(pdb_path))
            all_positions = sorted(self._ca_coords.keys())
            self._wt_contacts = EnergyDecomposer._find_contacts(
                self._ca_coords, all_positions, cutoff=8.0
            )
            logger.info("Contact map: %d residues, %d contacts loaded",
                        len(self._ca_coords),
                        sum(len(v) for v in self._wt_contacts.values()))
        except Exception as exc:
            logger.warning("Failed to load contact map: %s", exc)

    def _get_wt_hidden(self, device: torch.device) -> Optional[torch.Tensor]:
        """Cache the WT hidden states for contact preservation penalty."""
        if self._wt_hidden is not None and self._wt_hidden.device == device:
            return self._wt_hidden
        if self.embedder is None:
            return None
        emb_wt = self._get_cached_embedding(P53_WT)
        with torch.no_grad():
            h, _, _ = self.embedder.latent_forward_ascent(emb_wt)
            self._wt_hidden = h.detach()
        return self._wt_hidden

    def _get_dms_quality_tensor(self, device: torch.device) -> Optional[torch.Tensor]:
        """Build (seq_len, 20) tensor of raw DMS Z-scores for differentiable DMS penalty.

        Values: raw Nutlin-3 Z-scores with a compensatory dead zone applied:
        - Z < -0.5: kept as-is (strongly functional — safe rescue candidates)
        - -0.5 ≤ Z ≤ 2.0: zeroed out (neutral/compensatory zone — may rescue in context)
        - Z > 2.0: kept as-is (catastrophically LoF — penalize)

        Missing entries default to 0.0 (neutral — no guidance).

        The dead zone prevents penalizing mutations in the Z=+0.14 to +1.59 range
        where known intragenic suppressors (compensatory rescues) are found.
        """
        if not self._dms_lookup:
            return None
        if self._dms_quality_tensor is not None and self._dms_quality_tensor.device == device:
            return self._dms_quality_tensor

        assert self.embedder is not None
        # Map AA letters → index in AA_IDS
        aa_to_idx: Dict[str, int] = {}
        for i, aa_id in enumerate(AA_IDS):
            tokens = self.embedder.tokenizer.convert_ids_to_tokens([aa_id])
            if tokens:
                aa_to_idx[str(tokens[0]).upper()] = i

        tensor = torch.zeros(len(P53_WT), len(AA_IDS), device=device)
        filled = 0
        for (pos, aa), z_score in self._dms_lookup.items():
            if aa in aa_to_idx and 1 <= pos <= len(P53_WT):
                tensor[pos - 1, aa_to_idx[aa]] = z_score
                filled += 1

        # Zero out the compensatory zone where rescue mutations live
        # Z < -0.5: keep (strongly functional — safe rescue candidates)
        # -0.5 ≤ Z ≤ 2.0: zero (neutral — could be compensatory)
        # Z > 2.0: keep (catastrophically LoF — penalize)
        compensatory = (tensor >= -0.5) & (tensor <= 2.0)
        n_zeroed = int(compensatory.sum().item())
        tensor[compensatory] = 0.0

        logger.info(
            "DMS quality tensor: %d/%d entries filled, %d zeroed (compensatory dead zone) on %s",
            filled, len(self._dms_lookup), n_zeroed, device,
        )
        self._dms_quality_tensor = tensor
        return tensor

    def _load_models(self) -> None:
        if self.embedder is not None and self.oracle is not None and self.clinical is not None:
            return
        
        # Reset CUDA state for clean start
        if torch.cuda.is_available():
            try:
                torch.cuda.synchronize()
                torch.cuda.empty_cache()
            except Exception:
                pass
        
        self.embedder = ManifoldEmbedder(device=self.device)
        self.oracle = FunctionalOracle(model_path=self.oracle_path, device=self.device)
        self.clinical = ClinicalImpactEngine()

        device_obj = self.embedder.device if self.embedder else torch.device("cpu")
        is_windows_cuda = (device_obj.type == "cuda" and platform.system() == "Windows")
        enable_shallow = os.environ.get("P53CAD_ENABLE_SHALLOW", "0" if is_windows_cuda else "1") == "1"
        enable_sdpa = os.environ.get("P53CAD_ENABLE_SDPA", "0" if is_windows_cuda else "1") == "1"
        enable_batched = os.environ.get("P53CAD_ENABLE_BATCHED", "0") == "1"
        self._optimization_config = get_optimization_config(
            device_obj,
            enable_shallow_gradients=enable_shallow,
            enable_batched_trials=enable_batched,
            enable_sdpa=enable_sdpa,
        )
        if device_obj.type == "cuda":
            torch.backends.cuda.matmul.allow_tf32 = bool(self._optimization_config.get("tf32_enabled", True))
            torch.backends.cudnn.benchmark = bool(self._optimization_config.get("cudnn_benchmark", True))
        logger.info("Optimization config: %s", self._optimization_config)

    def _is_cuda_runtime_error(self, exc: Exception) -> bool:
        """Best-effort detection of CUDA runtime faults that poison context."""
        msg = str(exc).lower()
        if "cuda error" in msg:
            return True
        cuda_markers = (
            "cublas_status",
            "illegal memory access",
            "illegal instruction",
            "device-side assert",
            "invalid argument",
            "cudnn",
        )
        return any(m in msg for m in cuda_markers)

    def _recover_from_cuda_fault(self) -> None:
        """Drop model state and clear CUDA cache so next attempt starts fresh."""
        self.embedder = None
        self.oracle = None
        self._wt_hidden = None
        self._dms_quality_tensor = None
        self._emb_cache.clear()
        self._baseline_cache.clear()
        if torch.cuda.is_available():
            try:
                torch.cuda.synchronize()
            except Exception:
                pass
            try:
                torch.cuda.empty_cache()
            except Exception:
                pass

    def _run_scenario_with_recovery(
        self,
        scenario: ScenarioRuntime,
        *,
        pass_name: str,
        config: Dict[str, int],
        seed: int,
        allowed_profiles: Optional[List[str]] = None,
        budget: str = "high",
    ) -> Dict[str, Any]:
        retries = int(os.environ.get("P53CAD_CUDA_RECOVERY_RETRIES", "1"))
        for attempt in range(retries + 1):
            try:
                if self.embedder is None or self.oracle is None or self.clinical is None:
                    self._load_models()
                return self._run_scenario(
                    scenario,
                    pass_name=pass_name,
                    config=config,
                    seed=seed,
                    allowed_profiles=allowed_profiles,
                    budget=budget,
                )
            except Exception as exc:
                is_cuda_fault = self._is_cuda_runtime_error(exc)
                can_retry = (
                    is_cuda_fault
                    and torch.cuda.is_available()
                    and attempt < retries
                )
                if not can_retry:
                    raise
                logger.warning(
                    "CUDA fault in scenario %s (%s, attempt %d/%d): %s",
                    scenario.scenario_id,
                    pass_name,
                    attempt + 1,
                    retries + 1,
                    exc,
                )
                # Force conservative mode before retry.
                os.environ["P53CAD_FORCE_FP32"] = "1"
                os.environ["P53CAD_ENABLE_BF16_AUTOCAST"] = "0"
                os.environ["P53CAD_ENABLE_SDPA"] = "0"
                os.environ["P53CAD_ENABLE_SHALLOW"] = "0"
                os.environ["P53CAD_ENABLE_BATCHED"] = "0"
                os.environ["P53CAD_TRIAL_LR"] = os.environ.get("P53CAD_TRIAL_LR", "0.01")
                self._recover_from_cuda_fault()

    def _get_cached_embedding(self, seq: str) -> torch.Tensor:
        """Return detached embedding for *seq*, using cache to avoid redundant ESM-2 passes."""
        if seq not in self._emb_cache:
            assert self.embedder is not None
            self._emb_cache[seq] = self.embedder.get_embeddings(seq).detach()
        return self._emb_cache[seq]

    def _get_cached_baseline(
        self, target_seq: str, pooled_target_ref: torch.Tensor, mc_samples: int
    ) -> Dict[str, float]:
        """Return baseline metrics for *target_seq*, using cache for shared targets."""
        if target_seq not in self._baseline_cache:
            self._baseline_cache[target_seq] = self._compute_baseline_metrics(
                target_seq, pooled_target_ref, mc_samples=mc_samples
            )
        return self._baseline_cache[target_seq]

    def run(
        self,
        *,
        run_id: Optional[str] = None,
        seed: int = 42,
        resume: bool = True,
        include_pairs: bool = True,
        hotspots: Sequence[str] | None = None,
        delivery_methods: Sequence[str] | None = None,
        max_scenarios: Optional[int] = None,
        shortlist_n: int = 30,
        budget: str = "high",
        with_clinical: bool = True,
    ) -> Dict[str, Any]:
        with _prevent_sleep():
            return self._run_inner(
                run_id=run_id, seed=seed, resume=resume,
                include_pairs=include_pairs, hotspots=hotspots,
                delivery_methods=delivery_methods,
                max_scenarios=max_scenarios, shortlist_n=shortlist_n,
                budget=budget, with_clinical=with_clinical,
            )

    def _run_inner(
        self,
        *,
        run_id: Optional[str] = None,
        seed: int = 42,
        resume: bool = True,
        include_pairs: bool = True,
        hotspots: Sequence[str] | None = None,
        delivery_methods: Sequence[str] | None = None,
        max_scenarios: Optional[int] = None,
        shortlist_n: int = 30,
        budget: str = "high",
        with_clinical: bool = True,
    ) -> Dict[str, Any]:
        self._load_models()
        assert self.embedder is not None
        assert self.oracle is not None
        assert self.clinical is not None

        # Clear per-run caches
        self._emb_cache.clear()
        self._baseline_cache.clear()
        self._population_share.clear()
        self._async_logger.reset_counter()

        run_id = run_id or build_run_id("campaign")
        scenarios = build_scenario_matrix(
            hotspots=hotspots or BIG8_HOTSPOTS,
            delivery_methods=delivery_methods or DEFAULT_DELIVERY_METHODS,
            include_pairs=include_pairs,
        )
        if max_scenarios is not None:
            scenarios = scenarios[: int(max_scenarios)]

        config = {
            "seed": int(seed),
            "budget": budget,
            "include_pairs": bool(include_pairs),
            "hotspots": list(hotspots or BIG8_HOTSPOTS),
            "delivery_methods": list(delivery_methods or DEFAULT_DELIVERY_METHODS),
            "max_scenarios": int(max_scenarios) if max_scenarios is not None else None,
            "shortlist_n": int(shortlist_n),
            "with_clinical": bool(with_clinical),
            "pass_a": self._pass_a_config(budget),
            "pass_b": self._pass_b_config(budget),
        }

        runtime_caps = {
            "device": str(self.embedder.device),
            "oracle_path": str(self.oracle_path),
        }
        paths = self.store.init_run(run_id, config=config, runtime_caps=runtime_caps, resume=resume)
        self.store.update_manifest(run_id, {"status": "running"})

        ckpt_dir = paths.checkpoint_state_path.parent
        candidates_jsonl = ckpt_dir / "candidates.jsonl"
        trajectories_jsonl = ckpt_dir / "trajectories.jsonl"
        clinical_jsonl = ckpt_dir / "clinical.jsonl"
        scenarios_jsonl = ckpt_dir / "scenario_metrics.jsonl"

        done_df = self.store.load_scenario_checkpoints(run_id)
        done_pairs = set()
        if not done_df.empty:
            for _, row in done_df.iterrows():
                if str(row.get("status", "")) == "done":
                    done_pairs.add((str(row.get("scenario_id")), str(row.get("pass_name"))))

        pass_a_cfg = self._pass_a_config(budget)
        pass_b_cfg = self._pass_b_config(budget)

        scenario_objs = [
            ScenarioRuntime(
                scenario_id=s.scenario_id,
                target_label=s.target_label,
                targets=list(s.targets),
                delivery_method=s.delivery_method,
            )
            for s in scenarios
        ]

        # Pass A screening
        pass_a_best: Dict[str, float] = {}
        for idx, sc in enumerate(scenario_objs, start=1):
            key = (sc.scenario_id, "screen")
            if resume and key in done_pairs:
                continue
            logger.info("Pass A [%d/%d] %s", idx, len(scenario_objs), sc.scenario_id)
            out = self._run_scenario_with_recovery(
                sc, pass_name="screen", config=pass_a_cfg, seed=seed, budget=budget
            )
            self._append_rows(candidates_jsonl, out["candidates"])
            self._append_rows(trajectories_jsonl, out["trajectories"])
            self._append_rows(scenarios_jsonl, [out["scenario_metrics"]])
            pass_a_best[sc.scenario_id] = float(out["scenario_metrics"].get("best_score", -np.inf))
            self.store.append_scenario_checkpoint(
                run_id,
                {
                    "scenario_id": sc.scenario_id,
                    "pass_name": "screen",
                    "status": "done",
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                },
            )

        # Load pass A from checkpoint to select pass B scenarios.
        screen_df = self._read_jsonl(scenarios_jsonl)
        if screen_df.empty:
            raise RuntimeError("No screen metrics found; campaign cannot continue.")
        screen_df = screen_df[screen_df["pass_name"] == "screen"].copy()
        screen_df["best_score"] = pd.to_numeric(screen_df["best_score"], errors="coerce").fillna(-np.inf)
        screen_df = screen_df.sort_values("best_score", ascending=False)
        top_count = max(1, int(np.ceil(len(screen_df) * 0.33)))

        # Delivery diversity floor: ensure at least top-2 per delivery method
        # reach Pass B, then fill remaining with global top scores.
        delivery_floor = 2
        selected_pass_b: set[str] = set()
        if "scenario_id" in screen_df.columns:
            dm_col = screen_df["scenario_id"].astype(str).str.rsplit("__", n=1)
            screen_df["_delivery"] = dm_col.str[-1]
            for dm in DEFAULT_DELIVERY_METHODS:
                dm_rows = screen_df[screen_df["_delivery"] == dm]
                for sid in dm_rows.head(delivery_floor)["scenario_id"].astype(str):
                    selected_pass_b.add(sid)
            # Fill remaining slots from global ranking
            for sid in screen_df["scenario_id"].astype(str):
                if len(selected_pass_b) >= top_count:
                    break
                selected_pass_b.add(sid)
            screen_df.drop(columns=["_delivery"], inplace=True)
        else:
            selected_pass_b = set(screen_df.head(top_count)["scenario_id"].astype(str).tolist())

        # Profile pruning: pick top-3 profiles per scenario from Pass A candidates.
        screen_cands = self._read_jsonl(candidates_jsonl)
        profile_map: Dict[str, Optional[List[str]]] = {}
        if not screen_cands.empty and "profile" in screen_cands.columns and "score" in screen_cands.columns:
            sc_screen = screen_cands[screen_cands["pass_name"] == "screen"].copy() if "pass_name" in screen_cands.columns else screen_cands.copy()
            sc_screen["score"] = pd.to_numeric(sc_screen["score"], errors="coerce").fillna(-np.inf)
            for sid in selected_pass_b:
                sc_rows = sc_screen[sc_screen["scenario_id"] == sid]
                if sc_rows.empty:
                    profile_map[sid] = None  # fallback: all profiles
                else:
                    top_profiles = (
                        sc_rows.groupby("profile")["score"]
                        .max()
                        .sort_values(ascending=False)
                        .head(3)
                        .index.tolist()
                    )
                    profile_map[sid] = top_profiles if len(top_profiles) >= 1 else None
        logger.info(
            "Pass B: %d scenarios selected (%.0f%%), profiles pruned to top-3 per scenario",
            len(selected_pass_b),
            100.0 * len(selected_pass_b) / max(len(screen_df), 1),
        )

        # ESMFold quality gate: filter candidates before deep refinement
        if self._esmfold_gate_config.enabled and not screen_cands.empty:
            try:
                from p53cad.engine.physics_validation import LocalESMFoldPredictor

                gate_passed = []
                gate_failed = []
                screen_cols = screen_cands.columns.tolist()
                gate_skipped = 0
                gate_measured = 0
                gate_cache_dir = paths.run_dir / "esmfold_gate_cache"
                gate_device = str(self.embedder.device) if self.embedder is not None else "cpu"
                gate_predictor = LocalESMFoldPredictor(device=gate_device, cache_dir=gate_cache_dir)
                
                for sid in selected_pass_b:
                    sc_rows = screen_cands[screen_cands["scenario_id"] == sid]
                    if sc_rows.empty:
                        gate_passed.append(sid)
                        continue
                    
                    best_cand = sc_rows.nlargest(1, "score") if "score" in screen_cols else None
                    if best_cand is None or best_cand.empty:
                        gate_passed.append(sid)
                        continue
                    
                    seq = str(best_cand.iloc[0].get("sequence", ""))
                    if not seq or len(seq) != len(P53_WT):
                        gate_failed.append((sid, "invalid_sequence"))
                        continue

                    plddt_scores: List[float] = []
                    try:
                        pred = gate_predictor.predict(seq)
                        plddt_scores = list(pred.plddt_scores)
                    except Exception as exc:
                        gate_skipped += 1
                        logger.warning("ESMFold gate predictor failed for %s: %s", sid, exc)

                    passed, metrics = check_esmfold_gate(
                        plddt_scores=plddt_scores,
                        config=self._esmfold_gate_config,
                    )
                    if plddt_scores:
                        gate_measured += 1

                    if passed:
                        gate_passed.append(sid)
                    else:
                        gate_failed.append((sid, metrics.get("reason", "unknown")))
                
                n_filtered = len(selected_pass_b) - len(gate_passed)
                logger.info(
                    "ESMFold gate metrics available for %d/%d scenarios.",
                    gate_measured,
                    len(selected_pass_b),
                )
                if gate_skipped > 0:
                    logger.warning(
                        "ESMFold gate skipped for %d scenarios (predictor failure); not filtering those scenarios.",
                        gate_skipped,
                    )
                if n_filtered > 0:
                    logger.info(
                        "ESMFold gate: %d/%d passed, %d filtered (gate enabled=%s)",
                        len(gate_passed), len(selected_pass_b), n_filtered, self._esmfold_gate_config.enabled
                    )
                    selected_pass_b = set(gate_passed)
            except Exception as exc:
                logger.warning("ESMFold gate check failed: %s. Proceeding without gate.", exc)

        # Pass B deep refinement on top scenarios only.
        selected_objs = [sc for sc in scenario_objs if sc.scenario_id in selected_pass_b]
        for idx, sc in enumerate(selected_objs, start=1):
            key = (sc.scenario_id, "deep")
            if resume and key in done_pairs:
                continue
            prune = profile_map.get(sc.scenario_id)
            logger.info("Pass B [%d/%d] %s  profiles=%s", idx, len(selected_objs), sc.scenario_id, prune or "all")
            out = self._run_scenario_with_recovery(
                sc,
                pass_name="deep",
                config=pass_b_cfg,
                seed=seed,
                allowed_profiles=prune,
                budget=budget,
            )
            self._append_rows(candidates_jsonl, out["candidates"])
            self._append_rows(trajectories_jsonl, out["trajectories"])
            self._append_rows(scenarios_jsonl, [out["scenario_metrics"]])

            if with_clinical:
                clinical_rows = self._run_clinical_for_scenario(out["candidates"]) 
                self._append_rows(clinical_jsonl, clinical_rows)

            self.store.append_scenario_checkpoint(
                run_id,
                {
                    "scenario_id": sc.scenario_id,
                    "pass_name": "deep",
                    "status": "done",
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                },
            )

        # Build final parquet artifacts from checkpoints.
        scenario_df = self._read_jsonl(scenarios_jsonl)
        candidate_df = self._read_jsonl(candidates_jsonl)
        trajectory_df = self._read_jsonl(trajectories_jsonl)
        clinical_df = self._read_jsonl(clinical_jsonl)

        if not candidate_df.empty and not clinical_df.empty and "candidate_uid" in candidate_df.columns:
            clinical_agg = clinical_df.sort_values("clinical_score", ascending=False).drop_duplicates("candidate_uid")
            clinical_agg = clinical_agg[["candidate_uid", "clinical_score", "clinical_viability"]]
            candidate_df = candidate_df.merge(clinical_agg, on="candidate_uid", how="left")

        self.store.write_tables(
            run_id,
            scenarios_df=scenario_df,
            candidates_df=candidate_df,
            trajectories_df=trajectory_df,
            clinical_df=clinical_df,
        )

        top30 = select_presentation_shortlist(
            candidate_df[candidate_df["pass_name"] == "deep"].copy() if "pass_name" in candidate_df.columns else candidate_df,
            top_n=shortlist_n,
            max_per_target=4,
            delivery_methods=delivery_methods or DEFAULT_DELIVERY_METHODS,
        )
        self.store.write_top30(run_id, top30)

        summary = self._build_summary(run_id, scenario_df, candidate_df, top30)
        self.store.write_summary(run_id, summary)
        self.store.update_manifest(
            run_id,
            {
                "status": "completed",
                "n_scenarios": int(len(scenario_objs)),
                "n_candidates": int(len(candidate_df)),
                "n_trajectories": int(len(trajectory_df)),
                "n_shortlist": int(len(top30)),
            },
        )

        # Post-campaign analysis: ESMFold validation, MD scripts, PyMOL scripts
        post_results = self._post_campaign_analysis(run_id, top30, paths, budget=budget)

        return {
            "run_id": run_id,
            "run_dir": str(paths.run_dir),
            "n_scenarios": int(len(scenario_objs)),
            "n_candidates": int(len(candidate_df)),
            "n_shortlist": int(len(top30)),
            "selected_pass_b": int(len(selected_objs)),
            "post_analysis": post_results,
        }

    def _post_campaign_analysis(
        self,
        run_id: str,
        top_df: pd.DataFrame,
        paths: Any,
        budget: str = "high",
    ) -> Dict[str, Any]:
        """Run post-campaign analysis: ESMFold validation, MD scripts, PyMOL scripts."""
        results: Dict[str, Any] = {}
        run_dir = paths.run_dir

        # 1. ESMFold structural validation for all shortlisted candidates
        try:
            from p53cad.engine.md_validation import StructurePredictor, FreeEnergyCalculator

            predictor = StructurePredictor()
            fe_calc = FreeEnergyCalculator()
            validations = []

            for _, row in top_df.iterrows():
                seq = str(row.get("sequence", ""))
                target_label = str(row.get("target_label", "unknown"))
                muts = json.loads(str(row.get("mutations_json", "[]")))
                rescue_muts = [m for m in muts if m not in json.loads(str(row.get("targets_json", "[]")))]

                pred = predictor.predict_esmfold(seq, name=target_label)
                fe_results = [fe_calc.calculate_ddg(P53_WT, seq, m) for m in rescue_muts]

                validations.append({
                    "target_label": target_label,
                    "rescue_mutations": rescue_muts,
                    "mean_plddt": float(pred.mean_plddt),
                    "model": pred.model,
                    "n_residues": len(seq),
                    "mean_ddg": float(np.mean([fe.ddg_average for fe in fe_results])) if fe_results else 0.0,
                    "max_ddg": float(max((fe.ddg_average for fe in fe_results), default=0.0)),
                    "is_estimated": pred.pdb_string.startswith("HEADER    ESTIMATED"),
                })

            val_path = run_dir / "esmfold_validation.json"
            val_path.write_text(json.dumps(validations, indent=2), encoding="utf-8")
            results["esmfold_validations"] = len(validations)
            logger.info("ESMFold validation: %d candidates validated -> %s", len(validations), val_path)
        except Exception as exc:
            logger.warning("ESMFold validation failed: %s", exc)
            results["esmfold_error"] = str(exc)

        # 2. Generate MD validation scripts for top 5 candidates
        try:
            from p53cad.engine.md_validation import MDSimulationGenerator

            md_gen = MDSimulationGenerator()
            md_dir = run_dir / "md_scripts"
            md_dir.mkdir(exist_ok=True)
            n_md = 0

            for _, row in top_df.head(5).iterrows():
                seq = str(row.get("sequence", ""))
                target_label = str(row.get("target_label", "unknown"))
                muts = json.loads(str(row.get("mutations_json", "[]")))
                safe_name = target_label.replace("+", "_").replace(" ", "_")
                config = md_gen.generate_config(safe_name, seq, muts)

                for plat in ["openmm", "gromacs"]:
                    script = md_gen.generate_script(config, plat)
                    ext = ".py" if plat == "openmm" else ".sh"
                    script_path = md_dir / f"{safe_name}_{plat}{ext}"
                    script_path.write_text(script, encoding="utf-8")
                    n_md += 1

            results["md_scripts"] = n_md
            logger.info("MD scripts: %d files generated -> %s", n_md, md_dir)
        except Exception as exc:
            logger.warning("MD script generation failed: %s", exc)
            results["md_error"] = str(exc)

        # 3. Generate PyMOL visualization scripts for top 5 candidates
        try:
            from p53cad.viz.pymol import PyMolGenerator

            pymol_gen = PyMolGenerator()
            pymol_paths = pymol_gen.generate_campaign_pml(top_df, run_dir / "pymol", top_n=5)
            results["pymol_scripts"] = len(pymol_paths)
            logger.info("PyMOL scripts: %d files generated", len(pymol_paths))
        except Exception as exc:
            logger.warning("PyMOL script generation failed: %s", exc)
            results["pymol_error"] = str(exc)

        # 4. Physics-based validation (Tier 1 only — skip MD in auto mode)
        # Budget-aware: quick skips physics entirely, fast validates top 10 only
        if budget == "quick":
            logger.info("Skipping physics validation (budget=quick)")
        else:
            # Determine how many candidates to validate based on budget
            tier1_top_n = 10 if budget == "fast" else None  # None = all candidates
            try:
                from p53cad.engine.physics_validation import PhysicsValidationPipeline

                cancer_sequences = {}
                if "target_label" in top_df.columns:
                    for tl in top_df["target_label"].unique():
                        target_mut = str(tl).split("+")[0].strip()
                        cancer_seq = apply_mutation(P53_WT, target_mut)
                        if cancer_seq:
                            cancer_sequences[str(tl)] = cancer_seq

                pipeline = PhysicsValidationPipeline(
                    device=str(self.embedder.device) if self.embedder else "cpu",
                    cache_dir=run_dir / "esmfold_cache",
                )
                report = pipeline.validate_campaign(
                    run_id=run_id,
                    top_df=top_df,
                    output_dir=run_dir,
                    wt_sequence=P53_WT,
                    cancer_sequences=cancer_sequences,
                    tier1_top_n=tier1_top_n,
                    tier2_top_n=0,   # Skip MD in auto post-campaign mode
                    skip_md=True,
                    skip_esmfold=False,
                    skip_energy=False,
                    skip_dna=False,
                    skip_binding_sim=True,  # Expensive; opt-in via `p53cad validate`
                )
                val_path = run_dir / "physics_validation.json"
                val_path.write_text(json.dumps(report.to_dict(), indent=2), encoding="utf-8")
                results["physics_validation"] = report.n_candidates
                logger.info("Physics validation: %d candidates -> %s", report.n_candidates, val_path)
            except Exception as exc:
                logger.warning("Physics validation failed (non-fatal): %s", exc)
                results["physics_validation_error"] = str(exc)

        return results

    def rescore_candidates(
        self,
        candidates_df: pd.DataFrame,
        *,
        oracle_checkpoint: "str | Path",
        batch_log_interval: int = 25,
    ) -> pd.DataFrame:
        """Re-score all candidates with a different oracle checkpoint.

        This re-embeds every unique sequence through ESM-2 and evaluates the new
        oracle, updating the ``score``, ``score_raw``, and ``score_calibrated``
        columns in-place (preserving all other columns including uncertainty,
        OOD distance, and DMS fields which were computed during the original run).

        Parameters
        ----------
        candidates_df:
            The full candidates DataFrame from a previous campaign run.
        oracle_checkpoint:
            Path to the new oracle ``.pt`` checkpoint file to load.
        batch_log_interval:
            How often (every N unique sequences) to log progress.

        Returns
        -------
        pd.DataFrame
            A copy of *candidates_df* with updated score columns.
        """
        oracle_checkpoint = Path(oracle_checkpoint)
        logger.info("Re-scoring %d candidates with oracle: %s", len(candidates_df), oracle_checkpoint)

        # Load new oracle (temporary — does not replace self.oracle)
        new_oracle = FunctionalOracle(model_path=oracle_checkpoint, device=self.device)
        new_oracle.model.eval()

        # Ensure embedder is ready (lazy init; only loads ESM-2 once)
        if self.embedder is None:
            logger.info("Loading ESM-2 embedder for re-scoring...")
            self.embedder = ManifoldEmbedder(device=self.device)

        # Build per-sequence score cache (avoid re-embedding duplicates)
        score_cache: Dict[str, Dict[str, float]] = {}
        unique_seqs = candidates_df["sequence"].unique() if "sequence" in candidates_df.columns else []
        n_unique = len(unique_seqs)
        logger.info("Re-embedding %d unique sequences through ESM-2...", n_unique)

        for idx, seq in enumerate(unique_seqs):
            if (idx + 1) % batch_log_interval == 0 or idx == 0 or idx == n_unique - 1:
                logger.info("  Re-scoring sequence %d/%d", idx + 1, n_unique)
            try:
                with torch.no_grad():
                    word_emb = self.embedder.get_embeddings(seq)  # (1, L, D)
                    h, _logits, _probs = self.embedder.latent_forward_ascent(word_emb)[:3]
                    # h: (1, L, D) — last hidden state, same input the oracle saw during optimization
                    raw_score = float(new_oracle.model(h).squeeze(-1).mean().item())
                calibrated = self._calibrate_score(raw_score)
                score_cache[seq] = {"score_raw": raw_score, "score_calibrated": calibrated, "score": calibrated}
            except Exception as exc:
                logger.warning("Re-scoring failed for sequence [%.30s...]: %s", seq, exc)
                score_cache[seq] = {}

        # Apply updated scores to a copy of the DataFrame
        result_df = candidates_df.copy()
        for col in ("score", "score_raw", "score_calibrated"):
            if col not in result_df.columns:
                result_df[col] = float("nan")

        for seq, scores in score_cache.items():
            if not scores:
                continue
            mask = result_df["sequence"] == seq
            for col, val in scores.items():
                result_df.loc[mask, col] = val

        logger.info(
            "Re-scoring complete. score range: [%.3f, %.3f]",
            result_df["score"].min(),
            result_df["score"].max(),
        )
        return result_df

    def report_run(
        self,
        run_id: str,
        shortlist_n: int = 30,
        oracle_checkpoint: "str | Path | None" = None,
    ) -> Dict[str, Any]:
        bundle = self.store.load_run_bundle(run_id)
        candidate_df = bundle["candidates"]
        if candidate_df.empty:
            raise RuntimeError(f"No candidates found for run {run_id}")

        if oracle_checkpoint is not None:
            logger.info("Oracle re-scoring requested for report_run (checkpoint: %s)", oracle_checkpoint)
            candidate_df = self.rescore_candidates(candidate_df, oracle_checkpoint=oracle_checkpoint)

        top30 = select_presentation_shortlist(
            candidate_df[candidate_df["pass_name"] == "deep"].copy() if "pass_name" in candidate_df.columns else candidate_df,
            top_n=shortlist_n,
            max_per_target=4,
            delivery_methods=DEFAULT_DELIVERY_METHODS,
        )
        self.store.write_top30(run_id, top30)
        summary = self._build_summary(run_id, bundle["scenarios"], candidate_df, top30)
        self.store.write_summary(run_id, summary)
        self.store.update_manifest(run_id, {"status": "reported", "n_shortlist": int(len(top30))})
        return {
            "run_id": run_id,
            "n_candidates": int(len(candidate_df)),
            "n_shortlist": int(len(top30)),
            "run_dir": str(bundle["run_dir"]),
        }

    def _pass_a_config(self, budget: str) -> Dict[str, Any]:
        if budget == "quick":
            return {"steps": 25, "restarts": 1, "repeats": 1, "mc_samples": 1}
        if budget == "fast":
            return {"steps": 40, "restarts": 1, "repeats": 1, "mc_samples": 2}
        if budget == "medium":
            return {"steps": 60, "restarts": 1, "repeats": 1, "mc_samples": 3}
        return {"steps": 80, "restarts": 1, "repeats": 1, "mc_samples": 4}

    def _pass_b_config(self, budget: str) -> Dict[str, Any]:
        if budget == "quick":
            return {"steps": 60, "restarts": 1, "repeats": 1, "mc_samples": 2}
        if budget == "fast":
            return {"steps": 120, "restarts": 1, "repeats": 1, "mc_samples": 4}
        if budget == "medium":
            return {"steps": 200, "restarts": 2, "repeats": 2, "mc_samples": 8}
        # High: reduced from 3×3=9 to 2×2=4 trials per profile (55% fewer trials)
        return {"steps": 280, "restarts": 2, "repeats": 2, "mc_samples": 12}

    def _early_stop_patience(self, budget: str) -> int:
        """Number of check-intervals (each = 5 steps) without improvement before early stop."""
        if budget == "quick":
            return 3   # 15 gradient steps
        if budget == "fast":
            return 4   # 20 gradient steps
        return 6       # 30 gradient steps (medium/high)

    def _run_scenario(
        self,
        scenario: ScenarioRuntime,
        *,
        pass_name: str,
        config: Dict[str, int],
        seed: int,
        allowed_profiles: Optional[List[str]] = None,
        budget: str = "high",
    ) -> Dict[str, Any]:
        assert self.embedder is not None
        assert self.oracle is not None

        target_seq = self._build_target_sequence(scenario.targets)
        emb_target_ref = self._get_cached_embedding(target_seq)
        emb_wt = self._get_cached_embedding(P53_WT)
        with torch.no_grad():
            z_target_ref, _, _ = self.embedder.latent_forward_ascent(emb_target_ref)
            pooled_target_ref = z_target_ref.mean(dim=1)
            pooled_target_ref = self._align_oracle_dim(pooled_target_ref)

        wt_aa_tensor = self._wt_aa_tensor(device=emb_target_ref.device)

        min_identity = self._delivery_identity_floor(scenario.delivery_method)
        min_stability = -0.2
        min_binding = 5.0
        lock_positions: set[int] = set()
        for target in scenario.targets:
            pos = self._mutation_pos(target)
            if pos is not None:
                lock_positions.add(int(pos))
        lock_positions.update([175, 248, 273])
        lock_positions_sorted = sorted(lock_positions)
        max_mutations = int(len(P53_WT) * (100 - min_identity) / 100)
        locked_indices = [p - 1 for p in lock_positions_sorted]

        baseline = self._get_cached_baseline(target_seq, pooled_target_ref, mc_samples=config["mc_samples"])
        dms_quality = self._get_dms_quality_tensor(emb_target_ref.device)

        # Build pre-computed constants once per scenario for efficiency
        precomputed = build_precomputed_constants(
            device=emb_target_ref.device,
            wt_contacts=self._wt_contacts,
            ca_coords=self._ca_coords,
            locked_indices=locked_indices,
        )

        # Conditional rescue scores disabled - CUDA compatibility issues
        # The ESM-2 forward pass causes illegal instruction errors on some GPUs
        # Can re-enable once GPU/driver compatibility is resolved
        candidate_positions = [p for p in range(94, 293) if (p - 1) not in locked_indices]
        conditional_scores: Optional[Dict[int, Dict[str, float]]] = None

        candidates: List[Dict[str, Any]] = []
        trajectories: List[Dict[str, Any]] = []

        profiles = WEIGHT_PROFILES
        if allowed_profiles is not None:
            allowed_set = set(allowed_profiles)
            profiles = [p for p in WEIGHT_PROFILES if p["name"] in allowed_set]

        trial_counter = 0
        for repeat_idx in range(int(config["repeats"])):
            for profile in profiles:
                for restart_idx in range(int(config["restarts"])):
                    trial_counter += 1
                    trial_seed = int(seed + (repeat_idx * 1000) + (restart_idx * 100) + trial_counter)
                    cand, traj_rows = self._run_single_trial(
                        scenario=scenario,
                        target_seq=target_seq,
                        emb_target_ref=emb_target_ref,
                        emb_wt=emb_wt,
                        pooled_target_ref=pooled_target_ref,
                        wt_aa_tensor=wt_aa_tensor,
                        profile=profile,
                        pass_name=pass_name,
                        repeat_idx=repeat_idx + 1,
                        restart_idx=restart_idx + 1,
                        trial_idx=trial_counter,
                        trial_seed=trial_seed,
                        n_steps=int(config["steps"]),
                        mc_samples=int(config["mc_samples"]),
                        min_identity=min_identity,
                        min_stability=min_stability,
                        min_binding=min_binding,
                        max_mutations=max_mutations,
                        locked_indices=locked_indices,
                        baseline_score=float(baseline["score"]),
                        dms_quality=dms_quality,
                        conditional_scores=conditional_scores,
                        early_stop_patience=self._early_stop_patience(budget),
                        precomputed=precomputed,
                    )
                    candidates.append(cand)
                    trajectories.extend(traj_rows)

                    # Log progress every 5 trials
                    if trial_counter % 5 == 0:
                        logger.info(
                            "Progress: %d/%d trials, %d candidates, best=%.3f",
                            trial_counter,
                            int(config["repeats"]) * len(profiles) * int(config["restarts"]),
                            len(candidates),
                            max((c.get("score", -999) for c in candidates), default=-999),
                        )

        # Run one autoregressive trial per scenario for diversity
        try:
            ar_cand, ar_traj = self._run_autoregressive_trial(
                scenario=scenario,
                target_seq=target_seq,
                locked_indices=locked_indices,
                max_mutations=max_mutations,
                pass_name=pass_name,
                trial_idx=trial_counter + 1,
                trial_seed=seed + 9999,
                baseline_score=float(baseline["score"]),
            )
            candidates.append(ar_cand)
            trajectories.extend(ar_traj)
        except Exception as exc:
            logger.warning("Autoregressive trial failed: %s", exc)

        # Multi-objective Pareto ranking across all trial candidates
        try:
            from p53cad.engine.pareto import ParetoFront, ParetoSolution
            pareto = ParetoFront()
            for c in candidates:
                pareto.add(ParetoSolution(
                    candidate_uid=str(c.get("candidate_uid", "")),
                    objectives={
                        "oracle_score": float(c.get("score", 0.0)),
                        "pll": float(c.get("stability", 0.0)),
                        "dms_quality": -float(c.get("rescue_dms_mean", 0.0)),
                        "identity": float(c.get("identity", 100.0)),
                        "stability": float(c.get("stability", 0.0)),
                    },
                ))
            ranked = pareto.solutions
            uid_to_rank = {s.candidate_uid: s.pareto_rank for s in ranked}
            for c in candidates:
                c["pareto_rank"] = uid_to_rank.get(str(c.get("candidate_uid", "")), 999)
        except Exception as exc:
            logger.warning("Pareto ranking failed: %s", exc)
            for c in candidates:
                c["pareto_rank"] = 0

        # Population sharing: add high-scoring candidates to shared pool
        try:
            for c in candidates:
                muts_json = c.get("mutations_json", "[]")
                muts = json.loads(muts_json) if isinstance(muts_json, str) else muts_json
                mut_str = ",".join(sorted(muts)) if muts else ""
                score = c.get("score", -999)
                if score > -999:
                    self._population_share.add_candidate(mut_str, score)
            logger.debug("Population share pool size: %d", self._population_share.size)
        except Exception as exc:
            logger.warning("Population sharing failed: %s", exc)

        best_score = max((float(c.get("score", -np.inf)) for c in candidates), default=-np.inf)
        metrics = {
            "scenario_id": scenario.scenario_id,
            "target_label": scenario.target_label,
            "targets": json.dumps(scenario.targets),
            "delivery_method": scenario.delivery_method,
            "pass_name": pass_name,
            "n_candidates": int(len(candidates)),
            "best_score": float(best_score),
            "baseline_score": float(baseline["score"]),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }

        logger.info(
            "Scenario %s complete: %d candidates, best_score=%.3f, baseline=%.3f",
            scenario.scenario_id, len(candidates), best_score, baseline["score"],
        )

        return {"candidates": candidates, "trajectories": trajectories, "scenario_metrics": metrics}

    def _run_single_trial(
        self,
        *,
        scenario: ScenarioRuntime,
        target_seq: str,
        emb_target_ref: torch.Tensor,
        emb_wt: torch.Tensor,
        pooled_target_ref: torch.Tensor,
        wt_aa_tensor: torch.Tensor,
        profile: Dict[str, Any],
        pass_name: str,
        repeat_idx: int,
        restart_idx: int,
        trial_idx: int,
        trial_seed: int,
        n_steps: int,
        mc_samples: int,
        min_identity: float,
        min_stability: float,
        min_binding: float,
        max_mutations: int,
        locked_indices: List[int],
        baseline_score: float,
        dms_quality: Optional[torch.Tensor] = None,
        conditional_scores: Optional[Dict[int, Dict[str, float]]] = None,
        early_stop_patience: int = 6,
        precomputed: Optional[PrecomputedConstants] = None,
    ) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        assert self.embedder is not None
        assert self.oracle is not None

        torch.manual_seed(trial_seed)
        np.random.seed(trial_seed)

        emb = emb_target_ref.clone().detach().requires_grad_(True)
        with torch.no_grad():
            emb.data += torch.randn_like(emb) * 0.05

        is_windows_cuda = (emb.device.type == "cuda" and platform.system() == "Windows")
        default_lr = 0.015 if is_windows_cuda else 0.03
        trial_lr = float(os.environ.get("P53CAD_TRIAL_LR", str(default_lr)))
        optimizer = torch.optim.Adam([emb], lr=trial_lr)

        # Pre-compute WT hidden states for contact preservation penalty
        wt_hidden = self._get_wt_hidden(emb.device)

        # Pre-compute contact neighbor indices for locked (mutating) positions
        # These are the 0-indexed positions whose structural neighbours we monitor
        contact_neighbor_indices: List[List[int]] = []
        has_contacts = bool(self._wt_contacts and wt_hidden is not None)
        if has_contacts:
            for li in locked_indices:
                res_id = li + 1  # 1-indexed
                neighbors = self._wt_contacts.get(res_id, [])
                # Convert neighbor residue IDs to 0-indexed, filter to valid range
                ni = [r - 1 for r, _ in neighbors if 0 <= r - 1 < len(P53_WT)]
                contact_neighbor_indices.append(ni)

        # Pre-compute locked position set (used by epistasis, re-embedding, and autoregressive)
        locked_set = set(locked_indices)

        # Pre-compute pairwise Cα distances for epistasis proximity penalty
        # We track pairs of all mutatable (non-locked) positions that are close
        _epistasis_pairs: List[Tuple[int, int, float]] = []  # (idx_i, idx_j, distance)
        if self._ca_coords:
            all_pos = list(range(len(P53_WT)))
            for i_idx in range(len(all_pos)):
                if i_idx in locked_set:
                    continue
                ri = i_idx + 1  # 1-indexed
                if ri not in self._ca_coords:
                    continue
                for j_idx in range(i_idx + 1, len(all_pos)):
                    if j_idx in locked_set:
                        continue
                    rj = j_idx + 1
                    if rj not in self._ca_coords:
                        continue
                    diff = self._ca_coords[ri] - self._ca_coords[rj]
                    dist = float(np.sqrt(np.sum(diff ** 2)))
                    if dist < 10.0:
                        _epistasis_pairs.append((i_idx, j_idx, dist))
        cached_epistasis_loss = torch.zeros(1, device=emb.device)

        # Pre-compute epistasis pair indices as tensors to avoid .item() sync stalls in the loop
        if _epistasis_pairs:
            _ep_i = torch.tensor([p[0] for p in _epistasis_pairs], dtype=torch.long, device=emb.device)
            _ep_j = torch.tensor([p[1] for p in _epistasis_pairs], dtype=torch.long, device=emb.device)
            _ep_decay = torch.tensor([math.exp(-p[2] / 5.0) for p in _epistasis_pairs], dtype=torch.float32, device=emb.device)
        else:
            _ep_i = _ep_j = _ep_decay = None

        # Pre-compute position-weighted pooling kernel for mutation-neighborhood oracle.
        # Gaussian centered on each locked_index with σ=10 residues — upweights positions
        # near cancer sites so the oracle sees a stronger signal from rescue mutations.
        seq_len = emb.size(1)
        pool_weights = torch.ones(seq_len, device=emb.device)
        for li in locked_indices:
            dists = (torch.arange(seq_len, device=emb.device).float() - li)
            pool_weights += torch.exp(-dists**2 / 200.0)  # σ=10 → 2σ²=200
        pool_weights = pool_weights / pool_weights.sum()  # (L,)

        # Re-embedding interval: project embedding back to protein manifold every N steps.
        # Quicker budgets re-embed more often (tighter constraint) to prevent drift.
        reembed_interval = {25: 5, 60: 8, 80: 10, 120: 10, 200: 10, 280: 12}.get(n_steps, 10)

        trajectory: List[Dict[str, Any]] = []
        best_valid_state: Optional[Dict[str, Any]] = None
        best_valid_score = -float("inf")
        checks_without_improve = 0
        cached_uncertainty = 0.0

        for step_idx in range(1, int(n_steps) + 1):
            if not torch.isfinite(emb).all():
                logger.warning(
                    "Non-finite embedding detected in trial %d at step %d; sanitizing.",
                    trial_idx,
                    step_idx,
                )
                with torch.no_grad():
                    emb.data = torch.nan_to_num(emb.data, nan=0.0, posinf=1.0, neginf=-1.0)
            optimizer.zero_grad()
            z, logits, _ = self.embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            pooled = self._align_oracle_dim(pooled)

            # Attention oracle receives full per-position embeddings;
            # legacy MLP oracle receives mean-pooled embeddings.
            _uses_attention = isinstance(self.oracle.model, AttentionPoolingNet)
            if _uses_attention:
                z_oracle = self._align_oracle_dim(z)
                raw_score_t = self.oracle.model(z_oracle).squeeze(-1)
            else:
                raw_score_t = self.oracle.model(pooled).squeeze(-1)

            # Mutation-neighborhood oracle: position-weighted pooling that amplifies
            # the signal near cancer sites instead of averaging uniformly over 393 positions.
            pooled_local = (z.squeeze(0) * pool_weights.unsqueeze(1)).sum(dim=0, keepdim=True)  # (1, D)
            pooled_local = self._align_oracle_dim(pooled_local)
            if _uses_attention:
                # Attention oracle already does position weighting natively
                local_score_t = raw_score_t
            else:
                local_score_t = self.oracle.model(pooled_local).squeeze(-1)
            local_score_term = -2.0 * local_score_t

            probs_full = torch.softmax(logits, dim=-1)
            logits_aa = logits[:, :, AA_IDS]
            log_probs = F.log_softmax(logits_aa, dim=-1)
            stability_t = log_probs.max(dim=-1).values.mean(dim=-1)
            dna_force_t = self._batch_dna_force(z, probs_full)
            hydro_packing_t = self._batch_hydrophobic_packing(probs_full)
            ood_distance_t = torch.norm(pooled - pooled_target_ref, p=2, dim=-1)

            probs_aa = torch.softmax(logits_aa, dim=-1)

            # PLL: expected pseudo-log-likelihood — how natural this sequence looks to ESM-2
            expected_pll = (probs_aa * log_probs).sum(dim=-1).mean(dim=-1)  # [batch]

            wt_probs = probs_aa[:, torch.arange(len(P53_WT), device=emb.device), wt_aa_tensor]
            expected_mutations = (1.0 - wt_probs).sum(dim=-1)
            expected_identity = 100.0 * (1.0 - expected_mutations / float(len(P53_WT)))

            score_term = -raw_score_t * float(profile["function"])
            stability_term = -float(profile["stability"]) * stability_t
            binding_term = -float(profile["binding"]) * dna_force_t
            hydro_term = -3.0 * hydro_packing_t
            ood_penalty = self.ood_loss_weight * F.relu(ood_distance_t - self.ood_radius)
            mutation_penalty = 200.0 * F.relu(expected_mutations - max_mutations)
            identity_penalty = 1000.0 * F.relu(min_identity - expected_identity)
            stability_penalty = 100.0 * F.relu(min_stability - stability_t)
            binding_penalty = 80.0 * F.relu(min_binding - dna_force_t)
            if locked_indices:
                lock_penalty = 800.0 * (emb[:, locked_indices, :] - emb_target_ref[:, locked_indices, :]).pow(2).mean(dim=(1, 2))
            else:
                lock_penalty = torch.zeros_like(raw_score_t)
            if locked_indices:
                l1_mask = torch.ones(emb.size(1), device=emb.device, dtype=torch.bool)
                for li in locked_indices:
                    l1_mask[li] = False
                l1_penalty = 40.0 * (emb[:, l1_mask, :] - emb_wt[:, l1_mask, :]).abs().mean(dim=(1, 2))
            else:
                l1_penalty = 40.0 * (emb - emb_wt).abs().mean(dim=(1, 2))

            # DMS quality penalty: steer optimizer toward individually functional rescue mutations.
            # dms_quality has raw Nutlin-3 Z-scores: negative = functional (good), positive = LoF (bad).
            # Expected DMS = sum_over_AAs(P(aa|pos) * Z(pos, aa)), weighted by mutation probability.
            if dms_quality is not None:
                expected_dms = (probs_aa * dms_quality.unsqueeze(0)).sum(dim=-1)  # (batch, seq_len)
                mut_prob = 1.0 - wt_probs  # probability of mutation at each position
                weighted_dms = mut_prob * expected_dms
                # Only apply to non-locked positions (rescue sites, not cancer targets)
                dms_mask = l1_mask if locked_indices else torch.ones(emb.size(1), device=emb.device, dtype=torch.bool)
                dms_penalty = 10.0 * weighted_dms[:, dms_mask].mean(dim=-1)
            else:
                dms_penalty = torch.zeros_like(raw_score_t)

            # Epistasis penalty: computed every 10 steps to limit overhead.
            # Uses detached probs to avoid stale-graph errors on subsequent backward().
            if _ep_i is not None and step_idx % 10 == 1:
                probs_aa_d = probs_aa.detach()

                # A. Structural proximity: vectorized — no per-pair .item() sync stalls
                wt_aa_i = wt_aa_tensor[_ep_i]  # (n_pairs,)
                wt_aa_j = wt_aa_tensor[_ep_j]  # (n_pairs,)
                mp_i_t = 1.0 - probs_aa_d[0, _ep_i, wt_aa_i]  # (n_pairs,)
                mp_j_t = 1.0 - probs_aa_d[0, _ep_j, wt_aa_j]  # (n_pairs,)
                proximity_val = float((mp_i_t * mp_j_t * _ep_decay).mean().item())

                # B. Attention coupling: vectorized after single forward pass
                attn_val = 0.0
                try:
                    with torch.no_grad():
                        _, _, _, attns = self.embedder.latent_forward_ascent(
                            emb.detach(), return_attention=True
                        )
                    # attns: tuple of (1, heads, L, L) per layer — average all
                    attn_stack = torch.stack([a.mean(dim=1) for a in attns]).mean(dim=0)  # (1, L, L)
                    mutual_attn_t = (attn_stack[0, _ep_i, _ep_j] + attn_stack[0, _ep_j, _ep_i]) / 2.0  # (n_pairs,)
                    attn_val = float((mp_i_t * mp_j_t * mutual_attn_t).mean().item())
                except Exception:
                    pass

                cached_epistasis_loss = torch.tensor(
                    2.0 * (proximity_val + attn_val), device=emb.device
                )

            epistasis_penalty = cached_epistasis_loss.expand_as(raw_score_t)

            # PLL term: maximize sequence naturalness (negative because we minimize loss)
            pll_term = -3.0 * expected_pll

            # Cancer-site PLL: maximize log-probability specifically at locked (cancer) positions.
            # This directly measures "does the current sequence context make ESM-2 more confident
            # about the cancer residue?" — provides focused gradient that global PLL dilutes.
            if locked_indices:
                cancer_site_lp = log_probs[:, locked_indices, :].max(dim=-1).values.mean(dim=-1)
                cancer_pll_term = -5.0 * cancer_site_lp
            else:
                cancer_pll_term = torch.zeros_like(raw_score_t)

            # Contact preservation: penalize when mutated positions disrupt structural neighbors
            if has_contacts and wt_hidden is not None:
                cos_sims = []
                for ni_list in contact_neighbor_indices:
                    if not ni_list:
                        continue
                    ni_t = torch.tensor(ni_list, device=emb.device)
                    cur_h = z[:, ni_t, :]     # (batch, n_neighbors, D)
                    wt_h = wt_hidden[:, ni_t, :]  # (1, n_neighbors, D)
                    # Cosine similarity per neighbor, averaged
                    sim = F.cosine_similarity(cur_h, wt_h, dim=-1).mean(dim=-1)  # (batch,)
                    cos_sims.append(sim)
                if cos_sims:
                    mean_cos = torch.stack(cos_sims).mean(dim=0)  # (batch,)
                    contact_penalty = 5.0 * (1.0 - mean_cos)
                else:
                    contact_penalty = torch.zeros_like(raw_score_t)
            else:
                contact_penalty = torch.zeros_like(raw_score_t)

            # Conditional DMS term: reward mutations that ESM-2 predicts are favorable
            # in the cancer context (not just individually functional on WT).
            if conditional_scores:
                cond_log_prob = torch.zeros_like(raw_score_t)
                for pos, aa_scores in conditional_scores.items():
                    idx = pos - 1
                    if 0 <= idx < probs_aa.shape[1]:
                        for aa_char, lp in aa_scores.items():
                            tok_id = self.embedder.tokenizer.convert_tokens_to_ids(aa_char)
                            if tok_id in AA_IDS:
                                aa_idx = AA_IDS.index(tok_id)
                                cond_log_prob += probs_aa[:, idx, aa_idx] * lp
                cond_rescue_term = -2.0 * cond_log_prob / max(len(conditional_scores), 1)
            else:
                cond_rescue_term = torch.zeros_like(raw_score_t)

            # DBD structural confidence: use ESM-2's per-position confidence as a
            # fast proxy for structural integrity (no extra forward pass needed).
            dbd_range = list(range(93, 292))
            dbd_confidence = log_probs[:, dbd_range, :].max(dim=-1).values.mean(dim=-1)
            structure_term = -3.0 * dbd_confidence

            # NEW: Tier 1 structural loss terms using precomputed constants
            mut_conf_penalty = torch.zeros_like(raw_score_t)
            interface_floor_penalty = torch.zeros_like(raw_score_t)
            crater_penalty = torch.zeros_like(raw_score_t)
            contact_reg_penalty = torch.zeros_like(raw_score_t)

            if precomputed is not None:
                mut_conf_penalty = compute_mutation_confidence_penalty(
                    probs_aa=probs_aa,
                    wt_aa_tensor=precomputed.wt_aa_tensor,
                    helix_mask=precomputed.helix_mask,
                    sheet_mask=precomputed.sheet_mask,
                    locked_indices=locked_indices,
                    device=emb.device,
                )

                interface_floor_penalty = compute_interface_floor_penalty(
                    probs_aa=probs_aa,
                    wt_aa_tensor=precomputed.wt_aa_tensor,
                    interface_mask=precomputed.interface_mask,
                    device=emb.device,
                )

                crater_penalty = compute_crater_penalty(
                    probs_aa=probs_aa,
                    wt_aa_tensor=precomputed.wt_aa_tensor,
                    crater_kernel=precomputed.crater_kernel,
                    device=emb.device,
                )

                if step_idx % 10 == 1 and wt_hidden is not None:
                    contact_reg_penalty = compute_contact_map_regularization(
                        z=z,
                        wt_hidden=wt_hidden,
                        device=emb.device,
                    )

            loss_vec = (
                score_term
                + local_score_term
                + stability_term
                + binding_term
                + hydro_term
                + pll_term
                + cancer_pll_term
                + contact_penalty
                + epistasis_penalty
                + ood_penalty
                + mutation_penalty
                + identity_penalty
                + stability_penalty
                + binding_penalty
                + lock_penalty
                + l1_penalty
                + dms_penalty
                + cond_rescue_term
                + structure_term
                + mut_conf_penalty
                + interface_floor_penalty
                + crater_penalty
                + contact_reg_penalty
            )

            loss = loss_vec.mean()
            if not torch.isfinite(loss):
                logger.warning(
                    "Non-finite loss in trial %d at step %d; skipping step.",
                    trial_idx,
                    step_idx,
                )
                continue
            loss.backward()
            torch.nn.utils.clip_grad_norm_([emb], max_norm=1.0)
            optimizer.step()
            with torch.no_grad():
                emb.data = torch.nan_to_num(emb.data, nan=0.0, posinf=1.0, neginf=-1.0)

            # --- Periodic re-embedding: project back to protein manifold ---
            # After gradient steps the continuous embedding drifts off the manifold
            # of real ESM-2 token embeddings. Decode → hard-cap mutations → re-embed
            # keeps the optimizer in discrete protein space (projected gradient descent).
            if step_idx % reembed_interval == 0 and step_idx < n_steps:
                with torch.no_grad():
                    # 1. Decode current embedding to discrete sequence
                    projected_seq = self._decode_sequence(logits_aa)

                    # 2. Hard mutation cap: if too many mutations, keep only top-K
                    proj_muts = [
                        j for j in range(len(P53_WT))
                        if P53_WT[j] != projected_seq[j] and j not in locked_set
                    ]
                    n_total_muts = len(proj_muts) + len(locked_indices)
                    if n_total_muts > max_mutations:
                        # Rank non-locked mutations by ESM-2 confidence
                        mut_confidences = []
                        for j in proj_muts:
                            aa_idx = torch.argmax(logits_aa[0, j]).item()
                            conf = float(log_probs[0, j, aa_idx].item())
                            mut_confidences.append((j, conf))
                        mut_confidences.sort(key=lambda x: x[1], reverse=True)
                        budget = max(max_mutations - len(locked_indices), 0)
                        keep_positions = set(j for j, _ in mut_confidences[:budget])
                        # Reset low-confidence positions to WT
                        projected_list = list(projected_seq)
                        for j in proj_muts:
                            if j not in keep_positions:
                                projected_list[j] = P53_WT[j]
                        projected_seq = "".join(projected_list)

                    # 3. Re-embed the cleaned sequence through ESM-2
                    new_emb = self.embedder.get_embeddings(projected_seq).detach()

                    # 4. Preserve locked positions from original target embedding
                    for li in locked_indices:
                        new_emb[:, li, :] = emb_target_ref[:, li, :]

                    # 5. Replace optimized embedding and reset optimizer
                    n_after = sum(1 for j in range(len(P53_WT)) if P53_WT[j] != projected_seq[j])
                    seq_id = 100.0 * (1.0 - n_after / float(len(P53_WT)))
                    logger.debug(
                        "Re-embed step %d: %d muts -> %d after cap (identity %.1f%%)",
                        step_idx, n_total_muts, n_after, seq_id,
                    )
                    emb.data.copy_(new_emb)
                    emb.requires_grad_(True)
                    optimizer = torch.optim.Adam([emb], lr=trial_lr)

            if step_idx % 5 != 0 and step_idx != n_steps:
                continue

            with torch.no_grad():
                raw_score = float(raw_score_t[0].item())
                ood_distance = float(ood_distance_t[0].item())
                if step_idx == 5 or step_idx == n_steps or step_idx % 40 == 0:
                    # Use fewer MC samples for intermediate checks; full count only at final step
                    intermediate_mc = min(mc_samples, 4)
                    unc_samples = mc_samples if step_idx == n_steps else intermediate_mc
                    cached_uncertainty = self._estimate_uncertainty(pooled, mc_samples=unc_samples, z_full=z)

                score_bundle = self._build_trust_adjusted_score(
                    raw_score=raw_score,
                    pooled=pooled,
                    pooled_ref=pooled_target_ref,
                    cached_uncertainty=cached_uncertainty,
                )

                current_seq = self._decode_sequence(logits_aa)
                muts = [f"{P53_WT[j]}{j+1}{current_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != current_seq[j]]
                mut_positions = [self._mutation_pos(m) for m in muts]
                mut_positions = [int(p) for p in mut_positions if p is not None]
                seq_identity = 100.0 * (1.0 - len(muts) / float(len(P53_WT)))

                state = {
                    "step": int(step_idx),
                    "score": float(score_bundle["score_adjusted"]),
                    "score_raw": float(score_bundle["score_raw"]),
                    "score_calibrated": float(score_bundle["score_calibrated"]),
                    "stability": float(stability_t[0].item()),
                    "binding": float(dna_force_t[0].item()),
                    "identity": float(seq_identity),
                    "n_mutations": int(len(muts)),
                    "mutations_json": json.dumps(muts),
                    "mut_positions_json": json.dumps(mut_positions),
                    "sequence": current_seq,
                    "uncertainty": float(score_bundle["uncertainty"]),
                    "ood_distance": float(score_bundle["ood_distance"]),
                    "lx": float(pooled[0, 0].item()),
                    "ly": float(pooled[0, 1].item()),
                    "lz": float(pooled[0, 2].item()) if pooled.shape[-1] > 2 else float(score_bundle["score_adjusted"]),
                    "loss_total": float(loss_vec[0].item()),
                    "loss_score_term": float(score_term[0].item()),
                    "loss_stability_term": float(stability_term[0].item()),
                    "loss_binding_term": float(binding_term[0].item()),
                    "loss_hydrophobic_term": float(hydro_term[0].item()),
                    "loss_ood_penalty": float(ood_penalty[0].item()),
                    "loss_mutation_penalty": float(mutation_penalty[0].item()),
                    "loss_identity_penalty": float(identity_penalty[0].item()),
                    "loss_stability_penalty": float(stability_penalty[0].item()),
                    "loss_binding_penalty": float(binding_penalty[0].item()),
                    "loss_lock_penalty": float(lock_penalty[0].item()),
                    "loss_l1_penalty": float(l1_penalty[0].item()),
                    "loss_pll_term": float(pll_term[0].item()),
                    "loss_contact_penalty": float(contact_penalty[0].item()),
                    "loss_epistasis_penalty": float(epistasis_penalty[0].item()),
                    "loss_dms_penalty": float(dms_penalty[0].item()),
                    "loss_cancer_pll_term": float(cancer_pll_term[0].item()),
                    "loss_local_score_term": float(local_score_term[0].item()),
                    "loss_cond_rescue_term": float(cond_rescue_term[0].item()),
                    "loss_structure_term": float(structure_term[0].item()),
                    "loss_mut_conf_penalty": float(mut_conf_penalty[0].item()) if mut_conf_penalty.dim() > 0 else float(mut_conf_penalty.item()),
                    "loss_interface_floor_penalty": float(interface_floor_penalty[0].item()) if interface_floor_penalty.dim() > 0 else float(interface_floor_penalty.item()),
                    "loss_crater_penalty": float(crater_penalty[0].item()) if crater_penalty.dim() > 0 else float(crater_penalty.item()),
                    "loss_contact_reg_penalty": float(contact_reg_penalty[0].item()) if contact_reg_penalty.dim() > 0 else float(contact_reg_penalty.item()),
                }
                trajectory.append(state)

                is_valid = (
                    seq_identity >= min_identity
                    and float(stability_t[0].item()) >= min_stability
                    and float(dna_force_t[0].item()) >= min_binding
                )

                improved = False
                if is_valid and state["score"] > (best_valid_score + 1e-3):
                    best_valid_state = dict(state)
                    best_valid_score = float(state["score"])
                    improved = True

                if improved:
                    checks_without_improve = 0
                else:
                    checks_without_improve += 1
                # Early stop: min step threshold scales with patience (quick=15, fast=20, med/high=30)
                min_step_for_stop = early_stop_patience * 5
                if step_idx >= min_step_for_stop and checks_without_improve >= early_stop_patience:
                    logger.debug("Early stopping at step %d (no improvement for %d checks)", step_idx, early_stop_patience)
                    break

        final_state = best_valid_state if best_valid_state is not None else trajectory[-1]
        final_muts = json.loads(final_state.get("mutations_json", "[]"))
        final_positions = json.loads(final_state.get("mut_positions_json", "[]"))

        # Compute rescue DMS quality: mean raw Z-score of non-target mutations
        target_set = set(str(t).strip().upper() for t in scenario.targets)
        rescue_muts = [m for m in final_muts if m not in target_set]
        rescue_dms_scores = []
        n_functional_rescues = 0
        for rm in rescue_muts:
            parsed = parse_single_mutation(rm)
            if parsed is not None:
                _, pos, var_aa = parsed
                z = self._dms_lookup.get((pos, var_aa))
                if z is not None:
                    rescue_dms_scores.append(z)
                    if z < 0:
                        n_functional_rescues += 1
        rescue_dms_mean = float(np.mean(rescue_dms_scores)) if rescue_dms_scores else 0.0

        candidate_uid = (
            f"{scenario.scenario_id}|{pass_name}|rep{repeat_idx}|rst{restart_idx}|trial{trial_idx}"
        )
        candidate = {
            "candidate_uid": candidate_uid,
            "scenario_id": scenario.scenario_id,
            "target_label": scenario.target_label,
            "targets_json": json.dumps(scenario.targets),
            "delivery_method": scenario.delivery_method,
            "pass_name": pass_name,
            "repeat_idx": int(repeat_idx),
            "restart_idx": int(restart_idx),
            "trial_idx": int(trial_idx),
            "profile": profile["name"],
            "sequence": final_state["sequence"],
            "score": float(final_state["score"]),
            "score_raw": float(final_state.get("score_raw", final_state["score"])),
            "score_calibrated": float(final_state.get("score_calibrated", final_state["score"])),
            "score_gain_vs_target": float(final_state["score"] - baseline_score),
            "stability": float(final_state["stability"]),
            "binding": float(final_state["binding"]),
            "identity": float(final_state["identity"]),
            "n_mutations": int(final_state["n_mutations"]),
            "mutations_json": json.dumps(final_muts),
            "mut_positions_json": json.dumps(final_positions),
            "uncertainty": float(final_state.get("uncertainty", 0.0)),
            "ood_distance": float(final_state.get("ood_distance", 0.0)),
            "rescue_dms_mean": float(rescue_dms_mean),
            "n_functional_rescues": int(n_functional_rescues),
            "n_rescue_mutations": int(len(rescue_muts)),
            "meets_constraints": bool(
                float(final_state["identity"]) >= min_identity
                and float(final_state["stability"]) >= min_stability
                and float(final_state["binding"]) >= min_binding
            ),
            "selection_reason": "trust_adjusted_rank",
            "receptor_source": "n/a",
            "docking_backend": "n/a",
            "md_status": "n/a",
        }

        traj_rows = []
        for row in trajectory:
            payload = dict(row)
            payload["candidate_uid"] = candidate_uid
            payload["scenario_id"] = scenario.scenario_id
            payload["target_label"] = scenario.target_label
            payload["delivery_method"] = scenario.delivery_method
            payload["pass_name"] = pass_name
            payload["profile"] = profile["name"]
            traj_rows.append(payload)

        return candidate, traj_rows

    def _run_autoregressive_trial(
        self,
        *,
        scenario: ScenarioRuntime,
        target_seq: str,
        locked_indices: List[int],
        max_mutations: int,
        pass_name: str,
        trial_idx: int,
        trial_seed: int,
        baseline_score: float,
        n_candidates: int = 50,
        top_k: int = 3,
        temperature: float = 0.8,
    ) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        """Gibbs-like autoregressive sampling: propose mutations one position at a time.

        1. Start from the cancer sequence.
        2. Rank positions by ESM-2 attention (which positions the model focuses on).
        3. For top-N candidate positions (not locked):
           a. Mask the position in the current sequence.
           b. Get ESM-2's distribution over 20 AAs.
           c. Try top-K candidates; accept if oracle score improves.
        4. Iterate passes until convergence or max_mutations reached.
        """
        assert self.embedder is not None
        assert self.oracle is not None

        torch.manual_seed(trial_seed)
        np.random.seed(trial_seed)

        locked_set = set(locked_indices)
        current_seq = list(target_seq)
        best_seq = list(target_seq)
        best_score = baseline_score

        trajectory: List[Dict[str, Any]] = []
        mutations_applied: List[str] = []

        # Get attention weights to rank positions by importance
        emb = self.embedder.get_embeddings("".join(current_seq)).detach()
        with torch.no_grad():
            _, _, _, attns = self.embedder.latent_forward_ascent(emb, return_attention=True)
            # Average attention across layers and heads → (L, L)
            attn_avg = torch.stack([a.mean(dim=1) for a in attns]).mean(dim=0).squeeze(0)  # (L, L)
            # Per-position attention magnitude = sum of attention received from all other positions
            pos_importance = attn_avg.sum(dim=0)  # (L,)

        # Rank positions: highest attention, not locked, within DBD
        dbd_range = set(range(93, 292))
        candidate_positions = [
            i for i in range(len(P53_WT))
            if i not in locked_set and i in dbd_range
        ]
        candidate_positions.sort(key=lambda i: float(pos_importance[i].item()), reverse=True)
        candidate_positions = candidate_positions[:n_candidates]

        n_mutations_applied = 0
        for pass_num in range(3):  # Up to 3 passes over candidate positions
            improved_this_pass = False
            for pos_idx in candidate_positions:
                if n_mutations_applied >= max_mutations:
                    break

                # Mask this position and get ESM-2's distribution
                seq_str = "".join(current_seq)
                inputs = self.embedder.tokenizer(seq_str, return_tensors="pt", add_special_tokens=False).to(self.embedder.device)
                masked_ids = inputs.input_ids.clone()
                masked_ids[0, pos_idx] = self.embedder.tokenizer.mask_token_id

                with torch.no_grad():
                    outputs = self.embedder.model(input_ids=masked_ids, return_dict=True)
                    logits = outputs.logits[0, pos_idx]
                    probs = torch.softmax(logits / temperature, dim=-1)

                # Get top-K amino acid candidates at this position
                top_probs, top_ids = probs.topk(top_k)
                original_aa = current_seq[pos_idx]

                for k_idx in range(top_k):
                    candidate_aa_id = top_ids[k_idx].item()
                    candidate_tokens = self.embedder.tokenizer.convert_ids_to_tokens([candidate_aa_id])
                    if not candidate_tokens:
                        continue
                    candidate_aa = candidate_tokens[0]
                    if candidate_aa == original_aa or len(candidate_aa) != 1 or candidate_aa not in "ACDEFGHIKLMNPQRSTVWY":
                        continue

                    # Try this substitution
                    test_seq = list(current_seq)
                    test_seq[pos_idx] = candidate_aa
                    test_str = "".join(test_seq)

                    # Score with oracle
                    test_emb = self.embedder.get_embeddings(test_str).detach()
                    with torch.no_grad():
                        test_z, _, _ = self.embedder.latent_forward_ascent(test_emb)
                        if isinstance(self.oracle.model, AttentionPoolingNet):
                            test_oracle_input = self._align_oracle_dim(test_z)
                            test_score = float(self.oracle.model(test_oracle_input).squeeze(-1).item())
                        else:
                            test_pooled = test_z.mean(dim=1)
                            test_pooled = self._align_oracle_dim(test_pooled)
                            test_score = float(self.oracle.model(test_pooled).squeeze(-1).item())

                    if test_score > best_score:
                        current_seq = test_seq
                        best_score = test_score
                        best_seq = list(test_seq)
                        mut_label = f"{P53_WT[pos_idx]}{pos_idx+1}{candidate_aa}"
                        mutations_applied.append(mut_label)
                        n_mutations_applied += 1
                        improved_this_pass = True

                        trajectory.append({
                            "step": n_mutations_applied,
                            "score": best_score,
                            "mutation_applied": mut_label,
                            "n_mutations": n_mutations_applied,
                            "sequence": "".join(best_seq),
                        })
                        break  # Accept first improvement, move to next position

            if not improved_this_pass:
                break  # Converged

        # Build candidate dict
        final_seq = "".join(best_seq)
        all_muts = [f"{P53_WT[j]}{j+1}{best_seq[j]}" for j in range(len(P53_WT)) if P53_WT[j] != best_seq[j]]
        mut_positions = [j + 1 for j in range(len(P53_WT)) if P53_WT[j] != best_seq[j]]
        seq_identity = 100.0 * (1.0 - len(all_muts) / float(len(P53_WT)))

        # Compute stability / binding / uncertainty from final sequence
        final_stability = 0.0
        final_binding = 0.0
        final_uncertainty = 0.0
        final_ood_distance = 0.0
        try:
            final_emb = self.embedder.get_embeddings(final_seq).detach()
            with torch.no_grad():
                final_z, final_logits, _ = self.embedder.latent_forward_ascent(final_emb)
                final_logits_aa = final_logits[:, :, AA_IDS]
                final_stability = float(F.log_softmax(final_logits_aa, dim=-1).max(dim=-1).values.mean().item())
                final_probs_full = torch.softmax(final_logits, dim=-1)
                final_binding = float(self._batch_dna_force(final_z, final_probs_full)[0].item())
                final_pooled = final_z.mean(dim=1)
                final_pooled = self._align_oracle_dim(final_pooled)
                final_uncertainty = self._estimate_uncertainty(final_pooled, mc_samples=5, z_full=final_z)
                # OOD distance requires a reference; compute from WT embeddings
                wt_emb = self.embedder.get_embeddings(P53_WT).detach()
                wt_z, _, _ = self.embedder.latent_forward_ascent(wt_emb)
                wt_pooled = wt_z.mean(dim=1)
                wt_pooled = self._align_oracle_dim(wt_pooled)
                final_ood_distance = float(torch.norm(final_pooled - wt_pooled, p=2, dim=-1).item())
        except Exception as err:
            logger.warning("Autoregressive: failed to compute stability/binding for %s: %s", scenario.scenario_id, err)

        # DMS quality for rescue mutations
        target_set = set(str(t).strip().upper() for t in scenario.targets)
        rescue_muts = [m for m in all_muts if m not in target_set]
        rescue_dms_scores = []
        n_functional_rescues = 0
        for rm in rescue_muts:
            parsed = parse_single_mutation(rm)
            if parsed is not None:
                _, pos, var_aa = parsed
                z_val = self._dms_lookup.get((pos, var_aa))
                if z_val is not None:
                    rescue_dms_scores.append(z_val)
                    if z_val < 0:
                        n_functional_rescues += 1
        rescue_dms_mean = float(np.mean(rescue_dms_scores)) if rescue_dms_scores else 0.0

        candidate_uid = f"{scenario.scenario_id}|{pass_name}|autoregressive|trial{trial_idx}"
        candidate = {
            "candidate_uid": candidate_uid,
            "scenario_id": scenario.scenario_id,
            "target_label": scenario.target_label,
            "targets_json": json.dumps(scenario.targets),
            "delivery_method": scenario.delivery_method,
            "pass_name": pass_name,
            "repeat_idx": 1,
            "restart_idx": 1,
            "trial_idx": int(trial_idx),
            "profile": "Autoregressive",
            "sequence": final_seq,
            "score": float(best_score),
            "score_raw": float(best_score),
            "score_calibrated": float(self._calibrate_score(best_score)),
            "score_gain_vs_target": float(best_score - baseline_score),
            "stability": final_stability,
            "binding": final_binding,
            "identity": float(seq_identity),
            "n_mutations": int(len(all_muts)),
            "mutations_json": json.dumps(all_muts),
            "mut_positions_json": json.dumps(mut_positions),
            "uncertainty": final_uncertainty,
            "ood_distance": final_ood_distance,
            "rescue_dms_mean": float(rescue_dms_mean),
            "n_functional_rescues": int(n_functional_rescues),
            "n_rescue_mutations": int(len(rescue_muts)),
            "meets_constraints": bool(seq_identity >= 90.0),
            "selection_reason": "autoregressive_sampling",
            "receptor_source": "n/a",
            "docking_backend": "n/a",
            "md_status": "n/a",
        }

        traj_rows = []
        for row in trajectory:
            payload = dict(row)
            payload["candidate_uid"] = candidate_uid
            payload["scenario_id"] = scenario.scenario_id
            payload["target_label"] = scenario.target_label
            payload["delivery_method"] = scenario.delivery_method
            payload["pass_name"] = pass_name
            payload["profile"] = "Autoregressive"
            traj_rows.append(payload)

        return candidate, traj_rows

    def _run_clinical_for_scenario(self, candidates: List[Dict[str, Any]], top_k: int = 3) -> List[Dict[str, Any]]:
        assert self.clinical is not None

        if not candidates:
            return []
        df = pd.DataFrame(candidates)
        df = df.sort_values("score", ascending=False).head(int(max(top_k, 1)))

        rows: List[Dict[str, Any]] = []
        for _, row in df.iterrows():
            candidate_uid = str(row["candidate_uid"])
            targets = json.loads(str(row.get("targets_json", "[]")))
            cancer_mut = str(targets[0]) if targets else "R175H"
            rescue_muts = json.loads(str(row.get("mutations_json", "[]")))
            rescue_seq = str(row.get("sequence", P53_WT))
            report = self.clinical.generate_report(
                name=f"{candidate_uid.replace('|', '_')}",
                wt_sequence=P53_WT,
                rescue_sequence=rescue_seq,
                cancer_mutation=cancer_mut,
                rescue_mutations=rescue_muts,
            )
            rows.append(
                {
                    "candidate_uid": candidate_uid,
                    "scenario_id": str(row.get("scenario_id", "")),
                    "clinical_score": float(report.overall_clinical_score),
                    "clinical_viability": str(report.clinical_viability),
                    "us_annual_patients": int(report.patient_population.total_patients_per_year),
                    "global_estimate": int(report.patient_population.global_estimate),
                    "therapeutic_window": str(report.therapeutic_index.therapeutic_window),
                    "delivery_recommendation": str(report.delivery_options[0].method if report.delivery_options else "n/a"),
                }
            )
        return rows

    def _build_target_sequence(self, targets: Sequence[str]) -> str:
        seq = P53_WT
        for mut in targets:
            token = str(mut).strip().upper()
            if parse_single_mutation(token) is None:
                continue
            updated = apply_mutation(seq, token)
            if updated is not None:
                seq = updated
        return seq

    def _wt_aa_tensor(self, device: torch.device) -> torch.Tensor:
        assert self.embedder is not None
        wt_aa_indices = []
        for aa in P53_WT:
            aa_id = self.embedder.tokenizer.convert_tokens_to_ids(aa)
            if aa_id in AA_IDS:
                wt_aa_indices.append(AA_IDS.index(aa_id))
            else:
                wt_aa_indices.append(0)
        return torch.tensor(wt_aa_indices, device=device)

    def _align_oracle_dim(self, tensor: torch.Tensor) -> torch.Tensor:
        """Match embedding width to oracle input_dim by trimming or zero-padding."""
        assert self.oracle is not None
        target_dim = int(self.oracle.input_dim)
        current_dim = int(tensor.shape[-1])
        if current_dim == target_dim:
            return tensor
        if current_dim > target_dim:
            return tensor[..., :target_dim]
        return F.pad(tensor, (0, target_dim - current_dim), mode="constant", value=0.0)

    def _delivery_identity_floor(self, delivery_method: str) -> float:
        """Return minimum sequence identity floor based on delivery method.
        
        Tighter mutation budgets (94-95%) for surgical rescue edits instead of
        broad rewrites - forces biologically plausible, minimal-change solutions.
        """
        mode = str(delivery_method).strip().lower()
        if mode == "protein_therapy":
            return 95.0
        if mode == "mrna_therapy":
            return 94.0
        return 94.0

    def _batch_dna_force(self, z_batch: torch.Tensor, probs_full: torch.Tensor) -> torch.Tensor:
        hotspots = [119, 174, 240, 247, 272, 279]
        hotspots = [i for i in hotspots if i < z_batch.shape[1]]
        if not hotspots:
            return torch.zeros(z_batch.shape[0], device=z_batch.device)
        pos_charge_ids = [10, 15, 21]
        latent_force = z_batch[:, hotspots, :].norm(dim=-1).mean(dim=-1)
        charge_prob = probs_full[:, hotspots][:, :, pos_charge_ids].sum(dim=-1).mean(dim=-1)
        return 0.5 * latent_force + 5.0 * charge_prob

    def _batch_hydrophobic_packing(self, probs_full: torch.Tensor) -> torch.Tensor:
        hydro_ids = [4, 12, 7, 18, 22, 20]
        core_res = [i for i in range(93, 312) if i < probs_full.shape[1]]
        loops = set(range(111, 124)) | set(range(162, 195)) | set(range(235, 251))
        true_core = [i for i in core_res if i not in loops]
        if not true_core:
            return torch.zeros(probs_full.shape[0], device=probs_full.device)
        return probs_full[:, true_core][:, :, hydro_ids].sum(dim=-1).mean(dim=-1)

    def _decode_sequence(self, logits_aa: torch.Tensor) -> str:
        assert self.embedder is not None
        top_ids_aa = torch.argmax(logits_aa, dim=-1)[0]
        aa_local = [AA_IDS[idx] for idx in top_ids_aa.tolist()]
        tokens = self.embedder.tokenizer.convert_ids_to_tokens(aa_local)
        return "".join(tokens)[: len(P53_WT)]

    def _mutation_pos(self, mut: str) -> Optional[int]:
        digits = "".join(ch for ch in str(mut) if ch.isdigit())
        if not digits:
            return None
        try:
            return int(digits)
        except ValueError:
            return None

    def _compute_baseline_metrics(self, target_seq: str, pooled_target_ref: torch.Tensor, mc_samples: int) -> Dict[str, float]:
        assert self.embedder is not None
        assert self.oracle is not None

        emb = self.embedder.get_embeddings(target_seq).detach()
        with torch.no_grad():
            z, logits, _ = self.embedder.latent_forward_ascent(emb)
            pooled = z.mean(dim=1)
            pooled = self._align_oracle_dim(pooled)
            _uses_attention = isinstance(self.oracle.model, AttentionPoolingNet)
            if _uses_attention:
                z_oracle = self._align_oracle_dim(z)
                raw_score = float(self.oracle.model(z_oracle).item())
            else:
                raw_score = float(self.oracle.model(pooled).item())
            logits_aa = logits[:, :, AA_IDS]
            stability = float(F.log_softmax(logits_aa, dim=-1).max(dim=-1).values.mean().item())
            probs_full = torch.softmax(logits, dim=-1)
            binding = float(self._batch_dna_force(z, probs_full)[0].item())

        bundle = self._build_trust_adjusted_score(
            raw_score=raw_score,
            pooled=pooled,
            pooled_ref=pooled_target_ref,
            cached_uncertainty=self._estimate_uncertainty(pooled, mc_samples=mc_samples, z_full=z),
        )
        return {
            "score": float(bundle["score_adjusted"]),
            "score_raw": float(bundle["score_raw"]),
            "score_calibrated": float(bundle["score_calibrated"]),
            "stability": stability,
            "binding": binding,
            "uncertainty": float(bundle["uncertainty"]),
            "ood_distance": float(bundle["ood_distance"]),
        }

    def _score_calibration_profile(self) -> Dict[str, float]:
        fallback = {"clip_low": -3.0, "clip_high": 4.0, "center": 0.0, "scale": 1.0}
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

        return {"clip_low": clip_low, "clip_high": clip_high, "center": center, "scale": scale}

    def _calibrate_score(self, raw_score: float) -> float:
        center = float(self.calibration_profile.get("center", 0.0))
        scale = float(max(self.calibration_profile.get("scale", 1.0), 1e-6))
        clip_low = float(self.calibration_profile.get("clip_low", -3.0))
        clip_high = float(self.calibration_profile.get("clip_high", 4.0))
        squashed = center + scale * np.tanh((float(raw_score) - center) / scale)
        return float(np.clip(squashed, clip_low, clip_high))

    def _estimate_uncertainty(self, pooled: torch.Tensor, mc_samples: int = 8, z_full: Optional[torch.Tensor] = None) -> float:
        assert self.oracle is not None
        n_samples = max(int(mc_samples), 2)
        was_training = self.oracle.model.training
        _uses_attention = isinstance(self.oracle.model, AttentionPoolingNet)
        oracle_input = z_full if (_uses_attention and z_full is not None) else pooled
        if oracle_input is not None:
            oracle_input = self._align_oracle_dim(oracle_input)
        preds = []
        try:
            self.oracle.model.train()
            with torch.no_grad():
                for _ in range(n_samples):
                    preds.append(float(self.oracle.model(oracle_input).squeeze(-1).mean().item()))
        except Exception:
            return 0.0
        finally:
            if not was_training:
                self.oracle.model.eval()
        return float(np.std(np.asarray(preds, dtype=float)))

    def _build_trust_adjusted_score(
        self,
        *,
        raw_score: float,
        pooled: Optional[torch.Tensor],
        pooled_ref: Optional[torch.Tensor],
        cached_uncertainty: float,
    ) -> Dict[str, float]:
        calibrated = self._calibrate_score(raw_score)
        ood_distance = 0.0
        if pooled is not None and pooled_ref is not None:
            with torch.no_grad():
                ood_distance = float(torch.norm(pooled - pooled_ref, p=2, dim=-1).mean().item())

        adjusted = calibrated - self.uncertainty_weight * float(cached_uncertainty) - self.ood_rank_weight * max(
            0.0, ood_distance - self.ood_radius
        )
        adjusted = float(
            np.clip(
                adjusted,
                float(self.calibration_profile.get("clip_low", -3.0)),
                float(self.calibration_profile.get("clip_high", 4.0)),
            )
        )
        return {
            "score_raw": float(raw_score),
            "score_calibrated": float(calibrated),
            "score_adjusted": float(adjusted),
            "uncertainty": float(cached_uncertainty),
            "ood_distance": float(ood_distance),
        }

    def _append_rows(self, path: Path, rows: List[Dict[str, Any]]) -> None:
        if not rows:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")

    def _read_jsonl(self, path: Path) -> pd.DataFrame:
        if not path.exists():
            return pd.DataFrame()
        rows = []
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    logger.warning("Skipping malformed jsonl row from %s", path)
        return pd.DataFrame(rows)

    def _build_summary(
        self,
        run_id: str,
        scenario_df: pd.DataFrame,
        candidate_df: pd.DataFrame,
        top30_df: pd.DataFrame,
    ) -> str:
        lines: List[str] = []
        lines.append(f"# Campaign Summary: {run_id}")
        lines.append("")
        lines.append(f"Generated at: {datetime.now(timezone.utc).isoformat()}")
        lines.append("")
        lines.append("## Run Totals")
        unique_scenarios = 0
        pass_breakdown: Dict[str, int] = {}
        if scenario_df is not None and not scenario_df.empty:
            if "scenario_id" in scenario_df.columns:
                unique_scenarios = int(scenario_df["scenario_id"].astype(str).nunique())
            else:
                unique_scenarios = int(len(scenario_df))
            if "pass_name" in scenario_df.columns:
                pass_breakdown = (
                    scenario_df["pass_name"]
                    .astype(str)
                    .value_counts()
                    .sort_index()
                    .to_dict()
                )
        lines.append(f"- Scenarios (unique): {unique_scenarios}")
        if pass_breakdown:
            lines.append(f"- Scenario pass rows: {pass_breakdown}")
        lines.append(f"- Candidates: {len(candidate_df)}")
        lines.append(f"- Presentation shortlist: {len(top30_df)}")
        lines.append("")

        if not candidate_df.empty:
            best = candidate_df.sort_values("score", ascending=False).head(1).iloc[0]
            lines.append("## Best Candidate")
            lines.append(f"- Candidate UID: {best.get('candidate_uid', 'n/a')}")
            lines.append(f"- Target: {best.get('target_label', 'n/a')}")
            lines.append(f"- Delivery: {best.get('delivery_method', 'n/a')}")
            lines.append(f"- Score: {float(best.get('score', 0.0)):.3f}")
            lines.append(f"- Identity: {float(best.get('identity', 0.0)):.1f}%")
            lines.append(f"- Mutations: {int(best.get('n_mutations', 0))}")
            lines.append("")

        if not top30_df.empty:
            lines.append("## Top 10 from Shortlist")
            top10 = top30_df.head(10)
            for _, row in top10.iterrows():
                lines.append(
                    f"- #{int(row.get('presentation_rank', 0))}: "
                    f"{row.get('target_label', 'n/a')} | {row.get('delivery_method', 'n/a')} | "
                    f"score={float(row.get('score', 0.0)):.3f} | "
                    f"clinical={float(row.get('clinical_score', np.nan)) if pd.notna(row.get('clinical_score', np.nan)) else 'n/a'}"
                )

        return "\n".join(lines) + "\n"



def scenario_matrix_frame(
    *,
    include_pairs: bool = True,
    hotspots: Sequence[str] | None = None,
    delivery_methods: Sequence[str] | None = None,
) -> pd.DataFrame:
    scenarios = build_scenario_matrix(
        hotspots=hotspots or BIG8_HOTSPOTS,
        delivery_methods=delivery_methods or DEFAULT_DELIVERY_METHODS,
        include_pairs=include_pairs,
    )
    return scenarios_to_frame(scenarios)
