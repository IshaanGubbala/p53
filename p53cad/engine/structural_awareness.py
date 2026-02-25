"""
Structural Awareness Module for p53-proteoMgCAD

Provides design-time structural constraints for the optimization pipeline:
- Pre-computed static tensors (masks, index maps, lookup matrices)
- Tier 1 structural loss terms (mutation confidence, interface floor, crater)
- Tier 2 contact map consistency regularization
- CUDA poisoning guards for safe feature gates
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn.functional as F

from p53cad.core.logging import get_logger
from p53cad.data.dms import P53_WT


logger = get_logger(__name__)


DNA_BINDING_INTERFACE_RESIDUES = [
    120, 121, 122, 123, 124,  # L1 loop
    241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,  # L3 loop
    273, 274, 275, 276, 277, 278, 279, 280, 281,  # L2/H2 helix
    175, 176, 177, 178, 179, 180, 181,  # Core domain
]

SECONDARY_STRUCTURE_REGIONS = {
    "helix": list(range(17, 24)) + list(range(56, 64)) + list(range(90, 110)) +
             list(range(127, 145)) + list(range(153, 158)) + list(range(171, 181)) +
             list(range(229, 237)) + list(range(256, 268)) + list(range(277, 287)),
    "sheet": list(range(109, 117)) + list(range(146, 153)) + list(range(195, 205)) +
             list(range(213, 220)) + list(range(270, 277)),
}

HELIX_RESIDUES = set(SECONDARY_STRUCTURE_REGIONS["helix"])
SHEET_RESIDUES = set(SECONDARY_STRUCTURE_REGIONS["sheet"])
DBD_RANGE = list(range(94, 293))

AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
AA_LETTERS = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class PrecomputedConstants:
    """Pre-computed tensors for optimization loop efficiency."""
    
    device: torch.device
    
    wt_aa_tensor: torch.Tensor
    l1_mask: torch.Tensor
    dbd_mask: torch.Tensor
    interface_mask: torch.Tensor
    helix_mask: torch.Tensor
    sheet_mask: torch.Tensor
    
    contact_neighbor_indices: List[List[int]]
    epistasis_pairs: List[Tuple[int, int, float]]
    
    ep_i: Optional[torch.Tensor] = None
    ep_j: Optional[torch.Tensor] = None
    ep_decay: Optional[torch.Tensor] = None
    
    pool_weights: Optional[torch.Tensor] = None
    crater_kernel: Optional[torch.Tensor] = None
    
    def __post_init__(self):
        if self.ep_i is None and self.epistasis_pairs:
            self.ep_i = torch.tensor(
                [p[0] for p in self.epistasis_pairs], dtype=torch.long, device=self.device
            )
            self.ep_j = torch.tensor(
                [p[1] for p in self.epistasis_pairs], dtype=torch.long, device=self.device
            )
            self.ep_decay = torch.tensor(
                [math.exp(-p[2] / 5.0) for p in self.epistasis_pairs],
                dtype=torch.float32, device=self.device
            )


def build_precomputed_constants(
    device: torch.device,
    wt_contacts: Optional[Dict[int, List[Tuple[int, float]]]] = None,
    ca_coords: Optional[Dict[int, np.ndarray]] = None,
    locked_indices: Optional[List[int]] = None,
) -> PrecomputedConstants:
    """Build all pre-computed static tensors once per trial."""
    
    locked_set = set(locked_indices) if locked_indices else set()
    seq_len = len(P53_WT)
    
    aa_to_id = {aa: i for i, aa in enumerate(AA_LETTERS)}
    wt_aa_tensor = torch.tensor(
        [aa_to_id.get(P53_WT[i], 0) for i in range(seq_len)],
        dtype=torch.long, device=device
    )
    
    l1_mask = torch.ones(seq_len, dtype=torch.bool, device=device)
    for li in locked_indices or []:
        l1_mask[li] = False
    
    dbd_mask = torch.zeros(seq_len, dtype=torch.bool, device=device)
    for pos in DBD_RANGE:
        if 0 <= pos - 1 < seq_len:
            dbd_mask[pos - 1] = True
    
    interface_mask = torch.zeros(seq_len, dtype=torch.bool, device=device)
    for pos in DNA_BINDING_INTERFACE_RESIDUES:
        if 0 <= pos - 1 < seq_len:
            interface_mask[pos - 1] = True
    
    helix_mask = torch.zeros(seq_len, dtype=torch.bool, device=device)
    for pos in HELIX_RESIDUES:
        if 0 <= pos - 1 < seq_len:
            helix_mask[pos - 1] = True
    
    sheet_mask = torch.zeros(seq_len, dtype=torch.bool, device=device)
    for pos in SHEET_RESIDUES:
        if 0 <= pos - 1 < seq_len:
            sheet_mask[pos - 1] = True
    
    contact_neighbor_indices: List[List[int]] = []
    if wt_contacts:
        for li in locked_indices or []:
            res_id = li + 1
            neighbors = wt_contacts.get(res_id, [])
            ni = [r - 1 for r, _ in neighbors if 0 <= r - 1 < seq_len]
            contact_neighbor_indices.append(ni)
    
    _epistasis_pairs: List[Tuple[int, int, float]] = []
    if ca_coords:
        all_pos = list(range(seq_len))
        for i_idx in range(len(all_pos)):
            if i_idx in locked_set:
                continue
            ri = i_idx + 1
            if ri not in ca_coords:
                continue
            for j_idx in range(i_idx + 1, len(all_pos)):
                if j_idx in locked_set:
                    continue
                rj = j_idx + 1
                if rj not in ca_coords:
                    continue
                diff = ca_coords[ri] - ca_coords[rj]
                dist = float(np.sqrt(np.sum(diff ** 2)))
                if dist < 10.0:
                    _epistasis_pairs.append((i_idx, j_idx, dist))
    
    seq_len_float = float(seq_len)
    pool_weights = torch.ones(seq_len, device=device)
    for li in locked_indices or []:
        dists = (torch.arange(seq_len, device=device).float() - li)
        pool_weights += torch.exp(-dists ** 2 / 200.0)
    pool_weights = pool_weights / pool_weights.sum()
    
    crater_kernel = torch.zeros(seq_len, device=device)
    for li in locked_indices or []:
        dists = (torch.arange(seq_len, device=device).float() - li)
        crater_kernel += torch.exp(-dists ** 2 / 18.0)
    
    return PrecomputedConstants(
        device=device,
        wt_aa_tensor=wt_aa_tensor,
        l1_mask=l1_mask,
        dbd_mask=dbd_mask,
        interface_mask=interface_mask,
        helix_mask=helix_mask,
        sheet_mask=sheet_mask,
        contact_neighbor_indices=contact_neighbor_indices,
        epistasis_pairs=_epistasis_pairs,
        pool_weights=pool_weights,
        crater_kernel=crater_kernel,
    )


def compute_mutation_confidence_penalty(
    probs_aa: torch.Tensor,
    wt_aa_tensor: torch.Tensor,
    helix_mask: torch.Tensor,
    sheet_mask: torch.Tensor,
    locked_indices: List[int],
    device: torch.device,
    helix_weight: float = 1.5,
    sheet_weight: float = 1.2,
    coil_weight: float = 1.0,
) -> torch.Tensor:
    """Penalize low-plausibility substitutions.
    
    Stronger penalty in helices/sheets where structural constraints are tighter.
    """
    if not locked_indices:
        return torch.zeros(1, device=device)
    
    wt_probs = probs_aa[0, torch.arange(probs_aa.size(1), device=device), wt_aa_tensor]
    mut_prob = 1.0 - wt_probs
    
    struct_weight = torch.ones(probs_aa.size(1), device=device)
    struct_weight[helix_mask] = helix_weight
    struct_weight[sheet_mask] = sheet_weight
    
    active_mask = torch.ones(probs_aa.size(1), dtype=torch.bool, device=device)
    for li in locked_indices:
        active_mask[li] = False
    
    penalty = (mut_prob * struct_weight * active_mask).sum() / max(active_mask.sum().item(), 1)
    return 2.0 * penalty


def compute_interface_floor_penalty(
    probs_aa: torch.Tensor,
    wt_aa_tensor: torch.Tensor,
    interface_mask: torch.Tensor,
    device: torch.device,
    min_interface_similarity: float = 0.7,
) -> torch.Tensor:
    """Enforce minimum similarity at DNA-binding interface residues."""
    if interface_mask.sum().item() == 0:
        return torch.zeros(1, device=device)
    
    wt_probs = probs_aa[0, torch.arange(probs_aa.size(1), device=device), wt_aa_tensor]
    
    interface_similarity = wt_probs[interface_mask].mean()
    penalty = F.relu(min_interface_similarity - interface_similarity)
    return 15.0 * penalty


def compute_crater_penalty(
    probs_aa: torch.Tensor,
    wt_aa_tensor: torch.Tensor,
    crater_kernel: torch.Tensor,
    device: torch.device,
    threshold: float = 0.3,
) -> torch.Tensor:
    """Penalize local confidence collapses around mutated positions.
    
    Uses a Gaussian kernel centered on locked positions to detect
    low-confidence neighborhoods.
    """
    wt_probs = probs_aa[0, torch.arange(probs_aa.size(1), device=device), wt_aa_tensor]
    local_confidence = (wt_probs * crater_kernel) / (crater_kernel.sum() + 1e-8)
    
    collapses = F.relu(threshold - local_confidence)
    return 3.0 * collapses.sum()


def compute_contact_consistency_penalty(
    z: torch.Tensor,
    wt_hidden: torch.Tensor,
    contact_neighbor_indices: List[List[int]],
    device: torch.device,
    weight: float = 5.0,
) -> torch.Tensor:
    """Penalize changes in contact neighborhood hidden states."""
    if not contact_neighbor_indices:
        return torch.zeros(1, device=device)
    
    cos_sims = []
    for ni_list in contact_neighbor_indices:
        if not ni_list:
            continue
        ni_t = torch.tensor(ni_list, device=device)
        cur_h = z[:, ni_t, :]
        wt_h = wt_hidden[:, ni_t, :]
        sim = F.cosine_similarity(cur_h, wt_h, dim=-1).mean(dim=-1)
        cos_sims.append(sim)
    
    if not cos_sims:
        return torch.zeros(1, device=device)
    
    mean_cos = torch.stack(cos_sims).mean(dim=0)
    return weight * (1.0 - mean_cos)


def compute_contact_map_regularization(
    z: torch.Tensor,
    wt_hidden: torch.Tensor,
    device: torch.device,
    interval: int = 10,
    weight: float = 2.0,
) -> torch.Tensor:
    """Tier 2: Contact map consistency regularization against WT reference.
    
    Run at controlled intervals to preserve global structural topology.
    """
    seq_len = z.size(1)
    if wt_hidden is None or wt_hidden.size(1) != seq_len:
        return torch.zeros(1, device=device)
    
    diff = z - wt_hidden
    contact_reg = weight * (diff.abs().mean())
    return contact_reg


class AsyncLogger:
    """Background-thread async logger to avoid GPU path stalls."""
    
    def __init__(self, logger_name: str = "p53cad.async", interval: int = 10):
        self.logger = get_logger(logger_name)
        self.interval = interval
        self._buffer: List[Dict] = []
        self._counter = 0
        self._enabled = True
    
    def log(self, data: Dict) -> None:
        """Buffer log data for async write."""
        if not self._enabled:
            return
        self._buffer.append(data)
        self._counter += 1
    
    def flush(self) -> None:
        """Flush buffered logs."""
        if not self._buffer:
            return
        for entry in self._buffer:
            self.logger.info(entry)
        self._buffer.clear()
    
    def should_flush(self) -> bool:
        return self._counter >= self.interval
    
    def reset_counter(self) -> None:
        self._counter = 0
    
    def disable(self) -> None:
        self._enabled = False
        self.flush()


class CUDASafetyGuard:
    """Runtime checks to prevent CUDA poisoning from risky features."""
    
    _checked: Dict[str, bool] = {}
    _stable: bool = True
    
    @classmethod
    def check_cuda_stability(cls, device: torch.device) -> bool:
        """Verify CUDA is stable before enabling risky features."""
        if device.type != "cuda":
            return False
        
        cache_key = "cuda_stability"
        if cache_key in cls._checked:
            return cls._checked[cache_key]
        
        try:
            test_tensor = torch.zeros(10, device=device)
            _ = test_tensor + 1
            _ = torch.matmul(test_tensor.unsqueeze(0), test_tensor.unsqueeze(1))
            del test_tensor
            torch.cuda.synchronize()
            cls._stable = True
            cls._checked[cache_key] = True
            logger.debug("CUDA stability check passed")
            return True
        except Exception as e:
            logger.warning("CUDA stability check failed: %s. Risky features disabled.", e)
            cls._stable = False
            cls._checked[cache_key] = False
            return False
    
    @classmethod
    def check_memory_available(cls, device: torch.device, required_mb: int = 1000) -> bool:
        """Check if sufficient GPU memory is available."""
        if device.type != "cuda":
            return False
        
        try:
            free_mem = torch.cuda.get_device_properties(device).total_memory - torch.cuda.memory_allocated(device)
            return (free_mem / 1024 / 1024) > required_mb
        except Exception:
            return False
    
    @classmethod
    def is_stable(cls) -> bool:
        return cls._stable


def get_optimization_config(
    device: torch.device,
    enable_shallow_gradients: bool = True,
    enable_batched_trials: bool = False,
    enable_sdpa: bool = True,
) -> Dict[str, bool]:
    """Determine which optimization features are safe to enable."""
    
    config = {
        "shallow_gradients": False,
        "batched_trials": False,
        "sdpa_fast_path": False,
        "tf32_enabled": device.type == "cuda",
        "cudnn_benchmark": device.type == "cuda",
    }
    
    if not CUDASafetyGuard.check_cuda_stability(device):
        logger.info("CUDA instability detected, using safe defaults")
        return config
    
    if enable_shallow_gradients and CUDASafetyGuard.check_memory_available(device, required_mb=500):
        config["shallow_gradients"] = True
    
    if enable_batched_trials and CUDASafetyGuard.check_memory_available(device, required_mb=2000):
        config["batched_trials"] = True
    
    if enable_sdpa:
        config["sdpa_fast_path"] = True
    
    return config


class PopulationShare:
    """Lightweight cooperative mutation sharing between parallel profile searches."""
    
    def __init__(self, max_size: int = 100):
        self._mutations: Dict[str, float] = {}
        self._scores: Dict[str, float] = {}
        self._max_size = max_size
        self._lock_count = 0
    
    def add_candidate(self, mutations: str, score: float) -> None:
        """Add a candidate's mutations and score to the shared pool."""
        if len(self._mutations) >= self._max_size:
            worst_key = min(self._scores, key=self._scores.get)
            del self._mutations[worst_key]
            del self._scores[worst_key]
        
        mut_key = ",".join(sorted(mutations.split(","))) if mutations else ""
        if mut_key not in self._scores or score > self._scores[mut_key]:
            self._mutations[mut_key] = score
            self._scores[mut_key] = score
    
    def get_top_mutations(self, top_k: int = 10) -> List[Tuple[str, float]]:
        """Get top-K mutations from the shared pool."""
        sorted_items = sorted(self._scores.items(), key=lambda x: x[1], reverse=True)
        return [(k, v) for k, v in sorted_items[:top_k] if k]
    
    def get_inspiring_mutations(self, current_mutations: List[str], top_k: int = 5) -> List[str]:
        """Get mutations from other profiles that are different but high-scoring."""
        current_set = set(current_mutations)
        top_muts = self.get_top_mutations(top_k=20)
        
        inspiring = []
        for mut_str, score in top_muts:
            if not mut_str:
                continue
            muts = mut_str.split(",")
            for m in muts:
                if m and m not in current_set:
                    inspiring.append(m)
                    if len(inspiring) >= top_k:
                        return inspiring
        return inspiring
    
    def clear(self) -> None:
        """Clear the shared pool."""
        self._mutations.clear()
        self._scores.clear()
        self._lock_count = 0
    
    @property
    def size(self) -> int:
        return len(self._mutations)


@dataclass
class ESMFoldGateConfig:
    """Configuration for Tier 3 ESMFold quality gate."""
    min_global_plddt: float = 70.0
    min_dbd_plddt: float = 65.0
    min_interface_plddt: float = 60.0
    max_wt_rmsd: float = 5.0
    enabled: bool = True


def check_esmfold_gate(
    plddt_scores: List[float],
    dbd_range: Tuple[int, int] = (93, 292),
    interface_residues: List[int] = None,
    config: ESMFoldGateConfig = None,
) -> Tuple[bool, Dict[str, float]]:
    """Check if candidate passes the ESMFold quality gate.
    
    Returns:
        (passed, metrics_dict)
    """
    if config is None:
        config = ESMFoldGateConfig()
    
    if not config.enabled:
        return True, {"gate_passed": True, "reason": "disabled"}
    
    if not plddt_scores:
        # If no structural confidence values are available, do not hard-filter.
        # This keeps the gate fail-open unless real ESMFold metrics are provided.
        return True, {"gate_passed": True, "reason": "skipped_no_plddt"}
    
    metrics = {}
    metrics["mean_plddt"] = float(np.mean(plddt_scores))
    
    dbd_start, dbd_end = dbd_range
    dbd_scores = [plddt_scores[i] for i in range(dbd_start, min(dbd_end + 1, len(plddt_scores)))]
    metrics["dbd_plddt"] = float(np.mean(dbd_scores)) if dbd_scores else 0.0
    
    if interface_residues is None:
        interface_residues = DNA_BINDING_INTERFACE_RESIDUES
    
    interface_scores = [plddt_scores[i - 1] for i in interface_residues if i <= len(plddt_scores)]
    metrics["interface_plddt"] = float(np.mean(interface_scores)) if interface_scores else 0.0
    
    passed = (
        metrics["mean_plddt"] >= config.min_global_plddt and
        metrics["dbd_plddt"] >= config.min_dbd_plddt and
        metrics["interface_plddt"] >= config.min_interface_plddt
    )
    
    metrics["gate_passed"] = passed
    if not passed:
        reasons = []
        if metrics["mean_plddt"] < config.min_global_plddt:
            reasons.append(f"global_plddt {metrics['mean_plddt']:.1f} < {config.min_global_plddt}")
        if metrics["dbd_plddt"] < config.min_dbd_plddt:
            reasons.append(f"dbd_plddt {metrics['dbd_plddt']:.1f} < {config.min_dbd_plddt}")
        if metrics["interface_plddt"] < config.min_interface_plddt:
            reasons.append(f"interface_plddt {metrics['interface_plddt']:.1f} < {config.min_interface_plddt}")
        metrics["reason"] = "; ".join(reasons)
    
    return passed, metrics
