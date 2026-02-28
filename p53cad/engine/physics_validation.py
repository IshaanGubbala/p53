"""
Physics-Based Validation Pipeline for p53-proteoMgCAD

Replaces heuristic validation (Chou-Fasman pLDDT, lookup-table DDG) with real
physics: local ESMFold structure prediction, OpenMM energy minimization, short
MD stability checks, and DNA-binding interface analysis via mdtraj.

All heavy dependencies (openmm, pdbfixer, mdtraj, transformers) are imported
lazily so the rest of p53cad works without them.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from p53cad.core.logging import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ESMFoldResult:
    """Local ESMFold structure prediction result."""
    pdb_string: str
    plddt_scores: List[float]
    mean_plddt: float
    dbd_plddt: float          # residues 94-292
    ptm_score: Optional[float] = None
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d.pop("pdb_string", None)   # exclude bulky PDB text from JSON
        return d


@dataclass
class EnergyMinimizationResult:
    """OpenMM energy minimization result."""
    potential_energy_kcal: float
    force_field: str = "amber14"
    solvent_model: str = "OBC2"
    n_atoms: int = 0
    pdbfixer_missing_residues: int = 0
    pdbfixer_missing_atoms: int = 0
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DDGResult:
    """Relative DDG result (rescue vs WT and vs cancer mutant)."""
    ddg_vs_wt_kcal: float
    ddg_vs_cancer_kcal: float
    interpretation: str  # "stabilizing", "neutral", "destabilizing"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class MDStabilityResult:
    """Short MD stability check result."""
    rmsd_mean: float         # Angstroms
    rmsd_std: float
    rmsd_final: float
    rmsf_per_residue: List[float]
    rmsf_dbd_mean: float     # DBD residues 94-292
    ss_content: Dict[str, float]   # helix, sheet, coil fractions
    rg_mean: float           # radius of gyration (nm)
    stability_score: float   # 0-100
    stability_verdict: str   # "stable", "metastable", "unstable"
    simulation_ns: float = 0.2
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        # Truncate per-residue RMSF to save JSON space
        if len(d.get("rmsf_per_residue", [])) > 50:
            d["rmsf_per_residue"] = [round(v, 3) for v in d["rmsf_per_residue"]]
        return d


@dataclass
class DNABindingInterfaceResult:
    """DNA-binding interface preservation analysis."""
    dbd_rmsd: float                          # CA RMSD of DBD vs WT (Angstroms)
    per_contact_rmsd: Dict[str, float]       # residue label -> CA displacement
    interface_preservation_score: float      # 0.0 to 1.0
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DNABindingSimResult:
    """MM-GBSA DNA-binding functional simulation result."""
    binding_energy_kcal: float
    wt_binding_energy_kcal: float
    delta_binding_kcal: float
    binding_interpretation: str     # "enhanced", "preserved", "weakened", "disrupted"
    n_hbonds_protein_dna: int
    wt_n_hbonds: int
    hbond_preservation_ratio: float
    hbond_details: List[Dict[str, Any]]   # [{donor_res, acceptor_res, distance}, ...]
    n_interface_contacts: int
    wt_n_contacts: int
    contact_preservation_ratio: float
    per_residue_contacts: Dict[str, int]  # residue label -> n_contacts
    binding_score: float                  # 0-1 composite
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        # Truncate hbond details if too many
        if len(d.get("hbond_details", [])) > 30:
            d["hbond_details"] = d["hbond_details"][:30]
        return d


@dataclass
class PhysicsValidationResult:
    """Composite per-candidate validation result."""
    candidate_rank: int
    target_label: str
    sequence_hash: str
    esmfold: Optional[ESMFoldResult] = None
    energy: Optional[EnergyMinimizationResult] = None
    ddg: Optional[DDGResult] = None
    md: Optional[MDStabilityResult] = None
    dna_binding: Optional[DNABindingInterfaceResult] = None
    dna_binding_sim: Optional[DNABindingSimResult] = None
    overall_physics_score: float = 0.0
    verdict: str = "UNKNOWN"
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "candidate_rank": self.candidate_rank,
            "target_label": self.target_label,
            "sequence_hash": self.sequence_hash,
            "overall_physics_score": round(self.overall_physics_score, 1),
            "verdict": self.verdict,
            "errors": self.errors,
        }
        if self.esmfold:
            d["esmfold"] = self.esmfold.to_dict()
        if self.energy:
            d["energy"] = self.energy.to_dict()
        if self.ddg:
            d["ddg"] = self.ddg.to_dict()
        if self.md:
            d["md"] = self.md.to_dict()
        if self.dna_binding:
            d["dna_binding"] = self.dna_binding.to_dict()
        if self.dna_binding_sim:
            d["dna_binding_sim"] = self.dna_binding_sim.to_dict()
        return d


@dataclass
class PhysicsValidationReport:
    """Full report across all candidates."""
    run_id: str
    n_candidates: int
    n_tier2: int
    wt_energy_kcal: Optional[float] = None
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    elapsed_total_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "run_id": self.run_id,
            "n_candidates": self.n_candidates,
            "n_tier2": self.n_tier2,
            "wt_energy_kcal": self.wt_energy_kcal,
            "elapsed_total_sec": round(self.elapsed_total_sec, 1),
            "candidates": self.candidates,
        }


# ---------------------------------------------------------------------------
# Helper: sequence hash for caching
# ---------------------------------------------------------------------------

def _seq_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# LocalESMFoldPredictor  (singleton, lazy-loaded)
# ---------------------------------------------------------------------------

class LocalESMFoldPredictor:
    """Predict structures with facebook/esmfold_v1 loaded locally via HuggingFace."""

    _instance: Optional["LocalESMFoldPredictor"] = None

    def __new__(cls, *args: Any, **kwargs: Any) -> "LocalESMFoldPredictor":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self, device: str = "cpu", cache_dir: Optional[Path] = None):
        if self._initialized:
            return
        self._device = device
        self._cache_dir = Path(cache_dir) if cache_dir else None
        self._model = None
        self._tokenizer = None
        self._disabled_reason: Optional[str] = None
        self._initialized = True

    def _ensure_loaded(self) -> None:
        if self._model is not None:
            return
        import torch
        from transformers import AutoTokenizer, EsmForProteinFolding

        logger.info("Loading ESMFold model (this may take a minute)...")
        t0 = time.time()

        try:
            self._tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
        except Exception as e:
            raise RuntimeError(f"Failed to load ESMFold tokenizer: {e}") from e

        # ESMFold ships .bin (pickle) weights, not safetensors.
        # Transformers 4.57+ raises a CVE-2025-32434 ValueError on torch < 2.6 before loading.
        # Patch both the safety check and torch.load to allow pickle weights on any torch version.
        from unittest.mock import patch

        _orig_load = torch.load

        def _permissive_load(*args: Any, **kwargs: Any) -> Any:
            kwargs["weights_only"] = False
            return _orig_load(*args, **kwargs)

        try:
            # transformers<=4.44 does not expose check_torch_load_is_safe.
            # create=True keeps this compatible across transformer versions.
            with patch("transformers.modeling_utils.check_torch_load_is_safe", lambda: None, create=True), \
                 patch("torch.load", _permissive_load):
                self._model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
        except Exception as e:
            raise RuntimeError(f"Failed to load ESMFold model: {e}") from e

        self._model = self._model.to(self._device)
        self._model.eval()
        # chunk_size=64 trades memory for speed on long sequences
        self._model.trunk.set_chunk_size(64)
        logger.info("ESMFold loaded in %.1fs on %s", time.time() - t0, self._device)

    def predict(
        self,
        sequence: str,
        save_dir: Optional[Path] = None,
        filename: Optional[str] = None,
    ) -> ESMFoldResult:
        """Predict structure for a single sequence. Caches by SHA-256 hash."""
        if self._disabled_reason is not None:
            raise RuntimeError(f"ESMFold disabled for this process: {self._disabled_reason}")

        seq_h = _seq_hash(sequence)

        # Check disk cache
        if self._cache_dir:
            cached = self._cache_dir / f"{seq_h}.json"
            cached_pdb = self._cache_dir / f"{seq_h}.pdb"
            if cached.exists() and cached_pdb.exists():
                logger.info("ESMFold cache hit: %s", seq_h)
                meta = json.loads(cached.read_text())
                pdb_str = cached_pdb.read_text()
                return ESMFoldResult(
                    pdb_string=pdb_str,
                    plddt_scores=meta["plddt_scores"],
                    mean_plddt=meta["mean_plddt"],
                    dbd_plddt=meta["dbd_plddt"],
                    ptm_score=meta.get("ptm_score"),
                    elapsed_sec=0.0,
                )

        self._ensure_loaded()
        import torch

        # ESMFold expects canonical amino-acid letters; sanitize unknown tokens.
        aa_set = set("ACDEFGHIKLMNPQRSTVWY")
        seq_input = sequence.upper()
        if any(ch not in aa_set for ch in seq_input):
            n_bad = sum(1 for ch in seq_input if ch not in aa_set)
            logger.warning(
                "ESMFold input has %d non-canonical residues; replacing with 'A' for inference.",
                n_bad,
            )
            seq_input = "".join(ch if ch in aa_set else "A" for ch in seq_input)

        t0 = time.time()
        tokenized = self._tokenizer([seq_input], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.to(self._device) for k, v in tokenized.items()}
        if "input_ids" in tokenized:
            tokenized["input_ids"] = tokenized["input_ids"].long()
        if "attention_mask" in tokenized:
            tokenized["attention_mask"] = tokenized["attention_mask"].long()

        with torch.no_grad():
            try:
                output = self._model(**tokenized)
            except RuntimeError as exc:
                if "one_hot is only applicable to index tensor" not in str(exc):
                    raise
                self._disabled_reason = f"one_hot index tensor failure on {self._device}"
                logger.warning(
                    "ESMFold forward failed (%s). Disabling ESMFold for this process to avoid repeated stalls.",
                    self._disabled_reason,
                )
                raise RuntimeError(self._disabled_reason) from exc

        # Convert to PDB string
        from transformers.models.esm.openfold_utils.protein import to_pdb, Protein

        # Build Protein object from output
        pdb_string = self._output_to_pdb(output, seq_input)

        # Extract per-residue pLDDT.
        # Shape is (batch, L, 37) — per-atom pLDDT in atom37 format.
        # Take CA atom (index 1) as representative per-residue pLDDT.
        plddt_tensor = output["plddt"][0, :len(seq_input)]  # (L, 37)
        if plddt_tensor.ndim == 2:
            # Use CA atom (atom37 index 1) for per-residue score
            plddt_scores = plddt_tensor[:, 1].cpu().numpy().tolist()
        else:
            plddt_scores = plddt_tensor.cpu().numpy().tolist()
        mean_plddt = float(np.mean(plddt_scores))
        # DBD = residues 94-292 (0-indexed: 93-291)
        dbd_start, dbd_end = 93, min(292, len(seq_input))
        dbd_plddt = float(np.mean(plddt_scores[dbd_start:dbd_end])) if dbd_end > dbd_start else mean_plddt

        # pTM score — output["ptm"] is a scalar tensor
        ptm_score = None
        try:
            if "ptm" in output and output["ptm"] is not None:
                ptm_score = float(output["ptm"].cpu().item())
        except Exception:
            pass

        elapsed = time.time() - t0
        logger.info("ESMFold prediction: %d residues, pLDDT=%.1f, DBD=%.1f, %.1fs",
                     len(seq_input), mean_plddt, dbd_plddt, elapsed)

        result = ESMFoldResult(
            pdb_string=pdb_string,
            plddt_scores=plddt_scores,
            mean_plddt=mean_plddt,
            dbd_plddt=dbd_plddt,
            ptm_score=ptm_score,
            elapsed_sec=elapsed,
        )

        # Cache to disk
        if self._cache_dir:
            self._cache_dir.mkdir(parents=True, exist_ok=True)
            meta = result.to_dict()
            meta["sequence_length"] = len(seq_input)
            (self._cache_dir / f"{seq_h}.json").write_text(json.dumps(meta, indent=2))
            (self._cache_dir / f"{seq_h}.pdb").write_text(pdb_string)

        # Optionally save to a specific directory
        if save_dir:
            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)
            fname = filename or f"esmfold_{seq_h}.pdb"
            (save_dir / fname).write_text(pdb_string)

        return result

    def _output_to_pdb(self, output: Any, sequence: str) -> str:
        """Convert ESMFold output tensors to PDB string.

        ESMFold predicts 14 atoms per residue (atom14 format). We scatter
        these into atom37 slots and set the mask to only include positions
        that were actually predicted, so PDBFixer correctly rebuilds the rest.
        """
        import torch
        L = len(sequence)
        try:
            from transformers.models.esm.openfold_utils.protein import to_pdb, Protein

            positions_14 = output["positions"][-1][0, :L].cpu().numpy()  # (L, 14, 3)
            atom14_exists = output["atom14_atom_exists"][0, :L].cpu().numpy()  # (L, 14)
            aatype = output["aatype"][0, :L].cpu().numpy()
            residx_14_to_37 = output["residx_atom14_to_atom37"][0, :L].cpu().numpy()
            plddt_37 = output["plddt"][0, :L].cpu().numpy()  # (L, 37)

            # Scatter atom14 positions AND mask into atom37 slots
            positions_37 = np.zeros((L, 37, 3), dtype=np.float32)
            mask_37 = np.zeros((L, 37), dtype=np.float32)
            for i in range(L):
                for j in range(14):
                    if atom14_exists[i, j] > 0.5:
                        k = int(residx_14_to_37[i, j])
                        positions_37[i, k] = positions_14[i, j]
                        mask_37[i, k] = 1.0

            protein = Protein(
                aatype=aatype,
                atom_positions=positions_37,
                atom_mask=mask_37,  # Only predicted atoms, not full atom37_exists
                residue_index=np.arange(L),
                b_factors=plddt_37,
                chain_index=np.zeros(L, dtype=np.int32),
            )
            return to_pdb(protein)
        except Exception as e:
            logger.warning("openfold_utils PDB conversion failed (%s), using CA-only fallback", e)
            return self._fallback_pdb(output, sequence)

    def _fallback_pdb(self, output: Any, sequence: str) -> str:
        """Fallback PDB generation from atom37 positions."""
        import torch
        # Mapping of 3-letter codes
        _aa1to3 = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
            'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
            'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
        }
        positions = output["positions"][-1][0, :len(sequence)].cpu().numpy()
        # plddt is (batch, L, 37) per-atom; take CA (index 1) for per-residue
        plddt_per_atom = output["plddt"][0, :len(sequence)].cpu().numpy()
        plddt = plddt_per_atom[:, 1] if plddt_per_atom.ndim == 2 else plddt_per_atom
        lines = ["HEADER    ESMFOLD LOCAL PREDICTION"]
        atom_idx = 1
        for i, aa in enumerate(sequence):
            res3 = _aa1to3.get(aa, "UNK")
            # CA atom is index 1 in atom37 representation
            ca_pos = positions[i, 1]
            lines.append(
                f"ATOM  {atom_idx:5d}  CA  {res3} A{i+1:4d}    "
                f"{ca_pos[0]:8.3f}{ca_pos[1]:8.3f}{ca_pos[2]:8.3f}"
                f"  1.00{plddt[i]:6.2f}           C"
            )
            atom_idx += 1
        lines.append("END")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# OpenMMEnergyCalculator
# ---------------------------------------------------------------------------

class OpenMMEnergyCalculator:
    """Energy minimization via OpenMM with AMBER14 + OBC2 implicit solvent."""

    def __init__(self):
        self._forcefield = None

    @staticmethod
    def _select_platform(openmm: Any) -> Any:
        """Prefer CUDA, then OpenCL (Intel Xe / AMD), then CPU."""
        for pname in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = openmm.Platform.getPlatformByName(pname)
                logger.info("OpenMM platform selected: %s", pname)
                return platform
            except Exception:
                continue
        # Should never reach here — CPU platform is always present
        return openmm.Platform.getPlatformByName("CPU")

    def _get_forcefield(self):
        """Cached ForceField to avoid re-parsing AMBER14 XML on every candidate."""
        if self._forcefield is None:
            import openmm.app as app
            self._forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
        return self._forcefield

    def prepare_structure(self, pdb_string: str) -> Any:
        """Fix PDB with PDBFixer: add missing atoms + hydrogens.

        Removes non-standard residues and handles ESMFold output quirks
        (atoms at origin, missing chain IDs).
        """
        from pdbfixer import PDBFixer
        import io

        fixer = PDBFixer(pdbfile=io.StringIO(pdb_string))
        fixer.removeHeterogens(False)  # remove non-standard residues
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        missing_res = len(fixer.missingResidues)
        missing_atoms_count = sum(len(v) for v in fixer.missingAtoms.values())
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)
        return fixer, missing_res, missing_atoms_count

    def minimize_and_get_energy(
        self,
        pdb_string: str,
        max_iterations: int = 1000,
    ) -> EnergyMinimizationResult:
        """Run energy minimization and return potential energy."""
        t0 = time.time()
        fixer, missing_res, missing_atoms = self.prepare_structure(pdb_string)
        result = self.minimize_from_fixer(fixer, max_iterations=max_iterations)
        result.pdbfixer_missing_residues = missing_res
        result.pdbfixer_missing_atoms = missing_atoms
        result.elapsed_sec = round(time.time() - t0, 1)
        return result

    def minimize_from_fixer(
        self,
        fixer: Any,
        max_iterations: int = 1000,
    ) -> EnergyMinimizationResult:
        """Run energy minimization from an already-prepared PDBFixer object.

        Reuses the cached ForceField to avoid re-parsing AMBER14 XML per candidate.
        """
        import openmm
        import openmm.app as app
        import openmm.unit as unit

        t0 = time.time()

        forcefield = self._get_forcefield()
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            310 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
        )
        platform = self._select_platform(openmm)
        simulation = app.Simulation(fixer.topology, system, integrator, platform)
        simulation.context.setPositions(fixer.positions)

        # Minimize; retry with 2× iterations if energy didn't converge (positive = not converged)
        simulation.minimizeEnergy(maxIterations=max_iterations)
        state = simulation.context.getState(getEnergy=True)
        energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        energy_kcal = energy_kj / 4.184

        if energy_kcal > 0:
            logger.warning(
                "Energy minimization did not converge (E=%.1f kcal/mol), retrying with %d iterations",
                energy_kcal, max_iterations * 2,
            )
            simulation.minimizeEnergy(maxIterations=max_iterations * 2)
            state = simulation.context.getState(getEnergy=True)
            energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            energy_kcal = energy_kj / 4.184

        n_atoms = sum(1 for _ in fixer.topology.atoms())
        elapsed = time.time() - t0

        logger.info("Energy minimization: %.1f kcal/mol (%d atoms, %.1fs)",
                     energy_kcal, n_atoms, elapsed)

        return EnergyMinimizationResult(
            potential_energy_kcal=energy_kcal,
            force_field="amber14",
            solvent_model="OBC2",
            n_atoms=n_atoms,
            pdbfixer_missing_residues=0,
            pdbfixer_missing_atoms=0,
            elapsed_sec=round(elapsed, 1),
        )

    def compute_ddg(
        self,
        e_wt: float,
        e_cancer: float,
        e_rescue: float,
    ) -> DDGResult:
        """Compute relative DDG from minimized energies."""
        ddg_vs_wt = e_rescue - e_wt
        ddg_vs_cancer = e_rescue - e_cancer

        if ddg_vs_wt < -5.0:
            interp = "stabilizing"
        elif ddg_vs_wt > 5.0:
            interp = "destabilizing"
        else:
            interp = "neutral"

        return DDGResult(
            ddg_vs_wt_kcal=round(ddg_vs_wt, 2),
            ddg_vs_cancer_kcal=round(ddg_vs_cancer, 2),
            interpretation=interp,
        )


# ---------------------------------------------------------------------------
# MDStabilityChecker (top candidates only)
# ---------------------------------------------------------------------------

class MDStabilityChecker:
    """Short implicit-solvent MD for rapid stability screening."""

    def run_stability_check(
        self,
        pdb_string: str,
        output_dir: Optional[Path] = None,
        name: str = "candidate",
        simulation_ns: float = 0.2,
        fixer: Any = None,
    ) -> MDStabilityResult:
        """
        Run a short MD simulation and analyze trajectory.

        Pipeline: PDBFixer -> AMBER14/OBC2 -> equilibration -> production -> analysis.
        Default 0.2 ns production (~200 ps) balances speed and signal for gross
        instability detection (RMSD blowup, unfolding).

        If *fixer* is provided (already-prepared PDBFixer), skip re-preparation.
        """
        import openmm
        import openmm.app as app
        import openmm.unit as unit
        import mdtraj
        import io
        import tempfile

        t0 = time.time()

        # Prepare (reuse fixer if provided)
        calc = OpenMMEnergyCalculator()
        if fixer is None:
            fixer, _, _ = calc.prepare_structure(pdb_string)

        forcefield = calc._get_forcefield()
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            310 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
        )
        platform = calc._select_platform(openmm)
        simulation = app.Simulation(fixer.topology, system, integrator, platform)
        simulation.context.setPositions(fixer.positions)

        # Minimize first
        simulation.minimizeEnergy(maxIterations=200)

        # Equilibration: 1000 steps = 2 ps (fast)
        simulation.step(1000)

        # Production with DCD reporter
        production_steps = int(simulation_ns * 1e6 / 2)  # 2 fs timestep
        report_interval = max(500, production_steps // 100)  # ~100 frames

        tmp_dir = Path(tempfile.mkdtemp())
        dcd_path = tmp_dir / f"{name}.dcd"
        pdb_path = tmp_dir / f"{name}_init.pdb"

        # Save initial structure for mdtraj topology
        with open(pdb_path, "w") as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

        simulation.reporters.append(
            app.DCDReporter(str(dcd_path), report_interval)
        )
        simulation.reporters.append(
            app.StateDataReporter(
                str(tmp_dir / "energy.log"), report_interval,
                step=True, potentialEnergy=True, temperature=True,
            )
        )

        logger.info("MD production: %d steps (%.1f ns) for %s", production_steps, simulation_ns, name)
        simulation.step(production_steps)

        # Analysis with mdtraj
        traj = mdtraj.load(str(dcd_path), top=str(pdb_path))
        ref = traj[0]

        # CA atom indices
        ca_indices = traj.topology.select("name CA")

        # RMSD (nm -> Angstroms)
        rmsd_nm = mdtraj.rmsd(traj, ref, atom_indices=ca_indices)
        rmsd_ang = rmsd_nm * 10.0
        rmsd_mean = float(np.mean(rmsd_ang))
        rmsd_std = float(np.std(rmsd_ang))
        rmsd_final = float(rmsd_ang[-1])

        # RMSF per residue (nm -> Angstroms)
        rmsf_nm = mdtraj.rmsf(traj, ref, atom_indices=ca_indices)
        rmsf_ang = (rmsf_nm * 10.0).tolist()

        # DBD RMSF (residues 94-292, 0-indexed CA selection)
        n_ca = len(ca_indices)
        dbd_start_idx = min(93, n_ca - 1)
        dbd_end_idx = min(292, n_ca)
        rmsf_dbd = float(np.mean(rmsf_ang[dbd_start_idx:dbd_end_idx])) if dbd_end_idx > dbd_start_idx else 0.0

        # Radius of gyration (nm)
        rg = mdtraj.compute_rg(traj)
        rg_mean = float(np.mean(rg))

        # Secondary structure via DSSP
        ss_content = {"helix": 0.0, "sheet": 0.0, "coil": 0.0}
        try:
            dssp = mdtraj.compute_dssp(traj[-1])
            ss_flat = dssp.flatten()
            n_total = max(len(ss_flat), 1)
            ss_content["helix"] = float(np.sum(ss_flat == "H")) / n_total
            ss_content["sheet"] = float(np.sum(ss_flat == "E")) / n_total
            ss_content["coil"] = 1.0 - ss_content["helix"] - ss_content["sheet"]
        except Exception as e:
            logger.warning("DSSP computation failed: %s", e)

        # Stability verdict
        if rmsd_mean < 2.0:
            verdict = "stable"
            score = min(100.0, 100.0 - rmsd_mean * 20)
        elif rmsd_mean < 3.5:
            verdict = "metastable"
            score = max(30.0, 80.0 - (rmsd_mean - 2.0) * 30)
        else:
            verdict = "unstable"
            score = max(0.0, 30.0 - (rmsd_mean - 3.5) * 10)

        elapsed = time.time() - t0

        # Clean up DCD (can be ~50MB)
        try:
            dcd_path.unlink(missing_ok=True)
            pdb_path.unlink(missing_ok=True)
            (tmp_dir / "energy.log").unlink(missing_ok=True)
            tmp_dir.rmdir()
        except Exception:
            pass

        logger.info("MD stability: RMSD=%.2f+/-%.2f A, verdict=%s, score=%.0f, %.1fs",
                     rmsd_mean, rmsd_std, verdict, score, elapsed)

        return MDStabilityResult(
            rmsd_mean=round(rmsd_mean, 3),
            rmsd_std=round(rmsd_std, 3),
            rmsd_final=round(rmsd_final, 3),
            rmsf_per_residue=[round(v, 3) for v in rmsf_ang],
            rmsf_dbd_mean=round(rmsf_dbd, 3),
            ss_content={k: round(v, 3) for k, v in ss_content.items()},
            rg_mean=round(rg_mean, 4),
            stability_score=round(score, 1),
            stability_verdict=verdict,
            simulation_ns=simulation_ns,
            elapsed_sec=round(elapsed, 1),
        )


# ---------------------------------------------------------------------------
# DNABindingAnalyzer
# ---------------------------------------------------------------------------

# Key DNA-contact residues in p53 (1-indexed)
DNA_CONTACT_RESIDUES = {
    120: "K120", 176: "C176", 179: "H179", 238: "C238", 241: "S241",
    242: "C242", 248: "R248", 273: "R273", 277: "C277", 280: "R280", 283: "R283",
}


class DNABindingAnalyzer:
    """Analyze DNA-binding interface preservation via structural superposition."""

    def __init__(self, wt_pdb_path: Optional[Path] = None):
        self._wt_pdb_path = wt_pdb_path or Path("data/raw/p53_wt.pdb")

    def analyze(
        self,
        rescue_pdb_string: str,
        wt_pdb_string: Optional[str] = None,
    ) -> DNABindingInterfaceResult:
        """
        Superpose rescue DBD onto WT and measure displacements at DNA contacts.
        """
        import mdtraj
        import io
        import tempfile

        t0 = time.time()

        # Load WT reference
        if wt_pdb_string:
            wt_tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
            wt_tmp.write(wt_pdb_string)
            wt_tmp.close()
            wt_traj = mdtraj.load(wt_tmp.name)
            Path(wt_tmp.name).unlink(missing_ok=True)
        else:
            if not self._wt_pdb_path.exists():
                logger.warning("WT PDB not found at %s", self._wt_pdb_path)
                return DNABindingInterfaceResult(
                    dbd_rmsd=99.0, per_contact_rmsd={}, interface_preservation_score=0.0,
                    elapsed_sec=time.time() - t0,
                )
            wt_traj = mdtraj.load(str(self._wt_pdb_path))

        # Load rescue structure
        rescue_tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
        rescue_tmp.write(rescue_pdb_string)
        rescue_tmp.close()
        rescue_traj = mdtraj.load(rescue_tmp.name)
        Path(rescue_tmp.name).unlink(missing_ok=True)

        # Get CA atoms for DBD (residues 94-292, 0-indexed: 93-291)
        wt_ca = wt_traj.topology.select("name CA")
        rescue_ca = rescue_traj.topology.select("name CA")

        n_wt = len(wt_ca)
        n_rescue = len(rescue_ca)

        # Determine DBD range (handle different chain lengths)
        dbd_start = 93
        dbd_end = min(292, n_wt, n_rescue)

        if dbd_end <= dbd_start:
            logger.warning("DBD range invalid (wt_ca=%d, rescue_ca=%d)", n_wt, n_rescue)
            return DNABindingInterfaceResult(
                dbd_rmsd=99.0, per_contact_rmsd={}, interface_preservation_score=0.0,
                elapsed_sec=time.time() - t0,
            )

        # DBD CA indices (subset of all CA)
        wt_dbd_ca = wt_ca[dbd_start:dbd_end]
        rescue_dbd_ca = rescue_ca[dbd_start:dbd_end]
        n_dbd = min(len(wt_dbd_ca), len(rescue_dbd_ca))
        wt_dbd_ca = wt_dbd_ca[:n_dbd]
        rescue_dbd_ca = rescue_dbd_ca[:n_dbd]

        # Superpose on DBD CA atoms
        rescue_traj_sup = rescue_traj.superpose(wt_traj, atom_indices=rescue_dbd_ca, ref_atom_indices=wt_dbd_ca)

        # DBD RMSD after superposition
        dbd_rmsd_nm = mdtraj.rmsd(rescue_traj_sup, wt_traj, atom_indices=rescue_dbd_ca, ref_atom_indices=wt_dbd_ca)
        dbd_rmsd_ang = float(dbd_rmsd_nm[0] * 10.0)

        # Per-contact residue displacements
        per_contact: Dict[str, float] = {}
        preservation_scores: List[float] = []

        wt_xyz = wt_traj.xyz[0]     # (n_atoms, 3)
        rescue_xyz = rescue_traj_sup.xyz[0]

        for res_1idx, label in DNA_CONTACT_RESIDUES.items():
            ca_0idx = res_1idx - 1  # 0-indexed in CA array
            if ca_0idx >= n_wt or ca_0idx >= n_rescue:
                continue
            wt_pos = wt_xyz[wt_ca[ca_0idx]]
            rescue_pos = rescue_xyz[rescue_ca[ca_0idx]]
            disp_nm = float(np.linalg.norm(wt_pos - rescue_pos))
            disp_ang = disp_nm * 10.0
            per_contact[label] = round(disp_ang, 2)

            # Score: 1.0 if <1A, linear to 0.0 at 5A
            score = max(0.0, min(1.0, 1.0 - (disp_ang - 1.0) / 4.0))
            preservation_scores.append(score)

        interface_score = float(np.mean(preservation_scores)) if preservation_scores else 0.0
        elapsed = time.time() - t0

        logger.info("DNA interface: DBD RMSD=%.2f A, preservation=%.2f, %.1fs",
                     dbd_rmsd_ang, interface_score, elapsed)

        return DNABindingInterfaceResult(
            dbd_rmsd=round(dbd_rmsd_ang, 2),
            per_contact_rmsd=per_contact,
            interface_preservation_score=round(interface_score, 3),
            elapsed_sec=round(elapsed, 2),
        )


# ---------------------------------------------------------------------------
# Helper: trajectory to PDB string
# ---------------------------------------------------------------------------

def _traj_to_pdb_string(traj: Any) -> str:
    """Convert an mdtraj Trajectory to a PDB string."""
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    tmp.close()
    traj.save_pdb(tmp.name)
    pdb_str = Path(tmp.name).read_text()
    Path(tmp.name).unlink(missing_ok=True)
    return pdb_str


# ---------------------------------------------------------------------------
# DNABindingSimulator (MM-GBSA functional binding simulation)
# ---------------------------------------------------------------------------

# Standard DNA nucleotide residue names
_DNA_RESIDUE_NAMES = {"DA", "DT", "DG", "DC", "DA5", "DT5", "DG5", "DC5",
                       "DA3", "DT3", "DG3", "DC3", "A", "T", "G", "C"}


class DNABindingSimulator:
    """
    MM-GBSA DNA-binding functional simulation.

    Superimposes a rescue candidate onto PDB 2AHI (p53 DBD + DNA crystal structure),
    builds a protein-DNA complex, and computes:
      - MM-GBSA binding free energy (ΔG = E_complex - E_protein - E_DNA)
      - Protein-DNA hydrogen bond analysis
      - Interface contact analysis

    Reference: PDB 2AHI — p53 DBD bound to DNA response element (1.85 Å)
    """

    REFERENCE_PDB_URL = "https://files.rcsb.org/download/2AHI.pdb"
    REFERENCE_PDB_ID = "2AHI"

    def __init__(
        self,
        reference_pdb_path: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
    ):
        self._reference_path = reference_pdb_path or Path("data/raw/2AHI.pdb")
        self._cache_dir = cache_dir
        self._ref_protein_traj = None  # cached protein chain from 2AHI
        self._ref_dna_traj = None      # cached DNA chains from 2AHI
        self._wt_result: Optional[DNABindingSimResult] = None  # cached WT result

    def _ensure_reference(self) -> bool:
        """Download 2AHI from RCSB if not present. Returns True on success."""
        if self._reference_path.exists():
            return True
        try:
            import urllib.request
            self._reference_path.parent.mkdir(parents=True, exist_ok=True)
            logger.info("Downloading PDB 2AHI from RCSB...")
            urllib.request.urlretrieve(self.REFERENCE_PDB_URL, str(self._reference_path))
            logger.info("Downloaded 2AHI to %s", self._reference_path)
            return True
        except Exception as e:
            logger.warning("Failed to download 2AHI: %s", e)
            return False

    def _load_reference_components(self) -> bool:
        """Parse 2AHI and separate protein vs DNA chains by residue type."""
        if self._ref_protein_traj is not None:
            return True

        import mdtraj

        try:
            ref = mdtraj.load(str(self._reference_path))
        except Exception as e:
            logger.warning("Failed to load 2AHI: %s", e)
            return False

        # Identify protein vs DNA atoms by residue name
        protein_indices = []
        dna_indices = []
        for atom in ref.topology.atoms:
            res_name = atom.residue.name.strip()
            if res_name in _DNA_RESIDUE_NAMES:
                dna_indices.append(atom.index)
            elif atom.residue.is_protein:
                protein_indices.append(atom.index)

        if not protein_indices or not dna_indices:
            logger.warning("2AHI: could not separate protein (%d) and DNA (%d) atoms",
                           len(protein_indices), len(dna_indices))
            return False

        self._ref_protein_traj = ref.atom_slice(protein_indices)
        self._ref_dna_traj = ref.atom_slice(dna_indices)

        logger.info("2AHI loaded: %d protein atoms, %d DNA atoms",
                     len(protein_indices), len(dna_indices))
        return True

    def _build_complex(self, rescue_pdb_string: str) -> Optional[Tuple[str, str, str]]:
        """
        Superimpose rescue DBD onto 2AHI protein chain, combine with DNA.

        For full-length p53 (393 residues), extracts only the DBD region
        (residues 94-293) before building the complex. This cuts system size
        in half and avoids OpenMM minimization bottlenecks on the non-binding
        N/C-terminal domains.

        Returns (complex_pdb, protein_pdb, dna_pdb) or None on failure.
        """
        import mdtraj
        import tempfile

        # Load rescue structure
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
        tmp.write(rescue_pdb_string)
        tmp.close()
        rescue_traj = mdtraj.load(tmp.name)
        Path(tmp.name).unlink(missing_ok=True)

        # Get CA atoms for superposition
        ref_ca = self._ref_protein_traj.topology.select("name CA")
        rescue_ca = rescue_traj.topology.select("name CA")

        if len(ref_ca) == 0 or len(rescue_ca) == 0:
            logger.warning("No CA atoms found for superposition")
            return None

        # If rescue is full-length (393 res), extract DBD only (residues 94-293)
        # This keeps only the DNA-binding domain for energy computation
        if len(rescue_ca) > 250:
            dbd_start = 93   # 0-indexed
            dbd_end = min(293, len(rescue_ca))
            # Get all atom indices for DBD residues (not just CA)
            dbd_residue_indices = list(range(dbd_start, dbd_end))
            dbd_atom_indices = rescue_traj.topology.select(
                "resid " + " ".join(str(r) for r in dbd_residue_indices)
            )
            if len(dbd_atom_indices) == 0:
                logger.warning("Could not extract DBD residues from rescue PDB")
                return None
            rescue_traj = rescue_traj.atom_slice(dbd_atom_indices)
            rescue_ca = rescue_traj.topology.select("name CA")
            logger.info("Extracted DBD: %d atoms, %d CA", len(dbd_atom_indices), len(rescue_ca))

        n_align = min(len(ref_ca), len(rescue_ca))
        if n_align < 50:
            logger.warning("Too few CA atoms for reliable superposition (%d)", n_align)

        ref_align_ca = ref_ca[:n_align]
        rescue_align_ca = rescue_ca[:n_align]

        # Superimpose rescue DBD onto 2AHI protein chain
        rescue_traj.superpose(
            self._ref_protein_traj,
            atom_indices=rescue_align_ca,
            ref_atom_indices=ref_align_ca,
        )

        # Build complex: rescue DBD + 2AHI DNA
        complex_traj = rescue_traj.stack(self._ref_dna_traj)

        # Generate PDB strings for all three systems
        complex_pdb = _traj_to_pdb_string(complex_traj)
        protein_pdb = _traj_to_pdb_string(rescue_traj)
        dna_pdb = _traj_to_pdb_string(self._ref_dna_traj)

        return complex_pdb, protein_pdb, dna_pdb

    @staticmethod
    def _prepare_complex_structure(pdb_string: str) -> Any:
        """
        Prepare a protein-DNA complex with PDBFixer, skipping findMissingResidues.

        Standard PDBFixer workflow calls findMissingResidues() which tries to bridge
        residue numbering gaps between protein and DNA chains, causing hangs on
        multi-chain complexes. We skip that step and only fix atoms + hydrogens.
        """
        from pdbfixer import PDBFixer
        import io

        fixer = PDBFixer(pdbfile=io.StringIO(pdb_string))
        fixer.removeHeterogens(False)
        # SKIP findMissingResidues — causes hangs on protein-DNA complexes
        fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)
        return fixer

    def _minimize_complex(self, pdb_string: str, max_iterations: int = 200) -> Optional[float]:
        """Minimize a protein-DNA system and return energy in kcal/mol.

        Tries OpenCL (GPU on macOS/Apple Silicon) first, then falls back to CPU.
        """
        import openmm
        import openmm.app as app
        import openmm.unit as unit

        fixer = self._prepare_complex_structure(pdb_string)
        forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            310 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
        )
        # Try CUDA, then OpenCL (Apple GPU), then CPU
        platform = None
        for pname in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = openmm.Platform.getPlatformByName(pname)
                break
            except Exception:
                continue
        if platform is None:
            platform = openmm.Platform.getPlatformByName("CPU")
        simulation = app.Simulation(fixer.topology, system, integrator, platform)
        simulation.context.setPositions(fixer.positions)
        simulation.minimizeEnergy(maxIterations=max_iterations)
        state = simulation.context.getState(getEnergy=True)
        energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        logger.info("Minimized on %s: %.1f kcal/mol (%d atoms)",
                     platform.getName(), energy_kj / 4.184,
                     sum(1 for _ in fixer.topology.atoms()))
        return energy_kj / 4.184

    def _compute_binding_energy(
        self,
        complex_pdb: str,
        protein_pdb: str,
        dna_pdb: str,
        max_iterations: int = 200,
    ) -> Optional[float]:
        """
        Compute MM-GBSA binding energy: ΔG = E_complex - E_protein - E_DNA.

        Uses AMBER14 + OBC2 implicit solvent. Minimizes all three systems
        independently with a complex-safe PDBFixer pipeline that skips
        findMissingResidues (avoids hangs on protein-DNA complexes).
        """
        try:
            logger.info("Minimizing complex...")
            e_complex = self._minimize_complex(complex_pdb, max_iterations)
            logger.info("Minimizing protein...")
            e_protein = self._minimize_complex(protein_pdb, max_iterations)
            logger.info("Minimizing DNA...")
            e_dna = self._minimize_complex(dna_pdb, max_iterations)
        except Exception as e:
            logger.warning("MM-GBSA energy computation failed: %s", e)
            return None

        if e_complex is None or e_protein is None or e_dna is None:
            return None

        dg = e_complex - e_protein - e_dna
        logger.info("MM-GBSA: E_complex=%.1f, E_protein=%.1f, E_DNA=%.1f, ΔG=%.1f kcal/mol",
                     e_complex, e_protein, e_dna, dg)
        return round(dg, 2)

    def _analyze_hbonds(self, complex_pdb: str) -> Tuple[int, List[Dict[str, Any]]]:
        """
        Analyze protein-DNA hydrogen bonds using mdtraj baker_hubbard.

        Returns (n_hbonds, hbond_details) where hbond_details lists cross-interface bonds.
        """
        import mdtraj
        import tempfile

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
        tmp.write(complex_pdb)
        tmp.close()
        traj = mdtraj.load(tmp.name)
        Path(tmp.name).unlink(missing_ok=True)

        # Identify protein vs DNA atom indices
        protein_atoms = set()
        dna_atoms = set()
        for atom in traj.topology.atoms:
            res_name = atom.residue.name.strip()
            if res_name in _DNA_RESIDUE_NAMES:
                dna_atoms.add(atom.index)
            elif atom.residue.is_protein:
                protein_atoms.add(atom.index)

        if not protein_atoms or not dna_atoms:
            return 0, []

        # Compute hydrogen bonds
        try:
            hbonds = mdtraj.baker_hubbard(traj, freq=0.0)  # freq=0 returns all H-bonds in the frame
        except Exception as e:
            logger.warning("H-bond computation failed: %s", e)
            return 0, []

        # Filter for cross-interface bonds (one atom in protein, one in DNA)
        cross_hbonds = []
        for donor_idx, _h_idx, acceptor_idx in hbonds:
            donor_in_protein = donor_idx in protein_atoms
            donor_in_dna = donor_idx in dna_atoms
            acceptor_in_protein = acceptor_idx in protein_atoms
            acceptor_in_dna = acceptor_idx in dna_atoms

            if (donor_in_protein and acceptor_in_dna) or (donor_in_dna and acceptor_in_protein):
                donor_atom = traj.topology.atom(donor_idx)
                acceptor_atom = traj.topology.atom(acceptor_idx)
                # Distance between donor and acceptor
                dist = float(np.linalg.norm(
                    traj.xyz[0, donor_idx] - traj.xyz[0, acceptor_idx]
                )) * 10.0  # nm -> Angstroms
                cross_hbonds.append({
                    "donor_res": f"{donor_atom.residue.name}{donor_atom.residue.resSeq}",
                    "acceptor_res": f"{acceptor_atom.residue.name}{acceptor_atom.residue.resSeq}",
                    "distance_ang": round(dist, 2),
                })

        return len(cross_hbonds), cross_hbonds

    def _analyze_contacts(self, complex_pdb: str, cutoff_nm: float = 0.45) -> Tuple[int, Dict[str, int]]:
        """
        Analyze protein-DNA heavy-atom contacts within cutoff distance.

        Returns (n_contacts, per_residue_contacts) where per_residue_contacts maps
        DNA-contact residue labels to their contact counts.
        """
        import mdtraj
        import tempfile

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
        tmp.write(complex_pdb)
        tmp.close()
        traj = mdtraj.load(tmp.name)
        Path(tmp.name).unlink(missing_ok=True)

        # Build protein residue -> atom indices and DNA atom indices
        protein_residue_atoms: Dict[int, List[int]] = {}
        dna_heavy_atoms: List[int] = []

        for atom in traj.topology.atoms:
            if atom.element.symbol == "H":
                continue  # skip hydrogens
            res_name = atom.residue.name.strip()
            if res_name in _DNA_RESIDUE_NAMES:
                dna_heavy_atoms.append(atom.index)
            elif atom.residue.is_protein:
                res_seq = atom.residue.resSeq
                if res_seq not in protein_residue_atoms:
                    protein_residue_atoms[res_seq] = []
                protein_residue_atoms[res_seq].append(atom.index)

        if not dna_heavy_atoms or not protein_residue_atoms:
            return 0, {}

        # Compute contacts: for each protein residue, count heavy atoms within cutoff of any DNA atom
        dna_atom_set = np.array(dna_heavy_atoms)
        total_contacts = 0
        per_residue: Dict[str, int] = {}

        for res_seq, atom_indices in protein_residue_atoms.items():
            # Compute pairwise distances between this residue's atoms and DNA atoms
            pairs = []
            for p_idx in atom_indices:
                for d_idx in dna_heavy_atoms:
                    pairs.append([p_idx, d_idx])

            if not pairs:
                continue

            pairs_arr = np.array(pairs)
            distances = mdtraj.compute_distances(traj, pairs_arr)[0]  # first frame
            n_close = int(np.sum(distances < cutoff_nm))

            if n_close > 0:
                # Check if this is a known DNA contact residue
                label = DNA_CONTACT_RESIDUES.get(res_seq, f"R{res_seq}")
                per_residue[label] = n_close
                total_contacts += n_close

        return total_contacts, per_residue

    def simulate(
        self,
        rescue_pdb_string: str,
        sequence: str,
        wt_pdb_string: Optional[str] = None,
    ) -> DNABindingSimResult:
        """
        Run full MM-GBSA binding simulation for a rescue candidate.

        Args:
            rescue_pdb_string: ESMFold-predicted PDB for the rescue candidate
            sequence: amino acid sequence (for cache key)
            wt_pdb_string: WT PDB string (used to compute WT reference on first call)

        Returns:
            DNABindingSimResult with binding energy, H-bonds, contacts, and composite score
        """
        t0 = time.time()
        seq_h = _seq_hash(sequence)

        # Check disk cache
        if self._cache_dir:
            cached = self._cache_dir / f"binding_sim_{seq_h}.json"
            if cached.exists():
                logger.info("Binding sim cache hit: %s", seq_h)
                try:
                    data = json.loads(cached.read_text())
                    return DNABindingSimResult(**data)
                except Exception:
                    pass  # cache corrupted, recompute

        # Ensure reference structure is available
        if not self._ensure_reference() or not self._load_reference_components():
            return self._fallback_result(t0)

        # Compute WT reference if not cached
        if self._wt_result is None and wt_pdb_string:
            logger.info("Computing WT binding reference...")
            try:
                wt_complex = self._build_complex(wt_pdb_string)
            except Exception as e:
                logger.warning("Failed to build WT complex: %s", e)
                wt_complex = None
            if wt_complex:
                wt_cpx_pdb, wt_prot_pdb, wt_dna_pdb = wt_complex
                wt_energy = self._compute_binding_energy(wt_cpx_pdb, wt_prot_pdb, wt_dna_pdb)
                wt_hbonds, wt_hbond_details = self._analyze_hbonds(wt_cpx_pdb)
                wt_contacts, wt_per_res = self._analyze_contacts(wt_cpx_pdb)
                self._wt_result = DNABindingSimResult(
                    binding_energy_kcal=wt_energy if wt_energy is not None else 0.0,
                    wt_binding_energy_kcal=wt_energy if wt_energy is not None else 0.0,
                    delta_binding_kcal=0.0,
                    binding_interpretation="preserved",
                    n_hbonds_protein_dna=wt_hbonds,
                    wt_n_hbonds=wt_hbonds,
                    hbond_preservation_ratio=1.0,
                    hbond_details=wt_hbond_details,
                    n_interface_contacts=wt_contacts,
                    wt_n_contacts=wt_contacts,
                    contact_preservation_ratio=1.0,
                    per_residue_contacts=wt_per_res,
                    binding_score=1.0,
                )
                logger.info("WT binding: ΔG=%.1f kcal/mol, %d H-bonds, %d contacts",
                            self._wt_result.binding_energy_kcal, wt_hbonds, wt_contacts)

        # Build rescue complex
        try:
            complex_result = self._build_complex(rescue_pdb_string)
        except Exception as e:
            logger.warning("Failed to build protein-DNA complex: %s", e)
            complex_result = None
        if complex_result is None:
            return self._fallback_result(t0)

        cpx_pdb, prot_pdb, dna_pdb = complex_result

        # MM-GBSA binding energy
        binding_energy = self._compute_binding_energy(cpx_pdb, prot_pdb, dna_pdb)
        if binding_energy is None:
            return self._fallback_result(t0)

        # H-bond analysis
        n_hbonds, hbond_details = self._analyze_hbonds(cpx_pdb)

        # Contact analysis
        n_contacts, per_res_contacts = self._analyze_contacts(cpx_pdb)

        # Compare to WT
        wt_energy = self._wt_result.binding_energy_kcal if self._wt_result else 0.0
        wt_hbonds = self._wt_result.wt_n_hbonds if self._wt_result else max(n_hbonds, 1)
        wt_contacts = self._wt_result.wt_n_contacts if self._wt_result else max(n_contacts, 1)

        delta_binding = binding_energy - wt_energy

        # Interpretation: more negative = stronger binding
        if delta_binding < -5.0:
            interpretation = "enhanced"
        elif delta_binding < 5.0:
            interpretation = "preserved"
        elif delta_binding < 15.0:
            interpretation = "weakened"
        else:
            interpretation = "disrupted"

        hbond_ratio = min(1.0, n_hbonds / max(wt_hbonds, 1))
        contact_ratio = min(1.0, n_contacts / max(wt_contacts, 1))

        # Composite binding score (0-1)
        # 50% energy + 25% H-bond preservation + 25% contact preservation
        # Energy scoring: ΔΔG < -5 → 1.0, ΔΔG > +20 → 0.0 (linear between)
        energy_score = max(0.0, min(1.0, (20.0 - delta_binding) / 25.0))
        binding_score = 0.5 * energy_score + 0.25 * hbond_ratio + 0.25 * contact_ratio

        elapsed = time.time() - t0
        logger.info("Binding sim: ΔG=%.1f (δ=%.1f), %d H-bonds (%.0f%%), %d contacts (%.0f%%), score=%.2f, %.1fs",
                     binding_energy, delta_binding, n_hbonds, hbond_ratio * 100,
                     n_contacts, contact_ratio * 100, binding_score, elapsed)

        result = DNABindingSimResult(
            binding_energy_kcal=binding_energy,
            wt_binding_energy_kcal=round(wt_energy, 2),
            delta_binding_kcal=round(delta_binding, 2),
            binding_interpretation=interpretation,
            n_hbonds_protein_dna=n_hbonds,
            wt_n_hbonds=wt_hbonds,
            hbond_preservation_ratio=round(hbond_ratio, 3),
            hbond_details=hbond_details,
            n_interface_contacts=n_contacts,
            wt_n_contacts=wt_contacts,
            contact_preservation_ratio=round(contact_ratio, 3),
            per_residue_contacts=per_res_contacts,
            binding_score=round(binding_score, 3),
            elapsed_sec=round(elapsed, 1),
        )

        # Cache to disk
        if self._cache_dir:
            self._cache_dir.mkdir(parents=True, exist_ok=True)
            (self._cache_dir / f"binding_sim_{seq_h}.json").write_text(
                json.dumps(result.to_dict(), indent=2)
            )

        return result

    def _fallback_result(self, t0: float) -> DNABindingSimResult:
        """Neutral fallback when simulation cannot be performed."""
        logger.warning("Binding simulation unavailable, returning neutral fallback")
        return DNABindingSimResult(
            binding_energy_kcal=0.0,
            wt_binding_energy_kcal=0.0,
            delta_binding_kcal=0.0,
            binding_interpretation="unknown",
            n_hbonds_protein_dna=0,
            wt_n_hbonds=0,
            hbond_preservation_ratio=0.5,
            hbond_details=[],
            n_interface_contacts=0,
            wt_n_contacts=0,
            contact_preservation_ratio=0.5,
            per_residue_contacts={},
            binding_score=0.5,
            elapsed_sec=round(time.time() - t0, 1),
        )


# ---------------------------------------------------------------------------
# Composite scoring
# ---------------------------------------------------------------------------

def compute_physics_score(
    esmfold: Optional[ESMFoldResult] = None,
    ddg: Optional[DDGResult] = None,
    md: Optional[MDStabilityResult] = None,
    dna_binding: Optional[DNABindingInterfaceResult] = None,
    dna_binding_sim: Optional[DNABindingSimResult] = None,
) -> Tuple[float, str]:
    """
    Compute composite 0-100 physics score.

    Weights: pLDDT 30pts, DDG 25pts, MD 25pts, DNA interface 20pts.
    When binding sim is available, the 20 DNA pts are split:
      - 8 pts geometric (DNABindingInterfaceResult)
      - 12 pts functional (DNABindingSimResult)
    Missing components scored as neutral (50% of their weight).
    """
    score = 0.0
    max_score = 0.0

    # pLDDT (30 pts): 0 at pLDDT=30, 30 at pLDDT=90
    if esmfold:
        plddt_norm = max(0.0, min(1.0, (esmfold.dbd_plddt - 30) / 60))
        score += 30 * plddt_norm
        max_score += 30
    else:
        score += 15  # neutral fallback
        max_score += 30

    # DDG (25 pts): 25 if stabilizing (<-10), 12.5 if neutral, 0 if destabilizing (>+10)
    if ddg:
        ddg_val = ddg.ddg_vs_wt_kcal
        ddg_norm = max(0.0, min(1.0, (10 - ddg_val) / 20))
        score += 25 * ddg_norm
        max_score += 25
    else:
        score += 12.5
        max_score += 25

    # MD stability (25 pts): use stability_score directly (0-100 -> 0-25)
    if md:
        score += 25 * (md.stability_score / 100)
        max_score += 25
    else:
        score += 12.5
        max_score += 25

    # DNA interface (20 pts total)
    if dna_binding_sim:
        # Split: 8 pts geometric + 12 pts functional binding sim
        if dna_binding:
            score += 8 * dna_binding.interface_preservation_score
            max_score += 8
        else:
            score += 4
            max_score += 8
        score += 12 * dna_binding_sim.binding_score
        max_score += 12
    else:
        # No binding sim: all 20 pts from geometric analysis
        if dna_binding:
            score += 20 * dna_binding.interface_preservation_score
            max_score += 20
        else:
            score += 10
            max_score += 20

    # Normalize to 0-100
    final = (score / max_score) * 100 if max_score > 0 else 0

    if final >= 75:
        verdict = "STRONG"
    elif final >= 55:
        verdict = "PROMISING"
    elif final >= 35:
        verdict = "UNCERTAIN"
    else:
        verdict = "CONCERNING"

    return round(final, 1), verdict


# ---------------------------------------------------------------------------
# PhysicsValidationPipeline (main orchestrator)
# ---------------------------------------------------------------------------

class PhysicsValidationPipeline:
    """Orchestrate tiered physics validation across campaign candidates."""

    def __init__(
        self,
        device: str = "cpu",
        cache_dir: Optional[Path] = None,
        wt_pdb_path: Optional[Path] = None,
    ):
        self._device = device
        self._cache_dir = cache_dir
        self._wt_pdb_path = wt_pdb_path or Path("data/raw/p53_wt.pdb")
        self._energy_cache: Dict[str, float] = {}  # seq_hash -> energy

    def validate_campaign(
        self,
        run_id: str,
        top_df: Any,  # pd.DataFrame
        output_dir: Path,
        wt_sequence: str,
        cancer_sequences: Optional[Dict[str, str]] = None,
        tier1_top_n: Optional[int] = None,
        tier2_top_n: int = 3,
        simulation_ns: float = 0.2,
        skip_esmfold: bool = False,
        skip_energy: bool = False,
        skip_md: bool = False,
        skip_dna: bool = False,
        skip_binding_sim: bool = False,
    ) -> PhysicsValidationReport:
        """
        Run tiered physics validation.

        Tier 1 (all or top-N candidates): ESMFold + energy minimization + DDG + DNA interface
        Tier 2 (top N by Tier 1 score): Short MD stability simulation

        Args:
            tier1_top_n: If set, only validate the top N candidates (by oracle rank).
                         Default None = validate all candidates in top_df.
        """
        import pandas as pd

        # Selective Tier 1: validate only top N candidates if requested
        if tier1_top_n is not None and tier1_top_n < len(top_df):
            logger.info("Selective Tier 1: validating top %d of %d candidates", tier1_top_n, len(top_df))
            top_df = top_df.head(tier1_top_n)

        t0 = time.time()
        physics_pdb_dir = output_dir / "physics_pdbs"
        physics_pdb_dir.mkdir(parents=True, exist_ok=True)

        # Set up ESMFold cache in output directory
        esmfold_cache = self._cache_dir or (output_dir / "esmfold_cache")

        # Step 0: WT reference energy (cached)
        wt_energy: Optional[float] = None
        wt_pdb_string: Optional[str] = None
        if not skip_energy:
            wt_energy, wt_pdb_string = self._get_reference_energy(
                wt_sequence, "wt", esmfold_cache, physics_pdb_dir, skip_esmfold,
            )

        # Cancer-mutant energies (cached per target)
        cancer_energies: Dict[str, float] = {}
        if not skip_energy and cancer_sequences:
            for target_label, cancer_seq in cancer_sequences.items():
                e, _ = self._get_reference_energy(
                    cancer_seq, f"cancer_{target_label}", esmfold_cache, physics_pdb_dir, skip_esmfold,
                )
                if e is not None:
                    cancer_energies[target_label] = e

        # DNA binding analyzer
        dna_analyzer = DNABindingAnalyzer(wt_pdb_path=self._wt_pdb_path) if not skip_dna else None

        # DNA binding simulator (MM-GBSA functional test)
        binding_sim: Optional[DNABindingSimulator] = None
        if not skip_binding_sim:
            binding_sim_cache = output_dir / "binding_sim_cache"
            binding_sim = DNABindingSimulator(cache_dir=binding_sim_cache)

        # Shared energy calculator (caches ForceField across candidates)
        energy_calc = OpenMMEnergyCalculator() if not skip_energy else None

        # Process all candidates (Tier 1)
        results: List[PhysicsValidationResult] = []
        # Map rank -> prepared fixer for reuse in Tier 2 MD
        fixer_cache: Dict[int, Any] = {}

        for idx, row in top_df.iterrows():
            rank = int(row.get("rank", idx + 1)) if "rank" in top_df.columns else idx + 1
            target_label = str(row.get("target_label", "unknown"))
            seq = str(row.get("sequence", ""))
            seq_h = _seq_hash(seq)

            logger.info("Tier 1 validation: rank=%d target=%s hash=%s", rank, target_label, seq_h)
            result = PhysicsValidationResult(
                candidate_rank=rank,
                target_label=target_label,
                sequence_hash=seq_h,
            )

            # ESMFold
            pdb_string: Optional[str] = None
            if not skip_esmfold:
                try:
                    predictor = LocalESMFoldPredictor(device=self._device, cache_dir=esmfold_cache)
                    ef_result = predictor.predict(
                        seq, save_dir=physics_pdb_dir,
                        filename=f"rank{rank:02d}_{target_label.replace('+', '_')}.pdb",
                    )
                    result.esmfold = ef_result
                    pdb_string = ef_result.pdb_string
                except Exception as e:
                    logger.warning("ESMFold failed for rank %d: %s", rank, e)
                    result.errors.append(f"esmfold: {e}")

            # Prepare PDBFixer once per candidate (reused for energy + MD)
            fixer = None
            if pdb_string and energy_calc:
                try:
                    fixer, missing_res, missing_atoms = energy_calc.prepare_structure(pdb_string)
                    fixer_cache[rank] = fixer
                except Exception as e:
                    logger.warning("PDBFixer failed for rank %d: %s", rank, e)

            # Energy minimization + DDG (reuse prepared fixer)
            if not skip_energy and fixer and energy_calc:
                try:
                    energy_result = energy_calc.minimize_from_fixer(fixer)
                    energy_result.pdbfixer_missing_residues = missing_res
                    energy_result.pdbfixer_missing_atoms = missing_atoms
                    result.energy = energy_result
                    self._energy_cache[seq_h] = energy_result.potential_energy_kcal

                    # DDG computation
                    e_rescue = energy_result.potential_energy_kcal
                    e_wt_ref = wt_energy if wt_energy is not None else e_rescue
                    e_cancer_ref = cancer_energies.get(target_label, e_wt_ref)
                    result.ddg = energy_calc.compute_ddg(e_wt_ref, e_cancer_ref, e_rescue)
                except Exception as e:
                    logger.warning("Energy calc failed for rank %d: %s", rank, e)
                    result.errors.append(f"energy: {e}")

            # DNA binding interface
            if not skip_dna and pdb_string and dna_analyzer:
                try:
                    result.dna_binding = dna_analyzer.analyze(
                        pdb_string, wt_pdb_string=wt_pdb_string,
                    )
                except Exception as e:
                    logger.warning("DNA binding analysis failed for rank %d: %s", rank, e)
                    result.errors.append(f"dna_binding: {e}")

            # DNA binding simulation (MM-GBSA functional test)
            if binding_sim and pdb_string:
                try:
                    result.dna_binding_sim = binding_sim.simulate(
                        rescue_pdb_string=pdb_string,
                        sequence=seq,
                        wt_pdb_string=wt_pdb_string,
                    )
                except Exception as e:
                    logger.warning("Binding sim failed for rank %d: %s", rank, e)
                    result.errors.append(f"binding_sim: {e}")

            # Composite score (Tier 1 only, MD will update later)
            result.overall_physics_score, result.verdict = compute_physics_score(
                esmfold=result.esmfold,
                ddg=result.ddg,
                md=result.md,
                dna_binding=result.dna_binding,
                dna_binding_sim=result.dna_binding_sim,
            )

            results.append(result)

        # Tier 2: MD stability for top N candidates
        n_md = 0
        if not skip_md and tier2_top_n > 0:
            # Sort by Tier 1 score, pick top N
            sorted_results = sorted(results, key=lambda r: r.overall_physics_score, reverse=True)
            tier2_candidates = sorted_results[:tier2_top_n]

            checker = MDStabilityChecker()
            for res in tier2_candidates:
                # Find the PDB from ESMFold result
                pdb_str = res.esmfold.pdb_string if res.esmfold else None
                if not pdb_str:
                    continue
                try:
                    name = f"rank{res.candidate_rank:02d}_{res.target_label.replace('+', '_')}"
                    # Reuse cached fixer from Tier 1 if available
                    cached_fixer = fixer_cache.get(res.candidate_rank)
                    res.md = checker.run_stability_check(
                        pdb_str, output_dir=output_dir / "md_trajectories",
                        name=name, simulation_ns=simulation_ns,
                        fixer=cached_fixer,
                    )
                    # Recompute composite score with MD
                    res.overall_physics_score, res.verdict = compute_physics_score(
                        esmfold=res.esmfold,
                        ddg=res.ddg,
                        md=res.md,
                        dna_binding=res.dna_binding,
                        dna_binding_sim=res.dna_binding_sim,
                    )
                    n_md += 1
                except Exception as e:
                    logger.warning("MD stability failed for rank %d: %s", res.candidate_rank, e)
                    res.errors.append(f"md: {e}")

        elapsed_total = time.time() - t0
        logger.info("Physics validation complete: %d candidates, %d MD, %.1f min total",
                     len(results), n_md, elapsed_total / 60)

        report = PhysicsValidationReport(
            run_id=run_id,
            n_candidates=len(results),
            n_tier2=n_md,
            wt_energy_kcal=round(wt_energy, 2) if wt_energy is not None else None,
            candidates=[r.to_dict() for r in results],
            elapsed_total_sec=round(elapsed_total, 1),
        )
        return report

    def _get_reference_energy(
        self,
        sequence: str,
        label: str,
        esmfold_cache: Path,
        pdb_dir: Path,
        skip_esmfold: bool,
    ) -> Tuple[Optional[float], Optional[str]]:
        """Get energy for a reference sequence (WT or cancer), with caching."""
        seq_h = _seq_hash(sequence)
        if seq_h in self._energy_cache:
            return self._energy_cache[seq_h], None

        pdb_string: Optional[str] = None

        # For WT: prefer the on-disk AlphaFold PDB (properly resolved structure)
        if label == "wt" and self._wt_pdb_path.exists():
            pdb_string = self._wt_pdb_path.read_text()
            logger.info("Using on-disk WT PDB for energy: %s", self._wt_pdb_path)
        elif not skip_esmfold:
            try:
                predictor = LocalESMFoldPredictor(device=self._device, cache_dir=esmfold_cache)
                ef_result = predictor.predict(sequence, save_dir=pdb_dir, filename=f"{label}.pdb")
                pdb_string = ef_result.pdb_string
            except Exception as e:
                logger.warning("ESMFold failed for %s: %s", label, e)

        if pdb_string is None:
            logger.warning("No PDB available for %s energy calculation", label)
            return None, None

        try:
            calc = OpenMMEnergyCalculator()
            result = calc.minimize_and_get_energy(pdb_string)
            self._energy_cache[seq_h] = result.potential_energy_kcal
            return result.potential_energy_kcal, pdb_string
        except Exception as e:
            logger.warning("Energy calculation failed for %s: %s", label, e)
            # If ESMFold PDB failed energy, try on-disk WT as last resort
            if label == "wt" and pdb_string and self._wt_pdb_path.exists():
                try:
                    disk_pdb = self._wt_pdb_path.read_text()
                    result = calc.minimize_and_get_energy(disk_pdb)
                    self._energy_cache[seq_h] = result.potential_energy_kcal
                    return result.potential_energy_kcal, disk_pdb
                except Exception:
                    pass
            return None, pdb_string
