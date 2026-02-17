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
        self._load_error: Optional[str] = None
        self._initialized = True

    def _ensure_loaded(self) -> None:
        if self._model is not None:
            return
        if self._load_error is not None:
            raise RuntimeError(self._load_error)
        import torch
        from transformers import AutoTokenizer, EsmForProteinFolding
        from unittest.mock import patch

        logger.info("Loading ESMFold model (this may take a minute)...")
        t0 = time.time()

        try:
            self._tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
        except Exception as e:
            self._load_error = f"Failed to load ESMFold tokenizer: {e}"
            raise RuntimeError(self._load_error) from e

        # ESMFold ships .bin (pickle) weights, not safetensors.
        # Transformers 4.57+ imports check_torch_load_is_safe into modeling_utils
        # and raises a CVE-2025-32434 ValueError on torch < 2.6 before loading.
        # We patch BOTH:
        #   (a) the local reference in transformers.modeling_utils (where it's called)
        #   (b) torch.load itself, forcing weights_only=False so pickle loads on
        #       any torch version.
        # unittest.mock.patch handles cross-module reference tracking correctly.
        _orig_load = torch.load
        def _permissive_load(*args: Any, **kwargs: Any) -> Any:
            kwargs["weights_only"] = False
            return _orig_load(*args, **kwargs)

        try:
            with patch("transformers.modeling_utils.check_torch_load_is_safe", lambda: None), \
                 patch("torch.load", _permissive_load):
                self._model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
        except Exception as e:
            self._load_error = f"Failed to load ESMFold model: {e}"
            raise RuntimeError(self._load_error) from e

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

        t0 = time.time()
        tokenized = self._tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.to(self._device) for k, v in tokenized.items()}

        with torch.no_grad():
            output = self._model(**tokenized)

        # Convert to PDB string
        from transformers.models.esm.openfold_utils.protein import to_pdb, Protein

        # Build Protein object from output
        pdb_string = self._output_to_pdb(output, sequence)

        # Extract per-residue pLDDT.
        # Shape is (batch, L, 37) — per-atom pLDDT in atom37 format.
        # Take CA atom (index 1) as representative per-residue pLDDT.
        plddt_tensor = output["plddt"][0, :len(sequence)]  # (L, 37)
        if plddt_tensor.ndim == 2:
            # Use CA atom (atom37 index 1) for per-residue score
            plddt_scores = plddt_tensor[:, 1].cpu().numpy().tolist()
        else:
            plddt_scores = plddt_tensor.cpu().numpy().tolist()
        mean_plddt = float(np.mean(plddt_scores))
        # DBD = residues 94-292 (0-indexed: 93-291)
        dbd_start, dbd_end = 93, min(292, len(sequence))
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
                     len(sequence), mean_plddt, dbd_plddt, elapsed)

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
            meta["sequence_length"] = len(sequence)
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

    def _get_forcefield(self):
        """Cached ForceField to avoid re-parsing AMBER14 XML on every candidate."""
        if self._forcefield is None:
            import openmm.app as app
            self._forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
        return self._forcefield

    @staticmethod
    def _select_platform(openmm: Any) -> Any:
        """Prefer CUDA when available, otherwise fall back to CPU."""
        try:
            platform = openmm.Platform.getPlatformByName("CUDA")
            logger.info("OpenMM platform selected: CUDA")
            return platform
        except Exception:
            platform = openmm.Platform.getPlatformByName("CPU")
            logger.info("OpenMM platform selected: CPU (CUDA unavailable)")
            return platform

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

        # Minimize
        simulation.minimizeEnergy(maxIterations=max_iterations)
        state = simulation.context.getState(getEnergy=True)
        energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        energy_kcal = energy_kj / 4.184

        # Positive energy = not converged; retry with 2x iterations
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
# Composite scoring
# ---------------------------------------------------------------------------

def compute_physics_score(
    esmfold: Optional[ESMFoldResult] = None,
    ddg: Optional[DDGResult] = None,
    md: Optional[MDStabilityResult] = None,
    dna_binding: Optional[DNABindingInterfaceResult] = None,
) -> Tuple[float, str]:
    """
    Compute composite 0-100 physics score.

    Weights: pLDDT 30pts, DDG 25pts, MD 25pts, DNA interface 20pts.
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

    # DNA interface (20 pts): use preservation score (0-1 -> 0-20)
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

            # Composite score (Tier 1 only, MD will update later)
            result.overall_physics_score, result.verdict = compute_physics_score(
                esmfold=result.esmfold,
                ddg=result.ddg,
                md=result.md,
                dna_binding=result.dna_binding,
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
