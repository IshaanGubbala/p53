from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from src.core.hashing import file_sha256, json_sha256
from src.core.logging import get_logger

CACHE_PREFIX = "evoef2"
_VALIDATED_BINARIES: set[Path] = set()
_TOTAL_RE = re.compile(r"^Total\s*=\s*([-+]?\d+(?:\.\d+)?)", re.IGNORECASE)  # Fixed: single backslashes in raw string
EVOEF2_TIMEOUT = 300  # 5 minutes timeout for EvoEF2 operations


def canonicalize_mutation_set(mutations: list[str]) -> tuple[str, ...]:
    return tuple(sorted(mutations))


def _normalize_cfg(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _normalize_cfg(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_normalize_cfg(item) for item in value]
    return value


def mutation_set_id(pdb_hash: str, mutations: tuple[str, ...], engine_cfg: dict[str, Any]) -> str:
    payload = {"pdb": pdb_hash, "mutations": mutations, "engine_cfg": _normalize_cfg(engine_cfg)}
    return json_sha256(payload)


def _resolve_evoef2_binary(binary: str | None) -> Path:
    if not binary:
        raise RuntimeError("EvoEF2 binary path is required")
    candidate = Path(binary).expanduser()
    if candidate.is_file():
        if not os.access(candidate, os.X_OK):
            raise RuntimeError(f"EvoEF2 binary is not executable: {candidate}")
        return candidate
    resolved = shutil.which(binary)
    if resolved:
        return Path(resolved)
    raise RuntimeError(f"EvoEF2 binary not found: {binary}")


def _validate_evoef2_install(binary_path: Path) -> None:
    if binary_path in _VALIDATED_BINARIES:
        return
    base = binary_path.parent
    library_dir = base / "library"
    wread_path = base / "wread" / "weight_EvoEF2.txt"
    if not library_dir.is_dir():
        raise RuntimeError(f"EvoEF2 library directory missing: {library_dir}")
    if not wread_path.exists():
        raise RuntimeError(f"EvoEF2 weight file missing: {wread_path}")
    for path in library_dir.iterdir():
        if path.is_file() and path.stat().st_size > 0:
            _VALIDATED_BINARIES.add(binary_path)
            return
    raise RuntimeError("EvoEF2 library directory is empty; re-download EvoEF2")


def _infer_chain_id(pdb_path: Path) -> str:
    chain_residues: dict[str, set[int]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 26:
            continue
        chain_id = line[21].strip() or "A"
        resseq_text = line[22:26].strip()
        if not resseq_text.isdigit():
            continue
        resseq = int(resseq_text)
        chain_residues.setdefault(chain_id, set()).add(resseq)
    if not chain_residues:
        raise RuntimeError("Unable to infer chain ID from PDB")
    return max(chain_residues.items(), key=lambda item: len(item[1]))[0]


def _run_evoef2(command: list[str], workdir: Path) -> str:
    logger = get_logger(__name__)

    # Log the command being executed
    cmd_str = ' '.join(str(c) for c in command)
    logger.debug("Executing EvoEF2 command: %s", cmd_str)
    logger.debug("Working directory: %s", workdir)

    try:
        result = subprocess.run(
            command,
            cwd=workdir,
            check=False,
            capture_output=True,
            text=True,
            timeout=EVOEF2_TIMEOUT,
            stdin=subprocess.DEVNULL,  # Prevent hanging by redirecting stdin
        )
    except subprocess.TimeoutExpired as exc:
        logger.error("EvoEF2 command timed out after %d seconds", EVOEF2_TIMEOUT)
        logger.error("Command: %s", cmd_str)
        raise RuntimeError(
            f"EvoEF2 command timed out after {EVOEF2_TIMEOUT} seconds\n"
            f"Command: {cmd_str}"
        ) from exc

    output = (result.stdout or "") + "\n" + (result.stderr or "")

    # Log output for debugging
    if result.returncode == 0:
        logger.debug("EvoEF2 command succeeded")
        logger.debug("Output (first 500 chars): %s", output[:500])
    else:
        logger.error("EvoEF2 command failed with return code %d", result.returncode)
        logger.error("Full output: %s", output)

    if result.returncode != 0:
        message = (
            f"EvoEF2 command failed ({result.returncode})\n"
            f"Command: {' '.join(command)}\n"
            f"Output: {output}"
        )
        raise RuntimeError(message)
    return output


def _parse_total_energy(output: str) -> float:
    total = None
    for line in output.splitlines():
        match = _TOTAL_RE.search(line.strip())
        if match:
            total = float(match.group(1))
    if total is None:
        raise RuntimeError("EvoEF2 output missing Total energy line")
    return total


def compute_stability(pdb_path: Path, evoef2_cfg: dict[str, Any], workdir: Path) -> float:
    logger = get_logger(__name__)
    logger.info("Computing stability for PDB: %s", pdb_path)

    binary_path = _resolve_evoef2_binary(evoef2_cfg.get("binary"))
    logger.debug("Using EvoEF2 binary: %s", binary_path)
    _validate_evoef2_install(binary_path)

    # EvoEF2 has issues with absolute paths, so copy PDB to workdir and use relative path
    workdir.mkdir(parents=True, exist_ok=True)
    local_pdb = workdir / pdb_path.name
    if not local_pdb.exists() or local_pdb.stat().st_size != pdb_path.stat().st_size:
        shutil.copy2(pdb_path, local_pdb)
        logger.debug("Copied PDB to workdir: %s", local_pdb)

    cmd = [
        str(binary_path),
        "--command=ComputeStability",
        f"--pdb={pdb_path.name}",  # Use just the filename
    ]
    bbdep = evoef2_cfg.get("bbdep")
    if bbdep is True:
        cmd.append("--bbdep=enable")
    elif bbdep is False:
        cmd.append("--bbdep=disable")
    rotlib = evoef2_cfg.get("rotlib")
    if rotlib:
        cmd.append(f"--rotlib={rotlib}")

    output = _run_evoef2(cmd, workdir)
    energy = _parse_total_energy(output)
    logger.info("Computed stability energy: %.4f", energy)
    return energy


def compute_binding(pdb_path: Path, split_chains: str, evoef2_cfg: dict[str, Any], workdir: Path) -> float:
    """
    Compute binding free energy between protein chains using EvoEF2 ComputeBinding.

    Args:
        pdb_path: Path to complex PDB file (e.g., protein-DNA or protein-protein)
        split_chains: Chain splitting specification (e.g., "ABC,EF" for protein vs DNA)
        evoef2_cfg: EvoEF2 configuration dict
        workdir: Working directory for temporary files

    Returns:
        Binding free energy in kcal/mol (negative = favorable binding)

    Example:
        # DNA binding: protein chains A,B,C vs DNA chains E,F
        ddg_binding = compute_binding(pdb_path, "ABC,EF", cfg, workdir)

        # Protein-protein interface: chain A vs chain B
        ddg_interface = compute_binding(pdb_path, "A,B", cfg, workdir)
    """
    logger = get_logger(__name__)
    logger.info("Computing binding energy for PDB: %s (split: %s)", pdb_path, split_chains)

    binary_path = _resolve_evoef2_binary(evoef2_cfg.get("binary"))
    logger.debug("Using EvoEF2 binary: %s", binary_path)
    _validate_evoef2_install(binary_path)

    # Copy PDB to workdir
    workdir.mkdir(parents=True, exist_ok=True)
    local_pdb = workdir / pdb_path.name
    if not local_pdb.exists() or local_pdb.stat().st_size != pdb_path.stat().st_size:
        shutil.copy2(pdb_path, local_pdb)
        logger.debug("Copied PDB to workdir: %s", local_pdb)

    cmd = [
        str(binary_path),
        "--command=ComputeBinding",
        f"--pdb={pdb_path.name}",
        f"--split_chains={split_chains}",
    ]

    bbdep = evoef2_cfg.get("bbdep")
    if bbdep is True:
        cmd.append("--bbdep=enable")
    elif bbdep is False:
        cmd.append("--bbdep=disable")

    rotlib = evoef2_cfg.get("rotlib")
    if rotlib:
        cmd.append(f"--rotlib={rotlib}")

    output = _run_evoef2(cmd, workdir)
    binding_energy = _parse_total_energy(output)
    logger.info("Computed binding energy: %.4f kcal/mol", binding_energy)
    return binding_energy


def _find_mutant_model(workdir: Path) -> Path:
    matches = sorted(workdir.glob("*_Model_*.pdb"))
    if matches:
        return matches[0]
    raise RuntimeError("EvoEF2 BuildMutant did not produce a mutant model PDB")


def _format_mutation(token: str, chain_id: str) -> str:
    match = re.match(r"^([A-Za-z*])(\d+)([A-Za-z*])$", token.strip())
    if not match:
        raise ValueError(f"Invalid mutation token: {token}")
    ref, pos_text, alt = match.groups()
    if len(chain_id) != 1:
        raise ValueError("EvoEF2 chain_id must be a single character")
    return f"{ref.upper()}{chain_id}{int(pos_text)}{alt.upper()}"


def _write_mutant_file(mutations: tuple[str, ...], chain_id: str, out_file: Path) -> None:
    formatted = [_format_mutation(mut, chain_id) for mut in mutations]
    out_file.write_text(",".join(formatted) + ";\n", encoding="utf-8")


def score_mutation_set(
    mutations: list[str],
    pdb_path: Path,
    cache_dir: Path,
    evoef2_cfg: dict[str, Any],
    work_root: Path,
    base_energy: float | None = None,
    recompute: bool = False,
) -> float:
    logger = get_logger(__name__)
    mutations = list(mutations)
    canonical = canonicalize_mutation_set(mutations)
    logger.debug("Scoring mutations: %s", ', '.join(canonical))

    base_pdb = Path(evoef2_cfg.get("repaired_pdb") or pdb_path).expanduser()
    if not base_pdb.is_absolute():
        base_pdb = (Path.cwd() / base_pdb).resolve()
    if not base_pdb.exists():
        raise RuntimeError(f"EvoEF2 base PDB not found: {base_pdb}")

    pdb_hash = file_sha256(base_pdb)
    key = mutation_set_id(pdb_hash, canonical, evoef2_cfg)

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{CACHE_PREFIX}_{key}.json"
    cache_individual_energies = evoef2_cfg.get("cache", {}).get("cache_individual_energies", False)

    if cache_path.exists() and not recompute:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        if "ddg" in payload:
            logger.debug("Using cached ddG for mutations: %s", ', '.join(canonical))
            # Also use cached base_energy if available
            if base_energy is None and "base_energy" in payload:
                base_energy = float(payload["base_energy"])
            return float(payload["ddg"])

    workdir = work_root / key
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    chain_id = evoef2_cfg.get("chain_id") or _infer_chain_id(base_pdb)
    logger.debug("Using chain ID: %s", chain_id)
    mutant_file = workdir / "mutant.txt"
    _write_mutant_file(canonical, chain_id, mutant_file)
    logger.debug("Created mutant file: %s", mutant_file)

    # EvoEF2 has issues with absolute paths, so copy PDB to workdir and use relative path
    local_pdb = workdir / base_pdb.name
    if not local_pdb.exists() or local_pdb.stat().st_size != base_pdb.stat().st_size:
        shutil.copy2(base_pdb, local_pdb)
        logger.debug("Copied base PDB to workdir: %s", local_pdb)

    binary_path = _resolve_evoef2_binary(evoef2_cfg.get("binary"))
    _validate_evoef2_install(binary_path)

    cmd = [
        str(binary_path),
        "--command=BuildMutant",
        f"--pdb={base_pdb.name}",  # Use just the filename
        f"--mutant_file={mutant_file.name}",
    ]
    num_runs = evoef2_cfg.get("num_of_runs")
    if num_runs:
        cmd.append(f"--num_of_runs={int(num_runs)}")
    bbdep = evoef2_cfg.get("bbdep")
    if bbdep is True:
        cmd.append("--bbdep=enable")
    elif bbdep is False:
        cmd.append("--bbdep=disable")
    rotlib = evoef2_cfg.get("rotlib")
    if rotlib:
        cmd.append(f"--rotlib={rotlib}")

    logger.info("Building mutant model for: %s", ', '.join(canonical))
    _run_evoef2(cmd, workdir)
    mutant_pdb = _find_mutant_model(workdir)
    logger.debug("Found mutant model: %s", mutant_pdb)

    if base_energy is None:
        logger.info("Computing base energy (not cached)")
        base_energy = compute_stability(base_pdb, evoef2_cfg, workdir)
    logger.info("Computing mutant energy")
    mut_energy = compute_stability(mutant_pdb, evoef2_cfg, workdir)
    ddg = mut_energy - base_energy
    logger.info("Computed ddG = %.4f (base: %.4f, mutant: %.4f)", ddg, base_energy, mut_energy)

    payload = {
        "pdb_hash": pdb_hash,
        "mutations": canonical,
        "ddg": ddg,
    }

    # Store individual energies if requested (enables better cache reuse)
    if cache_individual_energies:
        payload["base_energy"] = base_energy
        payload["mutant_energy"] = mut_energy

    cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    logger.debug("Cached mutation score at %s", cache_path)
    return ddg


def build_mutant_model(
    mutations: list[str],
    pdb_path: Path,
    evoef2_cfg: dict[str, Any],
    work_root: Path,
    out_path: Path | None = None,
    recompute: bool = False,
) -> Path:
    logger = get_logger(__name__)
    canonical = canonicalize_mutation_set(list(mutations))

    base_pdb = Path(evoef2_cfg.get("repaired_pdb") or pdb_path).expanduser()
    if not base_pdb.is_absolute():
        base_pdb = (Path.cwd() / base_pdb).resolve()
    if not base_pdb.exists():
        raise RuntimeError(f"EvoEF2 base PDB not found: {base_pdb}")

    pdb_hash = file_sha256(base_pdb)
    key = mutation_set_id(pdb_hash, canonical, evoef2_cfg)
    workdir = work_root / f"model_{key}"

    if workdir.exists() and not recompute:
        try:
            existing = _find_mutant_model(workdir)
            if out_path:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(existing, out_path)
                return out_path
            return existing
        except RuntimeError:
            pass

    if workdir.exists() and recompute:
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    chain_id = evoef2_cfg.get("chain_id") or _infer_chain_id(base_pdb)
    mutant_file = workdir / "mutant.txt"
    _write_mutant_file(canonical, chain_id, mutant_file)

    local_pdb = workdir / base_pdb.name
    if not local_pdb.exists() or local_pdb.stat().st_size != base_pdb.stat().st_size:
        shutil.copy2(base_pdb, local_pdb)

    binary_path = _resolve_evoef2_binary(evoef2_cfg.get("binary"))
    _validate_evoef2_install(binary_path)

    cmd = [
        str(binary_path),
        "--command=BuildMutant",
        f"--pdb={base_pdb.name}",
        f"--mutant_file={mutant_file.name}",
    ]
    num_runs = evoef2_cfg.get("num_of_runs")
    if num_runs:
        cmd.append(f"--num_of_runs={int(num_runs)}")
    bbdep = evoef2_cfg.get("bbdep")
    if bbdep is True:
        cmd.append("--bbdep=enable")
    elif bbdep is False:
        cmd.append("--bbdep=disable")
    rotlib = evoef2_cfg.get("rotlib")
    if rotlib:
        cmd.append(f"--rotlib={rotlib}")

    logger.info("Building mutant model for animation: %s", ", ".join(canonical))
    _run_evoef2(cmd, workdir)
    mutant_pdb = _find_mutant_model(workdir)

    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(mutant_pdb, out_path)
        return out_path
    return mutant_pdb
