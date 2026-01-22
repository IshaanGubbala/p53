from __future__ import annotations

import json
import re
from pathlib import Path
from statistics import mean
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

USER_AGENT = "p53-stabilimut/0.1"
ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/files"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
DEFAULT_VERSION = "v4"
ALPHAFOLD_VERSIONS = [DEFAULT_VERSION, "v3", "v2", "v1"]


def _download(url: str, out_path: Path) -> None:
    if out_path.exists():
        return


    out_path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            payload = response.read()
    except HTTPError as exc:
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    out_path.write_bytes(payload)


def _alphafold_prefix(uniprot_id: str) -> str:
    return f"AF-{uniprot_id}-F1"


def _alphafold_model_filename(uniprot_id: str, version: str) -> str:
    return f"{_alphafold_prefix(uniprot_id)}-model_{version}.pdb"


def _alphafold_pae_filename(uniprot_id: str, version: str) -> str:
    return f"{_alphafold_prefix(uniprot_id)}-predicted_aligned_error_{version}.json"


def _iter_versions(version: str | None) -> list[str]:
    if version:
        return [version]
    return list(ALPHAFOLD_VERSIONS)


def _fetch_json(url: str) -> object:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            payload = response.read()
    except HTTPError as exc:
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    try:
        return json.loads(payload.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON from {url}") from exc


def _fetch_prediction_metadata(uniprot_id: str) -> list[dict[str, object]]:
    url = f"{ALPHAFOLD_API}/{uniprot_id}"
    data = _fetch_json(url)
    if not isinstance(data, list) or not data:
        raise RuntimeError(f"No AlphaFold predictions returned for {uniprot_id}")
    return data


def _select_prediction(predictions: list[dict[str, object]], uniprot_id: str) -> dict[str, object]:
    for prediction in predictions:
        if prediction.get("uniprotAccession") == uniprot_id:
            return prediction
    return predictions[0]


def _version_from_url(url: str) -> str | None:
    match = re.search(r"model_(v\\d+)", url)
    return match.group(1) if match else None


def _path_from_url(url: str) -> Path:
    parsed = urlparse(url)
    return Path(parsed.path).name


def download_alphafold_pdb(
    uniprot_id: str,
    out_dir: Path,
    version: str | None = None,
) -> tuple[Path, str, str | None]:
    """Download the AlphaFold PDB for a UniProt ID and return the file path, version, and PAE URL."""
    last_error: Exception | None = None
    if version is None:
        try:
            predictions = _fetch_prediction_metadata(uniprot_id)
            prediction = _select_prediction(predictions, uniprot_id)
            pdb_url = prediction.get("pdbUrl")
            pae_url = prediction.get("paeUrl")
            if isinstance(pdb_url, str) and pdb_url:
                filename = _path_from_url(pdb_url)
                out_path = out_dir / filename
                _download(pdb_url, out_path)
                resolved_version = prediction.get("modelVersion")
                if isinstance(resolved_version, str) and resolved_version:
                    version_name = resolved_version
                else:
                    version_name = _version_from_url(pdb_url) or DEFAULT_VERSION
                return out_path, version_name, pae_url if isinstance(pae_url, str) else None
        except FileNotFoundError as exc:
            last_error = exc
        except Exception as exc:  # noqa: BLE001
            last_error = exc

    for ver in _iter_versions(version):
        filename = _alphafold_model_filename(uniprot_id, ver)
        out_path = out_dir / filename
        url = f"{ALPHAFOLD_BASE}/{filename}"
        try:
            _download(url, out_path)
            pae_url = f"{ALPHAFOLD_BASE}/{_alphafold_pae_filename(uniprot_id, ver)}"
            return out_path, ver, pae_url
        except FileNotFoundError:
            continue
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            break

    raise RuntimeError(
        f"Failed to download AlphaFold PDB for {uniprot_id} "
        f"(tried versions: {', '.join(_iter_versions(version))})"
    ) from last_error


def download_alphafold_pae(
    uniprot_id: str,
    out_dir: Path,
    version: str,
    pae_url: str | None = None,
    optional: bool = True,
) -> Path | None:
    """Download the AlphaFold predicted aligned error JSON."""
    if pae_url:
        filename = _path_from_url(pae_url)
        url = pae_url
    else:
        filename = _alphafold_pae_filename(uniprot_id, version)
        url = f"{ALPHAFOLD_BASE}/{filename}"
    out_path = out_dir / filename
    try:
        _download(url, out_path)
    except FileNotFoundError:
        if optional:
            return None
        raise
    return out_path


def verify_model(out_pdb: Path) -> dict[str, object]:
    """Summarize residue counts and pLDDT statistics from an AlphaFold PDB."""
    residues: set[tuple[str, int]] = set()
    chains: set[str] = set()
    plddt_values: list[float] = []

    with out_pdb.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip() or "?"
            resseq_txt = line[22:26].strip()
            atom_name = line[12:16].strip()
            bfactor_txt = line[60:66].strip()

            if not resseq_txt.isdigit():
                continue

            resseq = int(resseq_txt)
            residues.add((chain_id, resseq))
            chains.add(chain_id)

            if atom_name == "CA" and bfactor_txt:
                try:
                    plddt_values.append(float(bfactor_txt))
                except ValueError:
                    continue

    residue_numbers = sorted({resseq for _, resseq in residues})
    missing = 0
    if residue_numbers:
        previous = residue_numbers[0]
        for resseq in residue_numbers[1:]:
            gap = resseq - previous
            if gap > 1:
                missing += gap - 1
            previous = resseq

    plddt_stats = {
        "mean": round(mean(plddt_values), 2) if plddt_values else None,
        "min": min(plddt_values) if plddt_values else None,
        "max": max(plddt_values) if plddt_values else None,
    }

    return {
        "chains": sorted(chains),
        "residues": len(residues),
        "min_residue": min(residue_numbers) if residue_numbers else None,
        "max_residue": max(residue_numbers) if residue_numbers else None,
        "missing_residues": missing,
        "plddt": plddt_stats,
    }


def save_summary(summary: dict[str, object], out_json: Path) -> None:
    """Write AlphaFold model summary to JSON."""
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
