from __future__ import annotations

from pathlib import Path
from typing import Any

from src.core.logging import get_logger
from src.data_fetch.fetch_clinvar import (
    download_clinvar_vcf,
    download_clinvar_xml_release,
    record_release_metadata,
)
from src.data_fetch.fetch_alphafold import (
    download_alphafold_pdb,
    download_alphafold_pae,
    save_summary,
    verify_model,
)
from src.data_fetch.fetch_uniprot import (
    download_uniprot_entry,
    extract_sequence,
    extract_features,
    save_annotations,
    save_fasta,
)
from src.data_fetch.fetch_nci_tp53 import (
    download_nci_tp53_tables,
    parse_hotspots,
    save_hotspots,
)


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
    }


def _fetch_uniprot(uniprot_id: str, raw_dir: Path, interim_dir: Path) -> None:
    uniprot_raw_dir = raw_dir / "uniprot"
    uniprot_interim_dir = interim_dir / "uniprot"
    annotations_dir = interim_dir / "annotations"

    _ensure_dir(uniprot_raw_dir)
    _ensure_dir(uniprot_interim_dir)
    _ensure_dir(annotations_dir)

    entry_path = uniprot_raw_dir / f"{uniprot_id}.txt"
    download_uniprot_entry(uniprot_id, entry_path)

    entry_text = entry_path.read_text(encoding="utf-8")
    sequence = extract_sequence(entry_text)
    if not sequence:
        raise RuntimeError("Failed to parse UniProt sequence")

    features = extract_features(entry_text)

    fasta_path = uniprot_interim_dir / f"{uniprot_id}.fasta"
    annotations_path = annotations_dir / f"{uniprot_id}.json"

    save_fasta(sequence, fasta_path)
    save_annotations(features, annotations_path)


def _fetch_alphafold(
    uniprot_id: str,
    raw_dir: Path,
    interim_dir: Path,
    alphafold_cfg: dict[str, Any],
) -> None:
    alphafold_raw_dir = raw_dir / "alphafold"
    alphafold_interim_dir = interim_dir / "alphafold"

    _ensure_dir(alphafold_raw_dir)
    _ensure_dir(alphafold_interim_dir)

    version = alphafold_cfg.get("version") if alphafold_cfg else None
    download_pae = True if alphafold_cfg is None else alphafold_cfg.get("download_pae", True)

    pdb_path, resolved_version, pae_url = download_alphafold_pdb(
        uniprot_id,
        alphafold_raw_dir,
        version=version,
    )
    if download_pae:
        download_alphafold_pae(
            uniprot_id,
            alphafold_raw_dir,
            version=resolved_version,
            pae_url=pae_url,
            optional=True,
        )

    summary = verify_model(pdb_path)
    summary_path = alphafold_interim_dir / f"{uniprot_id}_summary_{resolved_version}.json"
    save_summary(summary, summary_path)


def _fetch_clinvar(raw_dir: Path, clinvar_cfg: dict[str, Any]) -> None:
    clinvar_raw_dir = raw_dir / "clinvar"
    _ensure_dir(clinvar_raw_dir)

    download_xml = clinvar_cfg.get("download_xml", True)
    decompress_xml = clinvar_cfg.get("decompress_xml", True)
    allow_parts = clinvar_cfg.get("allow_parts", True)
    part_size_mb = clinvar_cfg.get("part_size_mb", 1024)

    if download_xml:
        download_clinvar_xml_release(
            clinvar_raw_dir,
            decompress=decompress_xml,
            allow_parts=allow_parts,
            part_size_mb=part_size_mb,
        )

    if clinvar_cfg.get("download_vcf", False):
        assembly = clinvar_cfg.get("vcf_assembly", "GRCh38")
        download_clinvar_vcf(
            clinvar_raw_dir,
            assembly=assembly,
            allow_parts=allow_parts,
            part_size_mb=part_size_mb,
        )

    record_release_metadata(clinvar_raw_dir)


def _select_hotspot_table(table_paths: list[Path]) -> Path:
    for keyword in ("hotspot", "hot_spot", "hot-spot"):
        for path in table_paths:
            if keyword in path.name.lower():
                return path
    return table_paths[0]


def _fetch_nci_tp53(raw_dir: Path, interim_dir: Path, nci_cfg: dict[str, Any]) -> None:
    nci_raw_dir = raw_dir / "nci_tp53"
    nci_interim_dir = interim_dir / "nci_tp53"

    _ensure_dir(nci_raw_dir)
    _ensure_dir(nci_interim_dir)

    urls = nci_cfg.get("urls") if nci_cfg else None
    table_paths = download_nci_tp53_tables(nci_raw_dir, urls=urls)
    hotspot_table = _select_hotspot_table(table_paths)
    hotspots = parse_hotspots(hotspot_table)
    save_hotspots(hotspots, nci_interim_dir / "hotspots.json")


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    requested = {
        "uniprot": args.all or args.uniprot,
        "alphafold": args.all or args.alphafold,
        "clinvar": args.all or args.clinvar,
        "nci": args.all or args.nci,
    }

    if not any(requested.values()):
        raise RuntimeError("No sources selected. Use --all or specify sources.")

    p53_cfg = configs.get("p53", {})
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    alphafold_cfg = p53_cfg.get("alphafold", {})
    clinvar_cfg = p53_cfg.get("clinvar", {})
    clinvar_enabled = bool(clinvar_cfg.get("enabled", True))
    nci_cfg = p53_cfg.get("nci_tp53", {})
    paths_cfg = configs.get("paths", {})
    paths = _paths_from_config(paths_cfg)

    logger.info("Fetching datasets for %s", uniprot_id)

    if requested["uniprot"]:
        logger.info("Fetching UniProt entry")
        _fetch_uniprot(uniprot_id, paths["raw"], paths["interim"])

    if requested["alphafold"]:
        logger.info("Fetching AlphaFold structure")
        _fetch_alphafold(uniprot_id, paths["raw"], paths["interim"], alphafold_cfg)

    if requested["clinvar"] and clinvar_enabled:
        logger.info("Fetching ClinVar release")
        _fetch_clinvar(paths["raw"], clinvar_cfg)
    elif requested["clinvar"] and not clinvar_enabled:
        logger.info("ClinVar disabled by config; skipping ClinVar fetch")

    if requested["nci"]:
        logger.info("Fetching NCI TP53 tables")
        _fetch_nci_tp53(paths["raw"], paths["interim"], nci_cfg)

    logger.info("Fetch complete")
    return 0
