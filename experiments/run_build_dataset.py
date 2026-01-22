from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.hashing import dict_sha256, file_sha256
from src.core.logging import get_logger
from src.variants.label_sets import export_label_sets, make_benign_set, make_pathogenic_set, quality_filter
from src.variants.map_to_protein import load_uniprot_sequence, map_and_filter, restrict_to_domain
from src.variants.normalize_variants import (
    drop_conflicted_labels,
    extract_transcript,
    filter_to_canonical_transcript,
    normalize_protein_change,
    dedupe_variants,
)
from src.variants.parse_clinvar_xml import parse_clinvar_xml_to_parquet


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    processed_dir = Path(data_cfg.get("processed", "data/processed"))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
        "processed": processed_dir,
    }


def _find_latest_file(paths: list[Path]) -> Path:
    return max(paths, key=lambda path: path.stat().st_mtime)


def _find_clinvar_xml(raw_dir: Path) -> Path | None:
    xml_candidates = list(raw_dir.glob("*.xml"))
    if xml_candidates:
        return _find_latest_file(xml_candidates)

    gz_candidates = list(raw_dir.glob("*.xml.gz"))
    if gz_candidates:
        return _find_latest_file(gz_candidates)

    parts_candidates = list(raw_dir.glob("*.parts.json"))
    if parts_candidates:
        return _find_latest_file(parts_candidates)

    logger = get_logger(__name__)
    logger.warning("No ClinVar XML files found, skipping ClinVar processing")
    return None


def _infer_transcripts(df: pd.DataFrame) -> pd.DataFrame:
    if "transcript" not in df.columns:
        df = df.copy()
        df["transcript"] = None

    def _infer(row) -> str | None:
        if row.get("transcript"):
            return row["transcript"]
        return extract_transcript([row.get("protein_hgvs"), row.get("coding_hgvs")])

    df = df.copy()
    df["transcript"] = df.apply(_infer, axis=1)
    return df


def _normalize_protein_changes(df: pd.DataFrame) -> pd.DataFrame:
    positions: list[int | None] = []
    refs: list[str | None] = []
    alts: list[str | None] = []

    for hgvs in df["protein_hgvs"].tolist():
        try:
            pos, ref, alt = normalize_protein_change(hgvs)
        except Exception:
            positions.append(None)
            refs.append(None)
            alts.append(None)
            continue
        positions.append(pos)
        refs.append(ref)
        alts.append(alt)

    df = df.copy()
    df["pos"] = positions
    df["ref"] = refs
    df["alt"] = alts
    df = df.dropna(subset=["pos", "ref", "alt"])
    df["pos"] = df["pos"].astype(int)
    return df


def run(args, configs: dict[str, Any]) -> int:
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))

    raw_dir = paths["raw"] / "clinvar"
    interim_dir = paths["interim"]
    processed_dir = paths["processed"]

    p53_cfg = configs.get("p53", {})
    clinvar_cfg = p53_cfg.get("clinvar", {})
    clinvar_enabled = bool(clinvar_cfg.get("enabled", True))
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    target_gene = clinvar_cfg.get("target_gene") or p53_cfg.get("gene_symbol")
    if target_gene:
        target_gene = str(target_gene).strip() or None
    canonical_transcript = p53_cfg.get("canonical_transcript")
    core_domain = p53_cfg.get("core_domain", {})
    domain_start = int(core_domain.get("start", 1))
    domain_end = int(core_domain.get("end", 9999))
    label_cfg = p53_cfg.get("labels", {})
    min_review_stars = int(label_cfg.get("min_review_stars", 0))

    if not clinvar_enabled:
        logger.info("ClinVar disabled by config; skipping dataset build")
        return 0

    xml_path = _find_clinvar_xml(raw_dir)
    if not xml_path:
        logger.warning("No ClinVar data available, skipping variant processing")
        return 0
        
    parsed_out = interim_dir / "clinvar.parquet"
    parser_mode = clinvar_cfg.get("parser", "auto")
    chunk_size = int(clinvar_cfg.get("chunk_size", 100000))
    num_workers = clinvar_cfg.get("num_workers")
    element_batch_size = int(clinvar_cfg.get("element_batch_size", 200))
    reparse = bool(clinvar_cfg.get("reparse", False))

    use_existing = False
    if parsed_out.exists() and not reparse:
        try:
            preview = pd.read_parquet(parsed_out, columns=["protein_hgvs"])
        except Exception:
            preview = None
        if preview is not None and preview["protein_hgvs"].notna().any():
            use_existing = True
        else:
            logger.warning("Existing ClinVar parse lacks protein HGVS; reparsing.")
            reparse = True

    if use_existing:
        logger.info("Using existing parsed variants at %s", parsed_out)
    else:
        logger.info("Parsing ClinVar XML: %s", xml_path)
        if parser_mode == "fast" or (parser_mode == "auto" and xml_path.name.endswith(".parts.json")):
            from src.variants.parse_clinvar_xml import parse_clinvar_xml_to_parquet_fast

            parse_clinvar_xml_to_parquet_fast(
                xml_path,
                parsed_out,
                chunk_size=chunk_size,
                element_batch_size=element_batch_size,
                num_workers=num_workers,
                target_gene=target_gene,
            )
        else:
            parse_clinvar_xml_to_parquet(
                xml_path,
                parsed_out,
                chunk_size=chunk_size,
                target_gene=target_gene,
            )

    fasta_path = paths["interim"] / "uniprot" / f"{uniprot_id}.fasta"
    normalized_out = interim_dir / "tp53_missense_normalized.parquet"
    signature_path = interim_dir / "tp53_missense_signature.json"
    labels_dir = processed_dir / "labels"

    signature = {
        "parsed_sha256": file_sha256(parsed_out),
        "uniprot_sha256": file_sha256(fasta_path),
        "canonical_transcript": canonical_transcript,
        "domain_start": domain_start,
        "domain_end": domain_end,
        "min_review_stars": min_review_stars,
        "target_gene": target_gene,
    }
    signature_hash = dict_sha256(signature)
    labels_ok = (labels_dir / "benign.parquet").exists() and (labels_dir / "pathogenic.parquet").exists()
    if normalized_out.exists() and signature_path.exists() and labels_ok:
        existing = json.loads(signature_path.read_text(encoding="utf-8"))
        if existing.get("signature_sha256") == signature_hash:
            logger.info("Normalized variants up to date; skipping rebuild")
            return 0

    logger.info("Loading parsed variants")
    df = pd.read_parquet(parsed_out)
    logger.info("Parsed %d raw TP53 missense rows", len(df))

    df = _infer_transcripts(df)
    df = _normalize_protein_changes(df)
    df = drop_conflicted_labels(df)
    df = filter_to_canonical_transcript(df, canonical_transcript)
    df = dedupe_variants(df)

    seq = load_uniprot_sequence(fasta_path)
    df = map_and_filter(df, seq)
    df = restrict_to_domain(df, domain_start, domain_end)

    normalized_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(normalized_out, index=False)
    logger.info("Saved normalized variants to %s", normalized_out)

    benign_df = make_benign_set(df)
    path_df = make_pathogenic_set(df)

    if min_review_stars > 0:
        benign_df = quality_filter(benign_df, min_review_stars=min_review_stars)
        path_df = quality_filter(path_df, min_review_stars=min_review_stars)

    export_label_sets(benign_df, path_df, processed_dir / "labels")
    logger.info("Exported benign/pathogenic label sets")

    # Generate train/test split
    from src.eval.train_test_split import create_split_indices, save_split

    validation_cfg = cfg.get("validation", {})
    if validation_cfg.get("enabled", False):
        logger.info("Generating train/test split...")

        split_ratio = validation_cfg.get("split_ratio", 0.8)
        split_seed = validation_cfg.get("split_seed", 42)
        stratify_by = validation_cfg.get("stratify_by", None)

        benign_path = processed_dir / "labels" / "benign.parquet"
        pathogenic_path = processed_dir / "labels" / "pathogenic.parquet"

        split_info = create_split_indices(
            benign_path,
            pathogenic_path,
            split_ratio=split_ratio,
            split_seed=split_seed,
            stratify_by=stratify_by,
        )

        split_path = processed_dir / "splits" / f"clinvar_split_seed{split_seed}.json"
        save_split(split_info, split_path)

        logger.info(f"Train/test split saved to {split_path}")
    else:
        logger.info("Validation splits disabled in config")

    signature_path.write_text(
        json.dumps({"signature_sha256": signature_hash, "inputs": signature}, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    return 0
