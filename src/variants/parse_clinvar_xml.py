from __future__ import annotations

import gzip
import io
import json
import multiprocessing as mp
import re
import xml.etree.ElementTree as ET
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, as_completed, wait
from pathlib import Path
from typing import Any, Callable, Iterator

import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

HGVS_PROTEIN_RE = re.compile(r"p\.[A-Za-z]{1,3}\d+[A-Za-z*]{1,3}(?=$|\s|\)|;|,)")
HGVS_CODING_RE = re.compile(r"c\.[^\s\)]+")
TRANSCRIPT_RE = re.compile(r"([A-Z]{2,4}_\d+\.\d+)")
GENE_TEXT_TAGS = {"Name", "PreferredName", "Expression", "Attribute", "ElementValue"}
_GENE_REGEX_CACHE: dict[str, re.Pattern[str]] = {}


def _strip_namespace(tag: str) -> str:
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def _unique(items: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def _gene_regex(target_gene: str) -> re.Pattern[str]:
    key = target_gene.upper()
    regex = _GENE_REGEX_CACHE.get(key)
    if regex is None:
        regex = re.compile(rf"(?<![A-Z0-9]){re.escape(key)}(?![A-Z0-9])", re.IGNORECASE)
        _GENE_REGEX_CACHE[key] = regex
    return regex


def _matches_gene_symbol(element: ET.Element, target_gene: str | None) -> bool:
    if not target_gene:
        return True

    gene_upper = target_gene.upper()
    token_re = _gene_regex(gene_upper)

    for node in element.iter():
        tag = _strip_namespace(node.tag)
        if tag == "GeneSymbol" and node.text:
            if node.text.strip().upper() == gene_upper:
                return True
        elif tag == "Gene":
            symbol = node.get("Symbol")
            if symbol and symbol.strip().upper() == gene_upper:
                return True
        elif tag in GENE_TEXT_TAGS and node.text:
            if token_re.search(node.text):
                return True

    return False


class _PartFileReader(io.RawIOBase):
    def __init__(self, paths: list[Path]) -> None:
        self._paths = list(paths)
        self._index = 0
        self._handle = None
        self._open_next()

    def _open_next(self) -> None:
        if self._handle is not None:
            self._handle.close()
        if self._index >= len(self._paths):
            self._handle = None
            return
        self._handle = self._paths[self._index].open("rb")
        self._index += 1

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> bytes:
        if self._handle is None:
            return b""
        if size == 0:
            return b""
        if size < 0:
            chunks: list[bytes] = []
            while self._handle is not None:
                data = self._handle.read()
                if data:
                    chunks.append(data)
                self._open_next()
            return b"".join(chunks)

        remaining = size
        chunks: list[bytes] = []
        while remaining > 0 and self._handle is not None:
            data = self._handle.read(remaining)
            if not data:
                self._open_next()
                continue
            chunks.append(data)
            remaining -= len(data)
        return b"".join(chunks)

    def readinto(self, b: bytearray) -> int:
        data = self.read(len(b))
        size = len(data)
        b[:size] = data
        return size

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        super().close()


def _open_xml_stream(path: Path) -> io.BufferedReader:
    if path.name.endswith(".parts.json"):
        manifest = json.loads(path.read_text(encoding="utf-8"))
        parts = [path.parent / part["path"] for part in manifest["parts"]]
        reader = _PartFileReader(parts)
        return io.BufferedReader(gzip.GzipFile(fileobj=reader))

    if path.suffix == ".gz":
        return io.BufferedReader(gzip.open(path, "rb"))

    return io.BufferedReader(path.open("rb"))


def _extract_transcripts(text: str) -> list[str]:
    return TRANSCRIPT_RE.findall(text)


def extract_hgvs(element: ET.Element) -> dict[str, list[str]]:
    protein_hgvs: list[str] = []
    coding_hgvs: list[str] = []
    transcripts: list[str] = []

    for node in element.iter():
        text = (node.text or "").strip()
        if not text:
            continue
        if "p." not in text and "c." not in text and "_" not in text:
            continue

        transcripts.extend(_extract_transcripts(text))
        protein_hgvs.extend(HGVS_PROTEIN_RE.findall(text))
        coding_hgvs.extend(HGVS_CODING_RE.findall(text))

    return {
        "protein_hgvs": _unique(protein_hgvs),
        "coding_hgvs": _unique(coding_hgvs),
        "transcripts": _unique(transcripts),
    }


def extract_significance(element: ET.Element) -> dict[str, str | None]:
    descriptions: list[str] = []
    review_status: list[str] = []
    last_evaluated: str | None = None

    for node in element.iter():
        if _strip_namespace(node.tag) != "ClinicalSignificance":
            continue
        if node.get("DateLastEvaluated"):
            last_evaluated = node.get("DateLastEvaluated")
        for child in node:
            child_tag = _strip_namespace(child.tag)
            if child_tag == "Description" and child.text:
                descriptions.append(child.text.strip())
            elif child_tag == "ReviewStatus" and child.text:
                review_status.append(child.text.strip())

    return {
        "clinical_significance": descriptions[0] if descriptions else None,
        "review_status": review_status[0] if review_status else None,
        "last_evaluated": last_evaluated,
    }


def _extract_gene_symbols(element: ET.Element) -> list[str]:
    symbols: list[str] = []
    for node in element.iter():
        tag = _strip_namespace(node.tag)
        if tag == "GeneSymbol" and node.text:
            symbols.append(node.text.strip())
        if tag == "Gene" and node.get("Symbol"):
            symbols.append(node.get("Symbol").strip())
    return _unique(symbols)


def _extract_molecular_consequences(element: ET.Element) -> list[str]:
    consequences: list[str] = []
    for node in element.iter():
        if _strip_namespace(node.tag) != "MolecularConsequence":
            continue
        if node.get("Type"):
            consequences.append(node.get("Type").strip())
        elif node.text:
            consequences.append(node.text.strip())
    return _unique([item.lower() for item in consequences if item])


def _extract_accession(element: ET.Element) -> str | None:
    for node in element.iter():
        if _strip_namespace(node.tag) != "ClinVarAccession":
            continue
        accession = node.get("Accession")
        if accession:
            return accession
    return None


def _extract_variation_id(element: ET.Element) -> str | None:
    variation_id = element.get("VariationID") or element.get("ID")
    if variation_id:
        return variation_id

    for node in element.iter():
        if _strip_namespace(node.tag) != "MeasureSet":
            continue
        if node.get("ID"):
            return node.get("ID")
    return None


def _build_record(element: ET.Element) -> dict[str, Any]:
    hgvs = extract_hgvs(element)
    sig = extract_significance(element)
    return {
        "variation_id": _extract_variation_id(element),
        "accession": _extract_accession(element),
        "gene_symbols": _extract_gene_symbols(element),
        "protein_hgvs": hgvs["protein_hgvs"],
        "coding_hgvs": hgvs["coding_hgvs"],
        "transcripts": hgvs["transcripts"],
        "molecular_consequences": _extract_molecular_consequences(element),
        **sig,
    }


def is_tp53(record: dict[str, Any]) -> bool:
    symbols = [symbol.upper() for symbol in record.get("gene_symbols", [])]
    if "TP53" in symbols:
        return True

    for entry in record.get("protein_hgvs", []) + record.get("coding_hgvs", []):
        if "TP53" in (entry or "").upper():
            return True
    return False


def is_missense(record: dict[str, Any]) -> bool:
    if any("missense" in consequence for consequence in record.get("molecular_consequences", [])):
        return True

    for hgvs in record.get("protein_hgvs", []):
        if HGVS_PROTEIN_RE.search(hgvs or ""):
            return True
    return False


def _choose_transcript(record: dict[str, Any]) -> str | None:
    transcripts = record.get("transcripts", [])
    if not transcripts:
        return None

    for prefix in ("NM_", "NR_"):
        for transcript in transcripts:
            if transcript.startswith(prefix):
                return transcript
    for prefix in ("NP_", "XP_", "XR_"):
        for transcript in transcripts:
            if transcript.startswith(prefix):
                return transcript
    return transcripts[0]


def to_flat_row(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "variation_id": record.get("variation_id"),
        "accession": record.get("accession"),
        "gene_symbol": (record.get("gene_symbols") or [None])[0],
        "protein_hgvs": record.get("protein_hgvs"),
        "coding_hgvs": record.get("coding_hgvs"),
        "transcript": record.get("transcript"),
        "molecular_consequence": (record.get("molecular_consequences") or [None])[0],
        "clinical_significance": record.get("clinical_significance"),
        "review_status": record.get("review_status"),
        "last_evaluated": record.get("last_evaluated"),
    }


def _explode_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    protein_hgvs = record.get("protein_hgvs", [])
    coding_hgvs = record.get("coding_hgvs", [])

    transcript = _choose_transcript(record)
    if not protein_hgvs:
        record_copy = {**record, "protein_hgvs": None, "coding_hgvs": coding_hgvs[0] if coding_hgvs else None}
        record_copy["transcript"] = transcript
        return [to_flat_row(record_copy)]

    rows: list[dict[str, Any]] = []
    for hgvs in protein_hgvs:
        record_copy = {**record, "protein_hgvs": hgvs}
        record_copy["coding_hgvs"] = coding_hgvs[0] if coding_hgvs else None
        record_copy["transcript"] = transcript
        rows.append(to_flat_row(record_copy))
    return rows


def iter_clinvar_variants(
    xml_path: Path,
    target_gene: str | None = None,
    on_record: Callable[[], None] | None = None,
) -> Iterator[dict[str, Any]]:
    with _open_xml_stream(xml_path) as handle:
        context = ET.iterparse(handle, events=("end",))
        for _, element in context:
            tag = _strip_namespace(element.tag)
            if tag not in {"ClinVarSet", "VariationArchive"}:
                continue

            if on_record:
                on_record()

            if target_gene and not _matches_gene_symbol(element, target_gene):
                element.clear()
                continue

            record = _build_record(element)
            if not is_missense(record):
                element.clear()
                continue
            if not target_gene and not is_tp53(record):
                element.clear()
                continue

            for row in _explode_record(record):
                yield row

            element.clear()


def _write_parquet(rows: list[dict[str, Any]], writer: pq.ParquetWriter | None, out_path: Path) -> pq.ParquetWriter:
    table = pa.Table.from_pylist(rows)
    if writer is None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(out_path, table.schema)
    writer.write_table(table)
    return writer


def _process_xml_batch(batch: list[bytes], target_gene: str | None = None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for xml_bytes in batch:
        element = ET.fromstring(xml_bytes)
        if target_gene and not _matches_gene_symbol(element, target_gene):
            continue
        record = _build_record(element)
        if not is_missense(record):
            continue
        if not target_gene and not is_tp53(record):
            continue
        rows.extend(_explode_record(record))
    return rows


def _consume_futures(
    futures: set,
    buffer: list[dict[str, Any]],
    writer: pq.ParquetWriter | None,
    out_path: Path,
    chunk_size: int,
) -> tuple[list[dict[str, Any]], pq.ParquetWriter | None]:
    for future in futures:
        rows = future.result()
        if rows:
            buffer.extend(rows)
        if len(buffer) >= chunk_size:
            writer = _write_parquet(buffer, writer, out_path)
            buffer.clear()
    return buffer, writer


def _parse_xml_with_multiprocessing(
    xml_path: Path,
    out_path: Path,
    chunk_size: int = 20000,
    element_batch_size: int = 200,
    num_workers: int | None = None,
    max_pending: int | None = None,
    target_gene: str | None = None,
) -> None:
    if num_workers is None:
        num_workers = min(mp.cpu_count(), 8)
    if num_workers < 1:
        num_workers = 1
    if max_pending is None:
        max_pending = max(2, num_workers * 2)

    if out_path.exists():
        out_path.unlink()

    buffer: list[dict[str, Any]] = []
    writer: pq.ParquetWriter | None = None
    pending: set = set()
    batch: list[bytes] = []
    record_counter = 0
    progress_every = max(1000, element_batch_size)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        with _open_xml_stream(xml_path) as handle:
            context = ET.iterparse(handle, events=("end",))
            with tqdm(desc="Processing ClinVar records", unit="records") as pbar:
                for _, element in context:
                    tag = _strip_namespace(element.tag)
                    if tag not in {"ClinVarSet", "VariationArchive"}:
                        continue

                    record_counter += 1
                    if record_counter >= progress_every:
                        pbar.update(record_counter)
                        record_counter = 0

                    if target_gene and not _matches_gene_symbol(element, target_gene):
                        element.clear()
                        continue

                    batch.append(ET.tostring(element, encoding="utf-8"))
                    element.clear()

                    if len(batch) >= element_batch_size:
                        future = executor.submit(_process_xml_batch, batch, target_gene)
                        pending.add(future)
                        batch = []

                        if len(pending) >= max_pending:
                            done, pending = wait(pending, return_when=FIRST_COMPLETED)
                            buffer, writer = _consume_futures(done, buffer, writer, out_path, chunk_size)

                if batch:
                    future = executor.submit(_process_xml_batch, batch, target_gene)
                    pending.add(future)

                if record_counter:
                    pbar.update(record_counter)

                for future in as_completed(pending):
                    buffer, writer = _consume_futures({future}, buffer, writer, out_path, chunk_size)

    if buffer:
        writer = _write_parquet(buffer, writer, out_path)
    if writer is not None:
        writer.close()


def parse_clinvar_xml_to_parquet(
    xml_path: Path,
    out_path: Path,
    chunk_size: int = 20000,
    progress_every: int = 1000,
    target_gene: str | None = None,
) -> None:
    if out_path.exists():
        out_path.unlink()

    buffer: list[dict[str, Any]] = []
    writer: pq.ParquetWriter | None = None
    record_counter = 0

    def _tick() -> None:
        nonlocal record_counter
        record_counter += 1

    with tqdm(desc="Processing ClinVar records", unit="records") as pbar:
        for row in iter_clinvar_variants(xml_path, target_gene=target_gene, on_record=_tick):
            if record_counter >= progress_every:
                pbar.update(record_counter)
                record_counter = 0
            buffer.append(row)
            if len(buffer) >= chunk_size:
                writer = _write_parquet(buffer, writer, out_path)
                buffer.clear()

        if record_counter:
            pbar.update(record_counter)

    if buffer:
        writer = _write_parquet(buffer, writer, out_path)

    if writer is not None:
        writer.close()


def parse_clinvar_xml_to_parquet_fast(
    xml_path: Path,
    out_path: Path,
    chunk_size: int = 50000,
    element_batch_size: int = 200,
    num_workers: int | None = None,
    target_gene: str | None = None,
) -> None:
    _parse_xml_with_multiprocessing(
        xml_path,
        out_path,
        chunk_size=chunk_size,
        element_batch_size=element_batch_size,
        num_workers=num_workers,
        target_gene=target_gene,
    )
