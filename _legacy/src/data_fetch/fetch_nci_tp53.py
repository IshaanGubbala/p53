from __future__ import annotations

import csv
import datetime as dt
import json
import re
import tarfile
import zipfile
from pathlib import Path
from typing import Iterable
from urllib.error import URLError
from urllib.request import Request, urlopen

from src.core.hashing import file_sha256

USER_AGENT = "p53-stabilimut/0.1"
DCEG_PAGE = "https://dceg.cancer.gov/tools/public-data/tp53-database"
DEFAULT_URLS = [
    "https://tp53.cancer.gov/view_data?bq_view_name=TumorVariantDownload",
    "https://tp53.cancer.gov/view_data?bq_view_name=CellLineDownload",
]


def _download(url: str, out_path: Path, chunk_size: int = 1024 * 1024) -> None:
    if out_path.exists():
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".partial")

    headers = {
        "User-Agent": USER_AGENT,
        "Accept": "text/csv,application/csv,text/plain,*/*",
        "Accept-Language": "en-US,en;q=0.9",
    }
    request = Request(url, headers=headers)
    try:
        with urlopen(request) as response, tmp_path.open("wb") as handle:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                handle.write(chunk)
    except URLError as exc:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    tmp_path.replace(out_path)


def _read_head(path: Path, size: int = 2048) -> bytes:
    with path.open("rb") as handle:
        return handle.read(size)


def _looks_like_html(path: Path) -> bool:
    head = _read_head(path).lower()
    return b"<html" in head or b"<!doctype" in head or b"<body" in head


def _extract_data_from_html(html_path: Path, out_dir: Path) -> list[Path] | None:
    """Try to extract CSV data from an HTML response."""
    try:
        import re
        import ast
        from bs4 import BeautifulSoup
        
        with html_path.open("r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        
        soup = BeautifulSoup(content, "html.parser")
        
        # Extract table headers first
        headers = []
        header_ths = soup.find_all("th")
        for th in header_ths:
            header_text = th.get_text(strip=True)
            if header_text and header_text not in headers:
                headers.append(header_text)
        
        # Look for data in script tags
        scripts = soup.find_all("script")
        for script in scripts:
            if script.string:
                script_content = script.string
                
                # Look for CSV data patterns
                csv_match = re.search(r'data:text/csv;base64,([A-Za-z0-9+/=]+)', script_content)
                if csv_match:
                    import base64
                    csv_data = base64.b64decode(csv_match.group(1))
                    csv_path = out_dir / "extracted_data.csv"
                    with csv_path.open("wb") as f:
                        f.write(csv_data)
                    return [csv_path]
                
                # Look for JavaScript data patterns with multiple regex attempts
                patterns = [
                    r'let\s+dataSet\s*=\s*assign_var\((\[[\s\S]*?\]\));',
                    r'dataSet\s*=\s*assign_var\((\[[\s\S]*?\]\));',
                    r'var\s+dataSet\s*=\s*assign_var\((\[[\s\S]*?\]\));'
                ]
                
                for pattern in patterns:
                    js_match = re.search(pattern, script_content)
                    if js_match:
                        try:
                            # Extract the JavaScript array and convert to Python
                            js_array_str = js_match.group(1)
                            # Replace JavaScript escape sequences and clean up
                            js_array_str = js_array_str.replace('\\u003e', '>').replace('\\u003c', '<').replace('\\"', '"')
                            
                            # Try to parse with a more robust method
                            try:
                                data = ast.literal_eval(js_array_str)
                            except SyntaxError:
                                # Fall back to json loads if ast fails
                                import json
                                data = json.loads(js_array_str.replace("'", '"'))
                            
                            # Convert to CSV with headers
                            import csv
                            csv_path = out_dir / "extracted_data.csv"
                            with csv_path.open("w", newline="", encoding="utf-8") as f:
                                writer = csv.writer(f)
                                
                                # Write headers if we have them
                                if headers:
                                    writer.writerow(headers)
                                
                                # Write data rows
                                if isinstance(data, list):
                                    for row in data:
                                        if isinstance(row, list):
                                            writer.writerow(row)
                                
                                if len(headers) > 0 and len(data) > 0:
                                    return [csv_path]
                                    
                        except (ValueError, SyntaxError, KeyError) as e:
                            print(f"JavaScript parsing error: {e}")
                            continue
        
        # Look for table data as fallback
        tables = soup.find_all("table")
        if tables:
            import csv
            table_path = out_dir / "table_data.csv"
            with table_path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                for table in tables:
                    rows = table.find_all("tr")
                    for row in rows:
                        cells = [cell.get_text(strip=True) for cell in row.find_all(["td", "th"])]
                        if cells and any(cell.strip() for cell in cells):  # Skip empty rows
                            writer.writerow(cells)
            
            # Check if we got more than just headers
            with table_path.open("r", encoding="utf-8") as f:
                lines = f.readlines()
                if len(lines) > 2:  # More than just header + one empty row
                    return [table_path]
                            
    except ImportError:
        # BeautifulSoup not available, skip HTML parsing
        pass
    except Exception as e:
        # Any other error, return None
        print(f"HTML parsing error: {e}")
        pass
    
    return None


def _is_gzip(path: Path) -> bool:
    head = _read_head(path, size=2)
    return head == b"\x1f\x8b"


def _classify_download(path: Path) -> str:
    if zipfile.is_zipfile(path):
        return "zip"
    if tarfile.is_tarfile(path):
        return "tar"
    if _is_gzip(path):
        return "gz"
    if _looks_like_html(path):
        return "html"
    return "unknown"


def _discover_urls() -> list[str]:
    request = Request(DCEG_PAGE, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            html = response.read().decode("utf-8", errors="ignore")
    except URLError:
        return []

    links = re.findall(r'href="([^"]+)"', html)
    candidates: list[str] = []
    for link in links:
        link_lower = link.lower()
        if not link_lower.endswith((".zip", ".gz", ".tar.gz", ".tgz", ".txt", ".tsv", ".csv", ".tab")):
            continue
        if "tp53" not in link_lower and "p53" not in link_lower:
            continue
        if link.startswith("http"):
            candidates.append(link)
        elif link.startswith("//"):
            candidates.append(f"https:{link}")
        else:
            candidates.append(f"https://dceg.cancer.gov{link}")
    return sorted(set(candidates))


def _extract_archive(archive_path: Path, out_dir: Path) -> list[Path]:
    extract_dir = out_dir / archive_path.stem
    if extract_dir.exists():
        return [p for p in extract_dir.rglob("*") if p.is_file()]

    extract_dir.mkdir(parents=True, exist_ok=True)
    if zipfile.is_zipfile(archive_path):
        with zipfile.ZipFile(archive_path) as zf:
            zf.extractall(extract_dir)
    elif tarfile.is_tarfile(archive_path):
        with tarfile.open(archive_path, "r:*") as tar:
            tar.extractall(extract_dir)
    else:
        raise RuntimeError(f"Unsupported archive format: {archive_path}")

    return [p for p in extract_dir.rglob("*") if p.is_file()]


def _normalize_table_paths(paths: Iterable[Path]) -> list[Path]:
    table_suffixes = {".tsv", ".txt", ".csv", ".tab"}
    tables = [path for path in paths if path.suffix.lower() in table_suffixes]
    return sorted(tables)


def _write_manifest(
    out_dir: Path,
    source_url: str,
    download_path: Path,
    table_paths: list[Path],
) -> None:
    manifest = {
        "generated_at": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "source_url": source_url,
        "archive": {
            "path": str(download_path.relative_to(out_dir)),
            "size": download_path.stat().st_size,
            "sha256": file_sha256(download_path),
        },
        "tables": [
            {
                "path": str(path.relative_to(out_dir)),
                "size": path.stat().st_size,
                "sha256": file_sha256(path),
            }
            for path in table_paths
        ],
    }

    out_path = out_dir / "nci_tp53_manifest.json"
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)


def download_nci_tp53_tables(out_dir: Path, urls: list[str] | None = None) -> list[Path]:
    """Download the NCI TP53 dataset tables and return extracted table paths."""
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate_urls = urls or DEFAULT_URLS
    last_error: Exception | None = None
    tables: list[Path] = []
    
    for url in candidate_urls:
        # Generate filename based on view name
        if "TumorVariantDownload" in url:
            filename = "tumor_variants.csv"
        elif "CellLineDownload" in url:
            filename = "cell_lines.csv"
        else:
            filename = url.split("/")[-1]
            if not filename.endswith((".csv", ".tsv", ".txt", ".tab")):
                filename += ".csv"
        
        download_path = out_dir / filename
        try:
            _download(url, download_path)
        except Exception as exc:
            last_error = exc
            continue

        if download_path.stat().st_size == 0:
            download_path.unlink(missing_ok=True)
            last_error = RuntimeError(f"Downloaded empty file for {url}")
            continue
            
        # Save the file regardless - HTML, CSV, whatever it is
        # You can process it later
        tables.append(download_path)
        _write_manifest(out_dir, url, download_path, [download_path])
        break  # Get the first successful download

    if not tables:
        raise RuntimeError(f"Failed to download NCI TP53 data: {last_error}")
    
    return tables


def _sniff_dialect(handle):
    sample = handle.read(4096)
    handle.seek(0)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;|")
    except csv.Error:
        class _Dialect(csv.Dialect):
            delimiter = "\t"
            quotechar = '"'
            doublequote = True
            escapechar = None
            lineterminator = "\n"
            quoting = csv.QUOTE_MINIMAL
            skipinitialspace = False

        return _Dialect()


def _parse_position(value: str | None) -> int | None:
    if value is None:
        return None
    match = re.search(r"(\d+)", value)
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def _preferred_keys(fieldnames) -> list[str]:
    priority = [
        "codon",
        "codon_number",
        "residue",
        "position",
        "pos",
        "aa_position",
        "protein_position",
        "amino_acid_position",
        "protein",
        "mutation",
    ]

    normalized = {name: name.strip().lower() for name in fieldnames}
    ordered: list[str] = []
    for key in priority:
        for original, norm in normalized.items():
            if norm == key:
                ordered.append(original)
    return ordered


def parse_hotspots(table_path: Path) -> list[dict]:
    """Parse hotspot residues from an NCI TP53 table."""
    hotspots: list[dict] = []
    with table_path.open("r", encoding="utf-8", errors="replace") as handle:
        dialect = _sniff_dialect(handle)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None:
            return hotspots

        keys = _preferred_keys(reader.fieldnames)
        fallback_keys = [
            name
            for name in reader.fieldnames
            if any(token in name.lower() for token in ["codon", "pos", "position", "mut", "protein", "aa"])
        ]

        for row in reader:
            pos = None
            for key in keys + fallback_keys:
                pos = _parse_position(row.get(key))
                if pos is not None:
                    break
            if pos is None:
                continue

            hotspots.append(
                {
                    "pos": pos,
                    "source": table_path.name,
                    "raw": row,
                }
            )

    return hotspots


def save_hotspots(hotspots: list[dict], out_path: Path) -> None:
    """Save hotspot list to JSON."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "count": len(hotspots),
        "hotspots": hotspots,
    }
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
