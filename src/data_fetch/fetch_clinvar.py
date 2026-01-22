from __future__ import annotations

import datetime as dt
import gzip
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from src.core.hashing import file_sha256

USER_AGENT = "p53-stabilimut/0.1"
CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
CLINVAR_XML_NAME = "ClinVarFullRelease_00-latest.xml.gz"


def _download(url: str, out_path: Path, chunk_size: int = 1024 * 1024) -> None:
    if out_path.exists():
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".partial")

    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response, tmp_path.open("wb") as handle:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                handle.write(chunk)
    except HTTPError as exc:
        tmp_path.unlink(missing_ok=True)
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    except URLError as exc:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    tmp_path.replace(out_path)


def _head_request(url: str) -> dict[str, str]:
    request = Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            return dict(response.headers)
    except HTTPError as exc:
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
        raise RuntimeError(f"Failed to fetch headers {url}: {exc}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to fetch headers {url}: {exc}") from exc


def _get_content_length(url: str) -> tuple[int | None, str | None]:
    try:
        headers = _head_request(url)
        length = headers.get("Content-Length")
        accept_ranges = headers.get("Accept-Ranges")
        if length and length.isdigit():
            return int(length), accept_ranges
    except FileNotFoundError:
        raise
    except Exception:
        pass

    request = Request(url, headers={"User-Agent": USER_AGENT, "Range": "bytes=0-0"})
    try:
        with urlopen(request) as response:
            content_range = response.headers.get("Content-Range")
            accept_ranges = response.headers.get("Accept-Ranges")
            if content_range and "/" in content_range:
                total = content_range.split("/")[-1]
                if total.isdigit():
                    return int(total), accept_ranges
            length = response.headers.get("Content-Length")
            if length and length.isdigit():
                return int(length), accept_ranges
    except HTTPError as exc:
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to fetch range headers {url}: {exc}") from exc

    return None, None


def _download_parts(
    url: str,
    out_dir: Path,
    base_name: str,
    part_size_bytes: int,
) -> Path:
    total_size, _ = _get_content_length(url)
    if total_size is None:
        raise RuntimeError(f"Unable to determine size for {url}")

    out_dir.mkdir(parents=True, exist_ok=True)
    parts: list[dict[str, Any]] = []
    total_parts = int(math.ceil(total_size / part_size_bytes))

    for idx in range(total_parts):
        start = idx * part_size_bytes
        end = min(total_size - 1, start + part_size_bytes - 1)
        expected_size = end - start + 1
        part_name = f"{base_name}.part{idx:03d}"
        part_path = out_dir / part_name

        if part_path.exists() and part_path.stat().st_size == expected_size:
            parts.append(
                {
                    "path": part_name,
                    "start": start,
                    "end": end,
                    "size": expected_size,
                }
            )
            continue

        tmp_path = part_path.with_suffix(part_path.suffix + ".partial")
        request = Request(
            url,
            headers={"User-Agent": USER_AGENT, "Range": f"bytes={start}-{end}"},
        )
        try:
            with urlopen(request) as response, tmp_path.open("wb") as handle:
                if response.status == 200 and total_size > part_size_bytes:
                    raise RuntimeError("Server does not honor range requests")
                if response.status == 200 and start != 0:
                    raise RuntimeError("Server does not honor range requests")
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

        if tmp_path.stat().st_size != expected_size:
            tmp_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"Incomplete download for {part_name} "
                f"(expected {expected_size} bytes)"
            )

        tmp_path.replace(part_path)
        parts.append(
            {
                "path": part_name,
                "start": start,
                "end": end,
                "size": expected_size,
            }
        )

    manifest = {
        "source_url": url,
        "total_size": total_size,
        "part_size": part_size_bytes,
        "parts": parts,
    }
    manifest_path = out_dir / f"{base_name}.parts.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    return manifest_path

def _fetch_listing(url: str) -> str:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request) as response:
            return response.read().decode("utf-8", errors="replace")
    except HTTPError as exc:
        if exc.code == 404:
            raise FileNotFoundError(f"Not found: {url}") from exc
        raise RuntimeError(f"Failed to fetch listing {url}: {exc}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to fetch listing {url}: {exc}") from exc


def _extract_filenames(listing: str, suffix: str) -> list[str]:
    names: set[str] = set()
    for match in re.findall(r'href="([^"]+)"', listing, flags=re.IGNORECASE):
        name = Path(match).name
        if name.endswith(suffix):
            names.add(name)

    for match in re.findall(rf"([A-Za-z0-9._-]+{re.escape(suffix)})", listing):
        names.add(match)

    return sorted(names)


def _latest_name(candidates: list[str], preferred: str | None = None) -> str:
    if preferred and preferred in candidates:
        return preferred

    def _sort_key(name: str) -> tuple[int, tuple[int, ...] | tuple[str, ...]]:
        lower = name.lower()
        if "latest" in lower:
            return (2, tuple())
        digits = tuple(int(value) for value in re.findall(r"\d+", name))
        if digits:
            return (1, digits)
        return (0, (name,))

    return max(candidates, key=_sort_key)


def _discover_latest_xml_name() -> str:
    listing = _fetch_listing(f"{CLINVAR_BASE}/xml/")
    names = _extract_filenames(listing, ".xml.gz")
    if not names:
        raise RuntimeError("No ClinVar XML files found in listing")
    return _latest_name(names, preferred=CLINVAR_XML_NAME)


def _discover_latest_vcf_name(assembly: str) -> str:
    listing = _fetch_listing(f"{CLINVAR_BASE}/vcf_{assembly}/")
    names = _extract_filenames(listing, ".vcf.gz")
    if not names:
        raise RuntimeError(f"No ClinVar VCF files found for {assembly}")
    return _latest_name(names, preferred="clinvar.vcf.gz")


def _decompress_gzip(gz_path: Path, out_path: Path) -> None:
    if out_path.exists():
        return

    try:
        with gzip.open(gz_path, "rb") as src, out_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    except OSError:
        out_path.unlink(missing_ok=True)
        raise


def download_clinvar_xml_release(
    out_dir: Path,
    decompress: bool = True,
    allow_parts: bool = True,
    part_size_mb: int = 1024,
) -> Path:
    """Download the ClinVar XML release and return the XML (or GZ) path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    xml_name = CLINVAR_XML_NAME
    url = f"{CLINVAR_BASE}/xml/{xml_name}"
    gz_path = out_dir / xml_name
    xml_path = out_dir / xml_name.replace(".gz", "")

    try:
        _download(url, gz_path)
    except FileNotFoundError:
        xml_name = _discover_latest_xml_name()
        url = f"{CLINVAR_BASE}/xml/{xml_name}"
        gz_path = out_dir / xml_name
        xml_path = out_dir / xml_name.replace(".gz", "")
        try:
            _download(url, gz_path)
        except OSError as exc:
            if exc.errno == 27 and allow_parts:
                manifest_path = _download_parts(
                    url,
                    out_dir,
                    xml_name,
                    part_size_mb * 1024 * 1024,
                )
                return manifest_path
            raise
    except OSError as exc:
        if exc.errno == 27 and allow_parts:
            manifest_path = _download_parts(
                url,
                out_dir,
                xml_name,
                part_size_mb * 1024 * 1024,
            )
            return manifest_path
        raise

    if decompress:
        _decompress_gzip(gz_path, xml_path)
        return xml_path
    return gz_path


def download_clinvar_vcf(
    out_dir: Path,
    assembly: str = "GRCh38",
    allow_parts: bool = True,
    part_size_mb: int = 1024,
) -> Path:
    """Download the ClinVar VCF for the requested assembly."""
    assembly_upper = assembly.upper()
    if assembly_upper not in {"GRCH37", "GRCH38"}:
        raise ValueError("assembly must be GRCh37 or GRCh38")
    assembly = "GRCh37" if assembly_upper == "GRCH37" else "GRCh38"

    vcf_dir = out_dir / f"vcf_{assembly}"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    vcf_name = "clinvar.vcf.gz"
    vcf_url = f"{CLINVAR_BASE}/vcf_{assembly}/{vcf_name}"
    tbi_url = f"{CLINVAR_BASE}/vcf_{assembly}/{vcf_name}.tbi"

    vcf_path = vcf_dir / vcf_name
    tbi_path = vcf_dir / f"{vcf_name}.tbi"

    try:
        _download(vcf_url, vcf_path)
    except FileNotFoundError:
        vcf_name = _discover_latest_vcf_name(assembly)
        vcf_url = f"{CLINVAR_BASE}/vcf_{assembly}/{vcf_name}"
        tbi_url = f"{CLINVAR_BASE}/vcf_{assembly}/{vcf_name}.tbi"
        vcf_path = vcf_dir / vcf_name
        tbi_path = vcf_dir / f"{vcf_name}.tbi"
        try:
            _download(vcf_url, vcf_path)
        except OSError as exc:
            if exc.errno == 27 and allow_parts:
                manifest_path = _download_parts(
                    vcf_url,
                    vcf_dir,
                    vcf_name,
                    part_size_mb * 1024 * 1024,
                )
                return manifest_path
            raise
    except OSError as exc:
        if exc.errno == 27 and allow_parts:
            manifest_path = _download_parts(
                vcf_url,
                vcf_dir,
                vcf_name,
                part_size_mb * 1024 * 1024,
            )
            return manifest_path
        raise

    try:
        _download(tbi_url, tbi_path)
    except FileNotFoundError:
        pass
    return vcf_path


def record_release_metadata(out_dir: Path) -> None:
    """Record hashes and sizes for ClinVar files in the output directory."""
    files: list[dict[str, Any]] = []
    for path in sorted(out_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.name.endswith(".partial"):
            continue
        if path.name.endswith(".parts.json"):
            files.append(
                {
                    "path": str(path.relative_to(out_dir)),
                    "size": path.stat().st_size,
                    "sha256": file_sha256(path),
                }
            )
            continue
        if ".part" in path.name:
            files.append(
                {
                    "path": str(path.relative_to(out_dir)),
                    "size": path.stat().st_size,
                    "sha256": file_sha256(path),
                }
            )
            continue
        if path.suffix.lower() not in {".xml", ".gz", ".tbi", ".vcf", ".json"}:
            continue
        files.append(
            {
                "path": str(path.relative_to(out_dir)),
                "size": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )

    metadata = {
        "generated_at": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "source_urls": {
            "xml": f"{CLINVAR_BASE}/xml/{CLINVAR_XML_NAME}",
            "vcf_GRCh38": f"{CLINVAR_BASE}/vcf_GRCh38/clinvar.vcf.gz",
            "vcf_GRCh37": f"{CLINVAR_BASE}/vcf_GRCh37/clinvar.vcf.gz",
        },
        "files": files,
    }

    out_path = out_dir / "clinvar_release.json"
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
