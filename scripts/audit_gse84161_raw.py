#!/usr/bin/env python3
"""Reproduce the metadata-only raw-CEL integrity audit for GSE84161.

The archive is authenticated against the pinned SHA-256 digest before the
tar members are opened.  Each member is checked as a safe regular tar member,
streamed through a gzip decoder to verify EOF and CRC, and hashed in both its
compressed and decompressed forms.  Only the fixed CEL binary header is
inspected for the expected format, dimensions, and chip name.  The intensity
payload is never interpreted and no expression value or score is computed.

The existing derived inventory and audit JSON are inputs.  They are compared
in memory and are never overwritten by this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import struct
import sys
import tarfile
import zlib
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterable
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_PATH = ROOT / "data/raw/gse84161/GSE84161_RAW.tar"
RAW_INVENTORY_PATH = ROOT / "data/derived/gse84161-raw-arrays.csv"
SAMPLE_MAP_PATH = ROOT / "data/derived/gse84161-samples.csv"
PRIOR_AUDIT_PATH = ROOT / "work/gse84161-study/raw-integrity-audit.json"
DOWNLOAD_RECEIPT_PATH = ROOT / "work/gse84161-study/raw-download-receipt.json"
CAPTURE_CONFIG_PATH = ROOT / "config/gse84161-raw-capture.json"
REPORT_PATH = ROOT / "reports/gse84161-raw-integrity.json"

EXPECTED_ARCHIVE_SHA256 = (
    "6d700c2e2cdcee9ef4eabb654b4c61f5d32952aa6f07bb00d11900f12ec78fff"
)
EXPECTED_SAMPLE_COUNT = 30
EXPECTED_CEL_MAGIC = 64
EXPECTED_CEL_VERSION = 4
EXPECTED_CEL_COLUMNS = 1164
EXPECTED_CEL_ROWS = 1164
EXPECTED_CEL_CELLS = EXPECTED_CEL_COLUMNS * EXPECTED_CEL_ROWS
EXPECTED_CHIP_HEADER = "HG-U133_Plus_2"
CEL_FIXED_HEADER_BYTES = 24
STREAM_CHUNK_BYTES = 1024 * 1024
HEADER_PREFIX_BYTES = 1024 * 1024

RAW_INVENTORY_FIELDS = [
    "gsm",
    "archive_member",
    "source_url",
    "compressed_bytes",
    "compressed_sha256",
    "decompressed_bytes",
    "decompressed_sha256",
    "gzip_integrity",
    "cel_magic",
    "cel_version",
    "columns",
    "rows",
    "intensity_cell_count",
    "chip_header_match",
]


class AuditError(RuntimeError):
    """Base class for a failed raw archive audit."""


class ArchiveHashMismatch(AuditError):
    """The archive did not match the pinned source digest."""


class InputMapError(AuditError):
    """The sample or existing inventory map is structurally invalid."""


class MemberSetError(AuditError):
    """The tar member set is unsafe or differs from the expected map."""


class GzipIntegrityError(AuditError):
    """A member did not decode to one complete gzip stream."""


class CelHeaderError(AuditError):
    """A decompressed member does not have the expected CEL header."""


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def sha256_file(path: Path, chunk_size: int = STREAM_CHUNK_BYTES) -> tuple[int, str]:
    """Return byte count and SHA-256 for a file without decoding its contents."""

    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
            total += len(chunk)
    return total, digest.hexdigest()


def _safe_regular_name(name: str) -> bool:
    """Accept only a plain POSIX basename with no traversal or separators."""

    if not name or "\x00" in name or "\\" in name:
        return False
    path = PurePosixPath(name)
    return (
        not path.is_absolute()
        and path.name == name
        and name not in {".", ".."}
        and ".." not in path.parts
    )


def _source_map(path: Path) -> dict[str, Any]:
    """Read the derived GSM-to-CEL URL map without touching expression data."""

    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise InputMapError(f"sample map has no header: {path}")
        rows = list(reader)

    required = {"sample_accession", "cel_url", "cel_available_route"}
    missing = sorted(required - set(reader.fieldnames))
    if missing:
        raise InputMapError(f"sample map missing columns: {missing}")
    if len(rows) != EXPECTED_SAMPLE_COUNT:
        raise InputMapError(
            f"sample map expected {EXPECTED_SAMPLE_COUNT} rows, found {len(rows)}"
        )

    by_member: dict[str, dict[str, str]] = {}
    by_gsm: dict[str, dict[str, str]] = {}
    order: list[str] = []
    errors: list[str] = []
    for row_number, row in enumerate(rows, start=2):
        gsm = row.get("sample_accession", "")
        url = row.get("cel_url", "")
        route = row.get("cel_available_route", "")
        if not gsm or not url:
            errors.append(f"row {row_number}: missing sample accession or CEL URL")
            continue
        if route.lower() != "true":
            errors.append(f"row {row_number}: CEL route is not marked true")
        parsed = urlparse(url)
        member = PurePosixPath(parsed.path).name
        if not _safe_regular_name(member) or not member.endswith(".CEL.gz"):
            errors.append(f"row {row_number}: unsafe or invalid CEL URL basename {member!r}")
            continue
        if not member.startswith(gsm + "_"):
            errors.append(f"row {row_number}: member {member!r} does not match {gsm!r}")
        if gsm in by_gsm:
            errors.append(f"row {row_number}: duplicate sample accession {gsm}")
        if member in by_member:
            errors.append(f"row {row_number}: duplicate archive member {member}")
        if gsm not in by_gsm and member not in by_member:
            # The source map preserves the FTP URL, while the existing
            # inventory canonicalizes that same public path to HTTPS.  URL
            # identity is checked by path below; output uses the stable HTTPS
            # spelling so the existing inventory can be reproduced exactly.
            record = {"gsm": gsm, "source_url": _canonical_url(url)}
            by_gsm[gsm] = record
            by_member[member] = record
            order.append(member)

    if errors:
        raise InputMapError("; ".join(errors))
    return {
        "rows": rows,
        "by_member": by_member,
        "by_gsm": by_gsm,
        "member_order": order,
        "sample_count": len(rows),
    }


def _canonical_url(url: str) -> str:
    """Canonicalize public FTP/HTTP sample URLs to the inventory spelling."""

    parsed = urlparse(url)
    if parsed.scheme.lower() == "ftp":
        parsed = parsed._replace(scheme="https")
    return parsed.geturl()


def _load_inventory(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Load the existing raw-array inventory for an in-memory comparison."""

    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise InputMapError(f"raw inventory has no header: {path}")
        fields = list(reader.fieldnames)
        if fields != RAW_INVENTORY_FIELDS:
            raise InputMapError(
                f"raw inventory columns differ from expected: {fields!r}"
            )
        rows = list(reader)

    if len(rows) != EXPECTED_SAMPLE_COUNT:
        raise InputMapError(
            f"raw inventory expected {EXPECTED_SAMPLE_COUNT} rows, found {len(rows)}"
        )
    members = [row["archive_member"] for row in rows]
    gsms = [row["gsm"] for row in rows]
    if len(set(members)) != len(members):
        raise InputMapError("raw inventory contains duplicate archive members")
    if len(set(gsms)) != len(gsms):
        raise InputMapError("raw inventory contains duplicate GSM identifiers")
    return fields, rows


def _member_table_entries(archive: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    """Return safe regular tar members and reject duplicate or unsafe names."""

    entries: dict[str, tarfile.TarInfo] = {}
    errors: list[str] = []
    for member in archive.getmembers():
        name = member.name
        if not _safe_regular_name(name):
            errors.append(f"unsafe member name {name!r}")
            continue
        if not member.isreg():
            errors.append(f"member is not a regular file {name!r}")
            continue
        if name in entries:
            errors.append(f"duplicate member name {name!r}")
            continue
        entries[name] = member
    if errors:
        raise MemberSetError("; ".join(errors))
    return entries


def _consume_decompressed(
    decompressed: bytes,
    digest: Any,
    prefix: bytearray,
    total: int,
) -> int:
    """Hash decompressed bytes and retain only a bounded header prefix."""

    if not decompressed:
        return total
    digest.update(decompressed)
    total += len(decompressed)
    if len(prefix) < HEADER_PREFIX_BYTES:
        prefix.extend(decompressed[: HEADER_PREFIX_BYTES - len(prefix)])
    return total


def _stream_gzip_member(stream: Any, name: str) -> dict[str, Any]:
    """Verify one gzip stream and return byte hashes plus a bounded prefix."""

    compressed_digest = hashlib.sha256()
    decompressed_digest = hashlib.sha256()
    compressed_total = 0
    decompressed_total = 0
    prefix = bytearray()
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)

    while True:
        chunk = stream.read(STREAM_CHUNK_BYTES)
        if not chunk:
            break
        compressed_digest.update(chunk)
        compressed_total += len(chunk)
        try:
            output = decoder.decompress(chunk)
        except zlib.error as error:
            raise GzipIntegrityError(f"{name}: gzip decode/CRC failure: {error}") from error
        decompressed_total = _consume_decompressed(
            output, decompressed_digest, prefix, decompressed_total
        )

    try:
        output = decoder.flush()
    except zlib.error as error:
        raise GzipIntegrityError(f"{name}: gzip flush failure: {error}") from error
    decompressed_total = _consume_decompressed(
        output, decompressed_digest, prefix, decompressed_total
    )

    if not decoder.eof:
        raise GzipIntegrityError(f"{name}: gzip stream did not reach EOF")
    if decoder.unused_data:
        raise GzipIntegrityError(f"{name}: trailing bytes after gzip stream")

    return {
        "compressed_bytes": str(compressed_total),
        "compressed_sha256": compressed_digest.hexdigest(),
        "decompressed_bytes": str(decompressed_total),
        "decompressed_sha256": decompressed_digest.hexdigest(),
        "gzip_integrity": "pass",
        "gzip_eof": True,
        "prefix": bytes(prefix),
    }


def _parse_cel_header(prefix: bytes, name: str) -> dict[str, str]:
    """Inspect only the fixed CEL header and its declared text header."""

    if len(prefix) < CEL_FIXED_HEADER_BYTES:
        raise CelHeaderError(f"{name}: truncated CEL fixed header")
    magic, version, columns, rows, cells, header_length = struct.unpack(
        "<6I", prefix[:CEL_FIXED_HEADER_BYTES]
    )
    if header_length > HEADER_PREFIX_BYTES - CEL_FIXED_HEADER_BYTES:
        raise CelHeaderError(
            f"{name}: declared header is larger than the bounded prefix ({header_length})"
        )
    end = CEL_FIXED_HEADER_BYTES + header_length
    if len(prefix) < end:
        raise CelHeaderError(f"{name}: truncated CEL text header")
    header_text = prefix[CEL_FIXED_HEADER_BYTES:end].decode("latin-1", errors="replace")

    if magic != EXPECTED_CEL_MAGIC:
        raise CelHeaderError(f"{name}: CEL magic {magic} != {EXPECTED_CEL_MAGIC}")
    if version != EXPECTED_CEL_VERSION:
        raise CelHeaderError(f"{name}: CEL version {version} != {EXPECTED_CEL_VERSION}")
    if columns != EXPECTED_CEL_COLUMNS or rows != EXPECTED_CEL_ROWS:
        raise CelHeaderError(
            f"{name}: dimensions {columns}x{rows} != "
            f"{EXPECTED_CEL_COLUMNS}x{EXPECTED_CEL_ROWS}"
        )
    if cells != EXPECTED_CEL_CELLS:
        raise CelHeaderError(
            f"{name}: intensity cell count {cells} != {EXPECTED_CEL_CELLS}"
        )
    if EXPECTED_CHIP_HEADER not in header_text:
        raise CelHeaderError(
            f"{name}: header does not contain {EXPECTED_CHIP_HEADER!r}"
        )

    return {
        "cel_magic": str(magic),
        "cel_version": str(version),
        "columns": str(columns),
        "rows": str(rows),
        "intensity_cell_count": str(cells),
        "chip_header_match": EXPECTED_CHIP_HEADER,
    }


def audit_archive(
    archive_path: Path,
    expected_archive_sha256: str,
    expected_members: Iterable[str],
    expected_by_member: dict[str, dict[str, str]] | None = None,
    member_order: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Authenticate, inventory, and header-check a raw GSE84161 tar archive."""

    archive_bytes, archive_sha256 = sha256_file(archive_path)
    if archive_sha256.lower() != expected_archive_sha256.lower():
        raise ArchiveHashMismatch(
            f"archive SHA-256 {archive_sha256} != pinned {expected_archive_sha256}"
        )

    expected_list_source = list(expected_members)
    expected_set = set(expected_list_source)
    expected_list = (
        list(member_order) if member_order is not None else list(expected_list_source)
    )
    if len(expected_set) != len(expected_list_source):
        raise MemberSetError("expected archive member map contains duplicates")
    if len(set(expected_list)) != len(expected_list) or set(expected_list) != expected_set:
        raise MemberSetError("member_order is not a unique ordering of expected_members")
    if not expected_set or any(not _safe_regular_name(name) for name in expected_set):
        raise MemberSetError("expected archive member map contains unsafe names")

    with tarfile.open(archive_path, mode="r:") as archive:
        entries = _member_table_entries(archive)
        actual_set = set(entries)
        if actual_set != expected_set:
            missing = sorted(expected_set - actual_set)
            unexpected = sorted(actual_set - expected_set)
            raise MemberSetError(
                f"archive member set mismatch; missing={missing!r}, "
                f"unexpected={unexpected!r}"
            )

        rows: list[dict[str, str]] = []
        dimensions: set[tuple[int, int, int]] = set()
        total_decompressed_bytes = 0
        for name in expected_list:
            member = entries[name]
            extracted = archive.extractfile(member)
            if extracted is None:
                raise MemberSetError(f"{name}: tar member could not be opened")
            with extracted:
                gzip_info = _stream_gzip_member(extracted, name)
            cel_info = _parse_cel_header(gzip_info.pop("prefix"), name)
            record: dict[str, str] = {
                "gsm": (expected_by_member or {}).get(name, {}).get("gsm", name.split("_", 1)[0]),
                "archive_member": name,
                "source_url": (expected_by_member or {}).get(name, {}).get("source_url", ""),
                "gzip_eof": str(bool(gzip_info.pop("gzip_eof"))).lower(),
            }
            record.update({key: value for key, value in gzip_info.items() if key in RAW_INVENTORY_FIELDS})
            record.update(cel_info)
            rows.append(record)
            dimensions.add(
                (
                    int(cel_info["columns"]),
                    int(cel_info["rows"]),
                    int(cel_info["intensity_cell_count"]),
                )
            )
            total_decompressed_bytes += int(record["decompressed_bytes"])

    return {
        "archive_bytes": archive_bytes,
        "archive_sha256": archive_sha256,
        "source_expected_members": len(expected_set),
        "archive_member_count": len(entries),
        "archive_regular_members": len(entries),
        "safe_regular_member_count": len(entries),
        "all_members_safe_regular": True,
        "all_names_match_SOFT_CEL_urls": actual_set == expected_set,
        "all_sample_map_members_present": expected_set <= actual_set,
        "all_gzip_eof_crc_checks_pass": all(
            row["gzip_integrity"] == "pass" and row["gzip_eof"] == "true"
            for row in rows
        ),
        "all_gzip_CRC_checks_pass": all(
            row["gzip_integrity"] == "pass" for row in rows
        ),
        "all_CEL_v4_headers_match_expected_array": all(
            row["cel_magic"] == str(EXPECTED_CEL_MAGIC)
            and row["cel_version"] == str(EXPECTED_CEL_VERSION)
            and row["columns"] == str(EXPECTED_CEL_COLUMNS)
            and row["rows"] == str(EXPECTED_CEL_ROWS)
            and row["intensity_cell_count"] == str(EXPECTED_CEL_CELLS)
            and row["chip_header_match"] == EXPECTED_CHIP_HEADER
            for row in rows
        ),
        "dimensions": [list(item) for item in sorted(dimensions)],
        "total_decompressed_bytes": total_decompressed_bytes,
        "members": rows,
    }


def _compare_inventory(
    computed_rows: list[dict[str, str]],
    existing_rows: list[dict[str, str]],
) -> dict[str, Any]:
    """Compare every inventory field while keeping the source CSV untouched."""

    mismatches: list[dict[str, Any]] = []
    computed_by_member = {row["archive_member"]: row for row in computed_rows}
    existing_by_member = {row["archive_member"]: row for row in existing_rows}
    for member in sorted(set(computed_by_member) | set(existing_by_member)):
        computed = computed_by_member.get(member)
        existing = existing_by_member.get(member)
        if computed is None or existing is None:
            mismatches.append(
                {"archive_member": member, "computed_present": computed is not None,
                 "existing_present": existing is not None}
            )
            continue
        for field in RAW_INVENTORY_FIELDS:
            if computed.get(field, "") != existing.get(field, ""):
                mismatches.append(
                    {
                        "archive_member": member,
                        "field": field,
                        "computed": computed.get(field, ""),
                        "existing": existing.get(field, ""),
                    }
                )
    existing_order = [row["archive_member"] for row in existing_rows]
    computed_order = [row["archive_member"] for row in computed_rows]
    if existing_order != computed_order:
        mismatches.append(
            {
                "field": "row_order",
                "computed": computed_order,
                "existing": existing_order,
            }
        )
    return {
        "inventory_matches_existing": not mismatches,
        "inventory_mismatches": mismatches,
        "existing_row_count": len(existing_rows),
        "computed_row_count": len(computed_rows),
    }


def _compare_prior_audit(
    computed: dict[str, Any],
    prior: dict[str, Any],
    inventory_sha256: str,
    inventory_path: Path,
    prior_path_label: str,
) -> dict[str, Any]:
    """Compare the fields of the earlier audit that this implementation reproduces."""

    expected = {
        "archive_bytes": computed["archive_bytes"],
        "archive_sha256": computed["archive_sha256"],
        "source_expected_members": computed["source_expected_members"],
        "archive_regular_members": computed["archive_regular_members"],
        "all_names_match_SOFT_CEL_urls": computed["all_names_match_SOFT_CEL_urls"],
        "all_gzip_CRC_checks_pass": computed["all_gzip_CRC_checks_pass"],
        "all_CEL_v4_headers_match_expected_array": computed[
            "all_CEL_v4_headers_match_expected_array"
        ],
        "dimensions": computed["dimensions"],
        "total_decompressed_bytes": computed["total_decompressed_bytes"],
        "inventory_path": _relative(inventory_path),
        "inventory_sha256": inventory_sha256,
        "assay_values_interpreted": False,
        "gene_expression_values_computed": False,
    }
    mismatches = [
        {"field": key, "computed": value, "prior": prior.get(key)}
        for key, value in expected.items()
        if prior.get(key) != value
    ]
    return {
        "prior_audit_matches": not mismatches,
        "prior_audit_mismatches": mismatches,
        "prior_audit_path": prior_path_label,
        "fields_compared": sorted(expected),
    }


def _source_hash(path: Path, role: str) -> dict[str, Any]:
    byte_count, digest = sha256_file(path)
    return {
        "path": _relative(path),
        "bytes": byte_count,
        "sha256": digest,
        "role": role,
    }


def build_report(
    archive_path: Path = ARCHIVE_PATH,
    raw_inventory_path: Path = RAW_INVENTORY_PATH,
    sample_map_path: Path = SAMPLE_MAP_PATH,
    prior_audit_path: Path | None = None,
    download_receipt_path: Path | None = None,
    capture_path: Path | None = CAPTURE_CONFIG_PATH,
    expected_archive_sha256: str = EXPECTED_ARCHIVE_SHA256,
    report_path: Path = REPORT_PATH,
) -> dict[str, Any]:
    """Build a reproducibility report without writing any input inventory."""

    source_map = _source_map(sample_map_path)
    inventory_fields, existing_rows = _load_inventory(raw_inventory_path)
    if capture_path is not None and (
        prior_audit_path is not None or download_receipt_path is not None
    ):
        raise InputMapError(
            "provide capture_path or both prior_audit_path and download_receipt_path"
        )
    if capture_path is not None:
        capture = json.loads(capture_path.read_text(encoding="utf-8"))
        try:
            prior = capture["initial_integrity_audit"]
            receipt = capture["raw_download_receipt"]
        except KeyError as error:
            raise InputMapError(
                f"capture config missing {error.args[0]!r}"
            ) from error
        prior_path_label = (
            f"{_relative(capture_path)}#initial_integrity_audit"
        )
        receipt_path_label = (
            f"{_relative(capture_path)}#raw_download_receipt"
        )
    else:
        if prior_audit_path is None or download_receipt_path is None:
            raise InputMapError(
                "both prior_audit_path and download_receipt_path are required "
                "when capture_path is omitted"
            )
        capture = None
        prior = json.loads(prior_audit_path.read_text(encoding="utf-8"))
        receipt = json.loads(download_receipt_path.read_text(encoding="utf-8"))
        prior_path_label = _relative(prior_audit_path)
        receipt_path_label = _relative(download_receipt_path)
    inventory_bytes, inventory_sha256 = sha256_file(raw_inventory_path)

    computed = audit_archive(
        archive_path=archive_path,
        expected_archive_sha256=expected_archive_sha256,
        expected_members=source_map["member_order"],
        expected_by_member=source_map["by_member"],
        member_order=source_map["member_order"],
    )
    computed["all_inventory_members_match_sample_map"] = set(
        row["archive_member"] for row in existing_rows
    ) == set(source_map["member_order"])
    if not computed["all_inventory_members_match_sample_map"]:
        raise InputMapError("existing raw inventory members differ from sample map")
    if any(
        existing.get("gsm") != source_map["by_member"][existing["archive_member"]]["gsm"]
        or existing.get("source_url")
        != source_map["by_member"][existing["archive_member"]]["source_url"]
        for existing in existing_rows
    ):
        raise InputMapError("existing raw inventory GSM or source URL differs from sample map")

    inventory_comparison = _compare_inventory(computed["members"], existing_rows)
    prior_comparison = _compare_prior_audit(
        computed,
        prior,
        inventory_sha256,
        raw_inventory_path,
        prior_path_label,
    )

    receipt_match = (
        receipt.get("status") == 200
        and receipt.get("bytes") == computed["archive_bytes"]
        and receipt.get("sha256", "").lower() == computed["archive_sha256"].lower()
    )
    if not receipt_match:
        raise InputMapError("raw download receipt does not match the authenticated archive")
    if not inventory_comparison["inventory_matches_existing"]:
        raise InputMapError(
            "computed raw member inventory differs from the existing inventory"
        )
    if not prior_comparison["prior_audit_matches"]:
        raise InputMapError(
            "computed raw audit differs from the existing integrity audit"
        )

    source_inputs = [
        _source_hash(archive_path, "authenticated raw archive"),
        _source_hash(sample_map_path, "derived SOFT/CEL sample map"),
        _source_hash(raw_inventory_path, "existing root-owned raw inventory"),
    ]
    if capture_path is not None:
        source_inputs.append(_source_hash(capture_path, "versioned raw audit capture"))
    else:
        source_inputs.extend(
            [
                _source_hash(prior_audit_path, "existing raw-integrity audit"),
                _source_hash(download_receipt_path, "raw download receipt"),
            ]
        )
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_scope": (
            "Raw tar member identity, gzip EOF/CRC, compressed/decompressed byte "
            "hashes, and fixed CEL header metadata for GSE84161."
        ),
        "archive": {
            "path": _relative(archive_path),
            "bytes": computed["archive_bytes"],
            "sha256": computed["archive_sha256"],
            "pinned_sha256": expected_archive_sha256,
            "hash_checked_before_tar_member_read": True,
            "download_receipt_matches": receipt_match,
        },
        "sample_map": {
            "path": _relative(sample_map_path),
            "sample_count": source_map["sample_count"],
            "expected_member_count": len(source_map["member_order"]),
            "all_cel_routes_marked_available": True,
            "map_is_used_as_the_SOFT_CEL_membership_source": True,
            "direct_SOFT_reparse_performed": False,
            "url_identity_rule": (
                "FTP sample URLs are canonicalized to HTTPS for comparison with "
                "the existing inventory; URL paths and basenames must match."
            ),
        },
        "member_accounting": {
            key: computed[key]
            for key in (
                "source_expected_members",
                "archive_member_count",
                "archive_regular_members",
                "safe_regular_member_count",
                "all_members_safe_regular",
                "all_names_match_SOFT_CEL_urls",
                "all_sample_map_members_present",
                "all_inventory_members_match_sample_map",
            )
        },
        "integrity": {
            key: computed[key]
            for key in (
                "all_gzip_eof_crc_checks_pass",
                "all_gzip_CRC_checks_pass",
                "all_CEL_v4_headers_match_expected_array",
                "dimensions",
                "total_decompressed_bytes",
            )
        },
        "members": computed["members"],
        "reproduction": {
            "existing_inventory": {
                "path": _relative(raw_inventory_path),
                "bytes": inventory_bytes,
                "sha256": inventory_sha256,
                "fields": inventory_fields,
                **inventory_comparison,
            },
            **prior_comparison,
        },
        "capture": {
            "path": _relative(capture_path) if capture_path is not None else None,
            "embedded_prior_audit": capture is not None,
            "embedded_download_receipt": capture is not None,
            "prior_audit_source": prior_path_label,
            "download_receipt_source": receipt_path_label,
        },
        "source_inputs": source_inputs,
        "output_hashes": {
            "preserved_raw_inventory_sha256": inventory_sha256,
            "member_compressed_and_decompressed_hashes": True,
            "report_path": _relative(report_path),
            "report_sha256_is_printed_after_write": True,
        },
        "limits": {
            "assay_values_interpreted": False,
            "gene_expression_values_computed": False,
            "probe_intensities_interpreted": False,
            "program_scores_computed": False,
            "raw_cel_payload_interpreted": False,
            "raw_cel_decompressed_for_crc_and_header_hashing": True,
            "note": (
                "Compressed members were hashed; gzip streams were decompressed "
                "to EOF/CRC and hashed; only the fixed CEL header and declared text "
                "header were inspected. Intensity payload bytes were not parsed."
            ),
        },
    }
    return report


def write_report(report: dict[str, Any], output_path: Path = REPORT_PATH) -> tuple[int, str]:
    """Write only the requested report and return its on-disk hash."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=False) + "\n"
    output_path.write_text(payload, encoding="utf-8")
    return sha256_file(output_path)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE_PATH)
    parser.add_argument("--raw-inventory", type=Path, default=RAW_INVENTORY_PATH)
    parser.add_argument("--sample-map", type=Path, default=SAMPLE_MAP_PATH)
    parser.add_argument("--capture-config", type=Path, default=CAPTURE_CONFIG_PATH)
    parser.add_argument("--prior-audit", type=Path, default=None)
    parser.add_argument("--download-receipt", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=REPORT_PATH)
    parser.add_argument("--expected-archive-sha256", default=EXPECTED_ARCHIVE_SHA256)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    try:
        capture_path = args.capture_config
        if args.prior_audit is not None or args.download_receipt is not None:
            capture_path = None
        report = build_report(
            archive_path=args.archive,
            raw_inventory_path=args.raw_inventory,
            sample_map_path=args.sample_map,
            prior_audit_path=args.prior_audit,
            download_receipt_path=args.download_receipt,
            capture_path=capture_path,
            expected_archive_sha256=args.expected_archive_sha256,
            report_path=args.output,
        )
    except (AuditError, OSError, json.JSONDecodeError) as error:
        print(f"raw GSE84161 audit failed: {error}", file=sys.stderr)
        return 2
    report_bytes, report_sha256 = write_report(report, args.output)
    print(
        json.dumps(
            {
                "report_path": _relative(args.output),
                "report_bytes": report_bytes,
                "report_sha256": report_sha256,
                "archive_sha256": report["archive"]["sha256"],
                "inventory_matches_existing": report["reproduction"][
                    "existing_inventory"
                ]["inventory_matches_existing"],
                "prior_audit_matches": report["reproduction"]["prior_audit_matches"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
