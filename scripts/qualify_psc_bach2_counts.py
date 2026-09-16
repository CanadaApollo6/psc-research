#!/usr/bin/env python3
"""Qualify one authorized BACH2 processed-count archive; never compute effects."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, HTTPRedirectHandler, build_opener

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/psc-bach2-qualification/counts"
WORK = ROOT / "work/psc-bach2-count-qualification"
FILENAME = "sample05_filtered_feature_bc_matrix.tar.gz"
EXPECTED_BYTES = 26_893_378
URL = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-14013/" + FILENAME
ALLOWED_HOSTS = {"www.ebi.ac.uk", "ftp.ebi.ac.uk"}
CAPS = {"outer_uncompressed_bytes": 700_000_000, "tar_members": 30,
        "tar_regular_member_bytes": 600_000_000, "nested_uncompressed_bytes": 600_000_000,
        "matrix_dimensions_each": 1_000_000, "matrix_coordinates": 50_000_000}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1_048_576):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def json_write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def check_url(url: str) -> None:
    parsed = urlsplit(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS or parsed.username or parsed.password:
        raise ValueError("Only the allowlisted ordinary public HTTPS file routes are permitted")
    if PurePosixPath(parsed.path).name != FILENAME or parsed.query:
        raise ValueError("Only the one authorized filtered archive is permitted")


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def acquire(raw_dir: Path = RAW) -> dict:
    """One exact file; no retries or arbitrary URLs; all retained bytes are charged."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = raw_dir / "acquisition-receipt.json"
    archive_path = raw_dir / FILENAME
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if (receipt.get("complete") and archive_path.exists()
                and archive_path.stat().st_size == EXPECTED_BYTES
                and file_sha256(archive_path) == receipt["sha256"]):
            return receipt
        raise ValueError("Prior attempt exists; no silent overwrite or fresh byte allowance")
    if archive_path.exists():
        raise ValueError("Archive path already exists without a verified receipt")
    receipt = {"url": URL, "started_at_utc": utc_now(), "requests": [], "complete": False,
               "body_bytes": 0, "error": None, "authorized_archive_bytes": EXPECTED_BYTES,
               "additional_metadata_body_bytes": 0}
    opener = build_opener(NoRedirect())
    response = None
    digest = hashlib.sha256()
    url = URL
    try:
        for _ in range(6):
            check_url(url)
            request = Request(url, method="GET", headers={"User-Agent": "PSC-public-count-qualification/1.0", "Accept-Encoding": "identity"})
            try:
                response = opener.open(request, timeout=30)
            except HTTPError as exc:
                response = exc
            receipt["requests"].append({"url": url, "method": "GET", "status": response.status,
                                        "headers": dict(response.headers.items())})
            if response.status in {301, 302, 303, 307, 308}:
                location = response.headers.get("Location")
                response.close()
                response = None
                if not location:
                    raise ValueError("Redirect lacks Location")
                url = urljoin(url, location)
                continue
            if response.status != 200:
                raise ValueError("Non-200 archive response; no body read")
            length = response.headers.get("Content-Length")
            if length is None or int(length) != EXPECTED_BYTES:
                raise ValueError("Declared length differs from exact authorized inventory; no body read")
            if response.headers.get("Content-Encoding", "identity") != "identity":
                raise ValueError("Unexpected transport encoding; no body read")
            # x mode never overwrites a source path. No tar path is extracted.
            with archive_path.open("xb") as out:
                while receipt["body_bytes"] < EXPECTED_BYTES:
                    block = response.read(min(65_536, EXPECTED_BYTES - receipt["body_bytes"]))
                    if not block:
                        break
                    out.write(block)
                    digest.update(block)
                    receipt["body_bytes"] += len(block)
            receipt["complete"] = receipt["body_bytes"] == EXPECTED_BYTES
            receipt["response_url"] = response.geturl()
            if not receipt["complete"]:
                raise ValueError("Short body; incomplete source retained and budget charged")
            break
        else:
            raise ValueError("Redirect limit exceeded")
    except Exception as exc:
        receipt["error"] = type(exc).__name__ + ": " + str(exc)
    finally:
        if response is not None:
            response.close()
    receipt["sha256"] = digest.hexdigest()
    receipt["completed_at_utc"] = utc_now()
    receipt["archive_path"] = str(archive_path.relative_to(ROOT)) if archive_path.is_relative_to(ROOT) else str(archive_path)
    json_write(receipt_path, receipt)
    if not receipt["complete"]:
        raise ValueError(receipt["error"] or "Incomplete download")
    return receipt




DEFINITION_URL = "https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.1/output/matrices"
DEFINITION_HOSTS = {"support.10xgenomics.com", "www.10xgenomics.com"}


def acquire_definition(raw_dir: Path = RAW, *, name: str = "count-definition", url: str = DEFINITION_URL) -> dict:
    """Bounded public format definitions share a separate cumulative 1 MB cap."""
    if name not in {"count-definition", "count-definition-current"}:
        raise ValueError("Unknown definition source name")
    receipt_path = raw_dir / (name + "-receipt.json")
    body_path = raw_dir / (name + ".body")
    if receipt_path.exists():
        previous = json.loads(receipt_path.read_text())
        if body_path.exists() and file_sha256(body_path) == previous["sha256"]:
            return previous
        raise ValueError("Definition receipt/cache mismatch; no overwrite")
    if body_path.exists():
        raise ValueError("Definition cache path already exists")
    opener = build_opener(NoRedirect())
    spent = sum(json.loads(p.read_text())["body_bytes"] for p in raw_dir.glob("count-definition*-receipt.json"))
    remaining = 1_000_000 - spent
    if remaining <= 0:
        raise ValueError("Cumulative definition metadata budget exhausted")
    receipt = {"url": url, "requests": [], "started_at_utc": utc_now(),
               "complete_response": False, "body_bytes": 0, "metadata_cap_bytes": 1_000_000, "error": None}
    response = None
    parts = []
    try:
        for _ in range(6):
            parsed = urlsplit(url)
            if parsed.scheme != "https" or parsed.hostname not in DEFINITION_HOSTS or parsed.username or parsed.password:
                raise ValueError("Definition redirect outside ordinary official HTTPS hosts")
            request = Request(url, headers={"User-Agent": "PSC-public-count-qualification/1.0", "Accept-Encoding": "identity"})
            try:
                response = opener.open(request, timeout=25)
            except HTTPError as exc:
                response = exc
            receipt["requests"].append({"url": url, "status": response.status, "headers": dict(response.headers.items())})
            if response.status in {301, 302, 303, 307, 308}:
                location = response.headers.get("Location")
                response.close()
                response = None
                if not location:
                    raise ValueError("Redirect lacks Location")
                url = urljoin(url, location)
                continue
            if response.status != 200:
                raise ValueError("Non-200 definition response; no body read")
            length = response.headers.get("Content-Length")
            if length and int(length) > remaining:
                raise ValueError("Definition body exceeds metadata allowance; no body read")
            if response.headers.get("Content-Encoding", "identity") != "identity":
                raise ValueError("Unexpected transport encoding")
            while receipt["body_bytes"] < remaining:
                block = response.read(min(65_536, remaining - receipt["body_bytes"]))
                if not block:
                    receipt["complete_response"] = True
                    break
                parts.append(block)
                receipt["body_bytes"] += len(block)
            if length and receipt["body_bytes"] == int(length):
                receipt["complete_response"] = True
            receipt["response_url"] = response.geturl()
            if not receipt["complete_response"]:
                raise ValueError("Definition cap reached without confirmed EOF")
            break
        else:
            raise ValueError("Definition redirect limit reached")
    except Exception as exc:
        receipt["error"] = type(exc).__name__ + ": " + str(exc)
    finally:
        if response is not None:
            response.close()
    body = b"".join(parts)
    with body_path.open("xb") as out:
        out.write(body)
    receipt["sha256"] = sha256(body)
    receipt["completed_at_utc"] = utc_now()
    json_write(receipt_path, receipt)
    return receipt


def decompress_outer(archive_path: Path, work_dir: Path = WORK) -> tuple[Path, dict]:
    """Fully consume gzip to validate CRC/ISIZE/EOF, with bounded disk output."""
    archive_hash = file_sha256(archive_path)
    destination = work_dir / "outer.tar"
    receipt_path = work_dir / "outer-gzip-integrity.json"
    work_dir.mkdir(parents=True, exist_ok=True)
    if destination.exists() or receipt_path.exists():
        if not (destination.exists() and receipt_path.exists()):
            raise ValueError("Incomplete work artifact exists; no overwrite")
        previous = json.loads(receipt_path.read_text())
        if (previous["archive_sha256"] != archive_hash or file_sha256(destination) != previous["tar_sha256"]
                or destination.stat().st_size != previous["uncompressed_bytes"]):
            raise ValueError("Outer-gzip work cache differs from receipt")
        return destination, previous
    digest = hashlib.sha256()
    total = 0
    with archive_path.open("rb") as source:
        if source.read(2) != b"\x1f\x8b":
            raise ValueError("Expected gzip magic")
        source.seek(0)
        with gzip.GzipFile(fileobj=source, mode="rb") as gz, destination.open("xb") as out:
            while True:
                block = gz.read(min(1_048_576, CAPS["outer_uncompressed_bytes"] - total + 1))
                if not block:
                    break
                total += len(block)
                if total > CAPS["outer_uncompressed_bytes"]:
                    raise ValueError("Outer gzip exceeds uncompressed-byte cap")
                digest.update(block)
                out.write(block)
    result = {"archive_sha256": archive_hash, "uncompressed_bytes": total,
              "tar_sha256": digest.hexdigest(), "gzip_CRC_ISIZE_and_EOF_passed": True}
    json_write(receipt_path, result)
    return destination, result


def canonical_member_name(name: str) -> str:
    if not name or "\x00" in name or "\\" in name or name.startswith("/") or re.match(r"^[A-Za-z]:", name):
        raise ValueError("Unsafe TAR member path")
    parts = name.split("/")
    if ".." in parts:
        raise ValueError("Parent traversal in TAR path")
    canonical = "/".join(x for x in parts if x not in {"", "."})
    if not canonical:
        raise ValueError("Empty canonical TAR member path")
    return canonical


def inventory_tar(tar_path: Path) -> list[dict]:
    """Validate every header/type/path, and require a complete zero TAR trailer."""
    result = []
    names = set()
    last_payload_end = 0
    with tarfile.open(tar_path, mode="r:") as archive:
        for member in archive:
            if len(result) >= CAPS["tar_members"]:
                raise ValueError("TAR member-count cap exceeded")
            canonical = canonical_member_name(member.name)
            if canonical in names:
                raise ValueError("Duplicate or aliased TAR member path")
            names.add(canonical)
            if not (member.isdir() or member.isfile()) or member.issparse():
                raise ValueError("TAR links, sparse and special member types are prohibited")
            if member.size < 0 or member.size > CAPS["tar_regular_member_bytes"]:
                raise ValueError("TAR member size exceeds cap")
            if member.isdir() and member.size:
                raise ValueError("Directory with nonzero payload")
            end = member.offset_data + ((member.size + 511) // 512) * 512
            if end > tar_path.stat().st_size:
                raise ValueError("TAR member extends beyond file")
            last_payload_end = max(last_payload_end, end)
            result.append({"source_path": member.name, "canonical_path": canonical,
                           "type": "directory" if member.isdir() else "regular_file",
                           "size_bytes": member.size, "offset_data": member.offset_data,
                           "pax_headers": member.pax_headers})
    regular = {x["canonical_path"] for x in result if x["type"] == "regular_file"}
    for item in result:
        parents = PurePosixPath(item["canonical_path"]).parents
        if any(str(parent) in regular for parent in parents):
            raise ValueError("A TAR regular file is also a path ancestor")
    trailing_bytes = tar_path.stat().st_size - last_payload_end
    if trailing_bytes < 1024 or tar_path.stat().st_size % 512:
        raise ValueError("Missing two-block TAR end marker or partial TAR block")
    with tar_path.open("rb") as stream:
        stream.seek(last_payload_end)
        while block := stream.read(1_048_576):
            if any(block):
                raise ValueError("Nonzero data after last TAR member; unparsed appendage")
    return result


def inventory(raw_dir: Path = RAW, work_dir: Path = WORK) -> dict:
    receipt = json.loads((raw_dir / "acquisition-receipt.json").read_text())
    archive_path = raw_dir / FILENAME
    if not receipt["complete"] or archive_path.stat().st_size != EXPECTED_BYTES or file_sha256(archive_path) != receipt["sha256"]:
        raise ValueError("Archive disagrees with acquisition receipt")
    tar_path, outer = decompress_outer(archive_path, work_dir)
    members = inventory_tar(tar_path)
    result = {"outer_gzip": outer, "members": members,
              "tar_member_headers_paths_types_and_zero_trailer_passed": True,
              "no_source_member_paths_extracted": True}
    json_write(work_dir / "archive-inventory.json", result)
    return result



def bounded_copy(source, destination: Path, byte_cap: int) -> tuple[int, str]:
    """Copy an entire stream without replacing any source or existing artifact."""
    digest = hashlib.sha256()
    total = 0
    with destination.open("xb") as output:
        while True:
            block = source.read(min(1_048_576, byte_cap - total + 1))
            if not block:
                break
            total += len(block)
            if total > byte_cap:
                raise ValueError("Decoded member exceeds byte cap")
            output.write(block)
            digest.update(block)
    return total, digest.hexdigest()


def prepare_members(tar_path: Path, members: list[dict], work_dir: Path = WORK) -> dict:
    roles = {"barcodes.tsv.gz": "barcodes", "features.tsv.gz": "features", "matrix.mtx.gz": "matrix"}
    regular = [x for x in members if x["type"] == "regular_file"]
    if Counter(PurePosixPath(x["source_path"]).name for x in regular) != Counter(roles.keys()):
        raise ValueError("Expected exactly the three explicit gzip-wrapped 10x matrix members")
    parents = {str(PurePosixPath(x["source_path"]).parent) for x in regular}
    if len(parents) != 1:
        raise ValueError("Feature, barcode and matrix members do not share one source directory")
    directory = work_dir / "members"
    directory.mkdir(parents=True, exist_ok=True)
    output = {}
    with tarfile.open(tar_path, mode="r:") as archive:
        for member in regular:
            basename = PurePosixPath(member["source_path"]).name
            role = roles[basename]
            compressed = directory / (role + ".gz")
            decoded = directory / (role + ".txt")
            pin_path = directory / (role + "-integrity.json")
            if pin_path.exists():
                pin = json.loads(pin_path.read_text())
                if (pin["source_path"] != member["source_path"]
                        or not compressed.exists() or not decoded.exists()
                        or file_sha256(compressed) != pin["compressed_sha256"]
                        or file_sha256(decoded) != pin["decoded_sha256"]):
                    raise ValueError("Member cache differs from provenance receipt")
                source = archive.extractfile(member["source_path"])
                if source is None:
                    raise ValueError("Cached regular member has no source stream")
                source_digest = hashlib.sha256()
                source_n = 0
                with source:
                    while block := source.read(1_048_576):
                        source_digest.update(block)
                        source_n += len(block)
                if source_n != pin["compressed_bytes"] or source_digest.hexdigest() != pin["compressed_sha256"]:
                    raise ValueError("Cached compressed member differs from the pinned TAR payload")
                pin["source_payload_matches_tar"] = True
                output[role] = pin
                continue
            if compressed.exists() or decoded.exists():
                raise ValueError("Unverified member artifact exists; no overwrite")
            source = archive.extractfile(member["source_path"])
            if source is None:
                raise ValueError("Regular member has no stream")
            with source:
                compressed_n, compressed_hash = bounded_copy(source, compressed, member["size_bytes"])
            if compressed_n != member["size_bytes"]:
                raise ValueError("Short TAR member read")
            with gzip.open(compressed, "rb") as uncompressed:
                decoded_n, decoded_hash = bounded_copy(uncompressed, decoded, CAPS["nested_uncompressed_bytes"])
            pin = {"source_path": member["source_path"], "compressed_bytes": compressed_n,
                   "compressed_sha256": compressed_hash, "decoded_bytes": decoded_n,
                   "decoded_sha256": decoded_hash, "nested_gzip_CRC_ISIZE_and_EOF_passed": True,
                   "decoded_path": str(decoded.relative_to(ROOT)) if decoded.is_relative_to(ROOT) else str(decoded)}
            pin["source_payload_matches_tar"] = True
            json_write(pin_path, pin)
            output[role] = pin
    return output


def parse_features(path: Path) -> tuple[list[list[str]], dict]:
    with path.open(newline="") as source:
        rows = list(csv.reader(source, delimiter="\t"))
    if not rows or {len(row) for row in rows} != {3}:
        raise ValueError("Expected a nonempty three-column feature axis; no inferred feature types")
    if any(not all(row) for row in rows):
        raise ValueError("Empty original feature field")
    ids = [row[0] for row in rows]
    symbols = [row[1] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate source feature IDs; do not silently collapse")
    return rows, {
        "rows": len(rows), "columns": 3,
        "column_convention": ["feature ID", "feature name", "feature type"],
        "unique_feature_ids": len(set(ids)), "unique_feature_names": len(set(symbols)),
        "duplicate_feature_name_rows": len(symbols) - len(set(symbols)),
        "feature_types_as_published": dict(sorted(Counter(r[2] for r in rows).items())),
        "ensembl_human_gene_shaped_IDs": sum(bool(re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", x)) for x in ids),
        "ensembl_human_gene_shaped_versioned_IDs": sum(bool(re.fullmatch(r"ENSG[0-9]+\.[0-9]+", x)) for x in ids),
        "full_deposited_axis_read": True, "independent_gene_identifier_reconciliation": False,
        "original_pre_gene_filter_universe_proven": False,
    }


def parse_barcodes(path: Path) -> list[str]:
    with path.open(newline="") as source:
        rows = list(csv.reader(source, delimiter="\t"))
    if not rows or any(len(row) != 1 or not row[0] for row in rows):
        raise ValueError("Expected one nonempty barcode per source row")
    barcodes = [r[0] for r in rows]
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("Duplicate source barcode")
    return barcodes


def exact_value(token: str, field: str):
    """Determine finite/nonnegative/integral status before any float conversion."""
    if len(token) > 100:
        raise ValueError("Numeric token exceeds safety limit")
    if field == "integer":
        if not re.fullmatch(r"[+-]?[0-9]+", token):
            raise ValueError("Noninteger token in integer MatrixMarket field")
        value = int(token)
    elif field == "real":
        if not re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", token):
            raise ValueError("Unsupported/nonfinite numeric token")
        value = Decimal(token)
        if not value.is_finite() or abs(value.adjusted()) > 308:
            raise ValueError("Nonfinite or extreme exponent")
        if value == value.to_integral_value():
            value = int(value)
        else:
            value = Fraction(value)  # Exact rational sums; no Decimal-context rounding.
    else:
        raise ValueError("Unsupported MatrixMarket field")
    if value < 0:
        raise ValueError("Negative matrix entry")
    return value


def parse_matrix(path: Path, feature_rows: list[list[str]], barcodes: list[str]) -> tuple[dict, dict]:
    from array import array
    types = [row[2] for row in feature_rows]
    sums = {kind: [0] * len(barcodes) for kind in sorted(set(types))}
    stored_by_type = Counter()
    nonzero_features = set()
    comment_lines = []
    zeros = fractions = nread = 0
    minimum = maximum = None
    keys = array("Q")
    previous_key = -1
    strictly_increasing = True
    with path.open() as source:
        header = source.readline().rstrip("\r\n")
        pieces = header.split()
        if len(pieces) != 5 or pieces[:3] != ["%%MatrixMarket", "matrix", "coordinate"] or pieces[3] not in {"integer", "real"} or pieces[4] != "general":
            raise ValueError("Unsupported MatrixMarket header")
        field = pieces[3]
        for line in source:
            if line.startswith("%"):
                comment_lines.append(line.rstrip("\r\n"))
                continue
            if line.strip():
                dimension_tokens = line.split()
                break
        else:
            raise ValueError("Missing MatrixMarket dimensions")
        if len(dimension_tokens) != 3 or any(not re.fullmatch(r"[0-9]+", x) for x in dimension_tokens):
            raise ValueError("Invalid dimension tokens")
        nfeatures, ncells, nnz = map(int, dimension_tokens)
        if not (0 < nfeatures <= CAPS["matrix_dimensions_each"] and 0 < ncells <= CAPS["matrix_dimensions_each"] and 0 <= nnz <= CAPS["matrix_coordinates"]):
            raise ValueError("Matrix dimensions exceed safety caps")
        if nfeatures != len(feature_rows) or ncells != len(barcodes):
            raise ValueError("Matrix/feature/barcode axes disagree")
        for line in source:
            if not line.strip() or line.startswith("%"):
                continue
            tokens = line.split()
            if len(tokens) != 3 or any(not re.fullmatch(r"[0-9]+", x) for x in tokens[:2]):
                raise ValueError("Malformed sparse coordinate")
            row, col = map(int, tokens[:2])
            if not (1 <= row <= nfeatures and 1 <= col <= ncells):
                raise ValueError("Out-of-bounds one-based coordinate")
            value = exact_value(tokens[2], field)
            nread += 1
            if nread > nnz:
                raise ValueError("More entries than declared")
            key = (col - 1) * nfeatures + row - 1
            keys.append(key)
            strictly_increasing = strictly_increasing and key > previous_key
            previous_key = key
            fractions += not isinstance(value, int)
            zeros += value == 0
            if value:
                nonzero_features.add(row)
            minimum = value if minimum is None or value < minimum else minimum
            maximum = value if maximum is None or value > maximum else maximum
            kind = types[row - 1]
            stored_by_type[kind] += 1
            sums[kind][col - 1] += value
    if nread != nnz:
        raise ValueError("Fewer entries than declared")
    if not strictly_increasing:
        ordered = sorted(keys)
        if any(a == b for a, b in zip(ordered, ordered[1:])):
            raise ValueError("Duplicate sparse coordinates; no summing or correction")
    export_metadata = [json.loads(line.split(":", 1)[1].strip())
                       for line in comment_lines if line.startswith("%metadata_json:")]
    if len(export_metadata) > 1:
        raise ValueError("Multiple exporter metadata comments; no implicit selection")
    result = {
        "header_as_published": header, "comment_lines_as_published": comment_lines,
        "export_metadata_as_published": export_metadata[0] if export_metadata else None,
        "dimensions": {"features": nfeatures, "barcodes": ncells, "stored_entries_declared": nnz},
        "stored_entries_observed": nread, "exact_numeric_before_float": True,
        "fractional_entries": fractions, "explicit_zero_entries": zeros,
        "minimum_stored_value_exact": str(minimum) if minimum is not None else None,
        "maximum_stored_value_exact": str(maximum) if maximum is not None else None,
        "finite_nonnegative_values": True,
        "duplicate_coordinates": 0,
        "strict_column_major_order": strictly_increasing,
        "features_with_any_positive_entry": len(nonzero_features),
        "all_zero_feature_rows": nfeatures - len(nonzero_features),
        "stored_entries_by_feature_type": dict(sorted(stored_by_type.items())),
        "sum_exact_by_feature_type": {kind: str(sum(values)) for kind, values in sums.items()},
        "zero_total_barcodes_by_feature_type": {kind: sum(value == 0 for value in values) for kind, values in sums.items()},
        "numeric_compatibility_is_not_original_UMI_provenance": True,
    }
    return result, sums


def fixed_metadata_join(barcodes: list[str], metadata_path: Path) -> tuple[dict, set[str]]:
    sample_name = FILENAME.removesuffix("_filtered_feature_bc_matrix.tar.gz")
    prefix = sample_name + "_"
    with metadata_path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        columns = reader.fieldnames
        if not columns or not {"sample_name", "barcode"}.issubset(columns):
            raise ValueError("Fixed metadata lacks exact join fields")
        selected = [r for r in reader if r["sample_name"] == sample_name]
    retained = set()
    for row in selected:
        if not row["barcode"].startswith(prefix):
            raise ValueError("Metadata barcode lacks fixed sample prefix")
        bare = row["barcode"][len(prefix):]
        if bare in retained:
            raise ValueError("Duplicate fixed metadata barcode")
        retained.add(bare)
    if not retained:
        raise ValueError("No retained metadata rows for explicit archive sample")
    matrix = set(barcodes)
    possible_margins = [col for col in columns if re.search(r"nCount|nFeature|RNA|SCT|UMI", col)]
    return {
        "rule": "Exact SDRF-linked archive basename; remove only the matching sample_name + '_' metadata prefix; no global barcode-only join",
        "matrix_barcodes": len(barcodes), "fixed_retained_metadata_rows": len(selected),
        "unique_retained_metadata_barcodes": len(retained),
        "matched_retained_barcodes": len(retained & matrix),
        "retained_barcodes_missing_from_matrix": len(retained - matrix),
        "matrix_barcodes_outside_fixed_retained_metadata": len(matrix - retained),
        "RNA_SCT_count_margin_columns_present": possible_margins,
        "RNA_SCT_margin_reproduction_available": bool(possible_margins),
        "retention_is_membership_not_an_expression_contrast": True,
    }, retained


def qualify(raw_dir: Path = RAW, work_dir: Path = WORK,
            output_path: Path = ROOT / "data/derived/psc-bach2-count-qualification.json") -> dict:
    plan = json.loads((raw_dir / "plan.json").read_text())
    for rel, pin in plan["preserve_prior_four_files"].items():
        if file_sha256(ROOT / rel) != pin["sha256"]:
            raise ValueError("Preceding metadata deliverable changed")
    inv = inventory(raw_dir, work_dir)
    decoded = prepare_members(work_dir / "outer.tar", inv["members"], work_dir)
    features, feature_summary = parse_features(ROOT / decoded["features"]["decoded_path"])
    barcodes = parse_barcodes(ROOT / decoded["barcodes"]["decoded_path"])
    matrix_summary, column_sums = parse_matrix(ROOT / decoded["matrix"]["decoded_path"], features, barcodes)
    metadata_dir = ROOT / "data/raw/psc-bach2-qualification"
    metadata_ledger = [json.loads(line) for line in (metadata_dir / "request-ledger.jsonl").read_text().splitlines()]
    metadata_entry = next(x for x in metadata_ledger if x["name"] == "scrna_meta.tsv")
    metadata_path = metadata_dir / metadata_entry["body_file"]
    if file_sha256(metadata_path) != metadata_entry["sha256"]:
        raise ValueError("Fixed retained metadata hash changed")
    join, retained = fixed_metadata_join(barcodes, metadata_path)
    if join["retained_barcodes_missing_from_matrix"]:
        raise ValueError("Retained metadata barcode missing from the selected matrix")
    # These per-cell totals are ignored artifacts. No per-gene values or state contrasts are computed.
    margin_rows = [{"barcode": barcode, "in_fixed_retained_metadata": barcode in retained,
                    **{kind: str(values[i]) for kind, values in column_sums.items()}}
                   for i, barcode in enumerate(barcodes)]
    json_write(work_dir / "single-sample-column-margins.json", margin_rows)
    receipt = json.loads((raw_dir / "acquisition-receipt.json").read_text())
    definition_attempts = [json.loads(p.read_text()) for p in sorted(raw_dir.glob("count-definition*-receipt.json"))]
    definition_public = [{"url": x["url"], "response_url": x.get("response_url"),
                          "http_statuses": [h["status"] for h in x["requests"]],
                          "body_bytes": x["body_bytes"], "sha256": x["sha256"],
                          "complete_response": x["complete_response"], "error": x["error"],
                          "retrieved_at_utc": x["completed_at_utc"]} for x in definition_attempts]
    metadata_bytes = sum(x["body_bytes"] for x in definition_attempts)
    if metadata_bytes > 1_000_000:
        raise ValueError("Additional metadata budget exceeded")
    result = {
        "scope": "One smallest filtered archive; structural, assay and fixed-barcode qualification only",
        "new_archive_body_bytes": receipt["body_bytes"],
        "additional_metadata_body_bytes": metadata_bytes,
        "additional_metadata_budget_bytes": 1_000_000,
        "official_format_definition_requests": definition_public,
        "ready_for_expression_analysis": False,
        "source": receipt,
        "archive_integrity": inv,
        "nested_member_integrity": decoded,
        "feature_axis": feature_summary,
        "matrix": matrix_summary,
        "retained_metadata_join": join,
        "fixed_metadata_sha256": metadata_entry["sha256"],
        "source_quantity_assessment": {
            "status": "Cell Ranger MEX-form RNA-feature count interface is supported; untouched original UMI versus later SCT/export lineage is not certified",
            "declared_export_metadata": matrix_summary["export_metadata_as_published"],
            "Methods_IDF_alignment_software_as_published": "Cellranger 3.0.2",
            "exporter_vs_Methods_version_disagreement_preserved": True,
            "Methods_IDF_alignment_reference_as_published": "GRCh38 human reference genome (release 93)",
            "actual_export_reference_manifest_in_archive": None,
            "normalized_or_SCT_layer_label_in_archive": None,
            "RNA_SCT_margins_in_fixed_metadata": "absent",
            "RNA_feature_rows_require_separation_from_ADT": True,
            "RNA_original_count_likelihood_gate_certified": False,
            "integer_values_or_filename_alone_are_not_sufficient": True,
            "no_normalization_or_historical_gene_filter_applied": True,
            "remaining_archives_can_check": ["all-eight source dimensions/types/exporter comments", "complete retained-barcode membership", "fixed-compartment measured feature coverage after separate authorization"],
            "remaining_archives_alone_cannot_supply": ["missing historical RNA/SCT margins", "an explicit original-RNA-count export record if the same source layout is repeated"],
        },
        "remaining_seven_archives": {"downloaded": False, "additional_published_bytes": 298_421_791,
                                      "authorization_received": False},
        "expression_effects_programs_or_model_fits": False,
        "prior_four_metadata_deliverable_hashes_unchanged": True,
    }
    json_write(output_path, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["acquire", "inventory", "qualify", "definition", "definition-current"])
    parser.add_argument("--output", type=Path, default=ROOT / "data/derived/psc-bach2-count-qualification.json")
    args = parser.parse_args()
    if args.command == "acquire":
        print(json.dumps(acquire(), sort_keys=True))
    elif args.command == "definition":
        print(json.dumps(acquire_definition(), sort_keys=True))
    elif args.command == "definition-current":
        print(json.dumps(acquire_definition(name="count-definition-current", url="https://www.10xgenomics.com/support/software/cell-ranger/6.1/analysis/outputs/cr-outputs-mex-matrices"), sort_keys=True))
    elif args.command == "inventory":
        print(json.dumps(inventory(), sort_keys=True))
    elif args.command == "qualify":
        result = qualify(output_path=args.output)
        print(json.dumps({"matrix": result["matrix"], "feature_axis": result["feature_axis"], "retained_metadata_join": result["retained_metadata_join"]}, sort_keys=True))


if __name__ == "__main__":
    main()
