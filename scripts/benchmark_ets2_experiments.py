"""Reproduce the prespecified ETS2 experimental summary benchmark.

The benchmark uses the publisher's Supplementary Tables S1 and S2 as the
primary source for the five editing/overexpression experiments.  The two MEK
inhibitor vectors are the released symbol/t-statistic rank files.  Those rank
files are only used directionally when a separate source audit supplies an
explicit, validated orientation and sign multiplier.

This is a summary-level reproduction.  It does not fit differential
expression models, treat genes as donor replicates, or infer a treatment
effect.  The default output directory is ``data/derived`` and the default
audit path is ``reports/ets2-experiment-benchmark.json``.  ``--output-directory``
and ``--audit-path`` make replay into an isolated directory possible without
changing the frozen inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "config/ets2-program-plan.json"
BENCHMARK_MANIFEST_PATH = ROOT / "config/ets2-benchmark-sources.json"
PROGRAM_MANIFEST_PATH = ROOT / "config/ets2-program-sources.json"
DEFAULT_OUTPUT_DIRECTORY = ROOT / "data/derived"
DEFAULT_AUDIT_PATH = ROOT / "reports/ets2-experiment-benchmark.json"

# The plan was frozen before this benchmark was implemented.  Keeping this
# value in the reader makes an accidental plan edit fail closed.
FROZEN_PLAN_SHA256 = "234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9"
EXCLUDED_ETS2_ID = "ENSG00000157557"
EXPECTED_SOURCE_SET_SIZE = 928
EXPECTED_EXPERIMENT_COMPARISONS = (
    ("g1", "g2"),
    ("g1", "chr21"),
    ("g1", "OE250"),
    ("g1", "OE500"),
    ("OE250", "OE500"),
    ("g1", "MEKi0.1"),
    ("g1", "MEKi0.5"),
)

EXPERIMENT_SPECS = {
    "g1": {"sheet": 1, "start_col": 3, "label": "ETS2 gRNA 1", "source_type": "publisher_table"},
    "g2": {"sheet": 1, "start_col": 9, "label": "ETS2 gRNA 2", "source_type": "publisher_table"},
    "chr21": {"sheet": 1, "start_col": 15, "label": "chr21q22 deletion", "source_type": "publisher_table"},
    "OE250": {"sheet": 2, "start_col": 3, "label": "ETS2 overexpression (250ng)", "source_type": "publisher_table"},
    "OE500": {"sheet": 2, "start_col": 9, "label": "ETS2 overexpression (500ng)", "source_type": "publisher_table"},
    "OE_combined": {"sheet": 2, "start_col": 15, "label": "ETS2 overexpression (250ng and 500ng, dose as covariate)", "source_type": "publisher_table"},
    "MEKi0.1": {"member": "RNA-seq/RNAseq_MEKi/PD_0.1_MEK.resSYMBOL.rank.csv", "label": "MEK inhibitor 0.1", "source_type": "released_rank_file"},
    "MEKi0.5": {"member": "RNA-seq/RNAseq_MEKi/PD_0.5_MEK.resSYMBOL.rank.csv", "label": "MEK inhibitor 0.5", "source_type": "released_rank_file"},
}

TABLE_FIELDS = ("EnsemblID", "geneID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
NUMERIC_TABLE_FIELDS = ("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
MEK_EXPERIMENTS = frozenset(("MEKi0.1", "MEKi0.5"))
_XLSX_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_XLSX_N = "{" + _XLSX_NS + "}"


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of *path* without loading it all at once."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_safe(value: Any) -> Any:
    if isinstance(value, float):
        return None if not math.isfinite(value) else value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _float_or_none(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _format_csv_value(value: Any) -> Any:
    """Keep missing values empty and floats stable enough for replay diffs."""

    if value is None:
        return ""
    if isinstance(value, float):
        return "" if not math.isfinite(value) else format(value, ".17g")
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (Mapping, list, tuple)):
        return json.dumps(_json_safe(value), sort_keys=True, separators=(",", ":"))
    return value


def _is_missing_identifier(value: Any) -> bool:
    return value is None or str(value).strip() in {"", "None", "NA", "NaN", "nan"}


def _relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _verify_manifest(manifest_path: Path, root: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_directory = root / str(manifest["source_directory"])
    verified: list[dict[str, Any]] = []
    for source in manifest.get("sources", []):
        relative_file = str(source["file"])
        path = source_directory / relative_file
        if not path.exists():
            raise FileNotFoundError(path)
        expected_bytes = int(source["bytes"])
        actual_bytes = path.stat().st_size
        if actual_bytes != expected_bytes:
            raise ValueError(f"Pinned source byte length changed: {relative_file}")
        expected_sha = str(source["sha256"])
        actual_sha = sha256_file(path)
        if actual_sha != expected_sha:
            raise ValueError(f"Pinned source SHA-256 changed: {relative_file}")
        verified.append({
            "id": str(source.get("id", relative_file)),
            "file": relative_file,
            "bytes": actual_bytes,
            "sha256": actual_sha,
            "url": source.get("url"),
        })
    archive_md5 = manifest.get("published_zip_md5")
    archive_md5_status = None
    if archive_md5:
        archive = source_directory / "ets2-public-release.zip"
        actual_md5 = md5_file(archive)
        if actual_md5 != str(archive_md5):
            raise ValueError("Publisher archive MD5 mismatch")
        archive_md5_status = {"expected": str(archive_md5), "actual": actual_md5, "verified": True}
    return {
        "path": _relative_or_absolute(manifest_path, root),
        "sha256": sha256_file(manifest_path),
        "source_directory": str(manifest.get("source_directory", "")),
        "sources": verified,
        "published_zip_md5": archive_md5_status,
    }


def verify_frozen_inputs(root: Path = ROOT, plan_path: Path = PLAN_PATH) -> dict[str, Any]:
    """Verify the frozen plan, its recorded input hashes, and both manifests."""

    plan_bytes = plan_path.read_bytes()
    plan_sha = hashlib.sha256(plan_bytes).hexdigest()
    if plan_sha != FROZEN_PLAN_SHA256:
        raise ValueError(f"Frozen ETS2 program plan changed: {plan_sha}")
    plan = json.loads(plan_bytes.decode("utf-8"))
    expected_comparisons = tuple(tuple(item) for item in plan["experimental_benchmark"]["planned_comparisons"])
    # The plan uses descriptive labels; this check ensures this script has not
    # silently broadened or reordered the prespecified comparison set.
    expected_labels = tuple((
        {"ETS2_g1": "g1", "ETS2_g2": "g2", "chr21_deletion": "chr21", "OE_250": "OE250", "OE_500": "OE500", "MEKi_0.1": "MEKi0.1", "MEKi_0.5": "MEKi0.5"}[left],
        {"ETS2_g1": "g1", "ETS2_g2": "g2", "chr21_deletion": "chr21", "OE_250": "OE250", "OE_500": "OE500", "MEKi_0.1": "MEKi0.1", "MEKi_0.5": "MEKi0.5"}[right],
    ) for left, right in expected_comparisons)
    if expected_labels != EXPECTED_EXPERIMENT_COMPARISONS:
        raise ValueError(f"Frozen experimental comparisons changed: {expected_labels}")

    input_hashes: dict[str, Any] = {}
    for relative_file, expected_sha in plan.get("input_sha256", {}).items():
        path = root / str(relative_file)
        if not path.exists():
            raise FileNotFoundError(path)
        actual_sha = sha256_file(path)
        input_hashes[str(relative_file)] = {
            "expected": str(expected_sha), "actual": actual_sha, "verified": actual_sha == str(expected_sha),
        }
        if actual_sha != str(expected_sha):
            raise ValueError(f"Frozen plan input changed: {relative_file}")

    benchmark_manifest = _verify_manifest(root / "config/ets2-benchmark-sources.json", root)
    program_manifest = _verify_manifest(root / "config/ets2-program-sources.json", root)
    return {
        "plan": {
            "file": _relative_or_absolute(plan_path, root),
            "sha256": plan_sha,
            "expected_sha256": FROZEN_PLAN_SHA256,
            "verified": True,
        },
        "plan_input_sha256": input_hashes,
        "source_manifests": {
            "config/ets2-benchmark-sources.json": benchmark_manifest,
            "config/ets2-program-sources.json": program_manifest,
        },
    }


def _column_number(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference)
    if not match:
        raise ValueError(f"Invalid worksheet cell reference: {reference}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - ord("A") + 1
    return result


def _read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(text.text or "" for text in item.iter(_XLSX_N + "t")) for item in root.findall(_XLSX_N + "si")]


def _read_worksheet_rows(archive: zipfile.ZipFile, sheet_number: int, shared_strings: Sequence[str]) -> list[dict[str, Any]]:
    """Read one worksheet only; callers deliberately never parse sheet 3."""

    root = ET.fromstring(archive.read(f"xl/worksheets/sheet{sheet_number}.xml"))
    sheet_data = root.find(_XLSX_N + "sheetData")
    if sheet_data is None:
        raise ValueError(f"Worksheet {sheet_number} has no sheetData")
    rows: list[dict[str, Any]] = []
    for row in sheet_data:
        values: dict[int, Any] = {}
        for cell in row:
            column = _column_number(str(cell.attrib["r"]))
            value_node = cell.find(_XLSX_N + "v")
            value = None if value_node is None else value_node.text
            cell_type = cell.attrib.get("t")
            if cell_type == "s" and value is not None:
                value = shared_strings[int(value)]
            elif cell_type == "inlineStr":
                value = "".join(text.text or "" for text in cell.iter(_XLSX_N + "t"))
            values[column] = value
        rows.append({"row_index": int(row.attrib.get("r", len(rows) + 1)), "cells": values})
    return rows


def _read_supplementary_tables(path: Path) -> tuple[dict[int, list[dict[str, Any]]], dict[str, Any]]:
    with zipfile.ZipFile(path) as archive:
        shared_strings = _read_shared_strings(archive)
        sheets = {1: _read_worksheet_rows(archive, 1, shared_strings), 2: _read_worksheet_rows(archive, 2, shared_strings)}
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet_titles = {
            int(item.attrib.get("sheetId", "0")): str(item.attrib.get("name", ""))
            for item in workbook.findall(_XLSX_N + "sheets/" + _XLSX_N + "sheet")
        }
        dimensions: dict[int, str | None] = {}
        for sheet_number in (1, 2):
            sheet_root = ET.fromstring(archive.read(f"xl/worksheets/sheet{sheet_number}.xml"))
            node = sheet_root.find(_XLSX_N + "dimension")
            dimensions[sheet_number] = node.attrib.get("ref") if node is not None else None
    return sheets, {
        "parsed_worksheets": [1, 2],
        "worksheet_titles": sheet_titles,
        "worksheet_dimensions": dimensions,
        "sheet3_oligo_sequences_read": False,
    }


def _parse_official_experiments(path: Path) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    sheets, parser_audit = _read_supplementary_tables(path)
    experiments: dict[str, dict[str, Any]] = {}
    audit: dict[str, dict[str, Any]] = {"parser": parser_audit}
    for experiment, spec in ((key, value) for key, value in EXPERIMENT_SPECS.items() if value["source_type"] == "publisher_table"):
        rows = sheets[int(spec["sheet"])]
        header = rows[1]["cells"]
        start = int(spec["start_col"])
        observed_header = (header.get(1), header.get(2)) + tuple(header.get(start + index) for index in range(6))
        if observed_header != TABLE_FIELDS:
            raise ValueError(f"Unexpected header for {experiment}: {observed_header}")
        records: list[dict[str, Any]] = []
        non_data_rows = 0
        non_data_row_details: list[dict[str, Any]] = []
        for row in rows[2:]:
            cells = row["cells"]
            ensembl_id = cells.get(1)
            if not isinstance(ensembl_id, str) or not ensembl_id.startswith("ENSG"):
                non_data_rows += 1
                note_text = " ".join(str(value).strip() for _, value in sorted(cells.items()) if value is not None).strip()
                non_data_row_details.append({
                    "row_index": int(row["row_index"]),
                    "reason": "publisher_note",
                    "text": note_text,
                })
                continue
            record: dict[str, Any] = {
                "experiment": experiment,
                "source_type": "publisher_table",
                "source_sheet": int(spec["sheet"]),
                "source_row": int(row["row_index"]),
                "ensembl_id": str(ensembl_id),
                "gene_id": cells.get(2),
            }
            for index, field in enumerate(NUMERIC_TABLE_FIELDS, start=start):
                raw = cells.get(index)
                record[f"{field}_raw"] = raw
                record[field] = _float_or_none(raw)
            records.append(record)
        counts = Counter(record["ensembl_id"] for record in records)
        duplicate_ids = {identifier: count for identifier, count in counts.items() if count > 1}
        gene_id_counts = Counter(
            str(record["gene_id"]).strip()
            for record in records
            if not _is_missing_identifier(record.get("gene_id"))
        )
        duplicate_gene_ids = {identifier: count for identifier, count in gene_id_counts.items() if count > 1}
        missing_gene_id_rows = sum(_is_missing_identifier(record.get("gene_id")) for record in records)
        literal_none_gene_id_rows = sum(str(record.get("gene_id")).strip() == "None" for record in records)
        unique_by_id = {
            record["ensembl_id"]: record
            for record in records
            if counts[record["ensembl_id"]] == 1
        }
        experiments[experiment] = {
            "experiment": experiment,
            "label": spec["label"],
            "source_type": "publisher_table",
            "source_sheet": int(spec["sheet"]),
            "records": records,
            "by_id": unique_by_id,
            "all_ids": [record["ensembl_id"] for record in records],
        }
        audit[experiment] = {
            "source_type": "publisher_table",
            "sheet": int(spec["sheet"]),
            "source_sheet": int(spec["sheet"]),
            "worksheet_title": parser_audit["worksheet_titles"].get(int(spec["sheet"])),
            "header_fields": list(TABLE_FIELDS),
            "label": spec["label"],
            "worksheet_rows": len(rows),
            "worksheet_row_elements": len(rows),
            "worksheet_max_row": parser_audit["worksheet_dimensions"].get(int(spec["sheet"])),
            "data_rows": len(records),
            "non_data_rows": non_data_rows,
            "non_data_row_details": non_data_row_details,
            "row_accounting": f"{len(rows)} stored row elements = 2 title/header rows + {len(records)} Ensembl data rows + {non_data_rows} publisher note rows; worksheet dimension is {parser_audit['worksheet_dimensions'].get(int(spec['sheet']))}",
            "publisher_contrast_description": (
                "Results shown with respect to expression in unedited cells"
                if int(spec["sheet"]) == 1 else
                "Results shown with respect to expression in cells transfected with equivalent amount of mRNA encoding reverse complement of ETS2"
            ),
            "unique_identifier_count": len(counts),
            "duplicate_identifier_rows": sum(count - 1 for count in duplicate_ids.values()),
            "duplicate_identifier_values": duplicate_ids,
            "missing_identifier_rows": 0,
            "unique_gene_id_count": len(gene_id_counts),
            "duplicate_gene_id_rows": sum(count - 1 for count in duplicate_gene_ids.values()),
            "duplicate_gene_id_values": duplicate_gene_ids,
            "missing_gene_id_rows": missing_gene_id_rows,
            "literal_none_gene_id_rows": literal_none_gene_id_rows,
            "numeric_counts": {field: sum(record[field] is not None for record in records) for field in NUMERIC_TABLE_FIELDS},
            "na_counts": {field: sum(record[field] is None for record in records) for field in NUMERIC_TABLE_FIELDS},
            "zero_counts": {field: sum(record[field] == 0 for record in records if record[field] is not None) for field in NUMERIC_TABLE_FIELDS},
        }
    return experiments, audit


class OneToOneXref:
    """The exact S1/S2 symbol cross-reference and its excluded ambiguities.

    A symbol is usable when it maps to one Ensembl ID.  Multiple symbols that
    happen to map to the same ID are retained as aliases in the map but are
    deduplicated by :func:`map_rank_rows` so a rank file cannot contribute that
    Ensembl gene twice.  A symbol mapping to multiple Ensembl IDs is excluded.
    """

    __slots__ = ("mapping", "ambiguous", "symbol_to_ids", "id_to_symbols", "reverse_alias_targets")

    def __init__(self, mapping: dict[str, str], ambiguous: dict[str, tuple[str, ...]],
                 symbol_to_ids: dict[str, tuple[str, ...]], id_to_symbols: dict[str, tuple[str, ...]],
                 reverse_alias_targets: dict[str, tuple[str, ...]]) -> None:
        self.mapping = mapping
        self.ambiguous = ambiguous
        self.symbol_to_ids = symbol_to_ids
        self.id_to_symbols = id_to_symbols
        self.reverse_alias_targets = reverse_alias_targets


def build_one_to_one_xref(records: Iterable[Any]) -> OneToOneXref:
    """Build a strict two-way symbol/Ensembl map from S1/S2 rows.

    ``records`` accepts dictionaries with ``ensembl_id``/``EnsemblID`` and
    ``gene_id``/``geneID``/``symbol`` keys, or two-tuples of ``(id, symbol)``.
    A symbol is usable only when it names exactly one Ensembl ID.  Multiple
    symbols for one target are recorded as aliases and deduplicated when rank
    rows are mapped; no rank row is silently counted twice.  A symbol mapping
    to multiple IDs is ambiguous and excluded.
    """

    symbol_to_ids: defaultdict[str, set[str]] = defaultdict(set)
    id_to_symbols: defaultdict[str, set[str]] = defaultdict(set)
    for item in records:
        if isinstance(item, Mapping):
            ensembl_id = item.get("ensembl_id", item.get("EnsemblID"))
            symbol = item.get("gene_id", item.get("geneID", item.get("symbol")))
        else:
            try:
                ensembl_id, symbol = item
            except (TypeError, ValueError) as exc:
                raise ValueError("xref records must be mappings or (EnsemblID, symbol) pairs") from exc
        if _is_missing_identifier(ensembl_id) or _is_missing_identifier(symbol):
            continue
        ensembl_id = str(ensembl_id)
        symbol = str(symbol)
        symbol_to_ids[symbol].add(ensembl_id)
        id_to_symbols[ensembl_id].add(symbol)
    usable: dict[str, str] = {}
    ambiguous: dict[str, tuple[str, ...]] = {}
    for symbol, ids in sorted(symbol_to_ids.items()):
        if len(ids) == 1:
            usable[symbol] = next(iter(ids))
        else:
            ambiguous[symbol] = tuple(sorted(ids))
    reverse_alias_targets = {
        identifier: tuple(sorted(symbols))
        for identifier, symbols in id_to_symbols.items()
        if len(symbols) > 1
    }
    return OneToOneXref(
        mapping=usable,
        ambiguous=ambiguous,
        symbol_to_ids={key: tuple(sorted(value)) for key, value in symbol_to_ids.items()},
        id_to_symbols={key: tuple(sorted(value)) for key, value in id_to_symbols.items()},
        reverse_alias_targets=reverse_alias_targets,
    )


def _load_inventory(root: Path) -> dict[str, dict[str, Any]]:
    path = root / "data/derived/ets2-benchmark-inventory.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as stream:
        return {str(row["member"]): row for row in csv.DictReader(stream)}


def _archive_member(root: Path, relative_member: str, inventory: Mapping[str, Mapping[str, Any]]) -> tuple[bytes, dict[str, Any]]:
    archive_path = root / "data/raw/ets2-benchmark/ets2-public-release.zip"
    with zipfile.ZipFile(archive_path) as archive:
        names = archive.namelist()
        if not names:
            raise ValueError("Publisher archive is empty")
        base = names[0]
        member_name = base + relative_member
        if member_name not in names:
            raise FileNotFoundError(f"Missing archive member: {relative_member}")
        payload = archive.read(member_name)
    actual_sha = hashlib.sha256(payload).hexdigest()
    inventory_row = inventory.get(relative_member)
    if inventory_row and actual_sha != str(inventory_row["sha256"]):
        raise ValueError(f"Archive member SHA-256 changed: {relative_member}")
    return payload, {
        "member": relative_member,
        "bytes": len(payload),
        "sha256": actual_sha,
        "inventory_sha256_verified": bool(inventory_row),
    }


def _read_rank_rows(payload: bytes, experiment: str, member: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    text = payload.decode("utf-8")
    for row_index, row in enumerate(csv.reader(io.StringIO(text)), start=1):
        identifier = row[0] if row else None
        raw_value = row[1] if len(row) > 1 else None
        rows.append({
            "experiment": experiment,
            "source_type": "released_rank_file",
            "source_row": row_index,
            "source_identifier": identifier,
            "raw_value": raw_value,
            "value": _float_or_none(raw_value),
            "member": member,
        })
    return rows


def map_rank_rows(rows: Sequence[Mapping[str, Any]], xref: OneToOneXref) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Map rank-file symbols once, retaining every input row's disposition."""

    identifiers = [row.get("source_identifier") for row in rows]
    identifier_counts = Counter(identifier for identifier in identifiers if not _is_missing_identifier(identifier))
    seen_targets: set[str] = set()
    mapped: dict[str, dict[str, Any]] = {}
    audited: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()
    for row in rows:
        symbol_value = row.get("source_identifier")
        symbol = "" if symbol_value is None else str(symbol_value)
        status = "mapped_one_to_one"
        ensembl_id: str | None = None
        if _is_missing_identifier(symbol_value):
            status = "missing_source_identifier"
        elif identifier_counts.get(symbol_value, 0) > 1:
            status = "duplicate_source_identifier"
        elif symbol in xref.ambiguous:
            status = "ambiguous_xref"
        elif symbol not in xref.mapping:
            status = "missing_xref"
        else:
            ensembl_id = xref.mapping[symbol]
            if ensembl_id in seen_targets:
                status = "duplicate_target_ensembl"
                ensembl_id = None
            else:
                seen_targets.add(ensembl_id)
                mapped[ensembl_id] = dict(row, ensembl_id=ensembl_id, mapping_status=status)
        status_counts[status] += 1
        audited.append(dict(row, ensembl_id=ensembl_id, mapping_status=status))
    value_rows = [row for row in rows if row.get("value") is not None]
    accounting = {
        "input_rows": len(rows),
        "unique_source_identifier_count": len(identifier_counts),
        "duplicate_source_identifier_rows": sum(count - 1 for count in identifier_counts.values() if count > 1),
        "missing_source_identifier_rows": sum(_is_missing_identifier(value) for value in identifiers),
        "numeric_value_rows": len(value_rows),
        "na_value_rows": len(rows) - len(value_rows),
        "zero_value_rows": sum(row.get("value") == 0 for row in value_rows),
        "mapped_rows": len(mapped),
        "mapped_unique_ensembl_count": len(mapped),
        "mapping_status_counts": dict(sorted(status_counts.items())),
        "xref_ambiguous_symbol_count": len(xref.ambiguous),
    }
    return mapped, audited, accounting


def _read_source_set(root: Path, inventory: Mapping[str, Mapping[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    member = "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv"
    payload, member_audit = _archive_member(root, member, inventory)
    source_rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    values: list[dict[str, Any]] = []
    for row_index, row in enumerate(source_rows, start=2):
        value = row.get("ETS2g1_DN")
        if _is_missing_identifier(value):
            continue
        values.append({
            "source_member": member,
            "column": "ETS2g1_DN",
            "row_index": row_index,
            "ensembl_id": str(value),
            "is_excluded": str(value) == EXCLUDED_ETS2_ID,
            "exclusion_reason": "exclude_target_ETS2" if str(value) == EXCLUDED_ETS2_ID else "",
        })
    counts = Counter(row["ensembl_id"] for row in values)
    duplicate_ids = {key: count for key, count in counts.items() if count > 1}
    if len(values) != EXPECTED_SOURCE_SET_SIZE:
        raise ValueError(f"Expected {EXPECTED_SOURCE_SET_SIZE} ETS2g1_DN rows, observed {len(values)}")
    if duplicate_ids:
        raise ValueError(f"ETS2g1_DN contains duplicate Ensembl IDs: {duplicate_ids}")
    if sum(row["is_excluded"] for row in values) != 1:
        raise ValueError("ETS2g1_DN does not contain exactly one ETS2 target row")
    return values, {
        "source_member": member,
        "column": "ETS2g1_DN",
        "source_rows": len(source_rows),
        "source_set_n": len(values),
        "unique_n": len(counts),
        "excluded_gene_id": EXCLUDED_ETS2_ID,
        "excluded_n": sum(row["is_excluded"] for row in values),
        "analysis_set_n": sum(not row["is_excluded"] for row in values),
        "member": member_audit,
    }


def _rank_values_for_experiment(experiment: Mapping[str, Any], orientation: Mapping[str, Any]) -> dict[str, float]:
    result: dict[str, float] = {}
    multiplier = orientation.get("sign_multiplier")
    if multiplier not in (-1, 1):
        return result
    for ensembl_id, row in experiment["by_id"].items():
        value = row.get("value")
        if value is not None:
            result[ensembl_id] = float(value) * int(multiplier)
    return result


def default_mek_orientation() -> dict[str, dict[str, Any]]:
    evidence = (
        "The source audit found that the released MEKi R script defines MEK-minus-control "
        "contrasts but only reads pre-existing PD_* rank files; the matching export is "
        "commented or mismatched, and Source Data 5i provides points/EGA rather than a "
        "MEK differential-expression table. No validated sign multiplier is available."
    )
    return {
        experiment: {
            "status": "unresolved",
            "sign_multiplier": None,
            "evidence": evidence,
            "directional_metrics_available": False,
        }
        for experiment in sorted(MEK_EXPERIMENTS)
    }


def _extract_orientation_payload(payload: Mapping[str, Any], source_path: Path | None) -> dict[str, dict[str, Any]]:
    """Accept only explicit validated orientation records; otherwise fail closed."""

    block: Any = None
    for key in ("mek_rank_orientation", "mek_orientation", "rank_orientation"):
        if key in payload:
            block = payload[key]
            break
    if isinstance(block, Mapping) and isinstance(block.get("experiments"), Mapping):
        block = block["experiments"]
    if not isinstance(block, Mapping):
        raise ValueError("Orientation file lacks an explicit mek_rank_orientation.experiments object")
    result = default_mek_orientation()
    for experiment in MEK_EXPERIMENTS:
        record = block.get(experiment)
        if not isinstance(record, Mapping):
            continue
        status = str(record.get("status", "unresolved"))
        multiplier = record.get("sign_multiplier", record.get("orientation_multiplier", record.get("multiplier")))
        evidence = str(record.get("evidence", record.get("basis", ""))).strip()
        validated = bool(record.get("validated", False)) or status.lower() in {"validated", "confirmed", "numeric_agreement"}
        try:
            multiplier_int = int(multiplier) if multiplier is not None else None
        except (TypeError, ValueError):
            multiplier_int = None
        if validated and multiplier_int in (-1, 1) and evidence:
            result[experiment] = {
                "status": "validated",
                "sign_multiplier": multiplier_int,
                "evidence": evidence,
                "directional_metrics_available": True,
                "source_path": str(source_path) if source_path else None,
            }
    return result


def load_mek_orientation(root: Path = ROOT, orientation_path: Path | None = None) -> tuple[dict[str, dict[str, Any]], str | None]:
    candidates = [orientation_path] if orientation_path else [
        root / "reports/ets2-source-audit.json",
        root / "data/derived/ets2-source-orientation.json",
    ]
    for candidate in candidates:
        if candidate is None or not candidate.exists():
            continue
        payload = json.loads(candidate.read_text(encoding="utf-8"))
        try:
            return _extract_orientation_payload(payload, candidate), str(candidate)
        except ValueError:
            if orientation_path is not None:
                raise
    return default_mek_orientation(), None


def _rankdata(values: Sequence[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    position = 0
    while position < len(ordered):
        end = position + 1
        while end < len(ordered) and ordered[end][1] == ordered[position][1]:
            end += 1
        rank = (position + 1 + end) / 2.0
        for index in range(position, end):
            ranks[ordered[index][0]] = rank
        position = end
    return ranks


def spearman_signed(left: Sequence[Any], right: Sequence[Any]) -> float | None:
    """Calculate Spearman rho with average ranks and constant-vector guards."""

    pairs = [(float(a), float(b)) for a, b in zip(left, right) if _finite(a) and _finite(b)]
    if len(pairs) < 2:
        return None
    left_values = [pair[0] for pair in pairs]
    right_values = [pair[1] for pair in pairs]
    left_ranks = _rankdata(left_values)
    right_ranks = _rankdata(right_values)
    left_mean = sum(left_ranks) / len(left_ranks)
    right_mean = sum(right_ranks) / len(right_ranks)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left_ranks, right_ranks))
    left_ss = sum((a - left_mean) ** 2 for a in left_ranks)
    right_ss = sum((b - right_mean) ** 2 for b in right_ranks)
    denominator = math.sqrt(left_ss * right_ss)
    if denominator == 0:
        return None
    return numerator / denominator


def sign_concordance(left: Sequence[Any], right: Sequence[Any]) -> dict[str, Any]:
    """Summarize signed agreement, excluding any pair containing a zero.

    Zero is retained in the accounting but is not called positive, negative,
    concordant, or discordant.  This prevents a zero statistic from becoming
    a spurious directional match.
    """

    counts = {"finite_pairs": 0, "informative_pairs": 0, "same_sign_n": 0, "opposite_sign_n": 0,
              "left_zero_n": 0, "right_zero_n": 0, "both_zero_n": 0, "one_zero_n": 0}
    for left_value, right_value in zip(left, right):
        if not (_finite(left_value) and _finite(right_value)):
            continue
        counts["finite_pairs"] += 1
        left_number, right_number = float(left_value), float(right_value)
        if left_number == 0 or right_number == 0:
            if left_number == 0:
                counts["left_zero_n"] += 1
            if right_number == 0:
                counts["right_zero_n"] += 1
            if left_number == 0 and right_number == 0:
                counts["both_zero_n"] += 1
            else:
                counts["one_zero_n"] += 1
            continue
        counts["informative_pairs"] += 1
        if (left_number < 0) == (right_number < 0):
            counts["same_sign_n"] += 1
        else:
            counts["opposite_sign_n"] += 1
    counts["fraction"] = (
        counts["same_sign_n"] / counts["informative_pairs"] if counts["informative_pairs"] else None
    )
    return counts


def _orientation_for(experiment: str, orientations: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    if experiment not in MEK_EXPERIMENTS:
        return {"status": "publisher_table", "sign_multiplier": 1, "directional_metrics_available": True}
    return dict(orientations.get(experiment, default_mek_orientation()[experiment]))


def _used_t(experiment: str, value: Any, orientations: Mapping[str, Mapping[str, Any]]) -> float | None:
    number = _float_or_none(value)
    if number is None:
        return None
    orientation = _orientation_for(experiment, orientations)
    if not orientation.get("directional_metrics_available"):
        return None
    return number * int(orientation.get("sign_multiplier", 1))


def build_comparison_row(left_name: str, right_name: str, experiments: Mapping[str, Mapping[str, Any]], orientations: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    left = experiments[left_name]
    right = experiments[right_name]
    shared = sorted(set(left["by_id"]) & set(right["by_id"]))
    left_orientation = _orientation_for(left_name, orientations)
    right_orientation = _orientation_for(right_name, orientations)
    directional = bool(left_orientation.get("directional_metrics_available") and right_orientation.get("directional_metrics_available"))
    left_t_raw = [left["by_id"][identifier].get("t") for identifier in shared]
    right_t_raw = []
    for identifier in shared:
        value = right["by_id"][identifier].get("t")
        right_t_raw.append(value if value is not None else right["by_id"][identifier].get("value"))
    left_t_used = [_used_t(left_name, value, orientations) for value in left_t_raw]
    right_t_used = [_used_t(right_name, value, orientations) for value in right_t_raw]
    t_pairs = [(a, b) for a, b in zip(left_t_raw, right_t_raw) if _finite(a) and _finite(b)]
    used_t_pairs = [
        (a, b) for a, b in zip(left_t_used, right_t_used)
        if a is not None and b is not None
    ]
    logfc_left = [left["by_id"][identifier].get("logFC") for identifier in shared]
    logfc_right = [right["by_id"][identifier].get("logFC") for identifier in shared]
    logfc_pairs = [(a, b) for a, b in zip(logfc_left, logfc_right) if _finite(a) and _finite(b)]
    t_agreement = sign_concordance(
        [a for a, b in zip(left_t_used, right_t_used) if a is not None and b is not None],
        [b for a, b in zip(left_t_used, right_t_used) if a is not None and b is not None],
    ) if directional else None
    logfc_agreement = sign_concordance(*zip(*logfc_pairs)) if logfc_pairs and directional else None
    row: dict[str, Any] = {
        "comparison": f"{left_name}-v-{right_name}",
        "left_experiment": left_name,
        "right_experiment": right_name,
        "n_shared_unique": len(shared),
        "n_t": len(t_pairs),
        "n_logFC": len(logfc_pairs),
        "spearman_signed_t": spearman_signed(*zip(*used_t_pairs)) if used_t_pairs and directional else None,
        "spearman_logFC": spearman_signed(*zip(*logfc_pairs)) if logfc_pairs and directional else None,
        "directional_status": "available" if directional else "unavailable_unresolved_orientation",
        "orientation_multiplier_left": left_orientation.get("sign_multiplier"),
        "orientation_multiplier_right": right_orientation.get("sign_multiplier"),
        "directional_note": "Both vectors have validated orientation" if directional else "MEKi rank-file orientation is unresolved; overlap counts are retained and directional metrics are unavailable",
    }
    if t_agreement is None:
        row.update({
            "direction_agreement_t": None, "direction_agreement_t_n": None,
            "direction_agreement_t_same_n": None, "direction_agreement_t_opposite_n": None,
            "direction_agreement_t_left_zero_n": None, "direction_agreement_t_right_zero_n": None,
            "direction_agreement_t_both_zero_n": None, "direction_agreement_t_one_zero_n": None,
        })
    else:
        row.update({
            "direction_agreement_t": t_agreement["fraction"], "direction_agreement_t_n": t_agreement["informative_pairs"],
            "direction_agreement_t_same_n": t_agreement["same_sign_n"], "direction_agreement_t_opposite_n": t_agreement["opposite_sign_n"],
            "direction_agreement_t_left_zero_n": t_agreement["left_zero_n"], "direction_agreement_t_right_zero_n": t_agreement["right_zero_n"],
            "direction_agreement_t_both_zero_n": t_agreement["both_zero_n"], "direction_agreement_t_one_zero_n": t_agreement["one_zero_n"],
        })
    if logfc_agreement is None:
        row.update({"direction_agreement_logFC": None, "direction_agreement_logFC_n": None,
                    "direction_agreement_logFC_same_n": None, "direction_agreement_logFC_opposite_n": None,
                    "direction_agreement_logFC_left_zero_n": None, "direction_agreement_logFC_right_zero_n": None,
                    "direction_agreement_logFC_both_zero_n": None, "direction_agreement_logFC_one_zero_n": None})
    else:
        row.update({"direction_agreement_logFC": logfc_agreement["fraction"], "direction_agreement_logFC_n": logfc_agreement["informative_pairs"],
                    "direction_agreement_logFC_same_n": logfc_agreement["same_sign_n"], "direction_agreement_logFC_opposite_n": logfc_agreement["opposite_sign_n"],
                    "direction_agreement_logFC_left_zero_n": logfc_agreement["left_zero_n"], "direction_agreement_logFC_right_zero_n": logfc_agreement["right_zero_n"],
                    "direction_agreement_logFC_both_zero_n": logfc_agreement["both_zero_n"], "direction_agreement_logFC_one_zero_n": logfc_agreement["one_zero_n"]})
    return row


def build_program_summary(source_set: Sequence[Mapping[str, Any]], experiments: Mapping[str, Mapping[str, Any]], orientations: Mapping[str, Mapping[str, Any]], source_set_audit: Mapping[str, Any]) -> list[dict[str, Any]]:
    source_ids = [str(row["ensembl_id"]) for row in source_set if not row["is_excluded"]]
    output: list[dict[str, Any]] = []
    for experiment in ("g1", "g2", "chr21", "OE250", "OE500", "MEKi0.1", "MEKi0.5"):
        info = experiments[experiment]
        orientation = _orientation_for(experiment, orientations)
        raw_values: list[float] = []
        directional_values: list[float] = []
        logfc_values: list[float] = []
        for identifier in source_ids:
            record = info["by_id"].get(identifier)
            if not record:
                continue
            raw_value = record.get("t")
            if raw_value is None:
                raw_value = record.get("value")
            if _finite(raw_value):
                raw_values.append(float(raw_value))
            used = _used_t(experiment, raw_value, orientations)
            if used is not None:
                directional_values.append(used)
            logfc = _float_or_none(record.get("logFC"))
            if logfc is not None:
                logfc_values.append(logfc)
        negative_n = sum(value < 0 for value in directional_values)
        positive_n = sum(value > 0 for value in directional_values)
        zero_n = sum(value == 0 for value in directional_values)
        relationship = {
            "g1": "circular_source_reproduction",
            "g2": "related_guide_companion_same_study",
            "chr21": "enhancer_deletion_companion_same_study",
            "OE250": "overexpression_companion_same_study",
            "OE500": "overexpression_companion_same_study",
            "MEKi0.1": "drug_response_companion_same_study",
            "MEKi0.5": "drug_response_companion_same_study",
        }[experiment]
        output.append({
            "experiment": experiment,
            "source_program": "ETS2g1_DN",
            "source_set_n": int(source_set_audit["source_set_n"]),
            "excluded_ets2_n": int(source_set_audit["excluded_n"]),
            "analysis_set_n": len(source_ids),
            "numeric_available_n": len(raw_values),
            "numeric_available_fraction": len(raw_values) / len(source_ids) if source_ids else None,
            "directional_available_n": len(directional_values) if orientation.get("directional_metrics_available") else None,
            "directional_available_fraction": len(directional_values) / len(source_ids) if orientation.get("directional_metrics_available") and source_ids else None,
            # The stored statistic is still a numeric summary when its sign
            # cannot be certified.  Directional fractions remain unavailable
            # until an audited orientation is supplied.
            "median_t": median(raw_values) if raw_values else None,
            "median_logFC": median(logfc_values) if logfc_values else None,
            "negative_n": negative_n if directional_values else None,
            "positive_n": positive_n if directional_values else None,
            "zero_n": zero_n if directional_values else None,
            "fraction_negative": negative_n / len(directional_values) if directional_values else None,
            "fraction_positive": positive_n / len(directional_values) if directional_values else None,
            "directional_status": "available" if orientation.get("directional_metrics_available") else "unavailable_unresolved_orientation",
            "orientation_multiplier": orientation.get("sign_multiplier"),
            "relationship_to_source_program": relationship,
            "independence_note": "Gene-set membership is source-defined; genes are not donor replicates. g1 is circular source reproduction and same-study companions are not independent validation.",
        })
    return output


def _official_statistics_rows(experiments: Mapping[str, Mapping[str, Any]], root: Path, supplementary_path: Path) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    source_file = _relative_or_absolute(supplementary_path, root)
    for experiment in ("g1", "g2", "chr21", "OE250", "OE500", "OE_combined"):
        for record in experiments[experiment]["records"]:
            result.append({
                "experiment": experiment,
                "source_type": "publisher_table",
                "source_file": source_file,
                "source_member": "",
                "source_row": record["source_row"],
                "source_identifier": record["ensembl_id"],
                "ensembl_id": record["ensembl_id"],
                "gene_id": record.get("gene_id"),
                "mapping_status": "official_exact_ensembl_id",
                "logFC_raw": record.get("logFC_raw"),
                "logFC": record.get("logFC"),
                "t_raw": record.get("t_raw"),
                "t_numeric": record.get("t"),
                "t_for_directional": record.get("t"),
                "orientation_multiplier": 1,
            })
    return result


def _rank_statistics_rows(audited_rows: Sequence[Mapping[str, Any]], root: Path, orientations: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in audited_rows:
        experiment = str(row["experiment"])
        orientation = _orientation_for(experiment, orientations)
        numeric = row.get("value")
        directional_value = _used_t(experiment, numeric, orientations)
        result.append({
            "experiment": experiment,
            "source_type": "released_rank_file",
            "source_file": "data/raw/ets2-benchmark/ets2-public-release.zip",
            "source_member": row.get("member"),
            "source_row": row.get("source_row"),
            "source_identifier": row.get("source_identifier"),
            "ensembl_id": row.get("ensembl_id"),
            "gene_id": row.get("source_identifier"),
            "mapping_status": row.get("mapping_status"),
            "logFC_raw": "",
            "logFC": None,
            "t_raw": row.get("raw_value"),
            "t_numeric": numeric,
            "t_for_directional": directional_value,
            "orientation_multiplier": orientation.get("sign_multiplier"),
        })
    return result


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field)) for field in fieldnames})
    return sha256_file(path)


def run_benchmark(
    *,
    root: Path = ROOT,
    plan_path: Path = PLAN_PATH,
    output_directory: Path | None = None,
    audit_path: Path | None = None,
    mek_orientation_path: Path | None = None,
) -> dict[str, Any]:
    """Run and write the frozen summary benchmark."""

    root = root.resolve()
    plan_path = plan_path.resolve()
    verified = verify_frozen_inputs(root, plan_path)
    output_directory = (output_directory or root / "data/derived").resolve()
    audit_path = (audit_path or root / "reports/ets2-experiment-benchmark.json").resolve()
    orientation, orientation_path_used = load_mek_orientation(root, mek_orientation_path)

    supplementary_path = root / "data/raw/ets2-benchmark/ets2-supplementary-tables.xlsx"
    official, official_audit = _parse_official_experiments(supplementary_path)
    inventory = _load_inventory(root)
    source_set, source_set_audit = _read_source_set(root, inventory)
    xref_records = []
    for experiment in ("g1", "g2", "chr21", "OE250", "OE500", "OE_combined"):
        xref_records.extend({"ensembl_id": row["ensembl_id"], "gene_id": row.get("gene_id")} for row in official[experiment]["records"])
    xref = build_one_to_one_xref(xref_records)

    rank_audit: dict[str, Any] = {}
    for experiment in sorted(MEK_EXPERIMENTS):
        member = str(EXPERIMENT_SPECS[experiment]["member"])
        payload, member_audit = _archive_member(root, member, inventory)
        rank_rows = _read_rank_rows(payload, experiment, member)
        mapped, audited_rows, accounting = map_rank_rows(rank_rows, xref)
        official[experiment] = {
            "experiment": experiment,
            "label": EXPERIMENT_SPECS[experiment]["label"],
            "source_type": "released_rank_file",
            "records": audited_rows,
            "by_id": mapped,
            "all_ids": [row.get("ensembl_id") for row in audited_rows if row.get("ensembl_id")],
        }
        rank_audit[experiment] = {**accounting, "source_type": "released_rank_file", "member": member_audit,
                                  "label": EXPERIMENT_SPECS[experiment]["label"], "orientation": orientation[experiment]}

    comparisons = [build_comparison_row(left, right, official, orientation) for left, right in EXPECTED_EXPERIMENT_COMPARISONS]
    program_summary = build_program_summary(source_set, official, orientation, source_set_audit)
    accounting_rows: list[dict[str, Any]] = []
    for experiment in ("g1", "g2", "chr21", "OE250", "OE500", "OE_combined"):
        row = dict(official_audit[experiment])
        row["experiment"] = experiment
        row["source_file"] = _relative_or_absolute(supplementary_path, root)
        row["source_sheet"] = row.get("sheet")
        row["orientation_status"] = "publisher_table"
        row["mapping_status_counts"] = {"official_exact_ensembl_id": row["data_rows"]}
        accounting_rows.append(row)
    for experiment in sorted(MEK_EXPERIMENTS):
        row = dict(rank_audit[experiment])
        row["experiment"] = experiment
        row["source_file"] = "data/raw/ets2-benchmark/ets2-public-release.zip"
        row["source_member"] = row["member"].get("member") if isinstance(row.get("member"), Mapping) else None
        row["orientation_status"] = orientation[experiment]["status"]
        accounting_rows.append(row)

    comparison_fields = list(comparisons[0].keys()) if comparisons else []
    program_fields = list(program_summary[0].keys()) if program_summary else []
    accounting_fields = [
        "experiment", "source_type", "source_file", "source_sheet", "worksheet_title", "header_fields", "publisher_contrast_description", "source_member", "worksheet_rows", "worksheet_row_elements", "worksheet_max_row", "row_accounting", "data_rows", "non_data_rows", "non_data_row_details",
        "input_rows", "unique_identifier_count", "unique_source_identifier_count", "unique_gene_id_count", "duplicate_identifier_rows", "duplicate_identifier_values", "duplicate_gene_id_rows", "duplicate_gene_id_values", "duplicate_source_identifier_rows",
        "missing_identifier_rows", "missing_gene_id_rows", "literal_none_gene_id_rows", "missing_source_identifier_rows", "numeric_counts", "na_counts", "zero_counts", "numeric_value_rows", "na_value_rows", "zero_value_rows",
        "mapped_rows", "mapped_unique_ensembl_count", "mapping_status_counts", "xref_ambiguous_symbol_count", "orientation_status",
    ]
    source_set_fields = ["source_member", "column", "row_index", "ensembl_id", "is_excluded", "exclusion_reason"]
    output_records: dict[str, Sequence[Mapping[str, Any]]] = {
        "ets2-experiment-comparisons.csv": comparisons,
        "ets2-experiment-program-summary.csv": program_summary,
        "ets2-experiment-input-accounting.csv": accounting_rows,
        "ets2-experiment-source-set.csv": source_set,
    }
    output_fieldnames = {
        "ets2-experiment-comparisons.csv": comparison_fields,
        "ets2-experiment-program-summary.csv": program_fields,
        "ets2-experiment-input-accounting.csv": accounting_fields,
        "ets2-experiment-source-set.csv": source_set_fields,
    }
    outputs: dict[str, str] = {}
    for filename, rows in output_records.items():
        outputs[_relative_or_absolute(output_directory / filename, root)] = _write_csv(output_directory / filename, rows, output_fieldnames[filename])

    source_audit_path = root / "reports/ets2-source-audit.json"
    source_audit_reference = None
    if source_audit_path.exists():
        source_audit_payload = json.loads(source_audit_path.read_text(encoding="utf-8"))
        source_audit_reference = {
            "path": _relative_or_absolute(source_audit_path, root),
            "sha256": sha256_file(source_audit_path),
            "used_for_orientation": False,
            "derived_output_sha256": source_audit_payload.get("derived_output_sha256", {}),
            "orientation_finding": (
                "The source audit reports that the MEKi R script defines MEK-minus-control "
                "contrasts but only reads pre-existing PD_* rank files; its matching export "
                "is commented or mismatched, so the saved rank-file sign cannot be certified."
            ),
        }
    audit = {
        "analysis_id": "psc-ets2-experimental-summary-benchmark-v1",
        "status": "complete_with_unresolved_meki_orientation" if any(item["directional_status"] != "available" for item in comparisons if "MEKi" in item["comparison"]) else "complete",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "plan_verification": verified["plan"],
        "plan_input_sha256": verified["plan_input_sha256"],
        "source_manifests": verified["source_manifests"],
        "parsed_supplementary_tables": official_audit["parser"],
        "source_set": source_set_audit,
        "identity_mapping": {
            "rule": "A rank symbol is usable only when it maps to one exact S1/S2 EnsemblID; symbols mapping to multiple IDs are excluded. Multiple symbols mapping to one ID are retained as aliases in the cross-reference but duplicate target IDs are counted once per rank file and the excluded rows are recorded.",
            "xref_scope": "All non-missing geneID-to-EnsemblID pairs from the six parsed S1/S2 contrast blocks; the union map is used for the MEKi rank files.",
            "xref_record_count": len(xref_records),
            "usable_symbol_count": len(xref.mapping),
            "ambiguous_symbol_count": len(xref.ambiguous),
            "ambiguous_symbols": {key: list(value) for key, value in sorted(xref.ambiguous.items())},
            "reverse_alias_target_count": len(xref.reverse_alias_targets),
            "reverse_alias_targets": {key: list(value) for key, value in sorted(xref.reverse_alias_targets.items())},
        },
        "experiments": {**{key: official_audit[key] for key in ("g1", "g2", "chr21", "OE250", "OE500", "OE_combined")}, **rank_audit},
        "mek_rank_orientation": {
            "orientation_file": orientation_path_used,
            "experiments": orientation,
            "rule": "Unresolved signs are reported as unavailable; no sign was selected to match expected biology.",
        },
        "source_audit_reference": source_audit_reference,
        "comparisons": comparisons,
        "program_summary": program_summary,
        "methods": {
            "primary_metrics": "Publisher Supplementary Tables S1/S2, exact EnsemblID rows and signed t/logFC values.",
            "rank_metrics": "Released MEKi symbol/t-statistic rank files, mapped only through the strict S1/S2 cross-reference.",
            "spearman": "Descriptive Spearman correlation of signed statistics on exact shared unique genes; average ranks for ties; constant vectors unavailable.",
            "direction_agreement": "Fraction of finite nonzero pairs with the same sign. Zero statistics remain counted separately and are excluded from the directional denominator.",
            "program": "ETS2g1_DN source set (928 rows), with ETS2 ENSG00000157557 excluded before availability and directional summaries.",
            "program_statistic_availability": "Median t reports the finite statistic as stored, even when a MEKi sign is unresolved; negative/positive fractions and other directional metrics remain unavailable until orientation is validated.",
            "independence": "Gene counts are overlap counts, not donor counts. g1 program effect is a circular source reproduction; guide, deletion, dose and drug comparisons are same-study companions.",
            "p_values": "No p-values or gene-as-independent-replicate inference are reported.",
        },
        "limitations": [
            "The public release lacks the raw macrophage count matrices for these CRISPR, overexpression and MEKi experiments; this is a summary-level reproduction.",
            "MEKi directional metrics remain unavailable unless a validated source-audit orientation is supplied; raw overlap and numeric availability accounting are retained.",
            "The rank-file symbol map is intentionally conservative and excludes ambiguous or missing identities; no aliases are guessed.",
            "This benchmark reproduces published experimental summaries and does not estimate PSC treatment benefit or clinical response.",
        ],
        "outputs_sha256": outputs,
        "analysis_script_sha256": sha256_file(Path(__file__)),
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(_json_safe(audit), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return audit


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=PLAN_PATH)
    parser.add_argument("--output-directory", type=Path, default=None, help="Directory for ets2-experiment-*.csv outputs (default: data/derived)")
    parser.add_argument("--audit-path", type=Path, default=None, help="JSON provenance report path (default: reports/ets2-experiment-benchmark.json)")
    parser.add_argument("--mek-orientation-file", type=Path, default=None, help="Optional source-audit JSON with explicit validated MEKi sign multipliers")
    args = parser.parse_args(argv)
    root = args.plan.resolve().parents[1]
    audit = run_benchmark(root=root, plan_path=args.plan, output_directory=args.output_directory, audit_path=args.audit_path, mek_orientation_path=args.mek_orientation_file)
    print(json.dumps(_json_safe(audit), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
