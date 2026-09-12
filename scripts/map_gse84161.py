#!/usr/bin/env python3
"""Map GSE84161 GPL570 probes and the frozen ETS2 source programs.

The implementation is deliberately limited to identifiers and provenance. It
reads the GPL570 annotation table, the pinned Ensembl 116 BioMart export and
the eight source gene-set columns. It never reads array intensities, computes
program scores, or estimates donor or treatment effects.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
SOURCE_MANIFEST = ROOT / "config/gse84161-gene-map-sources.json"
METHOD_LOCK = ROOT / "config/gse84161-normalization-method-lock.json"
AMENDED_RULES = ROOT / "config/ets2-validation-amended-rules.json"
PROGRAM_SOURCE_MANIFEST = ROOT / "config/ets2-program-sources.json"
GPL_SOFT = ROOT / "work/ets2-independent-search/GSE84161-family-soft.txt"
CROSSWALK = ROOT / "data/raw/gse84161/gene-map/ensembl116-human-gene-crosswalk.tsv"
CROSSWALK_COUNT = ROOT / "data/raw/gse84161/gene-map/ensembl116-biomart-count.tsv"
PROGRAM_ARCHIVE = ROOT / "data/raw/ets2-benchmark/ets2-public-release.zip"

PROBE_MAP = ROOT / "data/derived/gse84161-probe-map.csv.gz"
GENE_MAP = ROOT / "data/derived/gse84161-gene-map.csv.gz"
ENTREZ_GENE_UNIVERSE = ROOT / "data/derived/gse84161-entrez-gene-universe.csv.gz"
PROGRAM_MAP = ROOT / "data/derived/gse84161-program-mapping.csv.gz"
PROGRAM_COVERAGE = ROOT / "data/derived/gse84161-program-coverage.csv"
REPORT = ROOT / "reports/gse84161-mapping-audit.json"

EXPECTED_PROGRAM_IDS = (
    "ets2_g1_dn",
    "ets2_g1_up",
    "ets2_g2_dn",
    "chr21_dn",
    "inflammation",
    "interferon_gamma",
    "oxidative_stress",
    "apoptotic_signaling",
)
DERIVED_PROGRAM_ID = "ets2_g1_dn_without_comparators"
COMPARATOR_PROGRAM_IDS = (
    "inflammation",
    "interferon_gamma",
    "oxidative_stress",
    "apoptotic_signaling",
)
ETS2_STABLE_ID = "ENSG00000157557"
ETS2_SYMBOL = "ETS2"
EXPECTED_DENOMINATORS = {
    "ets2_g1_dn": 927,
    "ets2_g1_up": 668,
    "ets2_g2_dn": 875,
    "chr21_dn": 141,
    "inflammation": 815,
    "interferon_gamma": 139,
    "oxidative_stress": 436,
    "apoptotic_signaling": 588,
}
ENTREZ_UNRESOLVED = {"", "---", "NA", "N/A"}
ENSEMBL_RE = re.compile(r"^ENSG\d+$")
ENSEMBL_VERSIONED_RE = re.compile(r"^(ENSG\d+)\.\d+$")
ENTREZ_RE = re.compile(r"^\d+$")


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _file_provenance(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return {"path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest.hexdigest()}


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _verify_source(path: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    actual = _file_provenance(path)
    expected_bytes = int(source["bytes"])
    expected_sha = str(source["sha256"])
    if actual["bytes"] != expected_bytes:
        raise ValueError(f"Pinned source byte count changed for {path}")
    if actual["sha256"] != expected_sha:
        raise ValueError(f"Pinned source SHA-256 changed for {path}")
    return {
        **actual,
        "id": str(source.get("id", "")),
        "url": str(source.get("url", "")),
        "retrieved_at_utc": str(source.get("retrieved_at_utc", "")),
        "values_read": bool(source.get("values_read", False)),
    }


def _verify_inputs(manifest: Mapping[str, Any]) -> dict[str, Any]:
    method_record = manifest["method_lock"]
    method_actual = _file_provenance(METHOD_LOCK)["sha256"]
    if method_actual != str(method_record["sha256"]):
        raise ValueError("Normalization method lock hash changed")
    method = _read_json(METHOD_LOCK)
    for relative in (
        "docs/gse84161-normalization-method.md",
        "docs/gse84161-array-protocol.md",
        "config/ets2-validation-amended-rules.json",
        "config/ets2-benchmark-sources.json",
        "data/raw/ets2-benchmark/ets2-public-release.zip",
        "work/ets2-independent-search/GSE84161-family-soft.txt",
    ):
        expected = method.get("input_sha256", {}).get(relative)
        if expected is None:
            raise ValueError(f"Method lock does not pin {relative}")
        actual = _file_provenance(ROOT / relative)["sha256"]
        if actual != str(expected):
            raise ValueError(f"Method input hash changed for {relative}")
    mapping_doc_hash = _file_provenance(ROOT / "docs/gse84161-normalization-method.md")["sha256"]
    if mapping_doc_hash != str(method_record["mapping_method_sha256"]):
        raise ValueError("Mapping method document hash changed")
    verified_sources: list[dict[str, Any]] = []
    for source in manifest["source_files"]:
        path = ROOT / str(source["path"])
        verified_sources.append(_verify_source(path, source))
    amendment_hash = _file_provenance(AMENDED_RULES)["sha256"]
    expected_amendment_hash = method["input_sha256"]["config/ets2-validation-amended-rules.json"]
    if amendment_hash != str(expected_amendment_hash):
        raise ValueError("Amended program rules hash changed")
    return {
        "method_lock": method,
        "method_lock_sha256": method_actual,
        "source_manifest_sha256": _file_provenance(SOURCE_MANIFEST)["sha256"],
        "verified_sources": verified_sources,
    }


def _normalise_ensembl(value: str) -> tuple[str, str]:
    text = value.strip()
    if not text:
        return "", "blank"
    if ENSEMBL_RE.fullmatch(text):
        return text, "valid"
    versioned = ENSEMBL_VERSIONED_RE.fullmatch(text)
    if versioned:
        return versioned.group(1), "valid_version_stripped"
    return "", "invalid"


def _normalise_entrez(value: str) -> tuple[str, str]:
    text = value.strip()
    if not text:
        return "", "blank"
    if ENTREZ_RE.fullmatch(text):
        return text, "valid"
    return "", "invalid"


def _split_entrez_field(value: str) -> tuple[list[str], list[str]]:
    """Return unique valid numeric IDs and preserved invalid tokens."""

    valid: list[str] = []
    invalid: list[str] = []
    seen_valid: set[str] = set()
    seen_invalid: set[str] = set()
    for raw_token in value.split("///"):
        token = raw_token.strip()
        if token in ENTREZ_UNRESOLVED:
            continue
        if ENTREZ_RE.fullmatch(token):
            if token not in seen_valid:
                seen_valid.add(token)
                valid.append(token)
        elif token not in seen_invalid:
            seen_invalid.add(token)
            invalid.append(token)
    return valid, invalid


def _load_gpl570(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Read only the GPL570 platform table from the cached SOFT record."""

    header: list[str] | None = None
    rows: list[dict[str, Any]] = []
    in_table = False
    platform_seen = False
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            lower = line.lower()
            if lower == "!platform_table_begin":
                if in_table:
                    raise ValueError("Nested GPL570 platform table")
                in_table = True
                platform_seen = True
                header = None
                continue
            if in_table and lower == "!platform_table_end":
                in_table = False
                break
            if not in_table:
                continue
            if header is None:
                if not line:
                    continue
                header = line.split("\t")
                if header[0] != "ID":
                    raise ValueError(f"Unexpected GPL570 header: {header[:2]}")
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"GPL570 row has {len(fields)} fields; expected {len(header)}"
                )
            rows.append(dict(zip(header, fields)))
    if in_table or not platform_seen or header is None:
        raise ValueError("GPL570 platform table is missing or truncated")
    if len(rows) != len({row["ID"] for row in rows}):
        raise ValueError("GPL570 probe IDs are duplicated")
    return {"header": header, "row_count": len(rows)}, rows


def _probe_mapping(rows: Iterable[Mapping[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    output: list[dict[str, Any]] = []
    probes_by_entrez: dict[str, list[str]] = defaultdict(list)
    status_counts: Counter[str] = Counter()
    for source_row_number, row in enumerate(rows, start=1):
        probe_id = str(row.get("ID", "")).strip()
        raw_entrez = str(row.get("ENTREZ_GENE_ID", "")).strip()
        tokens, invalid = _split_entrez_field(raw_entrez)
        if invalid:
            status = "unmapped_invalid_entrez_token"
            mapped_entrez = ""
        elif not tokens:
            status = "unmapped_blank_or_unresolved"
            mapped_entrez = ""
        elif len(tokens) > 1:
            status = "ambiguous_multiple_entrez"
            mapped_entrez = ""
        else:
            status = "single_unique_entrez"
            mapped_entrez = tokens[0]
            probes_by_entrez[mapped_entrez].append(probe_id)
        status_counts[status] += 1
        output.append(
            {
                "probe_id": probe_id,
                "source_row_number": source_row_number,
                "entrez_gene_id_raw": raw_entrez,
                "entrez_unique_tokens": " /// ".join(tokens),
                "entrez_invalid_tokens": " /// ".join(invalid),
                "mapping_status": status,
                "entrez_gene_id": mapped_entrez,
                "gene_symbol": str(row.get("Gene Symbol", "")).strip(),
                "annotation_date": str(row.get("Annotation Date", "")).strip(),
                "eligible_probe": str(status == "single_unique_entrez").lower(),
            }
        )
    stats = {
        "probe_rows": len(output),
        "unique_probe_ids": len({row["probe_id"] for row in output}),
        "status_counts": dict(status_counts),
        "eligible_probe_rows": status_counts.get("single_unique_entrez", 0),
        "eligible_unique_entrez_genes": len(probes_by_entrez),
        "eligible_probe_counts_by_entrez": dict(probes_by_entrez),
        "annotation_dates": dict(Counter(row["annotation_date"] for row in output)),
    }
    return output, stats


def _zip_member_name(stream: zipfile.ZipFile, logical: str) -> str:
    exact = [name for name in stream.namelist() if name == logical]
    if len(exact) == 1:
        return exact[0]
    suffix = [name for name in stream.namelist() if name.endswith("/" + logical)]
    if len(suffix) != 1:
        raise ValueError(f"Expected one ZIP member {logical!r}; found {len(suffix)}")
    return suffix[0]


def _read_source_column(
    stream: zipfile.ZipFile, logical_member: str, column: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    csv.field_size_limit(max(csv.field_size_limit(), 1024 * 1024 * 1024))
    actual_member = _zip_member_name(stream, logical_member)
    content = stream.read(actual_member)
    reader = csv.reader(io.TextIOWrapper(io.BytesIO(content), encoding="utf-8", newline=""))
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValueError(f"Empty source member {logical_member}") from exc
    if column not in header:
        raise ValueError(f"Source member {logical_member} lacks {column}")
    index = header.index(column)
    values: list[dict[str, Any]] = []
    for row_number, row in enumerate(reader, start=2):
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        if len(row) > len(header):
            raise ValueError(f"Source member {logical_member} row {row_number} is too wide")
        value = row[index].strip()
        if value and value not in ENTREZ_UNRESOLVED:
            values.append({"source_row_number": row_number, "value": value})
    return values, {
        "logical_member": logical_member,
        "actual_member": actual_member,
        "bytes": len(content),
        "sha256": _sha256_bytes(content),
        "header_columns": len(header),
        "column": column,
    }


def _classify_source_namespace(values: list[str]) -> str:
    if not values:
        return "empty"
    ensembl_flags = [
        _normalise_ensembl(value)[1] in {"valid", "valid_version_stripped"}
        for value in values
    ]
    if all(ensembl_flags):
        return "Ensembl_gene_id"
    if any(ensembl_flags):
        return "mixed_source_namespace"
    return "HGNC_symbol"


def _relation_status(
    ensembl_id: str,
    entrez_id: str,
    ensembl_to_entrez: Mapping[str, set[str]],
    entrez_to_ensembl: Mapping[str, set[str]],
) -> str:
    if not ensembl_id:
        return "missing_ensembl_id"
    if not ensembl_to_entrez.get(ensembl_id):
        return "missing_entrez_relation"
    if len(ensembl_to_entrez[ensembl_id]) != 1:
        return "ambiguous_ensembl_to_entrez"
    only_entrez = next(iter(ensembl_to_entrez[ensembl_id]))
    if not entrez_id or entrez_id != only_entrez:
        return "missing_entrez_relation"
    if len(entrez_to_ensembl.get(only_entrez, set())) != 1:
        return "ambiguous_entrez_to_ensembl"
    return "reciprocal_one_to_one"


def _load_crosswalk(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    content = path.read_text(encoding="utf-8")
    lines = content.splitlines()
    if len(lines) < 3 or lines[-1] != "[success]":
        raise ValueError("Ensembl BioMart export lacks terminal [success] stamp")
    reader = csv.reader(io.StringIO("\n".join(lines[:-1])), delimiter="\t")
    header = next(reader)
    expected_header = [
        "Gene stable ID",
        "NCBI gene (formerly Entrezgene) ID",
        "HGNC symbol",
    ]
    if header != expected_header:
        raise ValueError(f"Unexpected Ensembl crosswalk header: {header}")
    raw_rows: list[dict[str, Any]] = []
    ensembl_to_entrez: dict[str, set[str]] = defaultdict(set)
    entrez_to_ensembl: dict[str, set[str]] = defaultdict(set)
    symbol_to_ensembl: dict[str, set[str]] = defaultdict(set)
    symbol_to_entrez: dict[str, set[str]] = defaultdict(set)
    # Keep the complete valid Ensembl-ID universe separate from the relation
    # map, because a valid gene row can have no Entrez value in this export.
    ensembl_universe: set[str] = set()
    seen_relationships: Counter[tuple[str, str, str]] = Counter()
    for row_number, row in enumerate(reader, start=1):
        if len(row) != 3:
            raise ValueError(f"Crosswalk row {row_number} has {len(row)} fields")
        ensembl_original, entrez_original, symbol_original = (value.strip() for value in row)
        ensembl_id, ensembl_status = _normalise_ensembl(ensembl_original)
        entrez_id, entrez_status = _normalise_entrez(entrez_original)
        symbol = symbol_original
        if ensembl_id:
            ensembl_universe.add(ensembl_id)
        if ensembl_id and entrez_id:
            ensembl_to_entrez[ensembl_id].add(entrez_id)
            entrez_to_ensembl[entrez_id].add(ensembl_id)
            seen_relationships[(ensembl_id, entrez_id, symbol)] += 1
        if symbol and ensembl_id:
            symbol_to_ensembl[symbol].add(ensembl_id)
        if symbol and entrez_id:
            symbol_to_entrez[symbol].add(entrez_id)
        if ensembl_status == "invalid":
            row_status = "invalid_ensembl_id"
        elif entrez_status == "invalid":
            row_status = "invalid_entrez_id"
        elif not ensembl_id:
            row_status = "missing_ensembl_id"
        elif not entrez_id:
            row_status = "missing_entrez_id"
        else:
            row_status = "valid_relationship"
        raw_rows.append(
            {
                "crosswalk_row_number": row_number,
                "ensembl_gene_id_original": ensembl_original,
                "ensembl_gene_id": ensembl_id,
                "entrez_gene_id_original": entrez_original,
                "entrez_gene_id": entrez_id,
                "hgnc_symbol_original": symbol_original,
                "hgnc_symbol": symbol,
                "ensembl_status": ensembl_status,
                "entrez_status": entrez_status,
                "row_status": row_status,
            }
        )
    for row in raw_rows:
        ensembl_id = row["ensembl_gene_id"]
        entrez_id = row["entrez_gene_id"]
        symbol = row["hgnc_symbol"]
        ensembl_set = ensembl_to_entrez.get(ensembl_id, set())
        entrez_set = entrez_to_ensembl.get(entrez_id, set())
        symbol_set = symbol_to_ensembl.get(symbol, set())
        row.update(
            {
                "ensembl_entrez_distinct_count": len(ensembl_set),
                "entrez_ensembl_distinct_count": len(entrez_set),
                "ensembl_entrez_reciprocal_one_to_one": str(
                    bool(
                        ensembl_id
                        and entrez_id
                        and len(ensembl_set) == 1
                        and len(entrez_set) == 1
                    )
                ).lower(),
                "hgnc_symbol_distinct_ensembl_count": len(symbol_set),
                "hgnc_symbol_unique_ensembl": str(bool(symbol and len(symbol_set) == 1)).lower(),
                "relationship_status": _relation_status(
                    ensembl_id, entrez_id, ensembl_to_entrez, entrez_to_ensembl
                ),
            }
        )
    reciprocal = {
        ensembl: next(iter(entrez_ids))
        for ensembl, entrez_ids in ensembl_to_entrez.items()
        if len(entrez_ids) == 1
        and len(entrez_to_ensembl[next(iter(entrez_ids))]) == 1
    }
    stats = {
        "header": header,
        "response_relationship_rows": len(raw_rows),
        "unique_ensembl_gene_ids": len(ensembl_universe),
        "unique_entrez_gene_ids_with_any_relation": len(entrez_to_ensembl),
        "valid_relationship_rows": sum(row["row_status"] == "valid_relationship" for row in raw_rows),
        "ensembl_ids_with_any_entrez": len(ensembl_to_entrez),
        "ensembl_ids_without_entrez": len(ensembl_universe - set(ensembl_to_entrez)),
        "reciprocal_one_to_one_ensembl_genes": len(reciprocal),
        "reciprocal_one_to_one_entrez_genes": len({value for value in reciprocal.values()}),
        "hgnc_symbols_with_ensembl": len(symbol_to_ensembl),
        "unique_hgnc_symbols": sum(len(values) == 1 for values in symbol_to_ensembl.values()),
        "duplicate_relationship_rows": sum(count - 1 for count in seen_relationships.values() if count > 1),
        "relationship_rows_with_duplicate_identity": sum(count > 1 for count in seen_relationships.values()),
        "terminal_success_stamp": "[success]",
    }
    maps = {
        "ensembl_to_entrez": ensembl_to_entrez,
        "entrez_to_ensembl": entrez_to_ensembl,
        "symbol_to_ensembl": symbol_to_ensembl,
        "symbol_to_entrez": symbol_to_entrez,
        "ensembl_universe": ensembl_universe,
        "reciprocal_ensembl_to_entrez": reciprocal,
    }
    return raw_rows, {"stats": stats, "maps": maps}


def _array_entrez_gene_map(
    probes_by_entrez: Mapping[str, list[str]],
    crosswalk_maps: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Join every eligible array Entrez gene to the full crosswalk.

    This table is keyed by the complete single-unique-Entrez GPL universe,
    including genes with no valid Ensembl relation in the pinned crosswalk.
    It records all Ensembl IDs before assigning one deterministic relation
    status, so array availability is not mistaken for crosswalk coverage.
    """

    ensembl_to_entrez = crosswalk_maps["ensembl_to_entrez"]
    entrez_to_ensembl = crosswalk_maps["entrez_to_ensembl"]
    rows: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()

    def entrez_sort_key(value: str) -> tuple[int, str]:
        return (int(value), value)

    def ensembl_sort_key(value: str) -> tuple[int, str]:
        suffix = value[4:]
        return (int(suffix), value)

    for entrez_id in sorted(probes_by_entrez, key=entrez_sort_key):
        probe_ids = list(probes_by_entrez[entrez_id])
        ensembl_ids = sorted(entrez_to_ensembl.get(entrez_id, set()), key=ensembl_sort_key)
        relation_counts = [
            f"{ensembl_id}:{len(ensembl_to_entrez.get(ensembl_id, set()))}"
            for ensembl_id in ensembl_ids
        ]
        if not ensembl_ids:
            status = "no_ensembl_relation"
            reciprocal_ensembl = ""
        elif len(ensembl_ids) != 1:
            status = "ambiguous_entrez_to_ensembl"
            reciprocal_ensembl = ""
        elif len(ensembl_to_entrez.get(ensembl_ids[0], set())) != 1:
            status = "ambiguous_ensembl_to_entrez"
            reciprocal_ensembl = ""
        else:
            status = "reciprocal_one_to_one"
            reciprocal_ensembl = ensembl_ids[0]
        status_counts[status] += 1
        rows.append(
            {
                "entrez_gene_id": entrez_id,
                "eligible_array_probe_count": len(probe_ids),
                "eligible_array_probe_ids": " /// ".join(probe_ids),
                "full_crosswalk_ensembl_gene_ids": " /// ".join(ensembl_ids),
                "full_crosswalk_ensembl_gene_count": len(ensembl_ids),
                "full_crosswalk_entrez_count_by_ensembl": " /// ".join(relation_counts),
                "crosswalk_mapping_status": status,
                "reciprocal_ensembl_gene_id": reciprocal_ensembl,
            }
        )
    stats = {
        "array_entrez_gene_count": len(rows),
        "status_counts": dict(sorted(status_counts.items())),
        "all_genes_have_eligible_probe": all(
            int(row["eligible_array_probe_count"]) > 0 for row in rows
        ),
        "status_definitions": {
            "reciprocal_one_to_one": "Exactly one full-crosswalk Ensembl ID and exactly one Entrez relation for that Ensembl ID.",
            "ambiguous_entrez_to_ensembl": "The full crosswalk assigns this Entrez ID to multiple Ensembl IDs.",
            "ambiguous_ensembl_to_entrez": "The full crosswalk assigns the sole Ensembl ID to multiple Entrez IDs.",
            "no_ensembl_relation": "No valid Ensembl-to-Entrez relationship for this Entrez ID in the pinned full crosswalk; this does not assert retirement.",
        },
    }
    return rows, stats


def _program_row_mapping(
    *,
    program_id: str,
    source_member_index: int,
    source_row_number: int,
    source_value: str,
    source_namespace: str,
    duplicate_source_member: bool,
    crosswalk_maps: Mapping[str, Any],
    probes_by_entrez: Mapping[str, list[str]],
) -> dict[str, Any]:
    ensembl_to_entrez = crosswalk_maps["ensembl_to_entrez"]
    entrez_to_ensembl = crosswalk_maps["entrez_to_ensembl"]
    symbol_to_ensembl = crosswalk_maps["symbol_to_ensembl"]
    ensembl_universe = crosswalk_maps.get("ensembl_universe", set(ensembl_to_entrez))
    reciprocal = crosswalk_maps["reciprocal_ensembl_to_entrez"]
    source_normalized = source_value.strip()
    ensembl_id = ""
    entrez_id = ""
    base_status = ""
    symbol_ensembl_count = 0
    if source_namespace == "Ensembl_gene_id":
        ensembl_id, status = _normalise_ensembl(source_normalized)
        if status == "invalid":
            base_status = "invalid_ensembl_id"
        elif not ensembl_id:
            base_status = "missing_ensembl_id"
        elif ensembl_id not in ensembl_universe:
            base_status = "missing_ensembl_id_in_snapshot"
        elif ensembl_id not in ensembl_to_entrez:
            base_status = "missing_entrez_relation"
        elif ensembl_id not in reciprocal:
            base_status = "ambiguous_ensembl_entrez_relation"
        else:
            entrez_id = str(reciprocal[ensembl_id])
            base_status = "reciprocal_one_to_one"
    elif source_namespace == "HGNC_symbol":
        symbol_ensembl_count = len(symbol_to_ensembl.get(source_normalized, set()))
        if not source_normalized:
            base_status = "missing_hgnc_symbol"
        elif symbol_ensembl_count == 0:
            base_status = "missing_hgnc_symbol"
        elif symbol_ensembl_count != 1:
            base_status = "ambiguous_hgnc_symbol"
        else:
            ensembl_id = next(iter(symbol_to_ensembl[source_normalized]))
            if ensembl_id not in ensembl_universe:
                base_status = "missing_ensembl_id_in_snapshot"
            elif ensembl_id not in ensembl_to_entrez:
                base_status = "missing_entrez_relation"
            elif ensembl_id not in reciprocal:
                base_status = "ambiguous_ensembl_entrez_relation"
            else:
                entrez_id = str(reciprocal[ensembl_id])
                base_status = "reciprocal_one_to_one"
    else:
        base_status = "invalid_source_namespace"
    ets2_excluded = bool(
        source_normalized == ETS2_SYMBOL or ensembl_id == ETS2_STABLE_ID
    )
    probe_ids = probes_by_entrez.get(entrez_id, []) if base_status == "reciprocal_one_to_one" else []
    if base_status == "reciprocal_one_to_one" and probe_ids:
        array_status = "mapped_eligible_array_gene"
    elif base_status == "reciprocal_one_to_one":
        array_status = "mapped_crosswalk_no_eligible_array_probe"
    else:
        array_status = "not_array_mapped"
    program_status = "excluded_ets2" if ets2_excluded else array_status if not duplicate_source_member else "duplicate_exact_source_member"
    return {
        "program_id": program_id,
        "parent_program_id": "",
        "source_member_index": source_member_index,
        "source_row_number": source_row_number,
        "source_member_original": source_value,
        "source_member_normalized": source_normalized,
        "source_namespace": source_namespace,
        "duplicate_exact_source_member": str(duplicate_source_member).lower(),
        "ets2_excluded": str(ets2_excluded).lower(),
        "intended_member_for_coverage": str(not duplicate_source_member and not ets2_excluded).lower(),
        "crosswalk_mapping_status": base_status,
        "ensembl_gene_id": ensembl_id,
        "entrez_gene_id": entrez_id,
        "hgnc_symbol_exact": source_normalized if source_namespace == "HGNC_symbol" else "",
        "symbol_distinct_ensembl_count": symbol_ensembl_count,
        "array_mapping_status": array_status,
        "eligible_array_probe_count": len(probe_ids),
        "eligible_array_probe_ids": " /// ".join(probe_ids),
        "program_member_status": program_status,
        "derived_crosswalk_comparator_overlap": "",
        "derived_array_comparator_overlap": "",
        "derived_member_for_coverage": "",
    }


def _load_program_rows(
    *,
    plan: Mapping[str, Any],
    source_zip: Path,
    crosswalk_maps: Mapping[str, Any],
    probes_by_entrez: Mapping[str, list[str]],
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, dict[str, set[str]]]]:
    program_items = [
        item for item in plan["programs"] if str(item.get("id")) in EXPECTED_PROGRAM_IDS
    ]
    if tuple(str(item.get("id")) for item in program_items) != EXPECTED_PROGRAM_IDS:
        raise ValueError("The frozen eight source program definitions changed")
    direct_rows: list[dict[str, Any]] = []
    source_records: dict[str, dict[str, Any]] = {}
    with zipfile.ZipFile(source_zip) as stream:
        for item in program_items:
            program_id = str(item["id"])
            logical_member = str(item["source_member"])
            column = str(item["column"])
            values, member = _read_source_column(stream, logical_member, column)
            source_values = [entry["value"] for entry in values]
            namespace = _classify_source_namespace(source_values)
            seen: set[str] = set()
            rows_for_program: list[dict[str, Any]] = []
            for member_index, entry in enumerate(values, start=1):
                source_value = entry["value"]
                duplicate = source_value in seen
                seen.add(source_value)
                rows_for_program.append(
                    _program_row_mapping(
                        program_id=program_id,
                        source_member_index=member_index,
                        source_row_number=int(entry["source_row_number"]),
                        source_value=source_value,
                        source_namespace=namespace,
                        duplicate_source_member=duplicate,
                        crosswalk_maps=crosswalk_maps,
                        probes_by_entrez=probes_by_entrez,
                    )
                )
            direct_rows.extend(rows_for_program)
            source_records[program_id] = {
                "source_member": logical_member,
                "column": column,
                "actual_member": member["actual_member"],
                "member_bytes": member["bytes"],
                "member_sha256": member["sha256"],
                "source_namespace": namespace,
                "source_member_count_raw": len(source_values),
                "source_unique_member_count": len(set(source_values)),
                "duplicate_exact_source_member_count": len(source_values) - len(set(source_values)),
                "source_values_are_nonblank": all(bool(value) for value in source_values),
            }
    direct_by_program: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in direct_rows:
        direct_by_program[row["program_id"]].append(row)
    crosswalk_gene_sets: dict[str, set[str]] = {}
    array_gene_sets: dict[str, set[str]] = {}
    for program_id in EXPECTED_PROGRAM_IDS:
        rows = direct_by_program[program_id]
        crosswalk_gene_sets[program_id] = {
            row["entrez_gene_id"]
            for row in rows
            if row["intended_member_for_coverage"] == "true" and row["entrez_gene_id"]
        }
        array_gene_sets[program_id] = {
            row["entrez_gene_id"]
            for row in rows
            if row["intended_member_for_coverage"] == "true"
            and row["array_mapping_status"] == "mapped_eligible_array_gene"
            and row["entrez_gene_id"]
        }
    crosswalk_comparator_union = set().union(*(crosswalk_gene_sets[pid] for pid in COMPARATOR_PROGRAM_IDS))
    array_comparator_union = set().union(*(array_gene_sets[pid] for pid in COMPARATOR_PROGRAM_IDS))
    primary_rows = direct_by_program["ets2_g1_dn"]
    derived_rows: list[dict[str, Any]] = []
    for original in primary_rows:
        row = dict(original)
        row["program_id"] = DERIVED_PROGRAM_ID
        row["parent_program_id"] = "ets2_g1_dn"
        ensembl_id = row["ensembl_gene_id"]
        entrez_id = row["entrez_gene_id"]
        crosswalk_overlap = bool(entrez_id and entrez_id in crosswalk_comparator_union)
        array_overlap = bool(entrez_id and entrez_id in array_comparator_union)
        row["derived_crosswalk_comparator_overlap"] = str(crosswalk_overlap).lower()
        row["derived_array_comparator_overlap"] = str(array_overlap).lower()
        if row["ets2_excluded"] == "true":
            row["program_member_status"] = "excluded_ets2"
            row["derived_member_for_coverage"] = "false"
        elif row["array_mapping_status"] == "mapped_eligible_array_gene" and array_overlap:
            row["program_member_status"] = "mapped_comparator_overlap"
            row["derived_member_for_coverage"] = "false"
        elif row["array_mapping_status"] == "mapped_eligible_array_gene":
            row["program_member_status"] = "mapped_nonoverlap"
            row["derived_member_for_coverage"] = "true"
        else:
            row["derived_member_for_coverage"] = "false"
        derived_rows.append(row)
    all_rows = direct_rows + derived_rows
    derived_sets = {
        DERIVED_PROGRAM_ID: {
            "crosswalk": crosswalk_gene_sets["ets2_g1_dn"] - crosswalk_comparator_union,
            "array": array_gene_sets["ets2_g1_dn"] - array_comparator_union,
            "crosswalk_comparator_union": crosswalk_comparator_union,
            "array_comparator_union": array_comparator_union,
            "parent_crosswalk": crosswalk_gene_sets["ets2_g1_dn"],
            "parent_array": array_gene_sets["ets2_g1_dn"],
        }
    }
    source_records[DERIVED_PROGRAM_ID] = {
        "source_member": "derived from ets2_g1_dn in Entrez space",
        "column": "mapped primary minus mapped comparator union",
        "source_namespace": "derived_entrez_space",
        "source_member_count_raw": source_records["ets2_g1_dn"]["source_member_count_raw"],
        "source_unique_member_count": source_records["ets2_g1_dn"]["source_unique_member_count"],
        "duplicate_exact_source_member_count": source_records["ets2_g1_dn"]["duplicate_exact_source_member_count"],
        "parent_program_id": "ets2_g1_dn",
        "comparator_program_ids": list(COMPARATOR_PROGRAM_IDS),
    }
    return all_rows, source_records, {**{key: {"crosswalk": crosswalk_gene_sets[key], "array": array_gene_sets[key]} for key in EXPECTED_PROGRAM_IDS}, **derived_sets}


def _coverage_rows(
    rows: list[dict[str, Any]],
    source_records: Mapping[str, Mapping[str, Any]],
    gene_sets: Mapping[str, Mapping[str, set[str]]],
    method: Mapping[str, Any],
) -> list[dict[str, Any]]:
    mapping_gates = method["mapping_gates"]
    rows_by_program: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        rows_by_program[row["program_id"]].append(row)
    coverage: list[dict[str, Any]] = []
    for program_id in (*EXPECTED_PROGRAM_IDS, DERIVED_PROGRAM_ID):
        source = source_records[program_id]
        member_rows = rows_by_program[program_id]
        if program_id == DERIVED_PROGRAM_ID:
            parent = next(row for row in coverage if row["program_id"] == "ets2_g1_dn")
            intended = int(parent["intended_member_count_after_ets2_exclusion"])
            excluded = int(parent["ets2_excluded_member_count"])
            parent_gate = parent["coverage_gate_pass"] == "true"
            crosswalk_gene_count = len(gene_sets[program_id]["crosswalk"])
            array_gene_count = len(gene_sets[program_id]["array"])
            coverage_pass = parent_gate and array_gene_count >= int(mapping_gates["minimum_genes"])
            # The derived set is deliberately gated by the parent's array
            # coverage. Keep its retained non-overlap fraction visible too,
            # so a passing parent gate cannot be mistaken for 652/927.
            retained_fraction = array_gene_count / intended if intended else 0.0
            gate_array_gene_count = int(parent["array_mapped_unique_entrez_gene_count"])
            gate_fraction = float(parent["mapping_fraction"])
            gate_basis = "parent ets2_g1_dn coverage gate plus derived remaining-gene minimum"
            overlap_count = len(gene_sets[program_id]["parent_array"] & gene_sets[program_id]["array_comparator_union"])
            parent_mapping_denominator = intended
            parent_array_gene_count = gate_array_gene_count
            parent_mapping_fraction = gate_fraction
        else:
            expected = EXPECTED_DENOMINATORS[program_id]
            intended = sum(
                row["intended_member_for_coverage"] == "true" for row in member_rows
            )
            if intended != expected:
                raise ValueError(
                    f"Frozen denominator changed for {program_id}: {intended} != {expected}"
                )
            excluded = sum(row["ets2_excluded"] == "true" for row in member_rows)
            crosswalk_gene_count = len(gene_sets[program_id]["crosswalk"])
            array_gene_count = len(gene_sets[program_id]["array"])
            retained_fraction = array_gene_count / intended if intended else 0.0
            gate_array_gene_count = array_gene_count
            gate_fraction = retained_fraction
            coverage_pass = gate_fraction >= float(mapping_gates["fraction"]) and array_gene_count >= int(mapping_gates["minimum_genes"])
            gate_basis = "source member denominator after ETS2 exclusion"
            overlap_count = 0
            parent_mapping_denominator = ""
            parent_array_gene_count = ""
            parent_mapping_fraction = ""
        status_counts = Counter(row["program_member_status"] for row in member_rows)
        crosswalk_status_counts = Counter(
            row["crosswalk_mapping_status"]
            for row in member_rows
            if row["intended_member_for_coverage"] == "true"
        )
        array_status_counts = Counter(
            row["array_mapping_status"]
            for row in member_rows
            if row["intended_member_for_coverage"] == "true"
        )
        if program_id == DERIVED_PROGRAM_ID:
            mapped_source_member_count = sum(
                row["derived_member_for_coverage"] == "true" for row in member_rows
            )
        else:
            mapped_source_member_count = sum(
                row["intended_member_for_coverage"] == "true"
                and row["array_mapping_status"] == "mapped_eligible_array_gene"
                for row in member_rows
            )
        coverage.append(
            {
                "program_id": program_id,
                "parent_program_id": "ets2_g1_dn" if program_id == DERIVED_PROGRAM_ID else "",
                "source_namespace": str(source["source_namespace"]),
                "source_member_count_raw": source["source_member_count_raw"],
                "source_unique_member_count": source["source_unique_member_count"],
                "ets2_excluded_member_count": excluded,
                "intended_member_count_after_ets2_exclusion": intended,
                "crosswalk_mapped_unique_entrez_gene_count": crosswalk_gene_count,
                "array_mapped_unique_entrez_gene_count": array_gene_count,
                "array_mapped_source_member_count": mapped_source_member_count,
                # mapping_fraction is the fraction evaluated by the gate.
                # For the derived row this is the parent's 801/927 fraction;
                # retained_mapping_fraction is the derived 652/927 fraction.
                "mapping_fraction": f"{gate_fraction:.12g}",
                "retained_mapping_fraction": f"{retained_fraction:.12g}",
                "gate_evaluated_mapping_fraction": f"{gate_fraction:.12g}",
                "gate_evaluated_mapped_unique_entrez_gene_count": gate_array_gene_count,
                "retained_mapped_unique_entrez_gene_count": array_gene_count,
                "retained_mapping_denominator": intended,
                "parent_mapping_denominator": parent_mapping_denominator,
                "parent_array_mapped_unique_entrez_gene_count": parent_array_gene_count,
                "parent_mapping_fraction": (
                    f"{parent_mapping_fraction:.12g}"
                    if isinstance(parent_mapping_fraction, float)
                    else parent_mapping_fraction
                ),
                "coverage_gate_fraction": mapping_gates["fraction"],
                "coverage_gate_minimum_genes": mapping_gates["minimum_genes"],
                "coverage_gate_basis": gate_basis,
                "coverage_gate_pass": str(coverage_pass).lower(),
                "mapped_comparator_overlap_unique_entrez_gene_count": overlap_count,
                "member_status_counts": json.dumps(dict(sorted(status_counts.items())), sort_keys=True),
                "crosswalk_status_counts": json.dumps(dict(sorted(crosswalk_status_counts.items())), sort_keys=True),
                "array_status_counts": json.dumps(dict(sorted(array_status_counts.items())), sort_keys=True),
                "source_member_definition": source["column"],
                "source_member_file": source["source_member"],
            }
        )
    return coverage


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fieldnames: list[str], compressed: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        with path.open("wb") as binary:
            with gzip.GzipFile(fileobj=binary, mode="wb", mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=fieldnames, lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
    else:
        with path.open("w", encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(text, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def _relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def build_audit(*, write_outputs: bool = True) -> dict[str, Any]:
    source_manifest = _read_json(SOURCE_MANIFEST)
    input_audit = _verify_inputs(source_manifest)
    method = input_audit["method_lock"]
    amendment = _read_json(AMENDED_RULES)
    amendment_programs = [
        item for item in amendment["programs"] if str(item.get("id")) in EXPECTED_PROGRAM_IDS
    ]
    if tuple(str(item["id"]) for item in amendment_programs) != EXPECTED_PROGRAM_IDS:
        raise ValueError("Amended rules do not expose exactly the eight source programs")
    gpl_summary, gpl_rows = _load_gpl570(GPL_SOFT)
    probe_rows, probe_stats = _probe_mapping(gpl_rows)
    probes_by_entrez = {
        gene_id: [row["probe_id"] for row in probe_rows if row["entrez_gene_id"] == gene_id]
        for gene_id in probe_stats["eligible_probe_counts_by_entrez"]
    }
    crosswalk_rows, crosswalk = _load_crosswalk(CROSSWALK)
    crosswalk_stats = crosswalk["stats"]
    count_text = CROSSWALK_COUNT.read_text(encoding="utf-8")
    expected_crosswalk_count = int(count_text.strip())
    if expected_crosswalk_count != crosswalk_stats["unique_ensembl_gene_ids"]:
        raise ValueError("BioMart count response does not match unique crosswalk Ensembl IDs")
    entrez_gene_rows, entrez_gene_stats = _array_entrez_gene_map(
        probes_by_entrez,
        crosswalk["maps"],
    )
    rows, source_records, gene_sets = _load_program_rows(
        plan=amendment,
        source_zip=PROGRAM_ARCHIVE,
        crosswalk_maps=crosswalk["maps"],
        probes_by_entrez=probes_by_entrez,
    )
    coverage = _coverage_rows(rows, source_records, gene_sets, method)
    source_member_fingerprints = {
        program_id: {
            key: value
            for key, value in source_records[program_id].items()
            if key in {"source_member", "column", "actual_member", "member_bytes", "member_sha256", "source_namespace"}
        }
        for program_id in EXPECTED_PROGRAM_IDS
    }
    report = {
        "audit_id": "gse84161-mapping-audit-v1",
        "status": "identifier_mapping_and_program_coverage_only",
        "created_at_utc": "2026-09-12",
        "scope": {
            "expression_values_read": False,
            "array_intensities_read": False,
            "program_scores_computed": False,
            "target_values_or_differences_read": False,
            "donor_or_treatment_effects_computed": False,
            "purpose": "fixed probe/gene crosswalk and source-program coverage before array RMA aggregation",
        },
        "source_manifest": {
            "path": _relative(SOURCE_MANIFEST),
            "sha256": input_audit["source_manifest_sha256"],
            "verified_source_files": input_audit["verified_sources"],
            "query": source_manifest["biomart_query"],
            "completeness_check": source_manifest["completeness_check"],
        },
        "method_lock": {
            "path": _relative(METHOD_LOCK),
            "sha256": input_audit["method_lock_sha256"],
            "mapping_document_sha256": source_manifest["method_lock"]["mapping_method_sha256"],
            "source_namespace_clarification": method["source_namespace_clarification"],
        },
        "gpl570": {
            "source_path": _relative(GPL_SOFT),
            "table_header": gpl_summary["header"],
            "table_row_count": gpl_summary["row_count"],
            "table_unique_probe_count": probe_stats["unique_probe_ids"],
            "probe_status_counts": probe_stats["status_counts"],
            "eligible_probe_rows": probe_stats["eligible_probe_rows"],
            "eligible_unique_entrez_genes": probe_stats["eligible_unique_entrez_genes"],
            "annotation_dates": probe_stats["annotation_dates"],
            "mapping_rule": source_manifest["mapping_rules"]["probe_to_entrez"],
        },
        "ensembl116_crosswalk": {
            "source_path": _relative(CROSSWALK),
            "count_source_path": _relative(CROSSWALK_COUNT),
            "crosswalk_sha256": _file_provenance(CROSSWALK)["sha256"],
            "count_response_sha256": _file_provenance(CROSSWALK_COUNT)["sha256"],
            "count_response_value": expected_crosswalk_count,
            **crosswalk_stats,
            "mapping_rule": source_manifest["mapping_rules"]["ensembl_to_entrez"],
            "symbol_mapping_rule": source_manifest["mapping_rules"]["symbol_programs"],
        },
        "array_entrez_gene_universe": {
            "output_path": _relative(ENTREZ_GENE_UNIVERSE),
            **entrez_gene_stats,
            "mapping_rule": "One row for each of the 20,486 eligible single-unique-Entrez GPL570 genes. Join all valid Ensembl IDs from the complete crosswalk, preserve eligible probe IDs, and classify the full relation as reciprocal_one_to_one, ambiguous_entrez_to_ensembl, ambiguous_ensembl_to_entrez or no_ensembl_relation.",
        },
        "source_programs": {
            "source_archive_path": _relative(PROGRAM_ARCHIVE),
            "source_archive_sha256": _file_provenance(PROGRAM_ARCHIVE)["sha256"],
            "member_fingerprints": source_member_fingerprints,
            "definitions": [
                {
                    "program_id": program_id,
                    "source_namespace": source_records[program_id]["source_namespace"],
                    "source_member": source_records[program_id]["source_member"],
                    "column": source_records[program_id]["column"],
                    "source_member_count_raw": source_records[program_id]["source_member_count_raw"],
                    "source_unique_member_count": source_records[program_id]["source_unique_member_count"],
                    "ets2_excluded_member_count": sum(
                        row["program_id"] == program_id and row["ets2_excluded"] == "true" for row in rows
                    ),
                    "intended_member_count_after_ets2_exclusion": EXPECTED_DENOMINATORS[program_id],
                }
                for program_id in EXPECTED_PROGRAM_IDS
            ],
            "derived_program": {
                "program_id": DERIVED_PROGRAM_ID,
                "parent_program_id": "ets2_g1_dn",
                "comparator_program_ids": list(COMPARATOR_PROGRAM_IDS),
                "mapped_crosswalk_unique_entrez_count": len(gene_sets[DERIVED_PROGRAM_ID]["crosswalk"]),
                "mapped_array_unique_entrez_count": len(gene_sets[DERIVED_PROGRAM_ID]["array"]),
                "mapped_array_comparator_union_count": len(gene_sets[DERIVED_PROGRAM_ID]["array_comparator_union"]),
                "mapped_array_comparator_overlap_count": len(
                    gene_sets[DERIVED_PROGRAM_ID]["parent_array"]
                    & gene_sets[DERIVED_PROGRAM_ID]["array_comparator_union"]
                ),
                "parent_intended_member_count": int(
                    coverage[0]["intended_member_count_after_ets2_exclusion"]
                ),
                "parent_array_mapped_unique_entrez_count": len(gene_sets[DERIVED_PROGRAM_ID]["parent_array"]),
                "parent_mapping_fraction": coverage[0]["mapping_fraction"],
                "retained_array_mapped_unique_entrez_count": len(gene_sets[DERIVED_PROGRAM_ID]["array"]),
                "retained_mapping_fraction": coverage[-1]["retained_mapping_fraction"],
                "derivation": source_manifest["mapping_rules"]["coverage_gate"],
            },
        },
        "program_coverage": coverage,
        "output_paths": {
            "probe_map": _relative(PROBE_MAP),
            "gene_map": _relative(GENE_MAP),
            "entrez_gene_universe": _relative(ENTREZ_GENE_UNIVERSE),
            "program_mapping": _relative(PROGRAM_MAP),
            "program_coverage": _relative(PROGRAM_COVERAGE),
            "report": _relative(REPORT),
        },
        "limitations": [
            "The BioMart full export has a terminal [success] stamp and an independent unique-Ensembl count response; this supports the returned export boundary but cannot prove absence of an unobserved remote backend omission.",
            "Reciprocal one-to-one is evaluated over all 91,748 returned relationship rows before array or program restriction. Duplicate relationships are retained and reported.",
            "A source Ensembl ID absent from the pinned valid-ID universe is reported separately from an in-snapshot ID without an Entrez relation; neither status is interpreted as gene retirement.",
            "The four GOBP source columns are exact HGNC symbols; no aliases, previous names, case changes or GPL symbol fallback are used.",
            "Coverage is identifier coverage only. No expression values, array effects, program scores or donor contrasts were read or computed.",
        ],
        "execution": {
            "script_path": _relative(Path(__file__).resolve()),
            "script_sha256": _file_provenance(Path(__file__).resolve())["sha256"],
            "method_lock_verified": True,
            "source_manifest_verified": True,
        },
    }
    if write_outputs:
        probe_fields = [
            "probe_id", "source_row_number", "entrez_gene_id_raw", "entrez_unique_tokens",
            "entrez_invalid_tokens", "mapping_status", "entrez_gene_id", "gene_symbol",
            "annotation_date", "eligible_probe",
        ]
        gene_fields = [
            "crosswalk_row_number", "ensembl_gene_id_original", "ensembl_gene_id",
            "entrez_gene_id_original", "entrez_gene_id", "hgnc_symbol_original",
            "hgnc_symbol", "ensembl_status", "entrez_status", "row_status",
            "ensembl_entrez_distinct_count", "entrez_ensembl_distinct_count",
            "ensembl_entrez_reciprocal_one_to_one", "hgnc_symbol_distinct_ensembl_count",
            "hgnc_symbol_unique_ensembl", "relationship_status",
        ]
        entrez_gene_fields = [
            "entrez_gene_id", "eligible_array_probe_count", "eligible_array_probe_ids",
            "full_crosswalk_ensembl_gene_ids", "full_crosswalk_ensembl_gene_count",
            "full_crosswalk_entrez_count_by_ensembl", "crosswalk_mapping_status",
            "reciprocal_ensembl_gene_id",
        ]
        program_fields = [
            "program_id", "parent_program_id", "source_member_index", "source_row_number",
            "source_member_original", "source_member_normalized", "source_namespace",
            "duplicate_exact_source_member", "ets2_excluded", "intended_member_for_coverage",
            "crosswalk_mapping_status", "ensembl_gene_id", "entrez_gene_id",
            "hgnc_symbol_exact", "symbol_distinct_ensembl_count", "array_mapping_status",
            "eligible_array_probe_count", "eligible_array_probe_ids", "program_member_status",
            "derived_crosswalk_comparator_overlap", "derived_array_comparator_overlap",
            "derived_member_for_coverage",
        ]
        coverage_fields = list(coverage[0]) if coverage else []
        _write_csv(PROBE_MAP, probe_rows, probe_fields, compressed=True)
        _write_csv(GENE_MAP, crosswalk_rows, gene_fields, compressed=True)
        _write_csv(
            ENTREZ_GENE_UNIVERSE,
            entrez_gene_rows,
            entrez_gene_fields,
            compressed=True,
        )
        _write_csv(PROGRAM_MAP, rows, program_fields, compressed=True)
        _write_csv(PROGRAM_COVERAGE, coverage, coverage_fields)
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    report = build_audit(write_outputs=True)
    print(
        json.dumps(
            {
                "gpl570_probes": report["gpl570"]["table_row_count"],
                "eligible_probe_rows": report["gpl570"]["eligible_probe_rows"],
                "eligible_entrez_genes": report["gpl570"]["eligible_unique_entrez_genes"],
                "crosswalk_relationship_rows": report["ensembl116_crosswalk"]["response_relationship_rows"],
                "crosswalk_unique_ensembl": report["ensembl116_crosswalk"]["unique_ensembl_gene_ids"],
                "programs": len(report["program_coverage"]),
                "expression_values_read": report["scope"]["expression_values_read"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
