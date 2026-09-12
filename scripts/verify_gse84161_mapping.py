#!/usr/bin/env python3
"""Independent, expression-free GSE84161 probe and program mapping check.

The script authenticates every pinned input before parsing the GPL570
annotation, Ensembl 116 relationship export, or source-program columns.  It
does not open CEL files or read expression values.  Its mapping calculation
is intentionally separate from the production analysis implementation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
METHOD_LOCK_PATH = ROOT / "config/gse84161-normalization-method-lock.json"
GENE_MAP_MANIFEST_PATH = ROOT / "config/gse84161-gene-map-sources.json"
DEFAULT_OUTPUT = ROOT / "reports/gse84161-mapping-independent-check.json"

EXPECTED_SOURCE_SPECS = {
    "ets2_g1_dn": (
        "ensembl_gene_id",
        "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv",
        "ETS2g1_DN",
    ),
    "ets2_g1_up": (
        "ensembl_gene_id",
        "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv",
        "ETS2g1_UP",
    ),
    "ets2_g2_dn": (
        "ensembl_gene_id",
        "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv",
        "ETS2g2_DN",
    ),
    "chr21_dn": (
        "ensembl_gene_id",
        "RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv",
        "CHR21_DN",
    ),
    "inflammation": (
        "hgnc_symbol",
        "RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv",
        "GOBP_INFLAMMATORY_RESPONSE",
    ),
    "interferon_gamma": (
        "hgnc_symbol",
        "RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv",
        "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
    ),
    "oxidative_stress": (
        "hgnc_symbol",
        "RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv",
        "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
    ),
    "apoptotic_signaling": (
        "hgnc_symbol",
        "RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv",
        "GOBP_APOPTOTIC_SIGNALING_PATHWAY",
    ),
}

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

PRIMARY_MAPPING_CANDIDATES = (
    ROOT / "reports/gse84161-normalization-mapping.json",
    ROOT / "reports/gse84161-normalization-plan.json",
    ROOT / "data/derived/gse84161-normalization-mapping.csv",
    ROOT / "data/derived/gse84161-program-coverage.csv",
    ROOT / "work/gse84161-normalization/primary-mapping.json",
)


class MappingCheckError(RuntimeError):
    """A pinned input or structural mapping invariant failed."""


def sha256_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            total += len(chunk)
    return total, digest.hexdigest()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _verify_file(path: Path, expected_bytes: int | None, expected_sha256: str | None) -> dict[str, Any]:
    if not path.exists():
        raise MappingCheckError(f"missing pinned input: {path}")
    byte_count, digest = sha256_file(path)
    if expected_bytes is not None and byte_count != int(expected_bytes):
        raise MappingCheckError(
            f"byte count changed for {path}: {byte_count} != {expected_bytes}"
        )
    if expected_sha256 is not None and digest.lower() != str(expected_sha256).lower():
        raise MappingCheckError(
            f"SHA-256 changed for {path}: {digest} != {expected_sha256}"
        )
    return {"path": _relative(path), "bytes": byte_count, "sha256": digest}


def _load_and_verify_sources() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Authenticate all locked inputs before any source table is parsed."""

    if not METHOD_LOCK_PATH.exists() or not GENE_MAP_MANIFEST_PATH.exists():
        raise MappingCheckError("method lock or gene-map manifest is missing")
    method_lock = json.loads(METHOD_LOCK_PATH.read_text(encoding="utf-8"))
    gene_manifest = json.loads(GENE_MAP_MANIFEST_PATH.read_text(encoding="utf-8"))

    expected_method_lock_sha = gene_manifest.get("method_lock", {}).get("sha256")
    lock_hash_record = _verify_file(
        METHOD_LOCK_PATH, None, str(expected_method_lock_sha) if expected_method_lock_sha else None
    )
    verified: list[dict[str, Any]] = [lock_hash_record]

    for relative, expected in sorted(method_lock.get("input_sha256", {}).items()):
        path = Path(str(relative))
        if not path.is_absolute():
            path = ROOT / path
        verified.append(_verify_file(path, None, str(expected)))

    manifest_sources: dict[str, dict[str, Any]] = {}
    for source in gene_manifest.get("source_files", []):
        path = ROOT / str(source["path"])
        record = _verify_file(path, int(source["bytes"]), str(source["sha256"]))
        manifest_sources[str(source["id"])] = {
            **record,
            "manifest_id": str(source["id"]),
            "role": str(source.get("scope", source.get("purpose", ""))),
        }
        verified.append(record)

    method_lock_record = gene_manifest.get("method_lock", {})
    method_lock_path = ROOT / str(method_lock_record.get("path", METHOD_LOCK_PATH))
    if method_lock_path.resolve() != METHOD_LOCK_PATH.resolve():
        raise MappingCheckError("gene-map manifest points to a different method lock")

    # The count response is intentionally checked as content only after its
    # byte hash was authenticated above.
    count_record = manifest_sources.get("ensembl116_human_gene_crosswalk_count")
    if count_record is None:
        raise MappingCheckError("gene-map manifest lacks the Ensembl count response")
    count_text = (
        (ROOT / "data/raw/gse84161/gene-map/ensembl116-biomart-count.tsv")
        .read_text(encoding="utf-8")
        .strip()
    )
    if count_text != str(
        gene_manifest["completeness_check"]["count_response_unique_ensembl_gene_ids"]
    ):
        raise MappingCheckError("Ensembl count response content differs from manifest")

    crosswalk_path = ROOT / "data/raw/gse84161/gene-map/ensembl116-human-gene-crosswalk.tsv"
    crosswalk_text = crosswalk_path.read_text(encoding="utf-8")
    if not crosswalk_text.rstrip("\r\n").endswith("[success]"):
        raise MappingCheckError("Ensembl crosswalk is missing the terminal success stamp")

    return method_lock, gene_manifest, {
        "verified_inputs": verified,
        "manifest_sources": manifest_sources,
        "method_lock_sha256": lock_hash_record["sha256"],
    }


def _resolve_zip_member(stream: zipfile.ZipFile, logical_member: str) -> str:
    exact = [name for name in stream.namelist() if name == logical_member]
    if len(exact) == 1:
        return exact[0]
    suffix = [name for name in stream.namelist() if name.endswith("/" + logical_member)]
    if len(suffix) != 1:
        raise MappingCheckError(
            f"expected one ZIP member for {logical_member!r}, found {len(suffix)}"
        )
    return suffix[0]


def _read_program_column(
    zip_path: Path, logical_member: str, column: str
) -> tuple[list[str], dict[str, Any]]:
    with zipfile.ZipFile(zip_path) as stream:
        actual_member = _resolve_zip_member(stream, logical_member)
        content = stream.read(actual_member)
    digest = hashlib.sha256(content).hexdigest()
    rows: list[str] = []
    reader = csv.reader(io.StringIO(content.decode("utf-8", errors="strict")))
    try:
        header = next(reader)
    except StopIteration as error:
        raise MappingCheckError(f"empty source member: {logical_member}") from error
    if len(header) != len(set(header)) or column not in header:
        raise MappingCheckError(f"invalid or missing source column {column!r}")
    index = header.index(column)
    for row_number, row in enumerate(reader, start=2):
        if len(row) != len(header):
            raise MappingCheckError(
                f"{logical_member} row {row_number} has {len(row)} fields, expected {len(header)}"
            )
        value = row[index].strip()
        if value:
            rows.append(value)
    return rows, {
        "member": actual_member,
        "bytes": len(content),
        "sha256": digest,
        "header_columns": len(header),
        "nonblank_rows": len(rows),
    }


def _parse_gpl(path: Path, expectations: dict[str, Any]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Parse only the exact GPL570 platform table from the pinned SOFT."""

    rows: list[dict[str, str]] = []
    in_platform = False
    ended = False
    header: list[str] | None = None
    with path.open("r", encoding="utf-8", errors="strict", newline="") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if line == "!platform_table_begin":
                if in_platform:
                    raise MappingCheckError("nested GPL platform table")
                in_platform = True
                continue
            if not in_platform:
                continue
            if line == "!platform_table_end":
                ended = True
                break
            if header is None:
                if not line:
                    continue
                header = line.split("\t")
                continue
            fields = line.split("\t")
            if len(fields) != len(header):
                raise MappingCheckError(
                    f"GPL row has {len(fields)} fields, expected {len(header)}"
                )
            rows.append(
                {
                    "probe_id": fields[header.index("ID")],
                    "entrez_field": fields[header.index("ENTREZ_GENE_ID")],
                }
            )
    if not ended or header is None:
        raise MappingCheckError("unterminated or missing GPL platform table")
    expected_rows = int(expectations["probe_rows"])
    if len(rows) != expected_rows:
        raise MappingCheckError(f"GPL row count {len(rows)} != {expected_rows}")
    probe_ids = [row["probe_id"] for row in rows]
    if len(set(probe_ids)) != len(probe_ids):
        raise MappingCheckError("GPL probe IDs are not unique")
    return rows, {
        "header": header,
        "probe_rows": len(rows),
        "unique_probe_ids": len(set(probe_ids)),
        "expected_platform": expectations.get("platform"),
        "annotation_date": expectations.get("annotation_date"),
    }


def _parse_crosswalk(path: Path, manifest: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build full-release relationships before applying any array restriction."""

    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    if not lines or lines[0].split("\t") != [
        "Gene stable ID",
        "NCBI gene (formerly Entrezgene) ID",
        "HGNC symbol",
    ]:
        raise MappingCheckError("unexpected Ensembl crosswalk header")
    if lines[-1] != "[success]":
        raise MappingCheckError("unexpected Ensembl crosswalk footer")

    rows: list[tuple[str, str, str]] = []
    for row_number, line in enumerate(lines[1:-1], start=2):
        fields = line.split("\t")
        if len(fields) != 3:
            raise MappingCheckError(
                f"crosswalk row {row_number} has {len(fields)} fields"
            )
        ensembl, entrez, symbol = (field.strip() for field in fields)
        if not ensembl or not re.fullmatch(r"ENSG\d+(?:\.\d+)?", ensembl):
            raise MappingCheckError(f"invalid Ensembl ID at crosswalk row {row_number}")
        if entrez and not re.fullmatch(r"\d+", entrez):
            raise MappingCheckError(f"invalid Entrez ID at crosswalk row {row_number}")
        rows.append((ensembl.split(".", 1)[0], entrez, symbol))

    expected_rows = int(manifest["completeness_check"]["observed_relationship_rows"])
    expected_ensembl = int(manifest["completeness_check"]["observed_unique_ensembl_gene_ids"])
    if len(rows) != expected_rows:
        raise MappingCheckError(f"crosswalk row count {len(rows)} != {expected_rows}")
    all_ensembl = {row[0] for row in rows}
    if len(all_ensembl) != expected_ensembl:
        raise MappingCheckError(
            f"crosswalk unique Ensembl count {len(all_ensembl)} != {expected_ensembl}"
        )

    relation_rows = set(rows)
    ensembl_to_entrez: dict[str, set[str]] = defaultdict(set)
    entrez_to_ensembl: dict[str, set[str]] = defaultdict(set)
    symbol_to_ensembl: dict[str, set[str]] = defaultdict(set)
    for ensembl, entrez, symbol in rows:
        if entrez:
            ensembl_to_entrez[ensembl].add(entrez)
            entrez_to_ensembl[entrez].add(ensembl)
        if symbol:
            symbol_to_ensembl[symbol].add(ensembl)

    reciprocal_ensembl = {
        ensembl
        for ensembl in all_ensembl
        if len(ensembl_to_entrez.get(ensembl, set())) == 1
        and len(
            entrez_to_ensembl[
                next(iter(ensembl_to_entrez[ensembl]))
            ]
        )
        == 1
    }
    crosswalk_summary = {
        "relationship_rows": len(rows),
        "unique_exact_relationship_rows": len(relation_rows),
        "duplicate_relationship_rows": len(rows) - len(relation_rows),
        "unique_ensembl_gene_ids": len(all_ensembl),
        "ensembl_rows_with_entrez": sum(bool(value) for value in ensembl_to_entrez.values()),
        "ensembl_without_entrez": len(all_ensembl - set(ensembl_to_entrez)),
        "unique_entrez_ids": len(entrez_to_ensembl),
        "reciprocal_one_to_one_ensembl_ids": len(reciprocal_ensembl),
        "symbols_with_one_ensembl": sum(
            len(value) == 1 for value in symbol_to_ensembl.values()
        ),
        "symbols_with_multiple_ensembl": sum(
            len(value) > 1 for value in symbol_to_ensembl.values()
        ),
        "unique_hgnc_symbols": len(symbol_to_ensembl),
    }
    return {
        "ensembl_to_entrez": ensembl_to_entrez,
        "entrez_to_ensembl": entrez_to_ensembl,
        "symbol_to_ensembl": symbol_to_ensembl,
        "all_ensembl": all_ensembl,
        "reciprocal_ensembl": reciprocal_ensembl,
    }, crosswalk_summary


def _probe_entrez_mapping(
    gpl_rows: list[dict[str, str]],
) -> tuple[set[str], dict[str, Any]]:
    eligible_entrez: set[str] = set()
    status_counts: Counter[str] = Counter()
    probes_per_entrez: Counter[str] = Counter()
    examples: dict[str, list[str]] = defaultdict(list)
    for row in gpl_rows:
        tokens = [token.strip() for token in row["entrez_field"].split("///")]
        numeric = {token for token in tokens if re.fullmatch(r"\d+", token)}
        unresolved = [
            token for token in tokens if token and not re.fullmatch(r"\d+", token)
        ]
        if len(numeric) == 1:
            status = "single_unique_entrez"
            entrez = next(iter(numeric))
            eligible_entrez.add(entrez)
            probes_per_entrez[entrez] += 1
        elif len(numeric) == 0:
            status = "blank_or_unresolved_entrez"
            entrez = ""
        else:
            status = "multiple_unique_entrez"
            entrez = ""
        status_counts[status] += 1
        if unresolved:
            status_counts["unresolved_tokens_discarded"] += 1
        if len(examples[status]) < 5:
            examples[status].append(row["probe_id"])
    expected = {
        "probe_rows": 54675,
        "single_unique_entrez": 41834,
        "blank_or_unresolved_entrez": 10541,
        "multiple_unique_entrez": 2300,
        "unique_eligible_entrez": 20486,
    }
    observed = {
        "probe_rows": len(gpl_rows),
        "single_unique_entrez": status_counts["single_unique_entrez"],
        "blank_or_unresolved_entrez": status_counts["blank_or_unresolved_entrez"],
        "multiple_unique_entrez": status_counts["multiple_unique_entrez"],
        "unique_eligible_entrez": len(eligible_entrez),
        "probes_with_unresolved_tokens_discarded": status_counts[
            "unresolved_tokens_discarded"
        ],
        "eligible_probe_count": sum(probes_per_entrez.values()),
        "probe_count_distribution": dict(
            Counter(probes_per_entrez.values())
        ),
        "status_examples": dict(examples),
    }
    if any(observed[key] != value for key, value in expected.items()):
        raise MappingCheckError(f"GPL Entrez accounting differs: {observed} != {expected}")
    return eligible_entrez, {"expected": expected, "observed": observed}


def _id_digest(values: set[str] | list[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in sorted(set(values))).encode("utf-8")
    ).hexdigest()


def _map_program(
    program_id: str,
    namespace: str,
    source_values: list[str],
    eligible_array_entrez: set[str],
    crosswalk: dict[str, Any],
) -> tuple[dict[str, Any], set[str]]:
    """Map one source column and retain one reason for every intended member."""

    source_keys: list[str] = []
    seen: set[str] = set()
    for value in source_values:
        key = (
            value.split(".", 1)[0]
            if namespace == "ensembl_gene_id" and re.fullmatch(
                r"ENSG\d+(?:\.\d+)?", value
            )
            else value
        )
        if key not in seen:
            seen.add(key)
            source_keys.append(key)
    source_before_exclusion = len(source_keys)
    excluded: list[str] = []
    intended: list[str] = []
    for value in source_keys:
        if value == "ENSG00000157557" or value == "ETS2":
            excluded.append(value)
        else:
            intended.append(value)

    ensembl_to_entrez = crosswalk["ensembl_to_entrez"]
    entrez_to_ensembl = crosswalk["entrez_to_ensembl"]
    symbol_to_ensembl = crosswalk["symbol_to_ensembl"]
    all_ensembl = crosswalk["all_ensembl"]
    mapped: set[str] = set()
    statuses: Counter[str] = Counter()
    examples: dict[str, list[str]] = defaultdict(list)

    for value in intended:
        status = ""
        ensembl = ""
        if namespace == "ensembl_gene_id":
            if value not in all_ensembl:
                status = "missing_ensembl_id"
            else:
                ensembl = value
        else:
            candidates = symbol_to_ensembl.get(value, set())
            if not candidates:
                status = "missing_hgnc_symbol"
            elif len(candidates) != 1:
                status = "ambiguous_hgnc_symbol"
            else:
                ensembl = next(iter(candidates))

        if not status:
            entrez_candidates = ensembl_to_entrez.get(ensembl, set())
            if not entrez_candidates:
                status = "missing_entrez_id"
            elif len(entrez_candidates) != 1:
                status = "ambiguous_ensembl_entrez"
            else:
                entrez = next(iter(entrez_candidates))
                if len(entrez_to_ensembl.get(entrez, set())) != 1:
                    status = "ambiguous_reciprocal_entrez"
                elif entrez not in eligible_array_entrez:
                    status = "no_eligible_gpl_probe"
                else:
                    status = "mapped"
                    mapped.add(entrez)
        statuses[status] += 1
        if len(examples[status]) < 10:
            examples[status].append(value)

    intended_count = len(intended)
    mapped_count = len(mapped)
    fraction = mapped_count / intended_count if intended_count else 0.0
    record = {
        "source_namespace": namespace,
        "source_unique_before_exclusion": source_before_exclusion,
        "excluded_target_gene_count": len(excluded),
        "excluded_target_values": excluded,
        "intended_source_members": intended_count,
        "mapped_unique_entrez_genes": mapped_count,
        "mapping_fraction": fraction,
        "mapping_gate_pass": bool(fraction >= 0.8 and mapped_count >= 10),
        "status_counts": dict(sorted(statuses.items())),
        "status_examples": dict(sorted(examples.items())),
        "mapped_entrez_sha256": _id_digest(mapped),
    }
    return record, mapped


def _find_primary_output(explicit: Path | None) -> Path | None:
    if explicit is not None:
        return explicit if explicit.exists() else None
    for path in PRIMARY_MAPPING_CANDIDATES:
        if path.exists():
            return path
    return None


def _compare_primary(
    path: Path | None,
    independent: dict[str, dict[str, Any]],
    derived: dict[str, Any],
) -> dict[str, Any]:
    if path is None:
        return {
            "status": "primary_mapping_output_not_found",
            "searched_paths": [_relative(candidate) for candidate in PRIMARY_MAPPING_CANDIDATES],
        }
    byte_count, digest = sha256_file(path)
    result: dict[str, Any] = {
        "status": "candidate_found_needs_structural_read",
        "path": _relative(path),
        "bytes": byte_count,
        "sha256": digest,
    }
    # Keep this comparison deliberately conservative: the production output
    # format may be finalized separately.  Recognize a simple JSON map or CSV
    # table when it exposes the frozen program/count fields.
    primary_records: dict[str, dict[str, Any]] = {}
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        entries = value.get("programs") or value.get("program_coverage") or []
        if isinstance(entries, dict):
            entries = [
                {"program_id": key, **(item if isinstance(item, dict) else {})}
                for key, item in entries.items()
            ]
        if isinstance(entries, list):
            for item in entries:
                if isinstance(item, dict) and item.get("program_id"):
                    primary_records[str(item["program_id"])] = item
    elif path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                if row.get("program_id"):
                    primary_records[str(row["program_id"])] = row
    if not primary_records:
        return result
    mismatches: list[dict[str, Any]] = []
    for program_id, expected in independent.items():
        observed = primary_records.get(program_id)
        if observed is None:
            mismatches.append({"program_id": program_id, "field": "presence"})
            continue
        for field in ("intended_source_members", "mapped_unique_entrez_genes"):
            candidate = observed.get(field)
            if candidate is None:
                candidate = observed.get(
                    "intended_unique_source_genes"
                    if field == "intended_source_members"
                    else "mapped_unique_genes"
                )
            if candidate is not None and int(float(candidate)) != int(expected[field]):
                mismatches.append(
                    {
                        "program_id": program_id,
                        "field": field,
                        "independent": expected[field],
                        "primary": candidate,
                    }
                )
    derived_primary = primary_records.get("ets2_g1_dn_without_comparators")
    if derived_primary:
        candidate = derived_primary.get("mapped_unique_entrez_genes") or derived_primary.get(
            "mapped_unique_genes"
        )
        if candidate is not None and int(float(candidate)) != int(
            derived["derived_mapped_unique_entrez_genes"]
        ):
            mismatches.append(
                {
                    "program_id": "ets2_g1_dn_without_comparators",
                    "field": "mapped_unique_entrez_genes",
                    "independent": derived["derived_mapped_unique_entrez_genes"],
                    "primary": candidate,
                }
            )
    result.update(
        {
            "status": "compared",
            "programs_compared": len(primary_records),
            "mismatches": mismatches,
            "counts_match": not mismatches,
        }
    )
    return result


def build_check(primary_output: Path | None = None) -> dict[str, Any]:
    method_lock, gene_manifest, verification = _load_and_verify_sources()
    if method_lock.get("source_member_denominators_after_ETS2_exclusion") != EXPECTED_DENOMINATORS:
        raise MappingCheckError("method-lock source denominators differ from independent expectations")

    gpl_source = verification["manifest_sources"]["gse84161_gpl570_soft"]
    gpl_path = ROOT / gpl_source["path"]
    gpl_rows, gpl_summary = _parse_gpl(
        gpl_path,
        gene_manifest["source_files"][0]["platform_expectations"],
    )
    eligible_entrez, probe_summary = _probe_entrez_mapping(gpl_rows)

    crosswalk_source = verification["manifest_sources"]["ensembl116_human_gene_crosswalk"]
    crosswalk_path = ROOT / crosswalk_source["path"]
    crosswalk, crosswalk_summary = _parse_crosswalk(crosswalk_path, gene_manifest)

    zip_source = verification["manifest_sources"]["ets2_program_archive"]
    zip_path = ROOT / zip_source["path"]
    source_values: dict[str, list[str]] = {}
    source_members: dict[str, dict[str, Any]] = {}
    for program_id, (namespace, member, column) in EXPECTED_SOURCE_SPECS.items():
        values, member_info = _read_program_column(zip_path, member, column)
        source_values[program_id] = values
        source_members.setdefault(member, member_info)

    program_records: dict[str, dict[str, Any]] = {}
    mapped_sets: dict[str, set[str]] = {}
    for program_id, (namespace, _, _) in EXPECTED_SOURCE_SPECS.items():
        record, mapped = _map_program(
            program_id,
            namespace,
            source_values[program_id],
            eligible_entrez,
            crosswalk,
        )
        expected_denominator = EXPECTED_DENOMINATORS[program_id]
        if record["intended_source_members"] != expected_denominator:
            raise MappingCheckError(
                f"{program_id} denominator {record['intended_source_members']} "
                f"!= frozen {expected_denominator}"
            )
        program_records[program_id] = record
        mapped_sets[program_id] = mapped

    comparator_union: set[str] = set()
    for program_id in (
        "inflammation",
        "interferon_gamma",
        "oxidative_stress",
        "apoptotic_signaling",
    ):
        comparator_union.update(mapped_sets[program_id])
    primary_set = mapped_sets["ets2_g1_dn"]
    derived_set = primary_set - comparator_union
    derived = {
        "parent_program_id": "ets2_g1_dn",
        "parent_intended_source_members": EXPECTED_DENOMINATORS["ets2_g1_dn"],
        "parent_mapped_unique_entrez_genes": len(primary_set),
        "comparator_union_mapped_unique_entrez_genes": len(comparator_union),
        "parent_comparator_overlap_mapped_unique_entrez_genes": len(
            primary_set & comparator_union
        ),
        "derived_mapped_unique_entrez_genes": len(derived_set),
        "derived_mapping_fraction": len(derived_set) / EXPECTED_DENOMINATORS["ets2_g1_dn"],
        "parent_mapping_gate_pass": program_records["ets2_g1_dn"]["mapping_gate_pass"],
        "derived_mapping_gate_pass": bool(
            program_records["ets2_g1_dn"]["mapping_gate_pass"] and len(derived_set) >= 10
        ),
        "comparator_union_entrez_sha256": _id_digest(comparator_union),
        "derived_entrez_sha256": _id_digest(derived_set),
    }

    primary_path = _find_primary_output(primary_output)
    return {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "check_id": "gse84161-independent-mapping-check-v1",
        "method_lock": {
            "path": _relative(METHOD_LOCK_PATH),
            "sha256": verification["method_lock_sha256"],
            "method_id": method_lock.get("method_id"),
            "frozen_at_utc": method_lock.get("frozen_at_utc"),
        },
        "verified_input_hashes": verification["verified_inputs"],
        "source_members": source_members,
        "gpl570": {
            "source": gpl_summary,
            "probe_to_entrez": probe_summary,
        },
        "ensembl116": {
            "source": crosswalk_source,
            "relationships": crosswalk_summary,
        },
        "programs": program_records,
        "derived_nonoverlap": derived,
        "primary_mapping_comparison": _compare_primary(
            primary_path, program_records, derived
        ),
        "limits": {
            "expression_values_read": False,
            "cel_files_opened": False,
            "probe_intensities_read": False,
            "program_scores_computed": False,
            "treatment_or_donor_effects_computed": False,
            "mapping_only": True,
            "note": (
                "Only GPL annotation fields, complete Ensembl relationship rows, "
                "and source-program identifiers were parsed. Crosswalk and array "
                "eligibility were resolved without opening CEL or expression files."
            ),
        },
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--primary-output", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    try:
        result = build_check(args.primary_output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
    except (MappingCheckError, OSError, KeyError, json.JSONDecodeError) as error:
        print(f"independent GSE84161 mapping check failed: {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "output": _relative(args.output),
                "output_sha256": sha256_file(args.output)[1],
                "program_count": len(result["programs"]),
                "gpl_eligible_probes": result["gpl570"]["probe_to_entrez"]["observed"][
                    "single_unique_entrez"
                ],
                "gpl_unique_entrez": result["gpl570"]["probe_to_entrez"]["observed"][
                    "unique_eligible_entrez"
                ],
                "primary_mapping_comparison": result["primary_mapping_comparison"][
                    "status"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
