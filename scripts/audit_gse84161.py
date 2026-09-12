#!/usr/bin/env python3
"""Audit GSE84161/GSE84162 metadata and array feature provenance.

This module is deliberately metadata-only.  It streams the cached GEO SOFT
files, keeps sample labels and feature identifiers, and never parses, stores,
or prints expression quantities from the VALUE columns.  The processed table
is described as a source-provided post-selection table; the full GPL570 probe
universe is kept as a separate annotation source.  No program score, target
effect, or sample difference is calculated here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ets2-independent-search"
OUT_WORK = ROOT / "work/gse84161-metadata"
DERIVED = ROOT / "data/derived"
REPORTS = ROOT / "reports"

GSE84161_SOFT = WORK / "GSE84161-family-soft.txt"
GSE84162_SOFT = WORK / "GSE84162-family-soft.txt"
EXPECTED_SOFT = {
    "GSE84161": (99895372, "620c3c8b1cb864cd61605b24e06a9fdc67fb81e6288a2918d22beb792ebf8777"),
    "GSE84162": (113260887, "d22b804584f1ce5dc21728ffd896cb0ef5ea2b0e00d98d16c039e738dcd4ac9a"),
}

GSE84161_SOFT_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    "acc=GSE84161&targ=all&form=text&view=quick"
)
GSE84162_SOFT_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    "acc=GSE84162&targ=all&form=text&view=quick"
)
GSE84161_GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84161"
GSE84162_GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84162"
GSE84161_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84161/"
    "matrix/GSE84161_series_matrix.txt.gz"
)
GSE84161_RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84161/"
    "suppl/GSE84161_RAW.tar"
)
GSE84162_RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84162/"
    "suppl/GSE84162_RAW.tar"
)

SAMPLE_CSV = DERIVED / "gse84161-samples.csv"
PROCESSING_CSV = DERIVED / "gse84161-processing-evidence.csv"
REPORT_JSON = REPORTS / "gse84161-metadata-audit.json"

ARM_ORDER = [
    "untreated_unstim",
    "tpl2i_unstim",
    "meki_unstim",
    "vehicle_lps",
    "tpl2i_lps",
    "meki_lps",
]

ARM_SPECS = {
    "untreated_unstim": {
        "treatment_name": "Untreated",
        "stimulation": "Unstimulated",
        "compound": "Vehicle",
        "control_condition_key": "",
        "control_role": "vehicle control for tpl2i_unstim and meki_unstim",
    },
    "tpl2i_unstim": {
        "treatment_name": "Tpl-2 SMI",
        "stimulation": "Unstimulated",
        "compound": "TPL2i",
        "control_condition_key": "untreated_unstim",
        "control_role": "treated arm matched to untreated_unstim",
    },
    "meki_unstim": {
        "treatment_name": "MEK SMI",
        "stimulation": "Unstimulated",
        "compound": "MEKi",
        "control_condition_key": "untreated_unstim",
        "control_role": "treated arm matched to untreated_unstim",
    },
    "vehicle_lps": {
        "treatment_name": "LPS",
        "stimulation": "LPS",
        "compound": "Vehicle",
        "control_condition_key": "",
        "control_role": "vehicle control for tpl2i_lps and meki_lps",
    },
    "tpl2i_lps": {
        "treatment_name": "Tpl-2 SMI + LPS",
        "stimulation": "LPS",
        "compound": "TPL2i",
        "control_condition_key": "vehicle_lps",
        "control_role": "treated arm matched to vehicle_lps",
    },
    "meki_lps": {
        "treatment_name": "MEK SMI + LPS",
        "stimulation": "LPS",
        "compound": "MEKi",
        "control_condition_key": "vehicle_lps",
        "control_role": "treated arm matched to vehicle_lps",
    },
}

SAMPLE_METADATA_FIELDS = [
    "Sample_title",
    "Sample_source_name_ch1",
    "Sample_characteristics_ch1",
    "Sample_treatment_protocol_ch1",
    "Sample_growth_protocol_ch1",
    "Sample_molecule_ch1",
    "Sample_data_processing",
    "Sample_platform_id",
    "Sample_supplementary_file",
    "Sample_supplementary_file_1",
    "Sample_data_row_count",
    "Sample_type",
    "Sample_instrument_model",
    "Sample_library_strategy",
    "Sample_series_id",
]


def _add_field(record: dict[str, Any], key: str, value: str) -> None:
    fields = record.setdefault("fields", {})
    if key not in fields:
        fields[key] = value
    elif isinstance(fields[key], list):
        fields[key].append(value)
    else:
        fields[key] = [fields[key], value]


def _field_list(record: dict[str, Any], key: str) -> list[str]:
    value = record.get("fields", {}).get(key, [])
    if isinstance(value, list):
        return [str(item) for item in value]
    if value == "":
        return []
    return [str(value)]


def _field(record: dict[str, Any], key: str, default: str = "") -> str:
    values = _field_list(record, key)
    return values[-1] if values else default


def _parse_characteristics(record: dict[str, Any]) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in _field_list(record, "Sample_characteristics_ch1"):
        if ":" in item:
            key, value = item.split(":", 1)
            result[key.strip()] = value.strip()
    return result


def _finalize_table(record: dict[str, Any], state: dict[str, Any]) -> None:
    digest = state["digest"].hexdigest()
    feature_ids = state.get("feature_ids", [])
    first_ids = state.get("first_ids", feature_ids[:3])
    last_ids = state.get("last_ids", feature_ids[-3:] if feature_ids else [])
    info = {
        "kind": state["kind"],
        "header": state.get("header", []),
        "row_count": state["row_count"],
        "unique_id_count": len(state["unique_ids"]),
        "duplicate_id_count": len(state["duplicate_ids"]),
        "duplicate_ids": sorted(state["duplicate_ids"]),
        "id_digest": digest,
        "first_ids": first_ids,
        "last_ids": last_ids,
        "feature_ids": feature_ids,
        "annotation_by_id": state.get("annotation_by_id", {}),
    }
    record["table"] = info


def parse_soft(path: Path) -> dict[str, Any]:
    """Parse SOFT metadata, table headers and feature IDs without values."""

    series: list[dict[str, Any]] = []
    platforms: list[dict[str, Any]] = []
    samples: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    current_kind = ""
    table_state: dict[str, Any] | None = None
    source_order = 0

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            lower = line.lower()

            if table_state is not None:
                if lower == "!platform_table_end" or lower == "!sample_table_end":
                    _finalize_table(table_state["record"], table_state)
                    table_state = None
                    continue
                if table_state["header"] is None:
                    if line:
                        table_state["header"] = line.split("\t")
                    continue
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                feature_id = fields[0]
                table_state["row_count"] += 1
                if feature_id in table_state["unique_ids"]:
                    table_state["duplicate_ids"].add(feature_id)
                table_state["unique_ids"].add(feature_id)
                table_state["digest"].update((feature_id + "\n").encode("utf-8"))
                if len(table_state["feature_ids"]) < 25000:
                    table_state["feature_ids"].append(feature_id)
                if len(table_state["first_ids"]) < 3:
                    table_state["first_ids"].append(feature_id)
                table_state["last_ids"].append(feature_id)
                if len(table_state["last_ids"]) > 3:
                    del table_state["last_ids"][0]
                if table_state["kind"] == "platform":
                    header = table_state["header"]
                    annotation_index = {
                        name: index
                        for index, name in enumerate(header)
                        if name in {"ENTREZ_GENE_ID", "Annotation Date"}
                    }
                    table_state["annotation_by_id"][feature_id] = {
                        name: fields[index] if index < len(fields) else ""
                        for name, index in annotation_index.items()
                    }
                continue

            if lower == "!platform_table_begin":
                if current is None:
                    raise ValueError(f"platform table without record in {path}")
                table_state = {
                    "kind": "platform",
                    "record": current,
                    "header": None,
                    "row_count": 0,
                    "unique_ids": set(),
                    "duplicate_ids": set(),
                    "digest": hashlib.sha256(),
                    "feature_ids": [],
                    "first_ids": [],
                    "last_ids": [],
                    "annotation_by_id": {},
                }
                continue
            if lower == "!sample_table_begin":
                if current is None or current_kind != "sample":
                    raise ValueError(f"sample table without sample in {path}")
                table_state = {
                    "kind": "sample",
                    "record": current,
                    "header": None,
                    "row_count": 0,
                    "unique_ids": set(),
                    "duplicate_ids": set(),
                    "digest": hashlib.sha256(),
                    "feature_ids": [],
                    "first_ids": [],
                    "last_ids": [],
                    "annotation_by_id": {},
                }
                continue

            if line.startswith("^SERIES = "):
                current = {"kind": "series", "accession": line.split("=", 1)[1].strip(), "fields": {}}
                series.append(current)
                current_kind = "series"
                continue
            if line.startswith("^PLATFORM = "):
                current = {
                    "kind": "platform",
                    "accession": line.split("=", 1)[1].strip(),
                    "fields": {},
                }
                platforms.append(current)
                current_kind = "platform"
                continue
            if line.startswith("^SAMPLE = "):
                source_order += 1
                current = {
                    "kind": "sample",
                    "accession": line.split("=", 1)[1].strip(),
                    "fields": {},
                    "comments": {},
                    "source_order": source_order,
                }
                samples.append(current)
                current_kind = "sample"
                continue

            if current is None:
                continue
            if line.startswith("!") and " = " in line:
                key, value = line[1:].split(" = ", 1)
                if key.startswith("Series_") or key.startswith("Platform_") or key.startswith("Sample_"):
                    _add_field(current, key, value)
                continue
            if current_kind == "sample" and line.startswith("#") and " = " in line:
                key, value = line[1:].split(" = ", 1)
                current.setdefault("comments", {})[key.strip()] = value

    if table_state is not None:
        raise ValueError(f"unterminated {table_state['kind']} table in {path}")

    return {
        "path": str(path.relative_to(ROOT)),
        "series": series,
        "platforms": platforms,
        "samples": samples,
    }


def _file_sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _series_summary(record: dict[str, Any]) -> dict[str, Any]:
    fields = record.get("fields", {})
    selected = [
        "Series_title",
        "Series_status",
        "Series_submission_date",
        "Series_last_update_date",
        "Series_overall_design",
        "Series_type",
        "Series_supplementary_file",
        "Series_relation",
        "Series_platform_id",
        "Series_sample_taxid",
    ]
    result: dict[str, Any] = {"accession": record["accession"]}
    for key in selected:
        values = _field_list(record, key)
        if values:
            result[key] = values if len(values) > 1 else values[0]
    result["series_sample_id_count"] = len(_field_list(record, "Series_sample_id"))
    return result


def _condition_key(sample: dict[str, Any]) -> str:
    chars = _parse_characteristics(sample)
    treatment = chars.get("treatment name", "")
    stimulation = chars.get("stimulation", "")
    compound = chars.get("compound", "")
    for key, spec in ARM_SPECS.items():
        if (
            treatment == spec["treatment_name"]
            and stimulation == spec["stimulation"]
            and compound == spec["compound"]
        ):
            return key
    return "unclassified"


def _replicate_number(sample: dict[str, Any]) -> int | None:
    chars = _parse_characteristics(sample)
    value = chars.get("replicate number", "")
    match = re.fullmatch(r"\s*(\d+)\s*", value)
    if match:
        return int(match.group(1))
    match = re.search(r"(?:^|\s)(\d+)\s*$", _field(sample, "Sample_title"))
    return int(match.group(1)) if match else None


def _sample_metadata(sample: dict[str, Any]) -> dict[str, Any]:
    chars = _parse_characteristics(sample)
    table = sample.get("table") or {}
    series_ids = _field_list(sample, "Sample_series_id")
    return {
        "sample_accession": sample["accession"],
        "source_order": sample.get("source_order"),
        "source_title": _field(sample, "Sample_title"),
        "source_name": _field(sample, "Sample_source_name_ch1"),
        "characteristics": chars,
        "characteristics_original": _field_list(sample, "Sample_characteristics_ch1"),
        "series_memberships": series_ids,
        "primary_tissue": chars.get("primary tissue", ""),
        "treatment_name": chars.get("treatment name", ""),
        "stimulation": chars.get("stimulation", ""),
        "compound": chars.get("compound", ""),
        "replicate_number": _replicate_number(sample),
        "condition_key": _condition_key(sample),
        "platform_id": _field(sample, "Sample_platform_id"),
        "source_data_row_count": _field(sample, "Sample_data_row_count"),
        "sample_type": _field(sample, "Sample_type"),
        "instrument_model": _field(sample, "Sample_instrument_model"),
        "library_strategy": _field(sample, "Sample_library_strategy"),
        "molecule": _field(sample, "Sample_molecule_ch1"),
        "data_processing": _field(sample, "Sample_data_processing"),
        "treatment_protocol": _field(sample, "Sample_treatment_protocol_ch1"),
        "growth_protocol": _field(sample, "Sample_growth_protocol_ch1"),
        "supplementary_files": _field_list(sample, "Sample_supplementary_file")
        + _field_list(sample, "Sample_supplementary_file_1"),
        "comments": dict(sample.get("comments", {})),
        "table": {
            key: table.get(key)
            for key in [
                "kind",
                "header",
                "row_count",
                "unique_id_count",
                "duplicate_id_count",
                "duplicate_ids",
                "id_digest",
                "first_ids",
                "last_ids",
            ]
        },
        "feature_ids": list(table.get("feature_ids", [])),
    }


def _processed_consistency(samples: list[dict[str, Any]]) -> dict[str, Any]:
    metadata = [_sample_metadata(sample) for sample in samples]
    tables = [item["table"] for item in metadata]
    reference = tables[0] if tables else {}
    same_digest = all(table.get("id_digest") == reference.get("id_digest") for table in tables)
    same_order_examples = all(
        table.get("first_ids") == reference.get("first_ids")
        and table.get("last_ids") == reference.get("last_ids")
        for table in tables
    )
    metadata_counts = sorted({item["source_data_row_count"] for item in metadata})
    return {
        "sample_count": len(samples),
        "table_presence_count": sum(bool(item.get("header")) for item in tables),
        "processed_table_header": reference.get("header", []),
        "idref_namespace": (
            "Affymetrix probe-set ID_REF"
            if reference.get("header", [""])[0:1] == ["ID_REF"]
            else "unresolved"
        ),
        "value_label": samples[0].get("comments", {}).get("VALUE", "") if samples else "",
        "processed_row_count": reference.get("row_count"),
        "processed_unique_id_count": reference.get("unique_id_count"),
        "processed_duplicate_id_count": reference.get("duplicate_id_count"),
        "processed_duplicate_ids": reference.get("duplicate_ids", []),
        "processed_id_digest": reference.get("id_digest"),
        "all_samples_same_ordered_ids": bool(same_digest and same_order_examples),
        "all_samples_same_digest": same_digest,
        "all_sample_metadata_row_counts": metadata_counts,
        "all_sample_table_row_counts": sorted({table.get("row_count") for table in tables}),
        "first_ids": reference.get("first_ids", []),
        "last_ids": reference.get("last_ids", []),
    }


def _entrez_parts(value: str) -> list[str]:
    if value.strip() in {"", "---", "NA", "N/A"}:
        return []
    parts: list[str] = []
    seen: set[str] = set()
    for raw_part in value.split("///"):
        part = raw_part.strip()
        if not part or part in {"---", "NA", "N/A"} or part in seen:
            continue
        seen.add(part)
        parts.append(part)
    return parts


def _annotation_stats(feature_ids: Iterable[str], platform_table: dict[str, Any]) -> dict[str, Any]:
    annotation_by_id = platform_table.get("annotation_by_id", {})
    raw_values = [annotation_by_id.get(feature_id, {}).get("ENTREZ_GENE_ID", "") for feature_id in feature_ids]
    status_counts = Counter()
    entrez_counts: Counter[str] = Counter()
    single_id_genes: Counter[str] = Counter()
    for value in raw_values:
        parts = _entrez_parts(value)
        if not parts:
            status_counts["blank_or_unresolved"] += 1
        elif len(parts) == 1:
            status_counts["single_entrez_field"] += 1
            single_id_genes[parts[0]] += 1
        else:
            status_counts["multiple_entrez_field"] += 1
        for part in parts:
            entrez_counts[part] += 1
    return {
        "probe_rows_considered": len(raw_values),
        "blank_or_unresolved_entrez_fields": status_counts["blank_or_unresolved"],
        "single_entrez_fields": status_counts["single_entrez_field"],
        "multiple_entrez_fields": status_counts["multiple_entrez_field"],
        "unique_genes_from_single_id_probes": len(single_id_genes),
        "genes_with_multiple_single_id_probes": sum(count > 1 for count in single_id_genes.values()),
        "total_entrez_assignments": sum(entrez_counts.values()),
        "unique_entrez_assignments": len(entrez_counts),
        "entrez_ids_assigned_to_more_than_one_probe": sum(
            count > 1 for count in entrez_counts.values()
        ),
        "probe_assignments_in_multi_probe_entrez_groups": sum(
            count for count in entrez_counts.values() if count > 1
        ),
        "max_probe_assignments_per_entrez": max(entrez_counts.values(), default=0),
    }


def _platform_summary(platform: dict[str, Any]) -> dict[str, Any]:
    table = platform.get("table") or {}
    annotation_dates = Counter(
        value.get("Annotation Date", "") for value in table.get("annotation_by_id", {}).values()
    )
    result = {
        "accession": platform["accession"],
        "title": _field(platform, "Platform_title"),
        "status": _field(platform, "Platform_status"),
        "submission_date": _field(platform, "Platform_submission_date"),
        "last_update_date": _field(platform, "Platform_last_update_date"),
        "technology": _field(platform, "Platform_technology"),
        "organism": _field(platform, "Platform_organism"),
        "table_present": bool(table),
        "table_header": table.get("header", []),
        "table_row_count": table.get("row_count", 0),
        "table_unique_id_count": table.get("unique_id_count", 0),
        "table_duplicate_id_count": table.get("duplicate_id_count", 0),
        "table_first_ids": table.get("first_ids", []),
        "table_last_ids": table.get("last_ids", []),
        "annotation_date_counts": dict(annotation_dates),
        "annotation_columns": [
            column
            for column in table.get("header", [])
            if column
            in {
                "ENTREZ_GENE_ID",
                "Gene Symbol",
                "RefSeq Transcript ID",
                "Representative Public ID",
                "Annotation Date",
            }
        ],
    }
    return result


def _build_control_map(metadata: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    by_replicate: dict[int, dict[str, dict[str, Any]]] = defaultdict(dict)
    errors: list[str] = []
    for item in metadata:
        replicate = item.get("replicate_number")
        condition = item.get("condition_key")
        if replicate is None:
            errors.append(f"missing replicate number for {item['sample_accession']}")
            continue
        if condition not in ARM_SPECS:
            errors.append(f"unclassified arm {condition!r} for {item['sample_accession']}")
            continue
        if condition in by_replicate[replicate]:
            errors.append(f"duplicate arm {condition} in replicate {replicate}")
        by_replicate[replicate][condition] = item

    groups: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    for replicate in sorted(by_replicate):
        arms = by_replicate[replicate]
        missing = [key for key in ARM_ORDER if key not in arms]
        if missing:
            errors.append(f"replicate {replicate} missing arms: {','.join(missing)}")
        arm_rows: list[dict[str, Any]] = []
        for condition in ARM_ORDER:
            item = arms.get(condition)
            if item is None:
                continue
            spec = ARM_SPECS[condition]
            arm_rows.append(
                {
                    "condition_key": condition,
                    "sample_accession": item["sample_accession"],
                    "source_title": item["source_title"],
                    "treatment_name": item["treatment_name"],
                    "stimulation": item["stimulation"],
                    "compound": item["compound"],
                    "control_condition_key": spec["control_condition_key"],
                    "control_role": spec["control_role"],
                }
            )
            if spec["control_condition_key"]:
                control = arms.get(spec["control_condition_key"])
                if control is None:
                    errors.append(
                        f"missing control {spec['control_condition_key']} for {condition} in replicate {replicate}"
                    )
                else:
                    pairs.append(
                        {
                            "replicate_number": replicate,
                            "replicate_group_label": f"replicate_{replicate}",
                            "treated_condition_key": condition,
                            "treated_sample_accession": item["sample_accession"],
                            "control_condition_key": spec["control_condition_key"],
                            "control_sample_accession": control["sample_accession"],
                            "status": "planned_unscored",
                        }
                    )
        groups.append(
            {
                "replicate_number": replicate,
                "replicate_group_label": f"replicate_{replicate}",
                "sample_level_donor_id": None,
                "donor_identity_status": "replicate number only; no sample-level donor identifier",
                "arm_count": len(arm_rows),
                "arms": arm_rows,
            }
        )
    return groups, pairs, errors


def _source_provenance(path: Path, url: str) -> dict[str, Any]:
    size, sha256 = _file_sha256(path)
    return {
        "path": str(path.relative_to(ROOT)),
        "url": url,
        "bytes": size,
        "sha256": sha256,
        "cached": True,
        "values_read": False,
    }


def _neutrophil_context(metadata: list[dict[str, Any]]) -> dict[str, Any]:
    neutrophils = [
        item
        for item in metadata
        if "neutrophil" in item.get("primary_tissue", "").lower()
    ]
    by_donor: dict[str, list[dict[str, Any]]] = defaultdict(list)
    donor_ids: dict[str, str] = {}
    for item in neutrophils:
        chars = item.get("characteristics", {})
        donor_tag = chars.get("patient number", "")
        patient_id = chars.get("inventory patient id", "")
        by_donor[donor_tag].append(item)
        donor_ids[donor_tag] = patient_id
    arms_by_donor = {
        donor: sorted(item["treatment_name"] for item in items)
        for donor, items in sorted(by_donor.items())
    }
    return {
        "accession": "GSE84153",
        "superseries_accession": "GSE84162",
        "cell_context": "purified human neutrophils",
        "sample_count": len(neutrophils),
        "explicit_donor_tags": sorted(by_donor),
        "explicit_donor_count": len(by_donor),
        "pseudonymous_inventory_patient_ids": donor_ids,
        "samples_per_donor": {donor: len(items) for donor, items in sorted(by_donor.items())},
        "arms_per_donor": arms_by_donor,
        "assay_platforms": sorted({item["platform_id"] for item in neutrophils}),
        "sample_table_rows": sorted({item["table"].get("row_count") for item in neutrophils}),
        "processed_file_format": "per-sample tab-delimited files with Entrez Gene IDs, gene exonic counts, gene length and RPKM columns",
        "source_treatment_labels_preserved": True,
        "substitution_for_monocytes": False,
        "interpretation": "inventory context only; neutrophil records are not monocyte donor evidence",
    }


def _compare_series(
    subseries: dict[str, Any], superseries: dict[str, Any], sub_metadata: list[dict[str, Any]], super_metadata: list[dict[str, Any]]
) -> dict[str, Any]:
    sub_by_id = {item["sample_accession"]: item for item in sub_metadata}
    super_by_id = {item["sample_accession"]: item for item in super_metadata}
    common_ids = sorted(set(sub_by_id) & set(super_by_id))
    differences: list[dict[str, Any]] = []
    compare_keys = [
        "source_title",
        "source_name",
        "characteristics_original",
        "treatment_protocol",
        "growth_protocol",
        "molecule",
        "data_processing",
        "platform_id",
        "source_data_row_count",
        "sample_type",
        "instrument_model",
        "library_strategy",
        "series_memberships",
        "supplementary_files",
    ]
    for sample_id in common_ids:
        for key in compare_keys:
            if sub_by_id[sample_id].get(key) != super_by_id[sample_id].get(key):
                differences.append({"sample_accession": sample_id, "field": key})
        if sub_by_id[sample_id].get("table", {}).get("id_digest") != super_by_id[sample_id].get("table", {}).get("id_digest"):
            differences.append({"sample_accession": sample_id, "field": "processed_id_digest"})

    super_only = sorted(set(super_by_id) - set(sub_by_id))
    super_only_metadata = [super_by_id[sample_id] for sample_id in super_only]
    neutrophils = _neutrophil_context(super_only_metadata)
    super_order = [item["sample_accession"] for item in super_metadata]
    numeric_gsms = [int(item["sample_accession"].replace("GSM", "")) for item in super_metadata]
    gaps = sorted({right - left for left, right in zip(numeric_gsms, numeric_gsms[1:])})
    return {
        "subseries_accession": subseries["accession"],
        "superseries_accession": superseries["accession"],
        "subseries_sample_count": len(sub_metadata),
        "superseries_sample_count": len(super_metadata),
        "common_monocyte_sample_count": len(common_ids),
        "subseries_only_sample_count": len(set(sub_by_id) - set(super_by_id)),
        "superseries_only_sample_count": len(super_only),
        "common_sample_ids_unique": len(common_ids) == len(set(common_ids)),
        "common_sample_ids": common_ids,
        "monocyte_metadata_exact_match": not differences,
        "monocyte_metadata_difference_fields": differences,
        "superseries_only_neutrophil_context": neutrophils,
        "superseries_sample_order": super_order,
        "superseries_numeric_gsm_step_values": gaps,
        "superseries_contains_nonsequential_gsm_block": any(gap != 1 for gap in gaps),
        "superseries_relation": _field_list(superseries, "Series_relation"),
    }


def _processing_rows(
    series: dict[str, Any],
    samples: list[dict[str, Any]],
    processed: dict[str, Any],
    platform: dict[str, Any],
    platform_summary: dict[str, Any],
    processed_annotation: dict[str, Any],
    full_annotation: dict[str, Any],
    gse84161_raw_url: str,
    cel_urls: list[str],
) -> list[dict[str, str]]:
    processing_text = samples[0].get("data_processing", "") if samples else ""
    heatmap_text = ""
    if "For visualization using heat maps" in processing_text:
        heatmap_text = processing_text.split("For visualization using heat maps", 1)[1].split("P-values", 1)[0].strip()
    design = _field(series, "Series_overall_design")
    rows = [
        {
            "evidence_id": "series_overall_design",
            "source_accession": "GSE84161",
            "scope": "series metadata",
            "field": "Series_overall_design",
            "observed_value": design,
            "interpretation": "Source states five total buffy-coat donors and describes the six-arm experiment.",
            "status": "source_preserved",
        },
        {
            "evidence_id": "sample_processing_rma",
            "source_accession": "GSE84161",
            "scope": "sample metadata",
            "field": "Sample_data_processing",
            "observed_value": processing_text,
            "interpretation": "Source describes RMA, probe filtering, highest-variance single probe per Entrez and a donor covariate; it does not define the variance-estimation sample scope.",
            "status": "source_preserved_scope_ambiguous",
        },
        {
            "evidence_id": "sample_processing_heatmap",
            "source_accession": "GSE84161",
            "scope": "sample metadata",
            "field": "heatmap_residual_normalization_text",
            "observed_value": heatmap_text,
            "interpretation": "The source describes additional donor-residual normalization for heat-map visualization; its applicability to the deposited VALUE table is not assumed.",
            "status": "source_preserved_not_assumed",
        },
        {
            "evidence_id": "processed_header",
            "source_accession": "GSE84161",
            "scope": "all 30 sample tables",
            "field": "ID_REF/VALUE header",
            "observed_value": "|".join(processed["processed_table_header"]),
            "interpretation": "ID_REF is an Affymetrix probe-set identifier and VALUE is labeled RMA normalized signal intensity.",
            "status": "identifier_header_only",
        },
        {
            "evidence_id": "processed_subset_dimensions",
            "source_accession": "GSE84161",
            "scope": "all 30 sample tables",
            "field": "processed table rows",
            "observed_value": f"{processed['processed_row_count']} rows; {processed['processed_unique_id_count']} unique IDs; {processed['processed_duplicate_id_count']} duplicate IDs; all 30 ordered ID lists identical",
            "interpretation": "The deposited table is a shared 19,944-probe subset rather than the 54,675-probe GPL570 universe.",
            "status": "structural_identifier_check_pass",
        },
        {
            "evidence_id": "processed_annotation_snapshot",
            "source_accession": "GSE84161/GPL570",
            "scope": "processed 19,944 probes mapped through current GPL snapshot",
            "field": "ENTREZ_GENE_ID annotation categories",
            "observed_value": json.dumps(processed_annotation, sort_keys=True),
            "interpretation": "Multiple Entrez annotations use the source delimiter ///; whole-field uniqueness is not treated as one-to-one mapping.",
            "status": "current_snapshot_not_historical_proof",
        },
        {
            "evidence_id": "full_gpl_dimensions",
            "source_accession": "GPL570",
            "scope": "full platform table",
            "field": "GPL570 platform table",
            "observed_value": f"{platform_summary['table_row_count']} rows; {platform_summary['table_unique_id_count']} unique probe IDs; {platform_summary['table_duplicate_id_count']} duplicate IDs",
            "interpretation": "This is the full current GPL570 annotation snapshot and must remain separate from the processed post-selection table.",
            "status": "identifier_annotation_only",
        },
        {
            "evidence_id": "full_gpl_annotation_columns",
            "source_accession": "GPL570",
            "scope": "full platform table",
            "field": "mapping columns and dates",
            "observed_value": json.dumps(
                {
                    "columns": platform_summary["annotation_columns"],
                    "annotation_date_counts": platform_summary["annotation_date_counts"],
                    "platform_last_update_date": platform_summary["last_update_date"],
                    "full_annotation": full_annotation,
                },
                sort_keys=True,
            ),
            "interpretation": "Current GPL570 annotation is dated Oct 6, 2014 in the table and the platform record was last updated Dec 14, 2020; this may differ from the historical analysis snapshot.",
            "status": "source_preserved_reproducibility_limit",
        },
        {
            "evidence_id": "cel_availability",
            "source_accession": "GSE84161",
            "scope": "30 monocyte samples",
            "field": "raw CEL routes",
            "observed_value": f"series archive {gse84161_raw_url}; {len(cel_urls)}/30 sample-level CEL URLs advertised",
            "interpretation": "Public raw-array routes are advertised for a future array-specific freeze; this audit did not download or inspect CEL contents.",
            "status": "availability_metadata_only",
        },
        {
            "evidence_id": "timing_specificity_difference",
            "source_accession": "GSE84161",
            "scope": "series versus sample protocol",
            "field": "drug pretreatment duration",
            "observed_value": "Series design: 1 hr; sample treatment protocol: at least 30 min",
            "interpretation": "Preserve both source statements; the one-hour design statement and the at-least-30-minute sample protocol are compatible but differ in specificity.",
            "status": "source_preserved_specificity_difference",
        },
        {
            "evidence_id": "donor_mapping_ambiguity",
            "source_accession": "GSE84161",
            "scope": "30 monocyte samples",
            "field": "donor labels",
            "observed_value": "Series states total of 5 donors; sample records expose replicate number 1-5 and no donor identifier",
            "interpretation": "Replicate groups are retained as public pseudonymous grouping labels only; donor identity is not inferred from expression or numbering.",
            "status": "donor_pairing_not_fully_auditable",
        },
    ]
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sample_rows(
    metadata: list[dict[str, Any]],
    processed: dict[str, Any],
    series_raw_url: str,
) -> list[dict[str, Any]]:
    by_replicate_condition = {
        (item["replicate_number"], item["condition_key"]): item
        for item in metadata
    }
    rows: list[dict[str, Any]] = []
    for item in sorted(metadata, key=lambda value: value["source_order"]):
        replicate = item["replicate_number"]
        condition = item["condition_key"]
        spec = ARM_SPECS.get(condition, {})
        control_condition = spec.get("control_condition_key", "")
        control_item = by_replicate_condition.get((replicate, control_condition)) if control_condition else None
        control_for = ""
        if condition == "untreated_unstim":
            control_for = "tpl2i_unstim;meki_unstim"
        elif condition == "vehicle_lps":
            control_for = "tpl2i_lps;meki_lps"
        cel_urls = [url for url in item["supplementary_files"] if ".CEL" in url.upper()]
        rows.append(
            {
                "sample_accession": item["sample_accession"],
                "source_order": item["source_order"],
                "source_title": item["source_title"],
                "source_name": item["source_name"],
                "series_memberships": ";".join(item["series_memberships"]),
                "primary_tissue": item["primary_tissue"],
                "treatment_name": item["treatment_name"],
                "stimulation": item["stimulation"],
                "compound": item["compound"],
                "replicate_number": replicate,
                "replicate_group_label": f"replicate_{replicate}",
                "donor_identity_status": "not provided at sample level; replicate number is not donor ID",
                "condition_key": condition,
                "control_condition_key": control_condition,
                "control_for_condition_keys": control_for,
                "matched_control_accession": control_item["sample_accession"] if control_item else "",
                "control_role": spec.get("control_role", "unclassified"),
                "source_characteristics": "; ".join(item["characteristics_original"]),
                "assay_platform": item["platform_id"],
                "quantity_label": item["comments"].get("VALUE", ""),
                "processed_table_header": "|".join(item["table"].get("header", [])),
                "processed_table_row_count": item["table"].get("row_count", ""),
                "processed_unique_id_count": item["table"].get("unique_id_count", ""),
                "processed_duplicate_id_count": item["table"].get("duplicate_id_count", ""),
                "processed_id_digest": item["table"].get("id_digest", ""),
                "processed_ids_match_all_samples": str(processed["all_samples_same_ordered_ids"]).lower(),
                "source_data_row_count": item["source_data_row_count"],
                "cel_available_route": str(bool(cel_urls)).lower(),
                "cel_url": cel_urls[0] if cel_urls else "",
                "series_raw_archive_url": series_raw_url,
            }
        )
    return rows


def build_audit(
    *,
    write_outputs: bool = True,
    gse84161_path: Path = GSE84161_SOFT,
    gse84162_path: Path = GSE84162_SOFT,
) -> dict[str, Any]:
    """Build the metadata audit and optionally write the requested artifacts."""

    for accession, path in (("GSE84161", gse84161_path), ("GSE84162", gse84162_path)):
        if _file_sha256(path) != EXPECTED_SOFT[accession]:
            raise ValueError(f"Source bytes changed for {accession}; review a new snapshot separately")
    parsed_84161 = parse_soft(gse84161_path)
    parsed_84162 = parse_soft(gse84162_path)
    series_84161 = parsed_84161["series"][0]
    series_84162 = parsed_84162["series"][0]
    samples_84161 = [_sample_metadata(sample) for sample in parsed_84161["samples"]]
    samples_84162 = [_sample_metadata(sample) for sample in parsed_84162["samples"]]
    platform_84161 = next(
        platform for platform in parsed_84161["platforms"] if platform["accession"] == "GPL570"
    )
    platform_summary = _platform_summary(platform_84161)
    processed = _processed_consistency(parsed_84161["samples"])
    processed_ids = samples_84161[0]["feature_ids"] if samples_84161 else []
    # The parser retains the first 25,000 IDs for compact sample/header
    # inspection.  Platform annotation rows are retained in annotation_by_id,
    # so use that complete key set for the full GPL statistics.
    full_ids = list((platform_84161.get("table") or {}).get("annotation_by_id", {}).keys())
    processed_annotation = _annotation_stats(processed_ids, platform_84161.get("table") or {})
    full_annotation = _annotation_stats(full_ids, platform_84161.get("table") or {})
    groups, pairs, mapping_errors = _build_control_map(samples_84161)
    superseries_metadata = _compare_series(
        series_84161,
        series_84162,
        samples_84161,
        samples_84162,
    )
    cel_urls = [
        url
        for sample in samples_84161
        for url in sample["supplementary_files"]
        if ".CEL" in url.upper()
    ]
    series_raw_url = _field(series_84161, "Series_supplementary_file", GSE84161_RAW_URL)
    source_files = [
        _source_provenance(gse84161_path, GSE84161_SOFT_URL),
        _source_provenance(gse84162_path, GSE84162_SOFT_URL),
    ]
    donor_characteristics = [
        item for item in _field_list(series_84161, "Series_overall_design") if item
    ]
    design_text = _field(series_84161, "Series_overall_design")
    donor_design_match = re.search(r"total of\s+(\d+)\s+donors", design_text, re.IGNORECASE)
    sample_level_donor_fields = [
        item["characteristics"].get("donor", "")
        for item in samples_84161
        if item["characteristics"].get("donor", "")
    ]
    series_timing = "1 hr"
    sample_timing = "at least 30 min"
    if "1 hr" not in design_text:
        series_timing = "not found"
    if not all(sample_timing in item["treatment_protocol"] for item in samples_84161):
        sample_timing = "not found"

    processing_rows = _processing_rows(
        series_84161,
        samples_84161,
        processed,
        platform_84161,
        platform_summary,
        processed_annotation,
        full_annotation,
        series_raw_url,
        cel_urls,
    )
    sample_rows = _sample_rows(samples_84161, processed, series_raw_url)

    report: dict[str, Any] = {
        "audit_id": "gse84161-metadata-audit-v1",
        "created_at_utc": "2026-09-12",
        "status": "metadata_and_identifier_audit_only",
        "scope": {
            "expression_values_read": False,
            "program_scores_computed": False,
            "target_values_or_differences_read": False,
            "raw_cel_downloaded_by_this_audit": False,
            "raw_cel_contents_inspected_by_this_audit": False,
            "purpose": "provenance and array-method preparation; no effect estimation",
        },
        "source_files": source_files,
        "series": {
            "gse84161": _series_summary(series_84161),
            "gse84162": _series_summary(series_84162),
        },
        "gse84161": {
            "geo_url": GSE84161_GEO_URL,
            "series_matrix_url": GSE84161_MATRIX_URL,
            "raw_archive_url": series_raw_url,
            "sample_count": len(samples_84161),
            "cell_context": "primary human peripheral-blood monocytes from buffy coats",
            "explicit_series_donor_count": int(donor_design_match.group(1)) if donor_design_match else None,
            "sample_level_donor_field_present": bool(sample_level_donor_fields),
            "sample_level_donor_field_values": sample_level_donor_fields,
            "donor_pairing_status": "series-level five-donor statement; sample-level replicate number 1-5 only; donor identity not auditable from these records",
            "replicate_groups": groups,
            "planned_control_pairs": pairs,
            "planned_control_pair_count": len(pairs),
            "mapping_errors": mapping_errors,
            "six_arm_condition_keys": ARM_ORDER,
            "processed_table": processed,
            "platform_gpl570": platform_summary,
            "processed_annotation_via_current_gpl570": processed_annotation,
            "full_gpl570_annotation": full_annotation,
            "cel_availability": {
                "series_raw_archive_url": series_raw_url,
                "sample_cel_url_count": len(cel_urls),
                "sample_count": len(samples_84161),
                "all_sample_cel_routes_advertised": len(cel_urls) == len(samples_84161),
                "downloaded_by_this_audit": False,
                "availability_limit": "route advertised in SOFT; raw-array archive and CEL contents are outside this metadata-only audit",
            },
            "processing_evidence_ids": [row["evidence_id"] for row in processing_rows],
            "timing_evidence": {
                "series_design_pretreatment": series_timing,
                "sample_protocol_pretreatment": sample_timing,
            },
        },
        "gse84162": {
            "geo_url": GSE84162_GEO_URL,
            "raw_archive_url": _field(series_84162, "Series_supplementary_file", GSE84162_RAW_URL),
            "sample_count": len(samples_84162),
            "cell_contexts": sorted({item["primary_tissue"] for item in samples_84162}),
            "superseries_only_neutrophil_context": superseries_metadata["superseries_only_neutrophil_context"],
        },
        "series_comparison": superseries_metadata,
        "contradictions_and_ambiguities": [
            {
                "id": "pretreatment_duration",
                "kind": "differing_specificity",
                "series_statement": "1 hr inhibitor uptake in Series_overall_design",
                "sample_statement": "at least 30 min pretreatment in Sample_treatment_protocol_ch1",
                "resolution": "preserve both; treat the one-hour design statement and at-least-30-minute sample statement as compatible descriptions with different specificity, then freeze the operational duration before analysis",
            },
            {
                "id": "donor_vs_replicate",
                "kind": "ambiguity",
                "statement": "Series says five total donors; sample records provide replicate number 1-5 and no donor identifier",
                "resolution": "use replicate_N as a public grouping label only; do not call it a donor ID",
            },
            {
                "id": "variance_scope",
                "kind": "ambiguity",
                "statement": "Source says highest variance probe per Entrez but does not define the sample set or filtering scope used for variance",
                "resolution": "do not treat the processed 19,944 rows as an independent reprocessing rule; raw-CEL method must freeze a deterministic mapping",
            },
            {
                "id": "heatmap_residual_quantity",
                "kind": "ambiguity",
                "statement": "Source describes donor-residual normalization for heat-map visualization after RMA, while the deposited VALUE quantity is labeled RMA normalized signal intensity",
                "resolution": "retain the narrative; do not assume the heat-map residual operation generated the deposited table",
            },
            {
                "id": "historical_gpl_annotation",
                "kind": "availability_limit",
                "statement": "Current GPL570 annotation snapshot is separate from the historical probe-selection snapshot",
                "resolution": "report full 54,675-probe GPL and processed 19,944-probe subset separately; current mapping does not prove an author violation",
            },
        ],
        "maximum_claim": "GSE84161 can support source-provenance and descriptive array-specific calibration for primary human monocyte MEK/TPL2/LPS context using source-labeled arms and public replicate groups. It cannot support unchanged raw-integer-count analysis or independent direct ETS2 causality; missing sample-level donor keys limit donor-level inference, and unresolved compound identity limits compound-specific causal or clinical claims.",
        "processing_rule_boundary": {
            "source_described_selection": "RMA; remove probes without Entrez; retain highest-variance probe set per Entrez",
            "processed_table": "19,944 Affymetrix probe-set rows per sample",
            "full_gpl570_table": "54,675 probe rows in current SOFT platform snapshot",
            "future_array_plan_required": "raw CEL RMA from all probes, current/pinned GPL annotation, deterministic Entrez mapping and duplicate policy before effects",
        },
        "output_paths": {
            "samples_csv": str(SAMPLE_CSV.relative_to(ROOT)),
            "processing_csv": str(PROCESSING_CSV.relative_to(ROOT)),
            "report_json": str(REPORT_JSON.relative_to(ROOT)),
        },
    }

    if write_outputs:
        sample_fieldnames = [
            "sample_accession",
            "source_order",
            "source_title",
            "source_name",
            "series_memberships",
            "primary_tissue",
            "treatment_name",
            "stimulation",
            "compound",
            "replicate_number",
            "replicate_group_label",
            "donor_identity_status",
            "condition_key",
            "control_condition_key",
            "control_for_condition_keys",
            "matched_control_accession",
            "control_role",
            "source_characteristics",
            "assay_platform",
            "quantity_label",
            "processed_table_header",
            "processed_table_row_count",
            "processed_unique_id_count",
            "processed_duplicate_id_count",
            "processed_id_digest",
            "processed_ids_match_all_samples",
            "source_data_row_count",
            "cel_available_route",
            "cel_url",
            "series_raw_archive_url",
        ]
        processing_fieldnames = [
            "evidence_id",
            "source_accession",
            "scope",
            "field",
            "observed_value",
            "interpretation",
            "status",
        ]
        _write_csv(SAMPLE_CSV, sample_rows, sample_fieldnames)
        _write_csv(PROCESSING_CSV, processing_rows, processing_fieldnames)
        REPORTS.mkdir(parents=True, exist_ok=True)
        REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true", help="parse and validate without writing artifacts")
    args = parser.parse_args()
    report = build_audit(write_outputs=not args.no_write)
    print(
        json.dumps(
            {
                "gse84161_samples": report["gse84161"]["sample_count"],
                "processed_rows": report["gse84161"]["processed_table"]["processed_row_count"],
                "processed_ids_identical": report["gse84161"]["processed_table"]["all_samples_same_ordered_ids"],
                "gpl570_rows": report["gse84161"]["platform_gpl570"]["table_row_count"],
                "planned_control_pairs": report["gse84161"]["planned_control_pair_count"],
                "gse84162_samples": report["gse84162"]["sample_count"],
                "neutrophil_samples": report["gse84162"]["superseries_only_neutrophil_context"]["sample_count"],
                "mapping_errors": report["gse84161"]["mapping_errors"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
