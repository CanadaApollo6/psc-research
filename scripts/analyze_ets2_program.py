"""Run the frozen, donor-aware ETS2 program analysis.

This module is deliberately separate from the five-gene liver-atlas analysis.
It reuses that module's immutable metadata and raw-count validation helpers,
but keeps the all-feature donor aggregation, expression-matched null sets and
program summaries in this file.  The public outputs contain aggregate tables;
donor-level expression and random-set selections are kept in ignored ``work``
files or in memory only.

No program score is meaningful until the frozen plan and all of its pinned
inputs pass validation.  The source-defined sets are read from the original
release ZIP using the exact columns named by the plan.  Saved rank-statistic
orientation is intentionally outside this module: this task only scores the
source memberships and the pinned liver atlas.
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
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

try:  # Imported lazily by the source/aggregation paths when available.
    from scipy import sparse
except ImportError:  # pragma: no cover - requirements pin scipy
    sparse = None  # type: ignore[assignment]

try:
    from analyze_liver_atlas import (
        ROOT as ATLAS_ROOT,
        _matrix_format,
        _stable_unique,
        _is_sparse,
        iter_validated_count_chunks,
        prepare_metadata,
        resolve_raw_counts,
        sha256_file,
    )
except ModuleNotFoundError:  # Import as ``scripts.analyze_ets2_program``.
    from scripts.analyze_liver_atlas import (
        ROOT as ATLAS_ROOT,
        _matrix_format,
        _stable_unique,
        _is_sparse,
        iter_validated_count_chunks,
        prepare_metadata,
        resolve_raw_counts,
        sha256_file,
    )


ROOT = Path(__file__).resolve().parents[1]
ATLAS_SCRIPT = ROOT / "scripts/analyze_liver_atlas.py"

DATASET_KEYS: tuple[str, ...] = ("sc", "sn")
AUTHOR_CELL_TYPES: tuple[str, ...] = ("Monocyte", "Kupffer", "LAM-like", "ActMac")
THRESHOLDS: tuple[int, ...] = (20, 50)
ORIGINAL_PROGRAM_IDS: tuple[str, ...] = (
    "ets2_g1_dn",
    "ets2_g1_up",
    "ets2_g2_dn",
    "chr21_dn",
    "inflammation",
    "interferon_gamma",
    "oxidative_stress",
    "apoptotic_signaling",
)
COMPARATOR_PROGRAM_IDS: tuple[str, ...] = (
    "inflammation",
    "interferon_gamma",
    "oxidative_stress",
    "apoptotic_signaling",
)
PROGRAM_IDS: tuple[str, ...] = ORIGINAL_PROGRAM_IDS + ("ets2_g1_dn_without_comparators",)
PRIMARY_PROGRAM_ID = "ets2_g1_dn"
ETS2_STABLE_ID = "ENSG00000157557"
ETS2_SYMBOL = "ETS2"
DEFAULT_BINS = 20
DEFAULT_REPLICATES = 1000
DEFAULT_BASE_SEED = 20260912
DEFAULT_MIN_MAPPED_FRACTION = 0.8
DEFAULT_MIN_MAPPED_GENES = 10
DEFAULT_MIN_PAIRED_DONORS = 3
DEFAULT_MIN_CORRELATION_DONORS = 4
# The plan is a frozen input.  This value is intentionally outside the JSON
# itself: a file cannot contain a trustworthy hash of its own final bytes.
FROZEN_PLAN_SHA256 = "234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9"


def _clean(value: Any) -> str:
    """Convert a scalar to a stable, blank-aware string."""

    if value is None:
        return ""
    try:
        missing = pd.isna(value)
        if missing is pd.NA or (isinstance(missing, (bool, np.bool_)) and bool(missing)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _is_missing(value: Any) -> bool:
    text = _clean(value)
    return not text or text.casefold() in {"nan", "none", "na", "n/a", "null"}


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def _json_safe(value: Any) -> Any:
    """Convert audit values to strict JSON, preserving unavailable as null."""

    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if math.isfinite(number) else None
    if value is pd.NA:
        return None
    return value


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_ensembl_id(value: Any) -> str | None:
    """Strip only a numeric Ensembl version suffix.

    The matcher intentionally does not canonicalize aliases, symbols or
    retired IDs.  A value is treated as a stable Ensembl gene ID only when it
    has the exact ENSG plus numeric form.
    """

    text = _clean(value)
    if re.fullmatch(r"ENSG\d+(?:\.\d+)?", text):
        return text.split(".", 1)[0]
    return None


def _assert_plan_invariants(plan: Mapping[str, Any]) -> None:
    """Fail closed if a plan changes the fixed analysis contract."""

    if str(plan.get("analysis_id", "")) != "psc-ets2-program-v1":
        raise ValueError("Unexpected ETS2 program analysis_id")
    if str(plan.get("cohort", "")) != "PSC":
        raise ValueError("ETS2 program analysis is frozen to PSC cells")
    if tuple(plan.get("datasets", [{}])[i].get("key") for i in range(len(plan.get("datasets", [])))) != DATASET_KEYS:
        raise ValueError("Frozen ETS2 dataset keys must be sc then sn")
    if tuple(plan.get("author_cell_types", [])) != AUTHOR_CELL_TYPES:
        raise ValueError("Frozen author-cell-type list changed")
    if tuple(int(value) for value in plan.get("cell_thresholds", [])) != THRESHOLDS:
        raise ValueError("Frozen cell thresholds changed")
    if int(plan.get("minimum_paired_donors", -1)) != DEFAULT_MIN_PAIRED_DONORS:
        raise ValueError("Frozen paired-donor threshold changed")
    if int(plan.get("minimum_correlation_donors", -1)) != DEFAULT_MIN_CORRELATION_DONORS:
        raise ValueError("Frozen correlation-donor threshold changed")
    if str(plan.get("count_source", "")) != "raw/X":
        raise ValueError("ETS2 program must use raw/X counts")
    if tuple(str(value) for value in plan.get("exclude_gene_ids", [])) != (ETS2_STABLE_ID,):
        raise ValueError("Frozen ETS2 exclusion ID changed")
    if tuple(str(value) for value in plan.get("exclude_gene_symbols", [])) != (ETS2_SYMBOL,):
        raise ValueError("Frozen ETS2 exclusion symbol changed")
    if float(plan.get("minimum_mapped_fraction", -1)) != DEFAULT_MIN_MAPPED_FRACTION:
        raise ValueError("Frozen mapping fraction gate changed")
    if int(plan.get("minimum_mapped_genes", -1)) != DEFAULT_MIN_MAPPED_GENES:
        raise ValueError("Frozen minimum mapped gene gate changed")
    programs = plan.get("programs", [])
    if tuple(str(item.get("id")) for item in programs) != PROGRAM_IDS:
        raise ValueError("Frozen ETS2 program list changed")
    expected_columns = {
        "ets2_g1_dn": ("RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv", "ETS2g1_DN"),
        "ets2_g1_up": ("RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv", "ETS2g1_UP"),
        "ets2_g2_dn": ("RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv", "ETS2g2_DN"),
        "chr21_dn": ("RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv", "CHR21_DN"),
        "inflammation": ("RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv", "GOBP_INFLAMMATORY_RESPONSE"),
        "interferon_gamma": ("RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv", "GOBP_RESPONSE_TO_INTERFERON_GAMMA"),
        "oxidative_stress": ("RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"),
        "apoptotic_signaling": ("RNA-seq/RNAseq_CRISPR/GOBP_genesets.csv", "GOBP_APOPTOTIC_SIGNALING_PATHWAY"),
    }
    by_id = {str(item.get("id")): item for item in programs}
    for program_id, (member, column) in expected_columns.items():
        item = by_id.get(program_id)
        if item is None or str(item.get("source_member")) != member or str(item.get("column")) != column:
            raise ValueError(f"Frozen source definition changed for {program_id}")
    derived = by_id.get("ets2_g1_dn_without_comparators", {})
    if str(derived.get("derived_from")) != PRIMARY_PROGRAM_ID or tuple(derived.get("subtract_union", [])) != COMPARATOR_PROGRAM_IDS:
        raise ValueError("Frozen non-overlap derivation changed")
    matching = plan.get("matching", {})
    if int(matching.get("bins", -1)) != DEFAULT_BINS:
        raise ValueError("Frozen expression-bin count changed")
    if int(matching.get("replicates", -1)) != DEFAULT_REPLICATES:
        raise ValueError("Frozen random-set replicate count changed")
    if int(matching.get("base_seed", -1)) != DEFAULT_BASE_SEED:
        raise ValueError("Frozen random-set seed changed")
    if int(matching.get("bins", -1)) <= 0 or int(matching.get("replicates", -1)) <= 0:
        raise ValueError("Frozen matching parameters must be positive")
    if tuple(str(item.get("left")) for item in plan.get("contrasts", [])) != ("ActMac", "LAM-like", "Monocyte"):
        raise ValueError("Frozen paired-contrast left populations changed")
    if any(str(item.get("right")) != "Kupffer" for item in plan.get("contrasts", [])):
        raise ValueError("Frozen paired-contrast right population changed")


def _verify_plan_inputs(root: Path, plan: Mapping[str, Any], plan_path: Path) -> dict[str, str]:
    """Verify all immutable files whose hashes are recorded in the plan."""

    recorded = plan.get("input_sha256", {})
    if not isinstance(recorded, Mapping):
        raise ValueError("Frozen plan is missing input_sha256")
    checked: dict[str, str] = {}
    for relative, expected in sorted(recorded.items()):
        path = Path(str(relative))
        if not path.is_absolute():
            path = root / path
        if not path.exists():
            raise FileNotFoundError(path)
        actual = sha256_file(path)
        if actual != str(expected):
            raise ValueError(f"Frozen input hash changed for {relative}")
        checked[str(relative)] = actual
    return checked


def load_frozen_plan(plan_path: Path, root: Path = ROOT) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    """Read and validate the frozen plan and all plan-level input hashes."""

    actual_plan_hash = sha256_file(plan_path)
    if actual_plan_hash != FROZEN_PLAN_SHA256:
        raise ValueError("Frozen ETS2 program plan hash changed")
    plan = _read_json(plan_path)
    _assert_plan_invariants(plan)
    hashes = _verify_plan_inputs(root, plan, plan_path)
    source_manifest_path = root / "config/ets2-program-sources.json"
    manifest = _read_json(source_manifest_path)
    expected_manifest_hash = plan.get("input_sha256", {}).get("config/ets2-program-sources.json")
    if expected_manifest_hash and sha256_file(source_manifest_path) != str(expected_manifest_hash):
        raise ValueError("ETS2 program source manifest changed")
    return plan, manifest, hashes


def verify_pinned_source(path: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    """Verify byte count and SHA-256 for one manifest source."""

    if not path.exists():
        raise FileNotFoundError(path)
    expected_bytes = int(source["bytes"])
    if path.stat().st_size != expected_bytes:
        raise ValueError(f"Pinned source byte length changed for {path}")
    expected_sha = str(source["sha256"])
    actual_sha = sha256_file(path)
    if actual_sha != expected_sha:
        raise ValueError(f"Pinned source SHA-256 changed for {path}")
    return {"bytes": expected_bytes, "sha256": actual_sha}


def _source_zip_from_manifest(root: Path, manifest: Mapping[str, Any]) -> Path:
    directory = root / str(manifest.get("source_directory", "data/raw/ets2-benchmark"))
    records = [item for item in manifest.get("sources", []) if str(item.get("file")) == "ets2-public-release.zip"]
    if len(records) != 1:
        raise ValueError("ETS2 benchmark manifest must contain one public release ZIP")
    path = directory / str(records[0]["file"])
    verify_pinned_source(path, records[0])
    published_md5 = manifest.get("published_zip_md5")
    if published_md5:
        with path.open("rb") as stream:
            actual_md5 = hashlib.file_digest(stream, "md5").hexdigest()
        if actual_md5 != str(published_md5):
            raise ValueError("Publisher ETS2 release ZIP MD5 mismatch")
    return path


def _resolve_zip_member(stream: zipfile.ZipFile, logical_member: str) -> str:
    exact = [name for name in stream.namelist() if name == logical_member]
    if len(exact) == 1:
        return exact[0]
    suffix = [name for name in stream.namelist() if name.endswith("/" + logical_member)]
    if len(suffix) != 1:
        raise ValueError(f"Expected one ZIP member for {logical_member!r}; found {len(suffix)}")
    return suffix[0]


def read_source_column(zip_path: Path, member: str, column: str) -> tuple[list[str], dict[str, Any]]:
    """Read one exact source column while preserving source row order."""

    with zipfile.ZipFile(zip_path) as stream:
        actual_member = _resolve_zip_member(stream, member)
        info = stream.getinfo(actual_member)
        content = stream.read(actual_member)
        text = io.TextIOWrapper(io.BytesIO(content), encoding="utf-8", newline="")
        reader = csv.reader(text)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Source member {member} is empty") from exc
        if len(header) != len(set(header)):
            raise ValueError(f"Source member {member} has duplicate column names")
        if column not in header:
            raise ValueError(f"Source member {member} is missing exact column {column}")
        index = header.index(column)
        values: list[str] = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            if len(row) > len(header):
                raise ValueError(f"Source member {member} row {row_number} has too many fields")
            value = _clean(row[index])
            if not _is_missing(value):
                values.append(value)
        return values, {"member": actual_member, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def source_member_fingerprints(zip_path: Path, members: Iterable[str]) -> dict[str, dict[str, Any]]:
    """Return pinned byte and SHA-256 metadata for selected archive members."""

    wanted = sorted(set(str(member) for member in members))
    with zipfile.ZipFile(zip_path) as stream:
        result: dict[str, dict[str, Any]] = {}
        for logical in wanted:
            actual = _resolve_zip_member(stream, logical)
            content = stream.read(actual)
            result[logical] = {
                "member": actual,
                "bytes": len(content),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
        return result


def load_source_programs(plan: Mapping[str, Any], zip_path: Path) -> tuple[dict[str, list[str]], dict[str, dict[str, Any]]]:
    """Load exactly the eight source-defined columns listed by the plan."""

    programs: dict[str, list[str]] = {}
    records: dict[str, dict[str, Any]] = {}
    for item in plan["programs"]:
        program_id = str(item["id"])
        if program_id == "ets2_g1_dn_without_comparators":
            continue
        values, _ = read_source_column(zip_path, str(item["source_member"]), str(item["column"]))
        programs[program_id] = values
        records[program_id] = {
            "source_member": str(item["source_member"]),
            "column": str(item["column"]),
            "source_rows_nonblank": len(values),
            "source_unique_values": len(set(values)),
        }
    if tuple(programs) != ORIGINAL_PROGRAM_IDS:
        raise ValueError("Source program loading did not produce the eight frozen programs")
    return programs, records


def build_feature_map(raw_var: pd.DataFrame, n_features: int) -> pd.DataFrame:
    """Build stable-ID and exact-name lookup columns for an atlas matrix."""

    if len(raw_var) != int(n_features):
        raise ValueError("raw/var and raw/X feature dimensions disagree")
    if "feature_name" not in raw_var.columns:
        raise ValueError("raw/var is missing feature_name")
    feature_ids = [_clean(value) for value in raw_var.index]
    if any(not value for value in feature_ids) or len(feature_ids) != len(set(feature_ids)):
        raise ValueError("raw/var feature identifiers are missing or duplicated")
    names = [_clean(value) for value in raw_var["feature_name"]]
    return pd.DataFrame(
        {
            "feature_index": np.arange(int(n_features), dtype=np.int64),
            "feature_id": feature_ids,
            "stable_feature_id": [normalize_ensembl_id(value) or "" for value in feature_ids],
            "feature_name": names,
        }
    )


def _lookup_feature_indices(features: pd.DataFrame) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    stable: dict[str, list[int]] = {}
    names: dict[str, list[int]] = {}
    for row in features.itertuples(index=False):
        if row.stable_feature_id:
            stable.setdefault(str(row.stable_feature_id), []).append(int(row.feature_index))
        if row.feature_name:
            names.setdefault(str(row.feature_name), []).append(int(row.feature_index))
    return stable, names


def _excluded_source(value: str) -> bool:
    return normalize_ensembl_id(value) == ETS2_STABLE_ID or value == ETS2_SYMBOL


def map_source_program(
    program_id: str,
    source_values: Sequence[str],
    features: pd.DataFrame,
    *,
    dataset_key: str,
    role: str = "",
) -> tuple[pd.DataFrame, set[int], dict[str, Any]]:
    """Map a source membership list by exact ID/name rules.

    The returned audit has one row for every unique source value.  Excluded
    ETS2 membership is retained as an explicit status, but never enters the
    target feature set or its coverage denominator.
    """

    stable_lookup, name_lookup = _lookup_feature_indices(features)
    # Deduplicate by the exact matching identity.  Versioned Ensembl values
    # with the same stable ID represent one intended gene; mixed IDs/symbols
    # remain exact, case-sensitive feature-name identities.
    unique_values: list[str] = []
    seen_source_keys: set[str] = set()
    for raw_value in source_values:
        value = _clean(raw_value)
        if _is_missing(value):
            continue
        source_key = normalize_ensembl_id(value) or value
        if source_key not in seen_source_keys:
            seen_source_keys.add(source_key)
            unique_values.append(value)
    records: list[dict[str, Any]] = []
    mapped_indices: set[int] = set()
    for source_index, value in enumerate(unique_values):
        normalized = normalize_ensembl_id(value)
        status: str
        match_mode: str
        candidates: list[int]
        if _excluded_source(value):
            records.append(
                {
                    "dataset_key": dataset_key,
                    "program_id": program_id,
                    "role": role,
                    "source_index": source_index,
                    "source_value": value,
                    "normalized_source_value": normalized or value,
                    "mapping_mode": "excluded_target",
                    "mapping_status": "excluded_target_gene",
                    "feature_index": "",
                    "feature_id": "",
                    "feature_name": ETS2_SYMBOL,
                }
            )
            continue
        if normalized is not None:
            match_mode = "exact_stable_ensembl_id"
            candidates = stable_lookup.get(normalized, [])
            if len(candidates) == 1:
                status = "mapped"
            elif len(candidates) == 0:
                status = "missing_ensembl_id"
            else:
                status = "ambiguous_ensembl_id"
        else:
            match_mode = "exact_feature_name"
            candidates = name_lookup.get(value, [])
            if len(candidates) == 1:
                status = "mapped"
            elif len(candidates) == 0:
                status = "missing_feature_name"
            else:
                status = "ambiguous_feature_name"
        if status == "mapped":
            index = candidates[0]
            mapped_indices.add(index)
            feature = features.iloc[index]
            feature_index: Any = int(feature.feature_index)
            feature_id: Any = str(feature.feature_id)
            feature_name: Any = str(feature.feature_name)
        else:
            feature_index = ""
            feature_id = ";".join(str(features.iloc[index].feature_id) for index in candidates)
            feature_name = value
        records.append(
            {
                "dataset_key": dataset_key,
                "program_id": program_id,
                "role": role,
                "source_index": source_index,
                "source_value": value,
                "normalized_source_value": normalized or value,
                "mapping_mode": match_mode,
                "mapping_status": status,
                "feature_index": feature_index,
                "feature_id": feature_id,
                "feature_name": feature_name,
            }
        )
    table = pd.DataFrame(records)
    excluded = int((table.mapping_status == "excluded_target_gene").sum()) if len(table) else 0
    intended = int(len(unique_values) - excluded)
    mapped = int(len(mapped_indices))
    missing = int(table.mapping_status.astype(str).str.startswith("missing").sum()) if len(table) else 0
    ambiguous = int(table.mapping_status.astype(str).str.startswith("ambiguous").sum()) if len(table) else 0
    coverage = mapped / intended if intended else 0.0
    inventory = {
        "dataset_key": dataset_key,
        "program_id": program_id,
        "source_unique_before_exclusion": len(unique_values),
        "excluded_target_gene_count": excluded,
        "intended_unique_source_genes": intended,
        "mapped_unique_genes": mapped,
        "missing_source_genes": missing,
        "ambiguous_source_genes": ambiguous,
        "mapping_fraction": coverage,
        "mapping_gate_pass": bool(coverage >= DEFAULT_MIN_MAPPED_FRACTION and mapped >= DEFAULT_MIN_MAPPED_GENES),
    }
    return table, mapped_indices, inventory


def _program_roles(plan: Mapping[str, Any]) -> dict[str, str]:
    return {str(item["id"]): str(item.get("role", "")) for item in plan["programs"]}


def build_program_mappings(
    plan: Mapping[str, Any],
    source_values: Mapping[str, Sequence[str]],
    feature_maps: Mapping[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, set[int]]]]:
    """Map all source programs and derive the ninth non-overlap program."""

    roles = _program_roles(plan)
    member_records: list[pd.DataFrame] = []
    inventory_records: list[dict[str, Any]] = []
    mapped_by_dataset: dict[str, dict[str, set[int]]] = {key: {} for key in feature_maps}
    direct_tables: dict[tuple[str, str], pd.DataFrame] = {}
    direct_inventory: dict[tuple[str, str], dict[str, Any]] = {}
    for dataset_key in feature_maps:
        for program_id in ORIGINAL_PROGRAM_IDS:
            table, mapped, inventory = map_source_program(
                program_id,
                source_values[program_id],
                feature_maps[dataset_key],
                dataset_key=dataset_key,
                role=roles.get(program_id, ""),
            )
            direct_tables[(dataset_key, program_id)] = table
            direct_inventory[(dataset_key, program_id)] = inventory
            mapped_by_dataset[dataset_key][program_id] = mapped
            member_records.append(table)
            inventory_records.append({**inventory, "role": roles.get(program_id, "")})

    # The derived program is defined after atlas mapping.  Parent coverage is
    # preserved as its gate; comparator overlap is recorded on every parent
    # membership row so mapping accountability is complete.
    derived_id = "ets2_g1_dn_without_comparators"
    for dataset_key, features in feature_maps.items():
        parent_table = direct_tables[(dataset_key, PRIMARY_PROGRAM_ID)].copy()
        comparator_union: set[int] = set()
        for program_id in COMPARATOR_PROGRAM_IDS:
            comparator_union.update(mapped_by_dataset[dataset_key][program_id])
        parent_mapped = mapped_by_dataset[dataset_key][PRIMARY_PROGRAM_ID]
        derived_indices = parent_mapped - comparator_union
        mapped_by_dataset[dataset_key][derived_id] = set(derived_indices)
        if len(parent_table):
            parent_table["derived_program_id"] = derived_id
            parent_table["removed_by_comparator"] = parent_table.apply(
                lambda row: bool(
                    row.mapping_status == "mapped"
                    and int(row.feature_index) in comparator_union
                ),
                axis=1,
            )
            parent_table["derived_mapping_status"] = parent_table.apply(
                lambda row: (
                    "mapped_nonoverlap"
                    if row.mapping_status == "mapped" and not bool(row.removed_by_comparator)
                    else "mapped_comparator_overlap"
                    if row.mapping_status == "mapped"
                    else str(row.mapping_status)
                ),
                axis=1,
            )
            # Keep only the derived view as additional membership rows, while
            # retaining the original table untouched in the direct program.
            derived_table = parent_table.copy()
            derived_table["program_id"] = derived_id
            derived_table["role"] = roles.get(derived_id, "")
            derived_table["mapping_status"] = derived_table["derived_mapping_status"]
            derived_table["mapping_mode"] = "derived_after_exact_mapping"
            member_records.append(derived_table.drop(columns=["derived_program_id", "removed_by_comparator", "derived_mapping_status"]))
        parent_inventory = direct_inventory[(dataset_key, PRIMARY_PROGRAM_ID)]
        intended_parent = int(parent_inventory["intended_unique_source_genes"])
        parent_gate = bool(parent_inventory["mapping_gate_pass"])
        source_overlap = int(len(parent_mapped & comparator_union))
        derived_records = {
            "dataset_key": dataset_key,
            "program_id": derived_id,
            "role": roles.get(derived_id, ""),
            "source_unique_before_exclusion": int(parent_inventory["source_unique_before_exclusion"]),
            "excluded_target_gene_count": int(parent_inventory["excluded_target_gene_count"]),
            # The denominator is the parent's intended source set by freeze;
            # no post hoc denominator based on expression or overlap is used.
            "intended_unique_source_genes": intended_parent,
            "mapped_unique_genes": int(len(derived_indices)),
            "missing_source_genes": int(parent_inventory["missing_source_genes"]),
            "ambiguous_source_genes": int(parent_inventory["ambiguous_source_genes"]),
            "mapping_fraction": (len(derived_indices) / intended_parent if intended_parent else 0.0),
            "derived_mapping_fraction": (len(derived_indices) / intended_parent if intended_parent else 0.0),
            "parent_mapping_fraction": float(parent_inventory["mapping_fraction"]),
            "mapping_gate_pass": bool(parent_gate and len(derived_indices) >= DEFAULT_MIN_MAPPED_GENES),
            "parent_mapping_gate_pass": parent_gate,
            "parent_program_id": PRIMARY_PROGRAM_ID,
            "parent_mapped_unique_genes": int(len(parent_mapped)),
            "comparator_union_mapped_genes": int(len(comparator_union)),
            "parent_comparator_overlap_genes": source_overlap,
            "derived_nonoverlap_mapped_genes": int(len(derived_indices)),
        }
        inventory_records.append(derived_records)
    membership = pd.concat(member_records, ignore_index=True) if member_records else pd.DataFrame()
    inventory = pd.DataFrame(inventory_records)
    if len(membership):
        membership = membership.sort_values(["dataset_key", "program_id", "source_index", "source_value"], kind="mergesort").reset_index(drop=True)
    if len(inventory):
        inventory = inventory.sort_values(["dataset_key", "program_id"], kind="mergesort").reset_index(drop=True)
    return membership, inventory, mapped_by_dataset


@dataclass
class AggregatedAtlas:
    dataset_key: str
    pseudobulk: pd.DataFrame
    feature_map: pd.DataFrame
    log_values: np.ndarray
    count_values: np.ndarray
    qc: dict[str, Any]


def _group_columns() -> list[str]:
    return [
        "dataset_key",
        "dataset_id",
        "modality",
        "suspension_type",
        "assay",
        "donor_id",
        "cohort",
        "author_cell_type",
    ]


def _matrix_values(value: Any) -> np.ndarray:
    if _is_sparse(value):
        return np.asarray(value.toarray())
    return np.asarray(value)


def aggregate_all_features(
    adata: Any,
    dataset_key: str,
    *,
    dataset_id: str,
    version_id: str = "",
    source_file: str = "",
    expected_observations: int | None = None,
    author_cell_types: Sequence[str] = AUTHOR_CELL_TYPES,
    cohort: str = "PSC",
    thresholds: Sequence[int] = THRESHOLDS,
    chunk_size: int = 4096,
) -> AggregatedAtlas:
    """Validate raw counts and aggregate all features for PSC populations."""

    if dataset_key not in DATASET_KEYS:
        raise ValueError(f"Unexpected dataset key {dataset_key!r}")
    n_obs = int(getattr(adata, "n_obs", getattr(adata, "shape", (0,))[0]))
    n_vars = int(getattr(adata, "n_vars", getattr(adata, "shape", (0, 0))[1]))
    if expected_observations is not None and n_obs != int(expected_observations):
        raise ValueError(f"{dataset_key} has {n_obs} observations; expected {expected_observations}")
    obs = pd.DataFrame(getattr(adata, "obs")).copy()
    if len(obs) != n_obs:
        raise ValueError("AnnData obs row count differs from n_obs")
    obs_names = pd.Index([str(value) for value in getattr(adata, "obs_names", obs.index)])
    if obs_names.has_duplicates:
        raise ValueError(f"{dataset_key} has duplicate cell IDs")
    obs.index = obs_names
    metadata = prepare_metadata(obs, dataset_key, dataset_id)
    raw_matrix, raw_var, count_source = resolve_raw_counts(adata)
    if raw_matrix.shape[0] != n_obs:
        raise ValueError("raw/X and obs dimensions disagree")
    if raw_matrix.shape[1] != n_vars and getattr(getattr(adata, "raw", None), "n_vars", raw_matrix.shape[1]) != raw_matrix.shape[1]:
        raise ValueError("raw/X feature dimension is inconsistent")
    normalized_x = getattr(adata, "X", None)
    normalized_shape = list(getattr(normalized_x, "shape", ())) if normalized_x is not None else []
    if normalized_x is None or len(normalized_shape) != 2 or normalized_shape[0] != n_obs or normalized_shape[1] != int(raw_matrix.shape[1]):
        raise ValueError("normalized AnnData.X is missing or has the wrong dimensions")
    features = build_feature_map(raw_var, int(raw_matrix.shape[1]))

    selected = metadata[
        metadata["retained"]
        & metadata["cohort"].eq(cohort)
        & metadata["author_cell_type"].isin(list(author_cell_types))
    ].copy()
    group_columns = _group_columns()
    keys = [tuple(row[column] for column in group_columns) for row in selected.to_dict("records")]
    unique_keys = sorted(set(keys), key=lambda value: tuple(str(item) for item in value))
    code_by_key = {key: code for code, key in enumerate(unique_keys)}
    selected_codes = np.asarray([code_by_key[key] for key in keys], dtype=np.int64)
    selected_rows = selected.index.to_numpy(dtype=np.int64)
    n_groups = len(unique_keys)
    count_values = np.zeros((n_groups, int(raw_matrix.shape[1])), dtype=np.int64)
    cell_counts = np.zeros(n_groups, dtype=np.int64)
    qc: dict[str, Any] = {
        "dataset_key": dataset_key,
        "dataset_id": dataset_id,
        "version_id": version_id,
        "source_file": source_file,
        "count_source": count_source,
        "raw_rows": n_obs,
        "raw_features": int(raw_matrix.shape[1]),
        "matrix_format": _matrix_format(raw_matrix),
        "normalized_matrix": "AnnData.X checked only for shape; never used as counts",
        "normalized_matrix_shape": normalized_shape,
        "source_cohort": cohort,
        "author_cell_types": list(author_cell_types),
    }
    validation_chunks = 0
    stored_nonzero = 0
    minimums: list[float] = []
    maximums: list[float] = []
    for start, stop, chunk, stats in iter_validated_count_chunks(raw_matrix, chunk_size, count_source):
        validation_chunks += 1
        stored_nonzero += int(stats["nonzero"])
        if stats["min"] is not None:
            minimums.append(float(stats["min"]))
            maximums.append(float(stats["max"]))
        row_codes = np.full(stop - start, -1, dtype=np.int64)
        if len(selected_rows):
            left = int(np.searchsorted(selected_rows, start, side="left"))
            right = int(np.searchsorted(selected_rows, stop, side="left"))
            if right > left:
                row_codes[selected_rows[left:right] - start] = selected_codes[left:right]
        valid_rows = row_codes >= 0
        if not np.any(valid_rows):
            continue
        codes = row_codes[valid_rows]
        positions = np.flatnonzero(valid_rows)
        if _is_sparse(chunk):
            if sparse is None:  # pragma: no cover
                raise RuntimeError("scipy is required for sparse aggregation")
            selector = sparse.csr_matrix(
                (np.ones(len(codes), dtype=np.int8), (codes, positions)),
                shape=(n_groups, stop - start),
            )
            # The count validator accepts integer-valued floating storage, but
            # sparse multiplication would otherwise accumulate in float32 and
            # lose exact counts above the 2^24 precision boundary.
            chunk_for_sum = chunk.astype(np.int64)
            sums = selector @ chunk_for_sum
            sums = np.asarray(sums.toarray() if _is_sparse(sums) else sums, dtype=np.int64)
            count_values += sums
        else:
            values = np.asarray(chunk, dtype=np.int64)
            # np.add.at handles repeated group codes and avoids a dense
            # n_cells-by-n_features temporary beyond the current chunk.
            np.add.at(count_values, codes, values[valid_rows])
        np.add.at(cell_counts, codes, 1)
    if n_obs and validation_chunks == 0:
        raise ValueError("raw/X yielded no validation chunks")
    all_feature_umis = count_values.sum(axis=1, dtype=np.int64)
    qc.update(
        {
            "validation_chunks": validation_chunks,
            "stored_nonzero_values": stored_nonzero,
            "minimum_stored_value": min(minimums) if minimums else 0.0,
            "maximum_stored_value": max(maximums) if maximums else 0.0,
            "total_cells": n_obs,
            "retained_cells_all_metadata": int(metadata["retained"].sum()),
            "selected_psc_cells": int(len(selected)),
            "selected_donor_groups": n_groups,
            "selected_donors": int(selected["donor_id"].nunique()) if len(selected) else 0,
            "excluded_cells_all_metadata": int((~metadata["retained"]).sum()),
            "excluded_by_reason": {str(k): int(v) for k, v in metadata.loc[~metadata["retained"], "exclusion_reason"].value_counts().items()},
        }
    )

    records: list[dict[str, Any]] = []
    for group_index, key in enumerate(unique_keys):
        values = dict(zip(group_columns, key))
        group_meta = selected[selected_codes == group_index]
        total = int(all_feature_umis[group_index])
        cells = int(cell_counts[group_index])
        row = {
            **values,
            "dataset_version_id": version_id,
            "source_file": source_file,
            "group_id": "|".join(str(value) for value in key),
            "cell_count": cells,
            "all_feature_umis": total,
            "positive_all_feature_umi_denominator": bool(total > 0),
            "passes_primary_cell_gate": bool(cells >= int(thresholds[0])),
            "passes_sensitivity_cell_gate": bool(cells >= int(thresholds[1])),
            "library_count": len(_stable_unique(group_meta["library_uuid"])),
        }
        for col in ("library_uuid", "cell_type", "disease"):
            if col in group_meta:
                row[f"{col}_values"] = ";".join(_stable_unique(group_meta[col]))
        records.append(row)
    pseudobulk = pd.DataFrame(records)
    if len(pseudobulk):
        pseudobulk = pseudobulk.sort_values(group_columns, kind="mergesort").reset_index(drop=True)
        # Reorder counts to the sorted pseudobulk rows so metadata and matrix
        # rows stay in exact lockstep.
        original_order = [unique_keys.index(tuple(row[col] for col in group_columns)) for row in pseudobulk.to_dict("records")]
        count_values = count_values[np.asarray(original_order, dtype=np.int64)]
    denominators = pseudobulk["all_feature_umis"].to_numpy(dtype=float) if len(pseudobulk) else np.asarray([], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        cpm = count_values.astype(float) / denominators[:, None] * 1_000_000 if len(pseudobulk) else np.empty_like(count_values, dtype=float)
    if len(pseudobulk) and (denominators <= 0).any():
        cpm[denominators <= 0, :] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        log_values = np.log2(cpm + 1.0)
    return AggregatedAtlas(dataset_key, pseudobulk, features, log_values, count_values, qc)


@dataclass
class ReferenceResult:
    reference: np.ndarray
    donor_reference: pd.DataFrame
    eligible_group_count: int
    eligible_donor_count: int


def compute_reference(aggregate: AggregatedAtlas, *, min_cells: int = 20) -> ReferenceResult:
    """Compute equal-donor-weighted reference expression for all features."""

    pb = aggregate.pseudobulk
    if not len(pb):
        return ReferenceResult(np.asarray([], dtype=float), pd.DataFrame(), 0, 0)
    eligible = pb[(pb.cell_count >= int(min_cells)) & pb.positive_all_feature_umi_denominator]
    if not len(eligible):
        return ReferenceResult(np.full(aggregate.log_values.shape[1], np.nan), pd.DataFrame(), 0, 0)
    rows = eligible.index.to_numpy(dtype=np.int64)
    logs = aggregate.log_values[rows]
    temporary = pd.DataFrame({"donor_id": eligible.donor_id.astype(str).to_numpy(), "row_index": np.arange(len(rows), dtype=np.int64)})
    donor_records: list[dict[str, Any]] = []
    donor_means: list[np.ndarray] = []
    for donor, group in temporary.groupby("donor_id", sort=True):
        positions = group.row_index.to_numpy(dtype=np.int64)
        mean = logs[positions].mean(axis=0)
        donor_means.append(mean)
        donor_records.append({"donor_id": str(donor), "eligible_group_count": int(len(positions))})
    reference = np.vstack(donor_means).mean(axis=0) if donor_means else np.full(logs.shape[1], np.nan)
    return ReferenceResult(reference, pd.DataFrame(donor_records), len(eligible), len(donor_means))


def expression_bins(reference: Sequence[float], feature_ids: Sequence[str], bins: int = DEFAULT_BINS) -> tuple[np.ndarray, list[np.ndarray]]:
    """Rank all raw features and split into equal-count expression bins."""

    values = np.asarray(reference, dtype=float)
    ids = np.asarray([str(value) for value in feature_ids], dtype=object)
    stable_ids = np.asarray([normalize_ensembl_id(value) or str(value) for value in ids], dtype=object)
    if len(values) != len(ids):
        raise ValueError("Reference vector and feature IDs have different lengths")
    if len(values) == 0:
        return np.asarray([], dtype=np.int64), [np.asarray([], dtype=np.int64) for _ in range(int(bins))]
    if not np.isfinite(values).all():
        raise ValueError("Cannot make expression bins with unavailable reference abundance")
    # Lower reference abundance belongs to the first bin.  The direction is
    # part of the frozen matching contract; stable feature ID breaks ties.
    order = np.asarray(sorted(range(len(ids)), key=lambda index: (float(values[index]), str(stable_ids[index]), str(ids[index]))), dtype=np.int64)
    bin_arrays = [np.asarray(chunk, dtype=np.int64) for chunk in np.array_split(order, int(bins))]
    feature_bin = np.full(len(ids), -1, dtype=np.int64)
    for bin_number, chunk in enumerate(bin_arrays):
        feature_bin[chunk] = bin_number
    return feature_bin, bin_arrays


def matched_seed(base_seed: int, dataset_key: str, program_id: str) -> int:
    """Implement the frozen SHA-256 seed rule exactly."""

    text = f"{int(base_seed)}:{dataset_key}:{program_id}".encode("utf-8")
    return int(hashlib.sha256(text).hexdigest()[:16], 16)


@dataclass
class MatchedSets:
    available: bool
    target_indices: np.ndarray
    target_bin_counts: dict[int, int]
    pool_bin_counts: dict[int, int]
    selections: np.ndarray
    seed: int
    reason: str
    underflow_bins: tuple[int, ...]
    # Actual eligible control members by bin.  Keeping this alongside the
    # draw matrix lets public QC report the same pool used for sampling.
    pool_indices: dict[int, np.ndarray] = field(default_factory=dict)


def sample_expression_matched_sets(
    bin_arrays: Sequence[np.ndarray],
    target_indices: Sequence[int],
    excluded_indices: Iterable[int],
    *,
    replicates: int = DEFAULT_REPLICATES,
    seed: int,
) -> MatchedSets:
    """Draw fixed random sets matching target-bin counts without replacement.

    No fallback is attempted.  A bin underflow returns an unavailable result,
    which callers propagate as missing matched scores.
    """

    target = np.asarray(sorted(set(int(value) for value in target_indices)), dtype=np.int64)
    excluded = {int(value) for value in excluded_indices}
    target_bin_counts: Counter[int] = Counter()
    feature_to_bin: dict[int, int] = {}
    for bin_number, members in enumerate(bin_arrays):
        for index in np.asarray(members, dtype=np.int64):
            feature_to_bin[int(index)] = int(bin_number)
    for index in target:
        if int(index) not in feature_to_bin:
            return MatchedSets(False, target, {}, {}, np.empty((0, len(target)), dtype=np.int64), int(seed), "target_outside_bins", tuple())
        target_bin_counts[feature_to_bin[int(index)]] += 1
    if len(target) == 0:
        return MatchedSets(False, target, dict(target_bin_counts), {}, np.empty((0, 0), dtype=np.int64), int(seed), "empty_target", tuple())
    pools: dict[int, np.ndarray] = {}
    pool_bin_counts: dict[int, int] = {}
    underflow: list[int] = []
    for bin_number, members in enumerate(bin_arrays):
        pool = np.asarray([int(index) for index in np.asarray(members, dtype=np.int64) if int(index) not in excluded], dtype=np.int64)
        pools[bin_number] = pool
        pool_bin_counts[bin_number] = int(len(pool))
        required = int(target_bin_counts.get(bin_number, 0))
        if len(pool) < required:
            underflow.append(bin_number)
    if underflow:
        return MatchedSets(
            False,
            target,
            dict(target_bin_counts),
            pool_bin_counts,
            np.empty((0, len(target)), dtype=np.int64),
            int(seed),
            "insufficient_control_pool",
            tuple(underflow),
            pools,
        )
    selections = np.empty((int(replicates), len(target)), dtype=np.int64)
    target_positions_by_bin: dict[int, np.ndarray] = {}
    offset = 0
    for bin_number in range(len(bin_arrays)):
        required = int(target_bin_counts.get(bin_number, 0))
        if required:
            target_positions_by_bin[bin_number] = np.arange(offset, offset + required, dtype=np.int64)
            offset += required
    rng = np.random.default_rng(int(seed))
    for replicate in range(int(replicates)):
        for bin_number in range(len(bin_arrays)):
            required = int(target_bin_counts.get(bin_number, 0))
            if required:
                positions = target_positions_by_bin[bin_number]
                selections[replicate, positions] = np.sort(rng.choice(pools[bin_number], size=required, replace=False))
    return MatchedSets(True, target, dict(target_bin_counts), pool_bin_counts, selections, int(seed), "available", tuple(), pools)


def _feature_index_for_ets2(feature_map: pd.DataFrame) -> int | None:
    matches = feature_map.index[feature_map.feature_name.eq(ETS2_SYMBOL)].tolist()
    if len(matches) == 1:
        return int(feature_map.iloc[matches[0]].feature_index)
    stable = feature_map.index[feature_map.stable_feature_id.eq(ETS2_STABLE_ID)].tolist()
    return int(feature_map.iloc[stable[0]].feature_index) if len(stable) == 1 else None


def _score_with_matched_sets(log_values: np.ndarray, target_indices: np.ndarray, matched: MatchedSets) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return unadjusted target mean, matched mean and conditional percentile."""

    if len(target_indices):
        target_mean = log_values[:, target_indices].mean(axis=1)
    else:
        target_mean = np.full(log_values.shape[0], np.nan)
    if not matched.available:
        return target_mean, np.full(log_values.shape[0], np.nan), np.full(log_values.shape[0], np.nan)
    # The replicate-by-group temporary is bounded by groups x 1000.  Work in
    # replicate chunks so a large gene set never produces a groups x reps x k
    # dense object.
    background = np.empty((log_values.shape[0], matched.selections.shape[0]), dtype=float)
    for start in range(0, matched.selections.shape[0], 64):
        stop = min(start + 64, matched.selections.shape[0])
        for local, selection in enumerate(matched.selections[start:stop], start=start):
            background[:, local] = log_values[:, selection].mean(axis=1)
    mean_background = background.mean(axis=1)
    percentile = (background <= target_mean[:, None]).mean(axis=1)
    return target_mean, mean_background, percentile


def matching_qc_table(
    dataset_key: str,
    program_id: str,
    matched: MatchedSets,
    reference: np.ndarray,
    feature_bin: np.ndarray,
    *,
    program_gate_pass: bool,
    replicates: int,
    n_bins: int = DEFAULT_BINS,
) -> pd.DataFrame:
    """Create one matching-bin QC row per bin, including underflow."""

    rows: list[dict[str, Any]] = []
    for bin_number in range(int(n_bins)):
        target_count = int(matched.target_bin_counts.get(bin_number, 0))
        pool_count = int(matched.pool_bin_counts.get(bin_number, 0))
        members = np.flatnonzero(feature_bin == bin_number)
        pool_members = np.asarray(matched.pool_indices.get(bin_number, np.asarray([], dtype=np.int64)), dtype=np.int64)
        # ``pool_indices`` is the exact set passed to the sampler.  Keep the
        # all-bin abundance as a separate diagnostic so its name cannot be
        # mistaken for the eligible control-pool abundance.
        pool_reference_mean = float(np.mean(reference[pool_members])) if len(pool_members) else math.nan
        rows.append(
            {
                "dataset_key": dataset_key,
                "program_id": program_id,
                "bin": bin_number,
                "bin_feature_count": int(len(members)),
                "target_count": target_count,
                "control_pool_count": pool_count,
                "required_controls_sufficient": bool(pool_count >= target_count),
                "underflow": bool(bin_number in matched.underflow_bins),
                "program_mapping_gate_pass": bool(program_gate_pass),
                "matching_available": bool(matched.available and program_gate_pass),
                "replicates": int(replicates),
                # Some frozen 64-bit seeds exceed signed int64; keep the
                # public decimal representation exact instead of allowing
                # pandas to coerce the column to float in CSV output.
                "seed": str(int(matched.seed)),
                "target_reference_mean": float(np.mean(reference[matched.target_indices[feature_bin[matched.target_indices] == bin_number]])) if target_count else math.nan,
                "bin_reference_mean": float(np.mean(reference[members])) if len(members) else math.nan,
                "control_pool_reference_mean": pool_reference_mean,
                "reason": matched.reason if not matched.available else "available",
            }
        )
    return pd.DataFrame(rows)


def score_dataset(
    aggregate: AggregatedAtlas,
    program_inventory: pd.DataFrame,
    mapped_by_program: Mapping[str, set[int]],
    reference_result: ReferenceResult,
    *,
    thresholds: Sequence[int] = THRESHOLDS,
    bins: int = DEFAULT_BINS,
    replicates: int = DEFAULT_REPLICATES,
    base_seed: int = DEFAULT_BASE_SEED,
    all_mapped_original: set[int],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Score all programs for one modality and return donor rows plus QC."""

    pb = aggregate.pseudobulk
    if len(reference_result.reference):
        feature_bin, bin_arrays = expression_bins(reference_result.reference, aggregate.feature_map.feature_id, bins=bins)
    else:
        feature_bin, bin_arrays = np.full(aggregate.log_values.shape[1], -1, dtype=np.int64), [np.asarray([], dtype=np.int64) for _ in range(int(bins))]
    rows: list[pd.DataFrame] = []
    qc_tables: list[pd.DataFrame] = []
    score_qc: dict[str, Any] = {"dataset_key": aggregate.dataset_key, "reference": {
        "eligible_group_count": reference_result.eligible_group_count,
        "eligible_donor_count": reference_result.eligible_donor_count,
        "reference_available": bool(len(reference_result.reference) and np.isfinite(reference_result.reference).all()),
    }, "programs": {}}
    ets2_index = _feature_index_for_ets2(aggregate.feature_map)
    for program_id in PROGRAM_IDS:
        inventory_row = program_inventory[program_inventory.program_id.eq(program_id)]
        if len(inventory_row) != 1:
            raise ValueError(f"Program inventory is incomplete for {aggregate.dataset_key}/{program_id}")
        inventory = inventory_row.iloc[0].to_dict()
        target_indices = np.asarray(sorted(mapped_by_program.get(program_id, set())), dtype=np.int64)
        mapping_gate_pass = bool(inventory.get("mapping_gate_pass", False))
        excluded = set(all_mapped_original) | ({ets2_index} if ets2_index is not None else set())
        seed = matched_seed(base_seed, aggregate.dataset_key, program_id)
        if len(reference_result.reference) and mapping_gate_pass:
            matched = sample_expression_matched_sets(bin_arrays, target_indices, excluded, replicates=replicates, seed=seed)
        else:
            matched = MatchedSets(False, target_indices, {}, {}, np.empty((0, len(target_indices)), dtype=np.int64), seed, "reference_unavailable" if not len(reference_result.reference) else "mapping_gate_failed", tuple())
        qc = matching_qc_table(aggregate.dataset_key, program_id, matched, reference_result.reference if len(reference_result.reference) else np.zeros(aggregate.log_values.shape[1]), feature_bin, program_gate_pass=mapping_gate_pass, replicates=replicates, n_bins=bins)
        if len(qc):
            qc["intended_unique_source_genes"] = int(inventory.get("intended_unique_source_genes", 0))
            qc["mapped_unique_genes"] = int(inventory.get("mapped_unique_genes", 0))
            qc["mapping_fraction"] = float(inventory.get("mapping_fraction", 0.0))
        qc_tables.append(qc)
        target_mean, matched_mean, percentile = _score_with_matched_sets(aggregate.log_values, target_indices, matched)
        score = target_mean - matched_mean
        target_detected = (aggregate.count_values[:, target_indices] > 0).sum(axis=1) if len(target_indices) else np.zeros(len(pb), dtype=np.int64)
        if ets2_index is None:
            ets2_log = np.full(len(pb), np.nan)
        else:
            ets2_log = aggregate.log_values[:, ets2_index]
        baseline_target = float(np.mean(reference_result.reference[target_indices])) if len(target_indices) and len(reference_result.reference) else math.nan
        baseline_control = math.nan
        if matched.available and len(reference_result.reference):
            baseline_control = float(np.mean([reference_result.reference[selection].mean() for selection in matched.selections]))
        score_qc["programs"][program_id] = {
            "mapping_gate_pass": mapping_gate_pass,
            "target_count": int(len(target_indices)),
            "matching_available": bool(matched.available and mapping_gate_pass),
            "matching_reason": matched.reason,
            "underflow_bins": list(matched.underflow_bins),
            "seed": int(seed),
            "baseline_target_reference_mean": baseline_target if np.isfinite(baseline_target) else None,
            "baseline_control_reference_mean": baseline_control if np.isfinite(baseline_control) else None,
            "baseline_mean_mismatch": baseline_target - baseline_control if np.isfinite(baseline_target) and np.isfinite(baseline_control) else None,
            "target_bin_counts": {str(k): int(v) for k, v in sorted(matched.target_bin_counts.items())},
            "control_pool_bin_counts": {str(k): int(v) for k, v in sorted(matched.pool_bin_counts.items())},
            "selection_shape": [int(value) for value in matched.selections.shape],
            "selections_sha256": hashlib.sha256(matched.selections.tobytes(order="C")).hexdigest(),
        }
        if len(pb):
            frame = pb[[
                "dataset_key", "dataset_id", "modality", "suspension_type", "assay", "donor_id", "cohort", "author_cell_type", "group_id", "cell_count", "all_feature_umis", "positive_all_feature_umi_denominator"
            ]].copy()
            frame["program_id"] = program_id
            frame["threshold"] = np.nan
            # Threshold is added in the loop below; duplicate calculations are
            # intentional because the primary random sets are reused at 50.
            frame["program_mapping_status"] = "evaluable" if mapping_gate_pass else "mapping_gate_failed"
            frame["source_intended_genes"] = int(inventory.get("intended_unique_source_genes", 0))
            frame["mapped_program_genes"] = int(len(target_indices))
            frame["genes_detected_any_umi"] = target_detected.astype(int)
            frame["program_mean_log2_cpm1"] = target_mean
            frame["matched_background_mean_log2_cpm1"] = matched_mean
            frame["program_score"] = score
            frame["conditional_percentile"] = percentile
            frame["ets2_log2_cpm1"] = ets2_log
            frame["matching_available"] = bool(matched.available and mapping_gate_pass)
            frame["matching_reason"] = matched.reason
            for threshold in thresholds:
                current = frame.copy()
                current["threshold"] = int(threshold)
                current["zero_all_feature_denominator"] = ~current.positive_all_feature_umi_denominator
                current["cell_gate_pass"] = current.cell_count.ge(int(threshold)) & ~current.zero_all_feature_denominator
                current["raw_mean_available"] = bool(mapping_gate_pass) & current.cell_gate_pass
                current["score_available"] = current.raw_mean_available & bool(matched.available)
                current.loc[~current.cell_gate_pass, ["program_score", "conditional_percentile", "program_mean_log2_cpm1", "matched_background_mean_log2_cpm1"]] = np.nan
                current.loc[~current.raw_mean_available, ["program_mean_log2_cpm1"]] = np.nan
                current["score_status"] = np.where(
                    current.zero_all_feature_denominator,
                    "zero_all_feature_denominator",
                    np.where(
                        ~current.cell_count.ge(int(threshold)),
                        "insufficient_cell_gate",
                        np.where(
                            not mapping_gate_pass,
                            "mapping_gate_failed",
                            np.where(~current.matching_available, "matched_controls_unavailable", "evaluable_descriptive"),
                        ),
                    ),
                )
                rows.append(current)
    donor_scores = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(donor_scores):
        donor_scores = donor_scores.sort_values(["dataset_key", "modality", "assay", "author_cell_type", "program_id", "threshold", "donor_id"], kind="mergesort").reset_index(drop=True)
    matching_qc = pd.concat(qc_tables, ignore_index=True) if qc_tables else pd.DataFrame()
    return donor_scores, matching_qc, score_qc


def _finite_values(frame: pd.DataFrame, column: str) -> np.ndarray:
    if column not in frame:
        return np.asarray([], dtype=float)
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def summarize_contexts(donor_scores: pd.DataFrame, aggregate: AggregatedAtlas, program_inventory: pd.DataFrame, thresholds: Sequence[int] = THRESHOLDS) -> pd.DataFrame:
    """Summarize donor-level program rows without exposing donor IDs."""

    if len(aggregate.pseudobulk) == 0:
        return pd.DataFrame()
    context_columns = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay"]
    output_columns = context_columns + ["author_cell_type"]
    observed_contexts = aggregate.pseudobulk[context_columns].drop_duplicates().sort_values(context_columns, kind="mergesort")
    records: list[dict[str, Any]] = []
    for base in observed_contexts.to_dict("records"):
        # The output is a Cartesian product of every observed dataset/assay
        # context and all four frozen author labels.  Missing labels remain
        # explicit rows with ``no_observed_context``.
        for author_cell_type in AUTHOR_CELL_TYPES:
            context = {**base, "author_cell_type": author_cell_type}
            subset_base = donor_scores
            for col, value in context.items():
                if col in subset_base:
                    subset_base = subset_base[subset_base[col].eq(value)]
                else:
                    subset_base = subset_base.iloc[0:0]
            observed = aggregate.pseudobulk
            for col, value in context.items():
                observed = observed[observed[col].eq(value)]
            observed_donors = int(observed["donor_id"].nunique())
            for program_id in PROGRAM_IDS:
                inv = program_inventory[(program_inventory.dataset_key.eq(str(base["dataset_key"]))) & (program_inventory.program_id.eq(program_id))]
                inventory = inv.iloc[0].to_dict() if len(inv) else {}
                mapping_gate_pass = bool(inventory.get("mapping_gate_pass", False))
                for threshold in thresholds:
                    rows = subset_base[subset_base.program_id.eq(program_id) & subset_base.threshold.eq(int(threshold))]
                    if "raw_mean_available" in rows:
                        raw_eligible = rows[rows.raw_mean_available.eq(True)]
                    else:
                        raw_eligible = rows[rows.score_status.isin(["evaluable_descriptive", "matched_controls_unavailable"])]
                    score_eligible = rows[rows.score_status.eq("evaluable_descriptive")]
                    score_values = _finite_values(score_eligible, "program_score")
                    target_values = _finite_values(raw_eligible, "program_mean_log2_cpm1")
                    background_values = _finite_values(score_eligible, "matched_background_mean_log2_cpm1")
                    percentile_values = _finite_values(score_eligible, "conditional_percentile")
                    detected = pd.to_numeric(raw_eligible.get("genes_detected_any_umi", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
                    detected = detected[np.isfinite(detected)]
                    zero_rows = rows[rows.score_status.eq("zero_all_feature_denominator")] if len(rows) else rows
                    cell_failed_rows = rows[rows.score_status.eq("insufficient_cell_gate")] if len(rows) else rows
                    if not len(rows):
                        status = "no_observed_context"
                    elif not mapping_gate_pass:
                        status = "mapping_gate_failed"
                    elif not len(raw_eligible):
                        status = "zero_all_feature_denominator" if len(zero_rows) == len(rows) and len(rows) else "no_eligible_donors"
                    elif not len(score_eligible):
                        status = "matched_controls_unavailable"
                    else:
                        status = "evaluable_descriptive"
                    raw_status = "mapping_gate_failed" if not mapping_gate_pass else (
                        "evaluable_descriptive" if len(raw_eligible) else
                        "zero_all_feature_denominator" if len(zero_rows) == len(rows) and len(rows) else
                        "insufficient_cell_gate" if len(cell_failed_rows) == len(rows) and len(rows) else
                        "no_eligible_donors"
                    )
                    records.append(
                        {
                            **context,
                            "program_id": program_id,
                            "threshold": int(threshold),
                            "status": status,
                            "raw_summary_status": raw_status,
                            "observed_donors": observed_donors,
                            # ``eligible_donors`` is the cell-gate donor unit;
                            # score-eligible donors are reported separately.
                            "eligible_donors": int(raw_eligible.donor_id.nunique()),
                            "score_eligible_donors": int(score_eligible.donor_id.nunique()),
                            "observed_groups": int(len(rows)),
                            "zero_denominator_groups": int(len(zero_rows)),
                            "cell_gate_failed_groups": int(len(cell_failed_rows)),
                            "source_intended_genes": int(inventory.get("intended_unique_source_genes", 0)),
                            "mapped_program_genes": int(inventory.get("mapped_unique_genes", 0)),
                            "mapping_fraction": float(inventory.get("mapping_fraction", 0.0)),
                            "donors_with_any_program_umi": int((detected > 0).sum()) if len(detected) else 0,
                            "median_genes_detected": float(np.median(detected)) if len(detected) else math.nan,
                            "minimum_genes_detected": float(np.min(detected)) if len(detected) else math.nan,
                            "maximum_genes_detected": float(np.max(detected)) if len(detected) else math.nan,
                            "median_program_mean_log2_cpm1": float(np.median(target_values)) if len(target_values) else math.nan,
                            "minimum_program_mean_log2_cpm1": float(np.min(target_values)) if len(target_values) else math.nan,
                            "maximum_program_mean_log2_cpm1": float(np.max(target_values)) if len(target_values) else math.nan,
                            "median_matched_background_log2_cpm1": float(np.median(background_values)) if len(background_values) else math.nan,
                            "median_program_score": float(np.median(score_values)) if len(score_values) else math.nan,
                            "minimum_program_score": float(np.min(score_values)) if len(score_values) else math.nan,
                            "maximum_program_score": float(np.max(score_values)) if len(score_values) else math.nan,
                            "positive_score_count": int((score_values > 0).sum()) if len(score_values) else 0,
                            "median_conditional_percentile": float(np.median(percentile_values)) if len(percentile_values) else math.nan,
                            "minimum_conditional_percentile": float(np.min(percentile_values)) if len(percentile_values) else math.nan,
                            "maximum_conditional_percentile": float(np.max(percentile_values)) if len(percentile_values) else math.nan,
                        }
                    )
    return pd.DataFrame(records).sort_values(output_columns + ["program_id", "threshold"], kind="mergesort").reset_index(drop=True)


def _paired_summary(deltas: np.ndarray, min_pairs: int) -> dict[str, Any]:
    if len(deltas) < min_pairs:
        return {
            "status": "insufficient_paired_donors",
            "paired_donors": int(len(deltas)),
            "median_delta": math.nan,
            "minimum_delta": math.nan,
            "maximum_delta": math.nan,
            "positive_delta_count": 0,
            "loo_median_min": math.nan,
            "loo_median_max": math.nan,
            "loo_evaluable_count": 0,
        }
    loo = [float(np.median(np.delete(deltas, i))) for i in range(len(deltas)) if len(deltas) - 1 >= min_pairs]
    return {
        "status": "evaluable_descriptive",
        "paired_donors": int(len(deltas)),
        "median_delta": float(np.median(deltas)),
        "minimum_delta": float(np.min(deltas)),
        "maximum_delta": float(np.max(deltas)),
        "positive_delta_count": int((deltas > 0).sum()),
        "loo_median_min": float(np.min(loo)) if loo else math.nan,
        "loo_median_max": float(np.max(loo)) if loo else math.nan,
        "loo_evaluable_count": int(len(loo)),
    }


def build_paired_contrasts(
    donor_scores: pd.DataFrame,
    aggregate_frames: Mapping[str, AggregatedAtlas],
    plan: Mapping[str, Any],
) -> pd.DataFrame:
    """Build exact within-donor-and-assay paired population contrasts."""

    if len(donor_scores) == 0:
        return pd.DataFrame()
    contexts = donor_scores[["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "program_id", "threshold"]].drop_duplicates()
    records: list[dict[str, Any]] = []
    contrast_specs = [(str(item["left"]), str(item["right"]), str(item.get("role", ""))) for item in plan["contrasts"]]
    for base in contexts.to_dict("records"):
        subset = donor_scores
        for column, value in base.items():
            subset = subset[subset[column].eq(value)]
        for left, right, role in contrast_specs:
            left_observed = subset[subset.author_cell_type == left]
            right_observed = subset[subset.author_cell_type == right]
            left_rows = left_observed[left_observed.score_status.eq("evaluable_descriptive")]
            right_rows = right_observed[right_observed.score_status.eq("evaluable_descriptive")]
            left_keys = ["dataset_key", "dataset_id", "modality", "suspension_type", "assay", "program_id", "threshold", "donor_id"]
            if left_observed.duplicated(left_keys).any() or right_observed.duplicated(left_keys).any():
                raise ValueError("Duplicate donor/assay rows prevent paired contrast")
            merged = left_rows.merge(right_rows, on=left_keys, suffixes=("_left", "_right"), how="inner")
            deltas = (merged["program_score_left"].to_numpy(dtype=float) - merged["program_score_right"].to_numpy(dtype=float)) if len(merged) else np.asarray([], dtype=float)
            summary = _paired_summary(deltas, int(plan["minimum_paired_donors"]))
            row = {
                **base,
                "left_population": left,
                "right_population": right,
                "contrast_role": role,
                "comparison": f"{left}_minus_{right}",
                "left_observed_donors": int(left_observed.donor_id.nunique()),
                "right_observed_donors": int(right_observed.donor_id.nunique()),
                "left_score_eligible_donors": int(left_rows.donor_id.nunique()),
                "right_score_eligible_donors": int(right_rows.donor_id.nunique()),
                "score_available_pairs": int(len(merged)),
                **summary,
                "interpretation": "descriptive_within_donor_assay_only",
            }
            records.append(row)
    return pd.DataFrame(records).sort_values(["dataset_key", "modality", "assay", "program_id", "threshold", "left_population"], kind="mergesort").reset_index(drop=True)


def _spearman(values_x: np.ndarray, values_y: np.ndarray) -> float:
    # Rank-average ties without depending on scipy.stats API details.
    x_rank = pd.Series(values_x).rank(method="average").to_numpy(dtype=float)
    y_rank = pd.Series(values_y).rank(method="average").to_numpy(dtype=float)
    if np.std(x_rank) == 0 or np.std(y_rank) == 0:
        return math.nan
    return float(np.corrcoef(x_rank, y_rank)[0, 1])


def build_correlations(donor_scores: pd.DataFrame, plan: Mapping[str, Any]) -> pd.DataFrame:
    """Build descriptive donor-level Spearman correlations for each context."""

    if len(donor_scores) == 0:
        return pd.DataFrame()
    primary = donor_scores[(donor_scores.program_id == PRIMARY_PROGRAM_ID) & donor_scores.score_status.eq("evaluable_descriptive")].copy()
    ets2 = donor_scores[donor_scores.ets2_log2_cpm1.notna()].copy()
    # Collapse duplicate assays within donor before the donor gate.  This
    # retains donor units even when the atlas has more than one assay label.
    donor_keys = ["dataset_key", "dataset_id", "modality", "suspension_type", "author_cell_type", "threshold", "donor_id"]
    primary_donor = primary.groupby(donor_keys, sort=True, as_index=False)["program_score"].mean().rename(columns={"program_score": "primary_score"})
    ets2_donor = ets2.groupby(donor_keys, sort=True, as_index=False)["ets2_log2_cpm1"].mean()
    program_donor: dict[str, pd.DataFrame] = {}
    for program_id in COMPARATOR_PROGRAM_IDS:
        frame = donor_scores[(donor_scores.program_id == program_id) & donor_scores.score_status.eq("evaluable_descriptive")]
        program_donor[program_id] = frame.groupby(donor_keys, sort=True, as_index=False)["program_score"].mean().rename(columns={"program_score": "comparator_score"})
    # Correlations are defined within population and modality.  If a modality
    # has multiple assay labels, donor-level values are averaged across those
    # assays before the donor gate, so assay is deliberately not an output
    # stratum here.
    context_columns = ["dataset_key", "dataset_id", "modality", "suspension_type", "threshold"]
    contexts = donor_scores[context_columns].drop_duplicates().sort_values(context_columns, kind="mergesort")
    records: list[dict[str, Any]] = []
    min_donors = int(plan["minimum_correlation_donors"])
    for base_context in contexts.to_dict("records"):
        # Keep all four frozen populations for every observed dataset/assay
        # context.  A missing label (for example ActMac/LAM-like in nuclei)
        # is an explicit unavailable row rather than a silently omitted test.
        for author_cell_type in AUTHOR_CELL_TYPES:
            context = {**base_context, "author_cell_type": author_cell_type}
            base = primary_donor
            for column, value in context.items():
                if column in base:
                    base = base[base[column].eq(value)]
            base = base.merge(ets2_donor, on=donor_keys, how="inner")
            observed = donor_scores
            for column, value in context.items():
                observed = observed[observed[column].eq(value)]
            observed_donors = int(observed.donor_id.nunique())
            comparisons: list[tuple[str, pd.DataFrame, str]] = [("ets2_log2_cpm1", base, "primary_program_vs_ETS2_log2CPM1")]
            for program_id in COMPARATOR_PROGRAM_IDS:
                comparator = program_donor[program_id]
                joined = base.merge(comparator, on=donor_keys, how="inner")
                comparisons.append((program_id, joined, f"primary_program_vs_{program_id}"))
            for comparison_id, joined, comparison_label in comparisons:
                if comparison_id == "ets2_log2_cpm1":
                    x = joined.primary_score.to_numpy(dtype=float) if len(joined) else np.asarray([], dtype=float)
                    y = joined.ets2_log2_cpm1.to_numpy(dtype=float) if len(joined) else np.asarray([], dtype=float)
                else:
                    x = joined.primary_score.to_numpy(dtype=float) if len(joined) else np.asarray([], dtype=float)
                    y = joined.comparator_score.to_numpy(dtype=float) if len(joined) else np.asarray([], dtype=float)
                finite = np.isfinite(x) & np.isfinite(y)
                x, y = x[finite], y[finite]
                if observed_donors == 0:
                    status = "no_observed_context"
                    rho = math.nan
                elif len(x) < min_donors:
                    status = "insufficient_donors"
                    rho = math.nan
                elif np.std(x) == 0 or np.std(y) == 0:
                    status = "constant_vector"
                    rho = math.nan
                else:
                    status = "evaluable_descriptive"
                    rho = _spearman(x, y)
                records.append(
                    {
                        **context,
                        "comparison_id": comparison_id,
                        "comparison": comparison_label,
                        "status": status,
                        "observed_donors": observed_donors,
                        "donors": int(len(x)),
                        "spearman_rho": rho,
                        "assays_collapsed_within_donor": True,
                        "interpretation": "descriptive_small_n_no_p_value",
                    }
                )
    return pd.DataFrame(records).sort_values(["dataset_key", "modality", "author_cell_type", "threshold", "comparison_id"], kind="mergesort").reset_index(drop=True)


def _write_csv(frame: pd.DataFrame, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")
    return sha256_file(path)


def _relative_or_string(path: Path, root: Path) -> str:
    return str(path.relative_to(root)) if path.is_relative_to(root) else str(path)


def _write_ignored_count_cache(aggregate: AggregatedAtlas, root: Path) -> dict[str, str]:
    """Write inspectable ignored all-feature counts plus compact metadata."""

    cache_directory = root / "work"
    cache_directory.mkdir(parents=True, exist_ok=True)
    count_path = cache_directory / f"ets2-program-all-feature-counts-{aggregate.dataset_key}.npz"
    np.savez_compressed(
        count_path,
        counts=aggregate.count_values,
        feature_index=aggregate.feature_map.feature_index.to_numpy(dtype=np.int64),
        feature_id=aggregate.feature_map.feature_id.astype(str).to_numpy(dtype=str),
        donor_id=aggregate.pseudobulk.donor_id.astype(str).to_numpy(dtype=str) if len(aggregate.pseudobulk) else np.asarray([], dtype=str),
        group_id=aggregate.pseudobulk.group_id.astype(str).to_numpy(dtype=str) if len(aggregate.pseudobulk) else np.asarray([], dtype=str),
    )
    feature_path = cache_directory / f"ets2-program-feature-map-{aggregate.dataset_key}.csv"
    feature_hash = _write_csv(aggregate.feature_map, feature_path)
    group_columns = [
        "dataset_key", "dataset_id", "modality", "suspension_type", "assay", "donor_id", "cohort",
        "author_cell_type", "group_id", "cell_count", "all_feature_umis", "positive_all_feature_umi_denominator",
    ]
    group_path = cache_directory / f"ets2-program-donor-groups-{aggregate.dataset_key}.csv"
    group_hash = _write_csv(aggregate.pseudobulk[[column for column in group_columns if column in aggregate.pseudobulk.columns]], group_path)
    return {
        _relative_or_string(count_path, root): sha256_file(count_path),
        _relative_or_string(feature_path, root): feature_hash,
        _relative_or_string(group_path, root): group_hash,
    }


def run_analysis(
    *,
    root: Path = ROOT,
    plan_path: Path | None = None,
    dataset_paths: Mapping[str, Path] | None = None,
    source_zip: Path | None = None,
    output_directory: Path | None = None,
    audit_path: Path | None = None,
    chunk_size: int = 4096,
    write_donor_cache: bool = True,
) -> dict[str, Any]:
    """Run the frozen ETS2 program analysis and write aggregate outputs."""

    plan_path = plan_path or (root / "config/ets2-program-plan.json")
    plan, program_manifest, plan_hashes = load_frozen_plan(plan_path, root=root)
    output_directory = output_directory or (root / "data/derived")
    audit_path = audit_path or (root / "reports/ets2-program-audit.json")
    benchmark_manifest = _read_json(root / "config/ets2-benchmark-sources.json")
    benchmark_manifest_path = root / "config/ets2-benchmark-sources.json"
    expected_benchmark_hash = plan.get("input_sha256", {}).get("config/ets2-benchmark-sources.json")
    if expected_benchmark_hash and sha256_file(benchmark_manifest_path) != str(expected_benchmark_hash):
        raise ValueError("ETS2 benchmark source manifest changed")
    source_zip = source_zip or _source_zip_from_manifest(root, benchmark_manifest)
    if source_zip != (root / str(benchmark_manifest.get("source_directory", "data/raw/ets2-benchmark")) / "ets2-public-release.zip"):
        # CLI replays may point to a copy, but the bytes must still match the
        # frozen source manifest record.
        source_record = next(item for item in benchmark_manifest["sources"] if item["file"] == "ets2-public-release.zip")
        verify_pinned_source(source_zip, source_record)
    source_values, source_records = load_source_programs(plan, source_zip)
    selected_members = [str(item["source_member"]) for item in plan["programs"] if "source_member" in item]
    member_hashes = source_member_fingerprints(source_zip, selected_members)

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("anndata is required for the ETS2 atlas analysis") from exc

    aggregates: dict[str, AggregatedAtlas] = {}
    for spec in plan["datasets"]:
        dataset_key = str(spec["key"])
        if dataset_paths and dataset_key in dataset_paths:
            path = Path(dataset_paths[dataset_key])
        else:
            manifest = _read_json(root / "config/liver-atlas-sources.json")
            path = root / str(manifest["source_directory"]) / str(spec["file"])
            source_record = next((item for item in manifest["sources"] if str(item["file"]) == str(spec["file"])), None)
            if source_record is None:
                raise ValueError(f"Liver-atlas manifest has no record for {spec['file']}")
            verify_pinned_source(path, source_record)
        data = ad.read_h5ad(path, backed="r")
        try:
            aggregates[dataset_key] = aggregate_all_features(
                data,
                dataset_key,
                dataset_id=str(spec["dataset_id"]),
                version_id=str(spec.get("version_id", "")),
                source_file=_relative_or_string(path, root),
                expected_observations=int(spec.get("expected_observations", 0)) or None,
                author_cell_types=plan["author_cell_types"],
                cohort=str(plan["cohort"]),
                thresholds=plan["cell_thresholds"],
                chunk_size=chunk_size,
            )
        finally:
            manager = getattr(data, "file", None)
            if manager is not None and hasattr(manager, "close"):
                manager.close()

    feature_maps = {key: aggregate.feature_map for key, aggregate in aggregates.items()}
    membership, inventory, mapped_by_dataset = build_program_mappings(plan, source_values, feature_maps)
    donor_frames: list[pd.DataFrame] = []
    matching_frames: list[pd.DataFrame] = []
    reference_audit: dict[str, Any] = {}
    score_qc_by_dataset: dict[str, Any] = {}
    for dataset_key, aggregate in aggregates.items():
        reference = compute_reference(aggregate, min_cells=int(plan["cell_thresholds"][0]))
        reference_audit[dataset_key] = {
            "eligible_group_count": reference.eligible_group_count,
            "eligible_donor_count": reference.eligible_donor_count,
            "reference_available": bool(len(reference.reference) and np.isfinite(reference.reference).all()),
        }
        inv = inventory[inventory.dataset_key.eq(dataset_key)]
        all_mapped_original = set().union(*(mapped_by_dataset[dataset_key].get(program_id, set()) for program_id in ORIGINAL_PROGRAM_IDS))
        donor_scores, matching_qc, score_qc = score_dataset(
            aggregate,
            inv,
            mapped_by_dataset[dataset_key],
            reference,
            thresholds=plan["cell_thresholds"],
            bins=int(plan["matching"]["bins"]),
            replicates=int(plan["matching"]["replicates"]),
            base_seed=int(plan["matching"]["base_seed"]),
            all_mapped_original=all_mapped_original,
        )
        donor_frames.append(donor_scores)
        matching_frames.append(matching_qc)
        score_qc_by_dataset[dataset_key] = score_qc
    donor_scores = pd.concat(donor_frames, ignore_index=True) if donor_frames else pd.DataFrame()
    matching_qc = pd.concat(matching_frames, ignore_index=True) if matching_frames else pd.DataFrame()
    contexts = pd.concat([summarize_contexts(donor_scores[donor_scores.dataset_key.eq(key)], aggregates[key], inventory[inventory.dataset_key.eq(key)], plan["cell_thresholds"]) for key in aggregates], ignore_index=True)
    contrasts = build_paired_contrasts(donor_scores, aggregates, plan)
    correlations = build_correlations(donor_scores, plan)

    outputs: dict[str, pd.DataFrame] = {
        "ets2-program-mapping.csv": inventory,
        "ets2-program-membership.csv": membership,
        "ets2-program-matching-qc.csv": matching_qc,
        "ets2-program-context-summary.csv": contexts,
        "ets2-program-contrasts.csv": contrasts,
        "ets2-program-correlations.csv": correlations,
    }
    output_hashes: dict[str, str] = {}
    for filename, frame in outputs.items():
        path = output_directory / filename
        output_hashes[_relative_or_string(path, root)] = _write_csv(frame, path)
    donor_cache_path = root / "work/ets2-program-donor-scores.csv"
    donor_cache_hash = None
    ignored_count_hashes: dict[str, str] = {}
    for aggregate in aggregates.values():
        ignored_count_hashes.update(_write_ignored_count_cache(aggregate, root))
    if write_donor_cache:
        donor_cache_hash = _write_csv(donor_scores, donor_cache_path)

    dataset_audit: list[dict[str, Any]] = []
    for key, aggregate in aggregates.items():
        dataset_audit.append({**aggregate.qc, "feature_ids": int(len(aggregate.feature_map)), "donor_cache_rows": int(len(donor_scores[donor_scores.dataset_key.eq(key)]))})
    audit: dict[str, Any] = {
        "analysis_id": str(plan["analysis_id"]),
        "plan_sha256": sha256_file(plan_path),
        "source_manifest_sha256": sha256_file(root / "config/ets2-program-sources.json"),
        "benchmark_manifest_sha256": sha256_file(benchmark_manifest_path),
        "atlas_manifest_sha256": sha256_file(root / "config/liver-atlas-sources.json"),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "atlas_helper_script_sha256": sha256_file(ATLAS_SCRIPT),
        "requirements_sha256": sha256_file(root / "requirements-ets2-program.txt"),
        "plan_input_sha256": plan_hashes,
        "source_zip": {"path": _relative_or_string(source_zip, root), "bytes": source_zip.stat().st_size, "sha256": sha256_file(source_zip)},
        "source_members": member_hashes,
        "source_column_inventory": source_records,
        "datasets": dataset_audit,
        "reference": reference_audit,
        "program_inventory_rows": int(len(inventory)),
        "matching_qc_rows": int(len(matching_qc)),
        "score_qc": score_qc_by_dataset,
        "context_rows": int(len(contexts)),
        "contrast_status_counts": {str(key): int(value) for key, value in contrasts.status.value_counts().items()} if len(contrasts) else {},
        "correlation_status_counts": {str(key): int(value) for key, value in correlations.status.value_counts().items()} if len(correlations) else {},
        "outputs_sha256": output_hashes,
        "ignored_outputs_sha256": {
            **ignored_count_hashes,
            **({"work/ets2-program-donor-scores.csv": donor_cache_hash} if donor_cache_hash else {}),
        },
        "count_source": "raw/X",
        "cohort": "PSC",
        "author_cell_types": list(AUTHOR_CELL_TYPES),
        "thresholds": list(THRESHOLDS),
        "limitations": "Donor-level descriptive PSC context only. Matched-set percentiles are calibration summaries, not population P-values. Negative program scores mean relative background centering, not biological inhibition. No disease comparison, cell-level test, causal mediation, treatment effect or clinical recommendation is estimated.",
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(_json_safe(audit), indent=2, sort_keys=True, allow_nan=False, default=_json_default) + "\n", encoding="utf-8")
    return audit


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=ROOT / "config/ets2-program-plan.json")
    parser.add_argument("--dataset", action="append", metavar="KEY=PATH", help="Override a pinned H5AD path; may be repeated")
    parser.add_argument("--source-zip", type=Path, default=None, help="Replay from a byte-identical copy of the pinned release ZIP")
    parser.add_argument("--output-directory", type=Path, default=None)
    parser.add_argument("--audit-path", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--no-donor-cache", action="store_true", help="Do not write ignored donor-level score cache")
    args = parser.parse_args(argv)
    overrides: dict[str, Path] = {}
    for item in args.dataset or []:
        if "=" not in item:
            parser.error("--dataset must be KEY=PATH")
        key, value = item.split("=", 1)
        if key not in DATASET_KEYS or not value:
            parser.error("--dataset must use KEY=PATH for sc or sn")
        overrides[key] = Path(value)
    audit = run_analysis(
        plan_path=args.plan,
        dataset_paths=overrides or None,
        source_zip=args.source_zip,
        output_directory=args.output_directory,
        audit_path=args.audit_path,
        chunk_size=args.chunk_size,
        write_donor_cache=not args.no_donor_cache,
    )
    print(json.dumps(audit, indent=2, sort_keys=True, default=_json_default))


if __name__ == "__main__":
    main()
