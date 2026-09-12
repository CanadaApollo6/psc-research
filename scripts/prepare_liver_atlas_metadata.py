"""Prepare an obs-only donor/annotation audit for the Andrews PSC atlas.

The two CELLxGENE H5AD inputs are opened in backed mode and only ``.obs`` is
copied.  This script deliberately does not touch ``X`` or ``raw.X``.  GEO
fields are retained as supplied by the existing imported table; title-derived
fields are used only as a label-level donor crosswalk and never as proof of a
shared library.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GEO = ROOT / "data/derived/geo-library-metadata.csv"
DEFAULT_METADATA_CONFIG = ROOT / "config/liver-atlas-metadata.json"
DEFAULT_SOURCES = ROOT / "config/liver-atlas-sources.json"
DEFAULT_PLAN = ROOT / "config/liver-atlas-plan.json"
DEFAULT_CSV = ROOT / "data/derived/liver-atlas-libraries.csv"
DEFAULT_REPORT = ROOT / "reports/liver-atlas-metadata-audit.json"

REQUIRED_OBS = ("donor_id", "assay", "author_cell_type", "disease")
ANNOTATION_FIELDS = (
    "author_cell_type",
    "cell_type",
    "cell_type_ontology_term_id",
)
DONOR_FIELDS = (
    "disease",
    "assay",
    "assay_ontology_term_id",
    "suspension_type",
    "library_uuid",
    "sample_uuid",
    "mapped_reference_annotation",
    "alignment_software",
    "suspension_dissociation_reagent",
)
GEO_RAW_FIELDS = (
    "series",
    "sample",
    "title",
    "donor_label_from_title",
    "cohort_from_title",
    "assay_from_title",
    "chemistry_from_title",
    "source_name",
    "description",
    "flags",
    "source_url",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_source_file(path: Path, declaration: Mapping[str, Any], label: str | None = None) -> dict[str, Any]:
    """Require exact declared bytes and SHA-256 before a source is parsed."""

    source_label = label or str(declaration.get("id") or path)
    if "bytes" not in declaration or "sha256" not in declaration:
        raise ValueError(f"Pinned source declaration is incomplete: {source_label}")
    if not path.exists():
        raise FileNotFoundError(path)
    expected_bytes = int(declaration["bytes"])
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(
            f"Pinned byte count mismatch for {source_label}: expected {expected_bytes}, got {actual_bytes}"
        )
    expected_sha256 = str(declaration["sha256"])
    actual_sha256 = sha256_file(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"Pinned SHA-256 mismatch for {source_label}: expected {expected_sha256}, got {actual_sha256}"
        )
    return {
        "id": source_label,
        "file": str(path.relative_to(ROOT)),
        "bytes": actual_bytes,
        "sha256": actual_sha256,
        "declared_bytes": expected_bytes,
        "declared_sha256": expected_sha256,
        "status": "verified_exact_bytes_and_sha256",
    }


def json_compact(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def text_value(value: Any) -> str:
    """Return a stable printable value while keeping source strings intact."""

    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value)


def normalized_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("").str.strip()


def unique_values(series: pd.Series) -> list[str]:
    values = {text_value(value) for value in series.tolist()}
    return sorted(value for value in values if value)


def joined_values(values: Iterable[str]) -> str:
    return ";".join(value for value in values if value)


def value_counts(series: pd.Series) -> dict[str, int]:
    values = series.astype("string").fillna("<missing>")
    values = values.replace({"": "<missing>"})
    counts = values.value_counts(dropna=False)
    return {str(key): int(counts[key]) for key in sorted(counts.index, key=str)}


def missing_mask(series: pd.Series) -> pd.Series:
    values = series.astype("string")
    return values.isna() | values.str.strip().fillna("").eq("")


def unknown_ontology_mask(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.strip().str.lower()
    return values.isna() | values.isin(("", "unknown", "--"))


def doublet_or_hybrid_mask(series: pd.Series) -> pd.Series:
    return series.astype("string").str.contains(r"doublet|--", case=False, regex=True, na=False)


def bool_true_mask(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.lower().eq("true")


def parse_title_fields(title: Any) -> dict[str, str]:
    """Parse only explicit title tokens; this does not alter the raw title."""

    raw = text_value(title)
    match = re.match(r"^([A-Za-z]+\d+[A-Za-z]*)(?:[_-].*)?$", raw)
    donor = match.group(1) if match else ""
    upper = raw.upper()
    if re.search(r"(?:^|[_-])SN(?:[_-]|$)", upper):
        suspension = "single-nucleus"
    elif re.search(r"(?:^|[_-])SC(?:[_-]|$)", upper):
        suspension = "single-cell"
    else:
        suspension = ""
    if re.search(r"5PR|5['’]", upper):
        chemistry = "5-prime"
    elif re.search(r"3PR|3['’]", upper):
        chemistry = "3-prime"
    else:
        chemistry = ""
    return {"donor": donor, "suspension": suspension, "chemistry": chemistry}


def description_chemistry(description: Any) -> str:
    value = text_value(description).lower()
    if re.search(r"5\s*pr|5['’]\s*", value):
        return "5-prime"
    if re.search(r"3\s*pr|3['’]\s*", value):
        return "3-prime"
    return ""


def description_suspension(description: Any) -> str:
    value = text_value(description).lower()
    if "single-nucleus" in value or "nuclei" in value:
        return "single-nucleus"
    if "scrna" in value or "single-cell" in value:
        return "single-cell"
    return ""


def read_geo_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [field for field in GEO_RAW_FIELDS if field not in table.columns]
    if missing:
        raise ValueError(f"GEO metadata is missing columns: {missing}")
    table = table.copy()
    parsed = table["title"].map(parse_title_fields)
    table["parsed_title_donor"] = [item["donor"] for item in parsed]
    table["parsed_title_suspension"] = [item["suspension"] for item in parsed]
    table["parsed_title_chemistry"] = [item["chemistry"] for item in parsed]
    table["description_chemistry"] = table["description"].map(description_chemistry)
    table["description_suspension"] = table["description"].map(description_suspension)
    table["title_description_chemistry_conflict"] = (
        table["chemistry_from_title"].map(text_value).ne("")
        & table["description_chemistry"].ne("")
        & table["chemistry_from_title"].ne(table["description_chemistry"])
    )
    table["title_description_suspension_conflict"] = (
        table["assay_from_title"].map(text_value).ne("")
        & table["description_suspension"].ne("")
        & table["assay_from_title"].ne(table["description_suspension"])
    )
    return table


def build_geo_index(table: pd.DataFrame) -> dict[str, list[dict[str, str]]]:
    """Index raw GEO records by the provisional donor label in the source CSV."""

    index: dict[str, list[dict[str, str]]] = {}
    for _, row in table.sort_values(["series", "sample"]).iterrows():
        donor = text_value(row.get("donor_label_from_title")) or text_value(row.get("parsed_title_donor"))
        if not donor:
            continue
        record = {field: text_value(row.get(field)) for field in GEO_RAW_FIELDS}
        index.setdefault(donor, []).append(record)
    return index


def annotation_catalog(obs: pd.DataFrame) -> list[dict[str, Any]]:
    """List original author labels alongside the supplied ontology labels."""

    available = [field for field in ANNOTATION_FIELDS if field in obs.columns]
    if not available or "author_cell_type" not in available:
        return []
    work = pd.DataFrame({field: obs[field].astype("string").fillna("<missing>") for field in available})
    for field in available:
        work[field] = work[field].replace({"": "<missing>"})
    grouped = work.groupby(available, observed=True, sort=True, dropna=False).size().reset_index(name="cells")
    records: list[dict[str, Any]] = []
    for row in grouped.to_dict(orient="records"):
        records.append({field: text_value(row.get(field)) for field in available} | {"cells": int(row["cells"])})
    return records


def exclusion_summary(obs: pd.DataFrame) -> dict[str, int]:
    n = len(obs)
    author = obs["author_cell_type"] if "author_cell_type" in obs else pd.Series([pd.NA] * n, index=obs.index)
    ontology = obs["cell_type"] if "cell_type" in obs else pd.Series([pd.NA] * n, index=obs.index)
    labels = [obs[field] if field in obs else pd.Series([pd.NA] * n, index=obs.index) for field in REQUIRED_OBS]
    doublet = doublet_or_hybrid_mask(author)
    unknown = unknown_ontology_mask(ontology)
    required_missing = pd.concat([missing_mask(series) for series in labels], axis=1).any(axis=1)
    any_exclusion = doublet | unknown | required_missing
    return {
        "cells": int(n),
        "doublet_or_hybrid_cells": int(doublet.sum()),
        "unknown_or_missing_ontology_cells": int(unknown.sum()),
        "required_metadata_missing_cells": int(required_missing.sum()),
        "any_exclusion_cells": int(any_exclusion.sum()),
        "retained_by_fixed_metadata_rules": int((~any_exclusion).sum()),
    }


def profile_dataset(
    obs: pd.DataFrame,
    dataset_key: str,
    shape: tuple[int, int] | None = None,
    raw_shape: tuple[int, int] | None = None,
    expected_observations: int | None = None,
) -> dict[str, Any]:
    """Summarize obs structure, donor metadata, annotation labels, and rules."""

    n_obs = len(obs)
    n_vars = int(shape[1]) if shape else None
    profile: dict[str, Any] = {
        "dataset_key": dataset_key,
        "observations": int(n_obs),
        "features": n_vars,
        "raw_shape": list(raw_shape) if raw_shape else None,
        "expected_observations": expected_observations,
        "observation_count_matches_manifest": expected_observations is None or n_obs == expected_observations,
        "obs_columns": list(map(str, obs.columns)),
        "missing_required_columns": [field for field in REQUIRED_OBS if field not in obs.columns],
        "duplicate_observation_ids": int(obs.index.duplicated().sum()),
        "field_value_counts": {},
        "donor_count": int(obs["donor_id"].nunique(dropna=False)) if "donor_id" in obs else 0,
        "donors": sorted(unique_values(obs["donor_id"])) if "donor_id" in obs else [],
        "annotation_catalog": annotation_catalog(obs),
        "annotation_fields_available": [field for field in ANNOTATION_FIELDS if field in obs.columns],
        "exclusion_summary": exclusion_summary(obs),
    }
    for field in (
        "donor_id",
        "disease",
        "assay",
        "assay_ontology_term_id",
        "suspension_type",
        "is_primary_data",
        "mapped_reference_annotation",
        "alignment_software",
        "suspension_dissociation_reagent",
        "author_cell_type",
        "cell_type",
        "cell_type_ontology_term_id",
    ):
        if field in obs:
            profile["field_value_counts"][field] = value_counts(obs[field])
    if "disease" in obs and "assay" in obs:
        profile["disease_by_assay"] = (
            pd.crosstab(obs["disease"].astype("string"), obs["assay"].astype("string"))
            .sort_index()
            .sort_index(axis=1)
            .astype(int)
            .to_dict()
        )
    else:
        profile["disease_by_assay"] = {}
    if "disease" in obs and "suspension_dissociation_reagent" in obs:
        profile["disease_by_dissociation_reagent"] = (
            pd.crosstab(
                obs["disease"].astype("string"),
                obs["suspension_dissociation_reagent"].astype("string"),
            )
            .sort_index()
            .sort_index(axis=1)
            .astype(int)
            .to_dict()
        )
    else:
        profile["disease_by_dissociation_reagent"] = {}
    return profile


def _group_metadata_values(group: pd.DataFrame, field: str) -> list[str]:
    return unique_values(group[field]) if field in group else []


def build_library_rows(
    obs: pd.DataFrame,
    dataset_key: str,
    dataset_id: str,
    version_id: str,
    source_file: str,
    geo_index: Mapping[str, list[dict[str, str]]],
) -> list[dict[str, Any]]:
    """Build one donor row per dataset; repeated libraries remain listed."""

    if "donor_id" not in obs:
        return []
    donor_values = normalized_series(obs["donor_id"])
    work = obs.copy()
    work["__donor_key"] = donor_values
    rows: list[dict[str, Any]] = []
    for donor, group in work.groupby("__donor_key", observed=True, sort=True, dropna=False):
        donor = text_value(donor)
        if not donor:
            donor = "<missing>"
        geo_records = list(geo_index.get(donor, []))
        donor_disease = _group_metadata_values(group, "disease")
        author_labels = value_counts(group["author_cell_type"]) if "author_cell_type" in group else {}
        ontology_labels = value_counts(group["cell_type"]) if "cell_type" in group else {}
        author = group["author_cell_type"] if "author_cell_type" in group else pd.Series([pd.NA] * len(group), index=group.index)
        ontology = group["cell_type"] if "cell_type" in group else pd.Series([pd.NA] * len(group), index=group.index)
        required = [group[field] if field in group else pd.Series([pd.NA] * len(group), index=group.index) for field in REQUIRED_OBS]
        required_missing = pd.concat([missing_mask(series) for series in required], axis=1).any(axis=1)
        primary = bool_true_mask(group["is_primary_data"]) if "is_primary_data" in group else pd.Series(False, index=group.index)
        metadata_conflicts = [
            field for field in DONOR_FIELDS if len(_group_metadata_values(group, field)) > 1
        ]
        raw_geo = json_compact(geo_records)
        raw_geo_fields = {
            field: joined_values(sorted({record[field] for record in geo_records if record.get(field)}))
            for field in GEO_RAW_FIELDS
            if field not in {"source_url"}
        }
        rows.append(
            {
                "dataset_key": dataset_key,
                "dataset_id": dataset_id,
                "dataset_version_id": version_id,
                "source_file": source_file,
                "donor_id": donor,
                "disease": joined_values(donor_disease),
                "assay": joined_values(_group_metadata_values(group, "assay")),
                "assay_ontology_term_id": joined_values(_group_metadata_values(group, "assay_ontology_term_id")),
                "suspension_type": joined_values(_group_metadata_values(group, "suspension_type")),
                "library_uuid_values": joined_values(_group_metadata_values(group, "library_uuid")),
                "sample_uuid_values": joined_values(_group_metadata_values(group, "sample_uuid")),
                "library_count": int(group["library_uuid"].nunique(dropna=True)) if "library_uuid" in group else 0,
                "sample_count": int(group["sample_uuid"].nunique(dropna=True)) if "sample_uuid" in group else 0,
                "cell_count": int(len(group)),
                "primary_data_cells": int(primary.sum()),
                "nonprimary_data_cells": int((~primary).sum()),
                "mapped_reference_annotation": joined_values(_group_metadata_values(group, "mapped_reference_annotation")),
                "alignment_software": joined_values(_group_metadata_values(group, "alignment_software")),
                "suspension_dissociation_reagent": joined_values(_group_metadata_values(group, "suspension_dissociation_reagent")),
                "author_cell_type_count": int(group["author_cell_type"].nunique(dropna=False)) if "author_cell_type" in group else 0,
                "ontology_cell_type_count": int(group["cell_type"].nunique(dropna=False)) if "cell_type" in group else 0,
                "unknown_or_missing_ontology_cells": int(unknown_ontology_mask(ontology).sum()),
                "doublet_or_hybrid_cells": int(doublet_or_hybrid_mask(author).sum()),
                "required_metadata_missing_cells": int(required_missing.sum()),
                "metadata_conflict_fields": joined_values(metadata_conflicts),
                "author_cell_type_counts_json": json_compact(author_labels),
                "cell_type_counts_json": json_compact(ontology_labels),
                "geo_donor_label_match": bool(geo_records),
                "geo_source_record_count": int(len(geo_records)),
                "geo_sample_ids": joined_values(sorted({record["sample"] for record in geo_records if record.get("sample")})),
                "geo_series_ids": joined_values(sorted({record["series"] for record in geo_records if record.get("series")})),
                "geo_assay_title_values": raw_geo_fields.get("assay_from_title", ""),
                "geo_chemistry_title_values": raw_geo_fields.get("chemistry_from_title", ""),
                "geo_suspension_title_values": raw_geo_fields.get("assay_from_title", ""),
                "geo_source_name_values": raw_geo_fields.get("source_name", ""),
                "geo_description_values": raw_geo_fields.get("description", ""),
                "geo_flags_values": raw_geo_fields.get("flags", ""),
                "geo_source_records_json": raw_geo,
                "geo_library_identity_established": False,
            }
        )
    return rows


def geo_reconciliation(table: pd.DataFrame, metadata_config: Mapping[str, Any]) -> dict[str, Any]:
    """Compare GEO title metadata with the author's donor summary without rewriting it."""

    series_records: list[dict[str, Any]] = []
    for series, group in table.groupby("series", sort=True, observed=True):
        flags = group["flags"].map(text_value)
        flags = flags[flags.ne("")]
        series_records.append(
            {
                "series": text_value(series),
                "rows": int(len(group)),
                "donor_labels": sorted(unique_values(group["donor_label_from_title"])),
                "cohort_labels": value_counts(group["cohort_from_title"]),
                "flags_rows": int(flags.size),
                "flag_values": value_counts(flags) if len(flags) else {},
                "title_description_chemistry_conflict_rows": int(group["title_description_chemistry_conflict"].sum()),
                "title_description_suspension_conflict_rows": int(group["title_description_suspension_conflict"].sum()),
            }
        )
    author_controls = sorted(
        str(value) for value in metadata_config.get("author_supplement", {}).get("control_donors", [])
    )
    geo_controls = sorted(
        unique_values(table.loc[table["cohort_from_title"].eq("Control label"), "donor_label_from_title"])
    )
    return {
        "input_file": "data/derived/geo-library-metadata.csv",
        "parent_series": "GSE243981",
        "subseries": ["GSE243977", "GSE247128"],
        "input_rows": int(len(table)),
        "input_unique_donor_labels": int(table["donor_label_from_title"].nunique()),
        "series": series_records,
        "author_supplement_control_donors": author_controls,
        "geo_control_donor_labels": geo_controls,
        "author_control_donors_missing_from_geo_title_labels": sorted(set(author_controls) - set(geo_controls)),
        "geo_control_labels_not_in_author_supplement": sorted(set(geo_controls) - set(author_controls)),
        "crosswalk_rule": "A matching donor string is a label-level link only; no exact GEO library identity is asserted.",
        "raw_values_preserved": True,
    }


def donor_overlap(obs_by_dataset: Mapping[str, pd.DataFrame]) -> dict[str, Any]:
    donor_sets: dict[str, set[str]] = {}
    cohort_sets: dict[str, dict[str, set[str]]] = {}
    for key, obs in obs_by_dataset.items():
        donors = set(unique_values(obs["donor_id"])) if "donor_id" in obs else set()
        donor_sets[key] = donors
        by_cohort: dict[str, set[str]] = {}
        if "disease" in obs:
            for donor, group in obs.groupby("donor_id", observed=True, sort=True, dropna=False):
                donor_text = text_value(donor)
                if not donor_text:
                    continue
                diseases = unique_values(group["disease"])
                for disease in diseases:
                    by_cohort.setdefault(disease, set()).add(donor_text)
        cohort_sets[key] = by_cohort
    overlap = sorted(set.intersection(*donor_sets.values())) if donor_sets else []
    shared_index: list[str] = []
    if len(obs_by_dataset) >= 2:
        series = [set(map(text_value, obs.index.tolist())) for obs in obs_by_dataset.values()]
        shared_index = sorted(set.intersection(*series))
    return {
        "donors_by_dataset": {key: sorted(values) for key, values in donor_sets.items()},
        "shared_donors": overlap,
        "shared_donor_count": len(overlap),
        "shared_donors_by_disease": {
            disease: sorted(set.intersection(*(cohort_sets[key].get(disease, set()) for key in obs_by_dataset)))
            for disease in sorted({disease for values in cohort_sets.values() for disease in values})
        },
        "shared_observation_ids": shared_index,
        "shared_observation_id_count": len(shared_index),
        "independence_note": "Shared donor labels are not independent sc/sn biological replicates; cell IDs are release-specific.",
    }


def paper_release_comparison(profiles: Iterable[Mapping[str, Any]], metadata_config: Mapping[str, Any]) -> list[dict[str, Any]]:
    paper = metadata_config.get("paper_counts", {}) or metadata_config.get("paper_source", {}).get("abstract_counts", {})
    result: list[dict[str, Any]] = []
    for profile in profiles:
        key = str(profile["dataset_key"])
        release_counts = profile.get("field_value_counts", {}).get("disease", {})
        row: dict[str, Any] = {
            "dataset_key": key,
            "release_disease_counts": release_counts,
            "paper_reference": paper.get(key, {}),
        }
        if key == "sc":
            row["scope_interpretation"] = "PSC/PBC counts agree with the paper abstract; distributed normal-cell scope is a six-donor atlas release and is smaller than the paper's 24-donor healthy comparison."
        elif key == "sn":
            row["scope_interpretation"] = "Distributed nuclei counts exceed the paper abstract counts; the release and paper subset/QC scope cannot be reconciled from labels alone."
        else:
            row["scope_interpretation"] = "Unknown dataset scope."
        result.append(row)
    return result


def source_provenance(
    metadata_config: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    source_manifest_path: Path,
    geo_path: Path,
    metadata_config_path: Path,
    verified_files: list[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "integrity_verification": {
            "status": "passed",
            "method": "Exact local byte count and SHA-256 compared with the pinned declaration before parsing any collection, H5AD obs, or GEO metadata.",
            "verified_files": verified_files,
        },
        "source_manifest": {
            "file": str(source_manifest_path.relative_to(ROOT)),
            "sha256": sha256_file(source_manifest_path),
            "declared_sources": source_manifest.get("sources", []),
        },
        "metadata_config": {
            "file": str(metadata_config_path.relative_to(ROOT)),
            "sha256": sha256_file(metadata_config_path),
            "role": "Curated provenance and comparison facts consumed by the script; not an external source byte-verification record.",
        },
        "geo_metadata": {
            "file": str(geo_path.relative_to(ROOT)),
            "sha256": sha256_file(geo_path),
            "series_urls": metadata_config.get("geo_series_urls", []),
            "snapshots": metadata_config.get("geo_snapshots", []),
        },
        "paper": metadata_config.get("paper_source", {}),
        "preprint": metadata_config.get("preprint_source", {}),
        "cellxgene_collection": metadata_config.get("cellxgene_collection", {}),
        "author_code": metadata_config.get("author_code", {}),
        "supplementary_evidence": metadata_config.get("author_supplement", {}).get("source", {}),
    }


def collection_inventory(path: Path, selected_dataset_ids: Iterable[str] = ()) -> dict[str, Any]:
    """Summarize public collection entries without exposing curator contact fields."""

    collection = json.loads(path.read_text())
    selected = {str(value) for value in selected_dataset_ids}

    def labels(items: Any) -> list[str]:
        if not isinstance(items, list):
            return []
        return sorted({text_value(item.get("label")) for item in items if isinstance(item, dict) and text_value(item.get("label"))})

    def assays(items: Any) -> list[str]:
        if not isinstance(items, list):
            return []
        return sorted({text_value(item.get("label")) for item in items if isinstance(item, dict) and text_value(item.get("label"))})

    entries: list[dict[str, Any]] = []
    for dataset in collection.get("datasets", []):
        dataset_id = text_value(dataset.get("dataset_id"))
        assets = dataset.get("assets") if isinstance(dataset.get("assets"), list) else []
        entries.append(
            {
                "dataset_id": dataset_id,
                "dataset_version_id": text_value(dataset.get("dataset_version_id")),
                "title": text_value(dataset.get("title")),
                "cell_count": dataset.get("cell_count"),
                "primary_cell_count": dataset.get("primary_cell_count"),
                "feature_count": dataset.get("feature_count"),
                "raw_data_location": text_value(dataset.get("raw_data_location")),
                "schema_version": text_value(dataset.get("schema_version")),
                "suspension_type": sorted(text_value(value) for value in dataset.get("suspension_type", []) if text_value(value)),
                "assay": assays(dataset.get("assay")),
                "disease": labels(dataset.get("disease")),
                "donor_ids": sorted(text_value(value) for value in dataset.get("donor_id", []) if text_value(value)),
                "cell_type_labels": labels(dataset.get("cell_type")),
                "asset_urls": sorted(text_value(asset.get("url")) for asset in assets if isinstance(asset, dict) and text_value(asset.get("url"))),
            }
        )
    entries.sort(key=lambda item: (item["title"], item["dataset_id"]))
    return {
        "collection_id": text_value(collection.get("collection_id")),
        "collection_version_id": text_value(collection.get("collection_version_id")),
        "name": text_value(collection.get("name")),
        "published_at": text_value(collection.get("published_at")),
        "revised_at": text_value(collection.get("revised_at")),
        "dataset_count": len(entries),
        "selected_dataset_ids": sorted(selected),
        "datasets": entries,
    }


def load_h5ad_obs(path: Path) -> tuple[pd.DataFrame, tuple[int, int], tuple[int, int] | None]:
    """Load only H5AD obs and shape metadata; never materialize X/raw.X."""

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("anndata is required to read the CELLxGENE H5AD metadata") from exc
    atlas = ad.read_h5ad(path, backed="r")
    try:
        shape = tuple(map(int, atlas.shape))
        raw_shape = tuple(map(int, atlas.raw.shape)) if atlas.raw is not None else None
        obs = atlas.obs.copy()
    finally:
        if getattr(atlas, "file", None) is not None:
            atlas.file.close()
    return obs, shape, raw_shape


def write_outputs(
    root: Path,
    metadata_config: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    source_manifest_path: Path,
    plan: Mapping[str, Any],
    geo_path: Path,
    output_csv: Path,
    output_report: Path,
) -> dict[str, Any]:
    source_by_id = {str(source.get("id")): source for source in source_manifest.get("sources", [])}
    collection_source = source_by_id.get("collection-metadata")
    if collection_source is None:
        raise ValueError("No collection metadata entry in source manifest")
    collection_path = root / "data/raw/liver-atlas" / str(collection_source["file"])
    verified_files: list[Mapping[str, Any]] = [
        verify_source_file(collection_path, collection_source, "collection-metadata")
    ]
    selected_dataset_ids = [str(dataset["dataset_id"]) for dataset in plan.get("datasets", [])]
    dataset_sources: list[tuple[Mapping[str, Any], Path]] = []
    for dataset in plan.get("datasets", []):
        source = source_by_id.get(str(dataset["dataset_id"]))
        if source is None:
            raise ValueError(f"No source manifest entry for {dataset['dataset_id']}")
        path = root / "data/raw/liver-atlas" / str(source["file"])
        dataset_sources.append((source, path))
        verified_files.append(verify_source_file(path, source, str(dataset["key"])))
    geo_declaration = {
        "id": "geo-library-metadata",
        "file": str(geo_path.relative_to(ROOT)),
        "bytes": metadata_config.get("geo_input_bytes"),
        "sha256": metadata_config.get("geo_input_sha256"),
    }
    verified_files.append(verify_source_file(geo_path, geo_declaration, "geo-library-metadata"))

    # Only after the preflight above do we parse metadata or open H5AD files.
    geo = read_geo_table(geo_path)
    geo_index = build_geo_index(geo)
    obs_by_dataset: dict[str, pd.DataFrame] = {}
    profiles: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for dataset, (source, path) in zip(plan.get("datasets", []), dataset_sources):
        key = str(dataset["key"])
        obs, shape, raw_shape = load_h5ad_obs(path)
        obs_by_dataset[key] = obs
        profiles.append(
            profile_dataset(
                obs,
                key,
                shape=shape,
                raw_shape=raw_shape,
                expected_observations=int(dataset.get("expected_observations", source.get("expected_observations"))),
            )
        )
        rows.extend(
            build_library_rows(
                obs,
                key,
                str(dataset["dataset_id"]),
                str(dataset["version_id"]),
                str(source["file"]),
                geo_index,
            )
        )

    columns = [
        "dataset_key", "dataset_id", "dataset_version_id", "source_file", "donor_id", "disease", "assay",
        "assay_ontology_term_id", "suspension_type", "library_uuid_values", "sample_uuid_values", "library_count",
        "sample_count", "cell_count", "primary_data_cells", "nonprimary_data_cells", "mapped_reference_annotation",
        "alignment_software", "suspension_dissociation_reagent", "author_cell_type_count", "ontology_cell_type_count",
        "unknown_or_missing_ontology_cells", "doublet_or_hybrid_cells", "required_metadata_missing_cells",
        "metadata_conflict_fields", "author_cell_type_counts_json", "cell_type_counts_json", "geo_donor_label_match",
        "geo_source_record_count", "geo_sample_ids", "geo_series_ids", "geo_assay_title_values",
        "geo_chemistry_title_values", "geo_suspension_title_values", "geo_source_name_values", "geo_description_values",
        "geo_flags_values", "geo_source_records_json", "geo_library_identity_established",
    ]
    table = pd.DataFrame(rows, columns=columns)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_csv, index=False, lineterminator="\n")

    report: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_scope": {
            "metadata_only": True,
            "candidate_expression_values_inspected": False,
            "matrix_values_loaded": False,
            "purpose": "Audit public cell annotations and donor/assay provenance before the bounded descriptive expression analysis.",
        },
        "inputs": source_provenance(
            metadata_config,
            source_manifest,
            source_manifest_path,
            geo_path,
            root / "config/liver-atlas-metadata.json",
            verified_files,
        ),
        "reproduction": {
            "script_reads": [
                "config/liver-atlas-sources.json: source file paths, declared H5AD hashes, retrieval metadata, and expected observation counts; collection and both H5AD files are byte/hash verified before metadata reads.",
                "config/liver-atlas-plan.json: selected SC/SN dataset IDs, versions, and expected observation counts.",
                "data/raw/liver-atlas/cellxgene-collection.json: public collection inventory, dataset assets, modalities, donor IDs, and ontology cell-type labels; curator contact fields are omitted from the report.",
                "data/raw/liver-atlas/psc-sc-ff713455.h5ad and psc-sn-34702b40.h5ad: backed obs plus shape/raw-shape metadata only; matrix values are never materialized.",
                "data/derived/geo-library-metadata.csv: exact pinned bytes/SHA-256 are verified before parsing; all original GEO columns are retained and title/description conflict counts are parsed mechanically.",
            ],
            "author_and_paper_evidence_pinned_in_config": [
                "Author repository file contents, commit metadata, local snapshot hashes, declared SC donor vector, and annotation-process descriptions are pinned in config/liver-atlas-metadata.json; the script does not re-execute the R code.",
                "The final-paper citation and abstract counts are pinned in config/liver-atlas-metadata.json; the script does not derive paper counts from H5AD metadata.",
                "The 24 control donor labels and duplicate/quality notes come from a preprint supplementary workbook pinned in config/liver-atlas-metadata.json; they are not treated as final-paper supplementary values.",
            ],
        },
        "cellxgene_collection_inventory": collection_inventory(collection_path, selected_dataset_ids),
        "dataset_profiles": profiles,
        "library_table": {
            "file": str(output_csv.relative_to(root)),
            "row_grain": "one row per CELLxGENE dataset and donor; repeated library/sample UUIDs are listed and combined for donor summaries",
            "rows": int(len(table)),
            "columns": columns,
            "geo_identity_policy": "GEO links use matching donor labels only. Exact GEO sample-to-CELLxGENE library identity is not established.",
        },
        "cross_dataset_donor_overlap": donor_overlap(obs_by_dataset),
        "paper_vs_release_counts": paper_release_comparison(profiles, metadata_config),
        "geo_reconciliation": geo_reconciliation(geo, metadata_config),
        "author_annotation_provenance": {
            "annotation_fields": {
                "author_cell_type": "CELLxGENE obs author_cell_type, retained verbatim",
                "cell_type": "CELLxGENE ontology-normalized cell type, retained verbatim",
                "cell_type_ontology_term_id": "CELLxGENE ontology term, retained verbatim",
            },
            "catalog_policy": "No cell labels are inferred from candidate expression; label/ontology pairs are reported as supplied.",
            "code_evidence": metadata_config.get("author_code", {}).get("annotation_evidence", []),
        },
        "quality_findings": [
            {
                "id": "sc-assay-condition-confounding",
                "severity": "high",
                "confidence": "high",
                "claim": "SC disease and assay are completely confounded in the distributed object.",
                "evidence": "The SC disease-by-assay table has all normal cells in 10x 5' v1 and all PSC/PBC cells in 10x 5' v2.",
                "impact": "A PSC-versus-normal SC expression contrast cannot separate disease from the captured chemistry in this release.",
                "action": "Retain SC for within-assay descriptive summaries and block cross-assay PSC-normal contrasts.",
            },
            {
                "id": "sn-release-paper-scope",
                "severity": "high",
                "confidence": "high",
                "claim": "Distributed SN counts exceed the paper abstract counts for PSC and PBC.",
                "evidence": "The release has 57,511 PSC and 21,754 PBC nuclei; the paper reports approximately 23,000 PSC and 20,202 PBC nuclei.",
                "impact": "Release-level cells cannot be silently treated as the paper's reported QC subset.",
                "action": "Use pinned release counts, retain the paper comparison, and label scope/QC reconciliation unresolved.",
            },
            {
                "id": "primary-and-processing-provenance-confounding",
                "severity": "high",
                "confidence": "high",
                "claim": "Several CELLxGENE processing fields differ systematically between normal controls and disease donors.",
                "evidence": "Normal cells/nuclei are is_primary_data=False while disease cells/nuclei are True; SN normal nuclei use NP40 or Tween-20 while cases use Tween-20, and normal alignment software is unknown while cases are Cell Ranger count v3.1.0. The CELLxGENE is_primary_data field is treated as reuse/provenance metadata, not as a quality label.",
                "impact": "Reference, pipeline, and preparation effects remain coupled to cohort labels in release metadata.",
                "action": "Keep these fields visible in donor tables and avoid attributing uncontrolled differences to PSC biology.",
            },
            {
                "id": "author-label-exclusions",
                "severity": "medium",
                "confidence": "high",
                "claim": "The public objects retain author doublet/hybrid and unknown ontology groups.",
                "evidence": "The fixed metadata rules identify 4,623 SC and 7,927 SN doublet/hybrid-labelled observations, plus 5,527 SC and 8,765 SN unknown ontology observations.",
                "impact": "Cell-type denominators depend on the frozen metadata exclusions.",
                "action": "Use original labels and apply the frozen exclusions before candidate aggregation; report excluded counts.",
            },
            {
                "id": "geo-library-identity-unresolved",
                "severity": "high",
                "confidence": "high",
                "claim": "GEO title donor labels can be linked to CELLxGENE donor_id, but exact library identity is unresolved.",
                "evidence": "The author code declares the SC donor set used for disease analysis, while GEO exposes sample titles and generic source descriptions; no shared GEO accession or CELLxGENE library UUID mapping is present in the inspected code/metadata.",
                "impact": "GEO chemistry/description conflicts must remain source-specific and cannot be used to overwrite CELLxGENE assay fields.",
                "action": "Preserve GEO raw fields, report donor-label matches separately, and do not claim GEO-identical libraries.",
            },
            {
                "id": "author-donor-name-variation",
                "severity": "medium",
                "confidence": "high",
                "claim": "The author SC code and CELLxGENE release contain a donor-name variation for the PSC014 sample.",
                "evidence": "The pinned author SC input vector uses PSC014X while CELLxGENE donor_id uses PSC014; the source strings are retained independently.",
                "impact": "A donor-level provenance comparison must distinguish source naming from an asserted identity mapping.",
                "action": "Do not rewrite either source value; treat PSC014X/PSC014 as an unresolved naming variation unless an explicit primary-source crosswalk is found.",
            },
            {
                "id": "shared-donor-sc-sn",
                "severity": "medium",
                "confidence": "high",
                "claim": "Eight donor labels occur in both SC and SN releases.",
                "evidence": "The exact overlap is reported in cross_dataset_donor_overlap; no observation IDs are shared.",
                "impact": "SC and SN should be reported as separate modalities and not counted as independent donor replication.",
                "action": "Use donor-within-dataset units and disclose the overlap in downstream summaries.",
            },
        ],
        "next_evidence_needed": [
            "An explicit public sample/library crosswalk would be required to claim GEO-identical CELLxGENE libraries.",
            "The paper's QC/subset definition would be required to reconcile distributed SN counts with abstract counts.",
        ],
    }
    output_report.parent.mkdir(parents=True, exist_ok=True)
    output_report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--geo", type=Path, default=DEFAULT_GEO)
    parser.add_argument("--metadata-config", type=Path, default=DEFAULT_METADATA_CONFIG)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output-report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    metadata_config = json.loads(args.metadata_config.resolve().read_text())
    source_manifest = json.loads(args.sources.resolve().read_text())
    plan = json.loads(args.plan.resolve().read_text())
    report = write_outputs(
        root,
        metadata_config,
        source_manifest,
        args.sources.resolve(),
        plan,
        args.geo.resolve(),
        args.output_csv.resolve(),
        args.output_report.resolve(),
    )
    print(json.dumps({"status": "ok", "library_rows": report["library_table"]["rows"], "report": str(args.output_report.resolve())}))


if __name__ == "__main__":
    main()
