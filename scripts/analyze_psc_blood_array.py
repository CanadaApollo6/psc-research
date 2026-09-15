#!/usr/bin/env python3
"""Reanalyse the complete adult GSE119600 array cohort under a pinned plan.

Use the project environment, for example::

    .venv/bin/python scripts/analyze_psc_blood_array.py \
        --plan config/psc-blood-plan.json \
        --output-dir data/derived/psc-blood/array \
        --work-dir work/psc-blood/array

The JSON may be a standalone array plan or contain ``array_analysis``. All
source paths resolve against the repository root. Each source requires path,
sha256, bytes, url, and retrieved_at_utc. Inputs are series_matrix (GEO series
matrix), platform_annotation (GEO annotation/SOFT or explicitly declared TSV),
and sample_metadata (a list of GEO records, or a dict keyed by GSE119600).
``--validate-only`` performs source, sample, feature, and numeric checks but
never estimates group effects. Full execution requires status="frozen",
allow_effects=true, a UTC frozen_at_utc, and a matching code_sha256. The optional
--plan-sha256 additionally pins the complete outer JSON file.

Method: exact single-Entrez mapping from the declared GPL column, log2 without
an offset, median across all eligible probes, then complete-family two-sided
Welch tests. HC3 is a distinct sensitivity with asymptotic-normal inference; its SE
is not the Welch SE. No expression-based selection or sample removal occurs.
All matrix cells, including excluded pediatric/unmapped-probe cells, must be
finite and strictly positive before logarithms are allowed. Positivity alone
does not establish that a source was deposited on a linear scale: the plan
must explicitly assert source_normalized_linear_intensity from source review.

Public tables are gene/probe or aggregate records only. Complete sample joins,
source sample fields, sample QC, and the adult gene-by-sample matrix stay under
ignored work/. Source intensities are never overwritten. Gene labels remain
as published. Entrez IDs are not silently converted to RNA Ensembl IDs or
matched through symbols. This is observational reanalysis of published data,
not untouched validation, proof of PSC-intrinsic specificity, or a clinical
recommendation.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
STUDY = "GSE119600"
PLATFORM = "GPL10558"
ADULT_COUNTS = {"PSC": 45, "Control": 47, "UC": 45, "PBC": 90, "CD": 48}
CHILD_COUNTS = {"UC": 48, "CD": 47}
CONDITION_MAP = {
    "primary sclerosing cholangitis": "PSC",
    "control": "Control",
    "ulcerative colitis": "UC",
    "primary biliary cholangitis": "PBC",
    "Crohn's disease": "CD",
}
CONTRASTS = (("PSC", "Control"), ("PSC", "UC"), ("PSC", "PBC"), ("PSC", "CD"))
METHODS = {
    "transform": "log2_no_shift",
    "probe_mapping": "unique_exact_entrez",
    "gene_aggregation": "median_of_log2_probes",
    "primary_test": "welch",
    "sensitivity_test": "ols_hc3_normal",
    "multiplicity": "BH_all_mapped_genes_per_contrast",
    "nonfinite_or_nonpositive": "block",
}
SOURCE_KEYS = ("series_matrix", "platform_annotation", "sample_metadata")
MISSING_ENTREZ = {"", "---", "na", "n/a", "null", "none", "nan"}
SHA_RE = re.compile(r"[0-9a-f]{64}")
ENTREZ_RE = re.compile(r"[1-9][0-9]*")
LIMITATIONS = [
    "Unadjusted established-disease whole-blood associations, not intrinsic PSC effects.",
    "Exact adult age, sex, cell proportions, treatment, and batch covariates are unavailable.",
    "The paper describes all PBC participants as female; no individual sex is imputed.",
    "All deposited adults are retained; the paper's unnamed outliers are not inferred.",
    "Distinct public sample labels do not independently certify distinct persons.",
    "The cohort and gene expression have been published and used in earlier PSC analyses.",
    "Array log2-intensity differences are not RNA-seq count-model coefficients.",
    "Array Entrez identity alone does not establish an exact RNA Ensembl-gene match.",
    "Reported BH families include every mapped gene, not only discovery hits.",
]


@dataclass
class SeriesMatrix:
    probe_ids: list[str]
    sample_ids: list[str]
    values: np.ndarray
    sample_fields: dict[str, list[list[str]]]
    series_fields: dict[str, list[list[str]]]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def open_text(path: Path):
    """Detect gzip by magic bytes and read through EOF for CRC/truncation checks."""
    with path.open("rb") as raw:
        magic = raw.read(2)
        raw.seek(0)
        if magic == b"\x1f\x8b":
            with gzip.GzipFile(fileobj=raw) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                    yield text
        else:
            with io.TextIOWrapper(raw, encoding="utf-8", newline="") as text:
                yield text


def _exact_unique_ids(values: list[str], what: str) -> None:
    if not values or any(not x or x != x.strip() for x in values):
        raise ValueError(f"Blank, padded, or missing {what}")
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate {what}")


def load_series_matrix(path: Path, study: str = STUDY) -> SeriesMatrix:
    probes: list[str] = []
    samples: list[str] = []
    rows: list[np.ndarray] = []
    sample_fields: dict[str, list[list[str]]] = defaultdict(list)
    series_fields: dict[str, list[list[str]]] = defaultdict(list)
    begun = ended = False
    with open_text(path) as stream:
        for fields in csv.reader(stream, delimiter="\t", strict=True):
            if not fields or fields == [""]:
                continue
            tag = fields[0]
            if tag == "!series_matrix_table_begin":
                if begun or ended or len(fields) != 1:
                    raise ValueError("Duplicate or malformed series matrix table start")
                begun = True
                continue
            if tag == "!series_matrix_table_end":
                if not begun or ended or not samples or len(fields) != 1:
                    raise ValueError("Unexpected series matrix table end")
                ended = True
                continue
            if ended:
                raise ValueError("Unexpected content after series matrix table end")
            if not begun:
                if tag.startswith("!Sample_"):
                    sample_fields[tag.removeprefix("!Sample_")].append(fields[1:])
                elif tag.startswith("!Series_"):
                    series_fields[tag.removeprefix("!Series_")].append(fields[1:])
                else:
                    raise ValueError(f"Unexpected series matrix metadata: {tag!r}")
                continue
            if not samples:
                if tag != "ID_REF":
                    raise ValueError("Series matrix lacks ID_REF header")
                samples = fields[1:]
                _exact_unique_ids(samples, "matrix sample identifiers")
                continue
            if len(fields) != len(samples) + 1:
                raise ValueError(f"Matrix feature row {len(probes) + 1} has incorrect width")
            probes.append(tag)
            try:
                rows.append(np.asarray(fields[1:], dtype=np.float64))
            except ValueError as exc:
                raise ValueError(f"Nonnumeric matrix feature row {len(probes)}") from exc
    if not begun or not ended or not rows:
        raise ValueError("Series matrix is missing, empty, or truncated")
    _exact_unique_ids(probes, "matrix probe identifiers")
    if series_fields.get("geo_accession") != [[study]]:
        raise ValueError("Series accession does not match the plan")
    if sample_fields.get("geo_accession") != [samples]:
        raise ValueError("Matrix headers and Sample_geo_accession order differ")
    for key, field_rows in sample_fields.items():
        if any(len(row) != len(samples) for row in field_rows):
            raise ValueError(f"Sample metadata field {key!r} has incorrect width")
    for count_row in sample_fields.get("data_row_count", []):
        if any(str(len(probes)) != count for count in count_row):
            raise ValueError("Sample_data_row_count differs from complete matrix rows")
    declared_series_samples = [x for row in series_fields.get("sample_id", [])
                               for cell in row for x in cell.split()]
    if declared_series_samples:
        _exact_unique_ids(declared_series_samples, "Series_sample_id identifiers")
        if set(declared_series_samples) != set(samples):
            raise ValueError("Series_sample_id does not match the complete matrix")
    return SeriesMatrix(probes, samples, np.stack(rows), dict(sample_fields), dict(series_fields))


def load_platform_annotation(
    path: Path, config: Mapping[str, Any], platform: str = PLATFORM,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    form = config["table_format"]
    if form not in {"geo_annotation", "geo_soft", "tsv"}:
        raise ValueError("Annotation format must be geo_annotation, geo_soft, or tsv")
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    begun = form == "tsv"
    ended = False
    table_kind = ""
    platform_ids: list[str] = []
    metadata: list[str] = []
    with open_text(path) as stream:
        for raw in stream:
            line = raw.rstrip("\r\n")
            lower = line.lower()
            if not line:
                continue
            if lower in {"!platform_table_begin", "!annotation_table_begin"}:
                if begun or ended:
                    raise ValueError("Duplicate annotation table")
                begun = True
                table_kind = lower.removesuffix("_begin")
                continue
            if lower in {"!platform_table_end", "!annotation_table_end"}:
                if not begun or ended or lower != table_kind + "_end":
                    raise ValueError("Unexpected annotation table end")
                ended = True
                continue
            if not begun:
                metadata.append(line)
                if "=" in line and line.split("=", 1)[0].strip() in {"^PLATFORM", "!Annotation_platform"}:
                    platform_ids.append(line.split("=", 1)[1].strip())
                continue
            if ended:
                raise ValueError("Unexpected content after annotation table")
            if line.startswith("#"):
                continue
            # SOFT tables are tab-separated; standard CSV quoting preserves raw labels.
            fields = next(csv.reader([line], delimiter="\t", strict=True))
            if header is None:
                header = fields
                _exact_unique_ids(header, "annotation column names")
                required = (config["probe_column"], config["entrez_column"], config["symbol_column"])
                if any(name not in header for name in required):
                    raise ValueError(f"Annotation does not contain required columns {required}")
                continue
            if len(fields) != len(header):
                raise ValueError(f"Annotation row {len(rows) + 1} has incorrect width")
            rows.append(dict(zip(header, fields)))
    if not begun or (form != "tsv" and not ended) or not rows or header is None:
        raise ValueError("Annotation table is missing, empty, or truncated")
    if any(value != platform for value in platform_ids):
        raise ValueError("Annotation platform accession differs from the plan")
    if form != "tsv" and not platform_ids:
        raise ValueError("GEO annotation does not declare a platform accession")
    _exact_unique_ids([row[config["probe_column"]] for row in rows], "annotation probe identifiers")
    return rows, {"columns": header, "row_count": len(rows),
                  "platform_accessions_as_published": platform_ids,
                  "metadata_as_published": metadata, "table_format": form}


def parse_entrez_field(raw: str) -> dict[str, Any]:
    """Keep a single unique canonical decimal ID; never use symbols or aliases.

    The explicit separator grammar is ///, semicolon, comma, or vertical bar.
    Repeated identical IDs are one unique ID. Any unresolved component alongside
    an ID makes the whole probe ineligible. No integer coercion, leading-zero
    repair, retired-ID replacement, or truncation of multi-ID fields is allowed.
    """
    tokens = [part.strip() for part in re.split(r"///|;|,|\|", raw.strip())]
    valid = sorted({x for x in tokens if ENTREZ_RE.fullmatch(x)}, key=int)
    unresolved = [x for x in tokens if x.lower() in MISSING_ENTREZ]
    invalid = [x for x in tokens if not ENTREZ_RE.fullmatch(x)
               and x.lower() not in MISSING_ENTREZ]
    if invalid:
        status = "unresolved_invalid_entrez_token"
    elif len(valid) > 1:
        status = "ambiguous_multiple_entrez"
    elif valid and unresolved:
        status = "unresolved_entrez_component"
    elif not valid:
        status = "unmapped_no_entrez"
    else:
        status = "single_unique_entrez"
    return {"entrez_gene_id_raw": raw, "entrez_gene_id": valid[0] if status == "single_unique_entrez" else "",
            "mapping_status": status, "entrez_tokens_json": json.dumps(valid),
            "unresolved_tokens_json": json.dumps(unresolved),
            "invalid_tokens_json": json.dumps(invalid)}


def build_feature_join(
    probe_ids: list[str], annotation: list[dict[str, str]], config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    _exact_unique_ids(probe_ids, "matrix probe identifiers")
    annotation_ids = [row[config["probe_column"]] for row in annotation]
    _exact_unique_ids(annotation_ids, "annotation probe identifiers")
    matrix_index = {probe: i for i, probe in enumerate(probe_ids)}
    annotation_index = {probe: i for i, probe in enumerate(annotation_ids)}
    all_ids = probe_ids + [x for x in annotation_ids if x not in matrix_index]
    feature_rows: list[dict[str, Any]] = []
    genes: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for probe in all_ids:
        ai = annotation_index.get(probe)
        mi = matrix_index.get(probe)
        source = annotation[ai] if ai is not None else {}
        parsed = parse_entrez_field(source.get(config["entrez_column"], ""))
        if ai is None:
            parsed["mapping_status"] = "missing_platform_annotation"
        join_status = ("matrix_and_platform" if ai is not None and mi is not None
                       else "matrix_only" if mi is not None else "platform_only")
        eligible = mi is not None and parsed["mapping_status"] == "single_unique_entrez"
        row = {"probe_id": probe, "matrix_row_number_1based": mi + 1 if mi is not None else "",
               "annotation_row_number_1based": ai + 1 if ai is not None else "",
               "feature_join_status": join_status, **parsed,
               "gene_symbol_as_published": source.get(config["symbol_column"], ""),
               "annotation_row_as_published_json": json.dumps(source, ensure_ascii=False, sort_keys=True),
               "eligible_for_gene_aggregation": eligible}
        feature_rows.append(row)
        if eligible:
            genes[parsed["entrez_gene_id"]].append(row)
    universe = []
    for gene in sorted(genes, key=int):
        features = genes[gene]
        universe.append({"entrez_gene_id": gene, "gene_id_namespace": "NCBI_Entrez_Gene",
                         "eligible_probe_count": len(features),
                         "probe_ids_json": json.dumps([row["probe_id"] for row in features]),
                         "gene_symbols_as_published_json": json.dumps(sorted({row["gene_symbol_as_published"]
                                                                               for row in features}), ensure_ascii=False)})
    audit = {"matrix_probe_count": len(probe_ids), "platform_probe_count": len(annotation),
             "union_probe_count": len(feature_rows),
             "join_status_counts": dict(Counter(row["feature_join_status"] for row in feature_rows)),
             "mapping_status_counts_all_platform_rows": dict(Counter(row["mapping_status"]
                            for row in feature_rows if row["annotation_row_number_1based"] != "")),
             "mapping_status_counts_measured_rows": dict(Counter(row["mapping_status"]
                            for row in feature_rows if row["matrix_row_number_1based"] != "")),
             "eligible_measured_probes": sum(row["eligible_for_gene_aggregation"] for row in feature_rows),
             "mapped_gene_family_size": len(universe),
             "probe_selection_uses_expression_or_effects": False,
             "ensembl_crosswalk_applied": False}
    return feature_rows, universe, audit


def load_sample_metadata(path: Path, study: str = STUDY) -> Any:
    """Read independent GEO records, skipping any quick-SOFT expression table."""
    with open_text(path) as stream:
        first = stream.read(1)
    if first in {"{", "["}:
        with open_text(path) as stream:
            return json.load(stream)
    records: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    in_table = False
    with open_text(path) as stream:
        for raw in stream:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if in_table:
                    raise ValueError("Truncated sample table before the next GEO record")
                current = {"sample_accession": line.split(" = ", 1)[1]}
                records.append(current)
            elif line == "!sample_table_begin":
                if in_table or current is None:
                    raise ValueError("Unexpected sample table start")
                in_table = True
            elif line == "!sample_table_end":
                if not in_table:
                    raise ValueError("Unexpected sample table end")
                in_table = False
            elif not in_table and line.startswith("!Sample_"):
                if current is None or " = " not in line:
                    raise ValueError("GEO sample metadata is malformed")
                tag, value = line.split(" = ", 1)
                current.setdefault(tag.removeprefix("!Sample_"), []).append(value)
            elif line and not in_table and not line.startswith("#"):
                raise ValueError("Unexpected content in independent sample metadata")
    if in_table or not records:
        raise ValueError("Independent sample metadata is empty or truncated")
    for record in records:
        if record.get("geo_accession") != [record["sample_accession"]] or record.get("series_id") != [study]:
            raise ValueError("Independent GEO sample/series accessions do not match")
    return records


def _characteristics(values: list[str]) -> dict[str, str]:
    parsed = {}
    for value in values:
        if ":" not in value:
            raise ValueError("Sample characteristic lacks an explicit key:value pair")
        key, text = value.split(":", 1)
        if key in parsed:
            raise ValueError(f"Duplicated sample characteristic key {key!r}")
        parsed[key] = text.strip()
    return parsed


def join_samples(
    matrix: SeriesMatrix, metadata: Any, expected_total: int,
    expected_adult: Mapping[str, int], expected_child: Mapping[str, int],
    study: str = STUDY, platform: str = PLATFORM,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    source_records = metadata[study] if isinstance(metadata, dict) else metadata
    if not isinstance(source_records, list):
        raise ValueError("Sample metadata must be a list or a dict keyed by study")
    source_ids = [row["sample_accession"] for row in source_records]
    _exact_unique_ids(source_ids, "independent metadata sample identifiers")
    _exact_unique_ids(matrix.sample_ids, "matrix sample identifiers")
    if len(matrix.sample_ids) != expected_total or set(source_ids) != set(matrix.sample_ids):
        raise ValueError("Full matrix and independent metadata sample sets/counts differ")
    by_sample = {row["sample_accession"]: (i, row) for i, row in enumerate(source_records)}
    fields_to_join = ("title", "source_name_ch1", "characteristics_ch1", "platform_id", "data_processing")
    rows = []
    adult_counts: Counter[str] = Counter()
    child_counts: Counter[str] = Counter()
    titles = []
    for index, sample in enumerate(matrix.sample_ids):
        source_index, record = by_sample[sample]
        matrix_record = {key: [line[index] for line in matrix.sample_fields[key]]
                         for key in matrix.sample_fields}
        for key in fields_to_join:
            expected = record.get(key)
            observed = matrix_record.get(key)
            if not isinstance(expected, list) or not observed or Counter(expected) != Counter(observed):
                raise ValueError(f"Exact metadata join failed for {sample}, {key}")
        if len(record["title"]) != 1 or record["platform_id"] != [platform]:
            raise ValueError("Sample title or platform metadata is not unique/exact")
        title = record["title"][0]
        titles.append(title)
        chars = _characteristics(record["characteristics_ch1"])
        if chars.get("tissue") != "whole blood":
            raise ValueError("Unexpected sample tissue")
        if chars.get("condition") not in CONDITION_MAP:
            raise ValueError("Unknown source condition; no guessed group assignment is allowed")
        group = CONDITION_MAP[chars["condition"]]
        age_group = chars.get("age group")
        if age_group not in {"adult", "child"}:
            raise ValueError("Unknown age group; no guessed adult selection is allowed")
        selected = age_group == "adult"
        (adult_counts if selected else child_counts)[group] += 1
        rows.append({"sample_accession": sample, "matrix_column_number_1based": index + 2,
                     "metadata_row_number_1based": source_index + 1,
                     "title_as_published": title, "source_name_as_published": record["source_name_ch1"][0],
                     "condition_as_published": chars["condition"], "age_group_as_published": age_group,
                     "analysis_group": group, "selected_adult": selected,
                     "selection_status": "retained_all_deposited_adults" if selected else "excluded_pediatric_by_metadata",
                     "sample_join_status": "exact_accession_and_metadata", "platform_id": platform,
                     "characteristics_as_published_json": json.dumps(record["characteristics_ch1"], ensure_ascii=False),
                     "source_metadata_as_published_json": json.dumps(record, ensure_ascii=False, sort_keys=True),
                     "matrix_metadata_as_published_json": json.dumps(matrix_record, ensure_ascii=False, sort_keys=True),
                     "exact_age_available": "age" in chars, "sex_available": any(key.lower() == "sex" for key in chars)})
    # Exact duplicate titles are a source-identity failure, not an expression-QC exclusion.
    if any(not title for title in titles) or len(titles) != len(set(titles)):
        raise ValueError("Blank or duplicate exact sample titles")
    if dict(adult_counts) != dict(expected_adult) or dict(child_counts) != dict(expected_child):
        raise ValueError(f"Observed adult/child groups differ from the plan: {dict(adult_counts)}, {dict(child_counts)}")
    return rows, {"total_deposited_samples": len(rows), "adult_samples_retained": sum(adult_counts.values()),
                  "pediatric_samples_excluded": sum(child_counts.values()),
                  "adult_group_counts": dict(adult_counts), "pediatric_group_counts": dict(child_counts),
                  "exact_sample_accession_join": True, "exact_metadata_field_join": True,
                  "unique_matrix_sample_accessions": len(set(matrix.sample_ids)),
                  "unique_exact_sample_titles": len(set(titles)),
                  "adults_with_exact_age": sum(row["selected_adult"] and row["exact_age_available"] for row in rows),
                  "adults_with_sex": sum(row["selected_adult"] and row["sex_available"] for row in rows),
                  "numeric_covariates_imputed": False, "effect_based_sample_exclusions": 0,
                  "person_level_identity_independently_confirmed": False}


def linear_value_audit(values: np.ndarray) -> dict[str, Any]:
    if values.ndim != 2 or not values.size:
        raise ValueError("Source matrix must be nonempty and two-dimensional")
    finite = np.isfinite(values)
    valid_values = values[finite]
    zero = int(np.count_nonzero(values == 0))
    negative = int(np.count_nonzero(values < 0))
    nonfinite = int(values.size - np.count_nonzero(finite))
    return {"probe_rows": values.shape[0], "sample_columns": values.shape[1],
            "source_cell_count": int(values.size), "nonfinite_cells": nonfinite,
            "zero_cells": zero, "negative_cells": negative,
            "finite_min": float(np.min(valid_values)) if valid_values.size else None,
            "finite_max": float(np.max(valid_values)) if valid_values.size else None,
            "finite_quantiles_01_50_99": np.quantile(valid_values, [0.01, 0.5, 0.99]).tolist() if valid_values.size else [],
            "log2_without_shift_valid": nonfinite == 0 and zero == 0 and negative == 0,
            "log_validity_scope": "all_cells_all_deposited_samples_all_probes",
            "positivity_proves_linear_source_scale": False}


def log2_no_shift(values: np.ndarray) -> np.ndarray:
    if values.ndim != 2 or not values.size or not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("log2 requires finite strictly positive source intensities; offsets are forbidden")
    result = np.log2(values)
    if not np.all(np.isfinite(result)):
        raise ValueError("Nonfinite log2 output")
    return result


def aggregate_probes(
    log2_values: np.ndarray, probe_ids: list[str], feature_rows: list[dict[str, Any]],
) -> tuple[list[str], np.ndarray]:
    if log2_values.ndim != 2 or log2_values.shape[0] != len(probe_ids) or not np.all(np.isfinite(log2_values)):
        raise ValueError("Invalid transformed probe matrix")
    _exact_unique_ids(probe_ids, "aggregation probe identifiers")
    by_probe = {probe: i for i, probe in enumerate(probe_ids)}
    seen: set[str] = set()
    groups: dict[str, list[int]] = defaultdict(list)
    for row in feature_rows:
        if row["probe_id"] in seen:
            raise ValueError("Duplicate feature join row")
        seen.add(row["probe_id"])
        if row["eligible_for_gene_aggregation"]:
            if row["probe_id"] not in by_probe or not ENTREZ_RE.fullmatch(row["entrez_gene_id"]):
                raise ValueError("Eligible feature lacks exact matrix/Entrez identity")
            groups[row["entrez_gene_id"]].append(by_probe[row["probe_id"]])
    if not set(probe_ids).issubset(seen):
        raise ValueError("Incomplete matrix-feature join at aggregation")
    genes = sorted(groups, key=int)
    if not genes:
        raise ValueError("No mapped genes remain; no tests are possible")
    return genes, np.stack([np.median(log2_values[groups[gene], :], axis=0) for gene in genes])


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """BH over the entire supplied family; unavailable tests retain missing q.

    Nonestimable tests occupy the family as P=1 for the correction arithmetic
    only. Their source P and reported adjusted P remain NaN, never a null result.
    """
    p = np.asarray(pvalues, dtype=np.float64)
    if p.ndim != 1 or not p.size or np.any(np.isinf(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("BH requires a nonempty one-dimensional valid P family")
    missing = np.isnan(p)
    arithmetic_p = np.where(missing, 1.0, p)
    order = np.argsort(arithmetic_p, kind="stable")
    ordered_q = arithmetic_p[order] * p.size / np.arange(1, p.size + 1)
    ordered_q = np.minimum(1.0, np.minimum.accumulate(ordered_q[::-1])[::-1])
    q = np.empty(p.size)
    q[order] = ordered_q
    q[missing] = np.nan
    return q


def pairwise_statistics(case: np.ndarray, reference: np.ndarray) -> dict[str, np.ndarray]:
    """Gene-wise case-minus-reference Welch primary and two-group OLS-HC3.

    Both arrays have shape genes x samples. HC3 variance is s_case^2/(n_case-1)
    + s_reference^2/(n_reference-1), not the Welch variance. HC3 P/intervals use
    an asymptotic standard-normal reference, not an exact t distribution.
    No covariates are fitted.
    """
    case = np.asarray(case, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if (case.ndim != 2 or reference.ndim != 2 or case.shape[0] != reference.shape[0]
            or not case.shape[0] or min(case.shape[1], reference.shape[1]) < 2
            or not np.all(np.isfinite(case)) or not np.all(np.isfinite(reference))):
        raise ValueError("Pairwise tests need finite matched genes and at least two samples per group")
    n1, n0 = case.shape[1], reference.shape[1]
    difference = np.mean(case, axis=1) - np.mean(reference, axis=1)
    s1, s0 = np.var(case, axis=1, ddof=1), np.var(reference, axis=1, ddof=1)
    terms1, terms0 = s1 / n1, s0 / n0
    welch_var = terms1 + terms0
    hc3_var = s1 / (n1 - 1) + s0 / (n0 - 1)
    valid = (welch_var > 0) & np.isfinite(welch_var) & np.isfinite(difference)
    welch_df = np.full(case.shape[0], np.nan)
    welch_df[valid] = welch_var[valid] ** 2 / (terms1[valid] ** 2 / (n1 - 1)
                                                          + terms0[valid] ** 2 / (n0 - 1))
    result: dict[str, np.ndarray] = {"log2_intensity_difference": difference}
    for prefix, variance, dfs in (("welch", welch_var, welch_df),
                                 ("hc3", hc3_var, np.full(case.shape[0], n1 + n0 - 2.0))):
        se = np.sqrt(variance)
        statistic = np.full(case.shape[0], np.nan)
        p = np.full(case.shape[0], np.nan)
        lower = np.full(case.shape[0], np.nan)
        upper = np.full(case.shape[0], np.nan)
        good = valid & np.isfinite(variance) & np.isfinite(dfs) & (dfs > 0)
        statistic[good] = difference[good] / se[good]
        if prefix == "welch":
            p[good] = 2 * stats.t.sf(np.abs(statistic[good]), dfs[good])
            critical = stats.t.ppf(0.975, dfs[good])
        else:
            p[good] = 2 * stats.norm.sf(np.abs(statistic[good]))
            critical = stats.norm.ppf(0.975)
        lower[good] = difference[good] - critical * se[good]
        upper[good] = difference[good] + critical * se[good]
        result.update({prefix + "_status": np.where(good, "estimable", "not_estimable_zero_or_invalid_variance"),
                       prefix + "_se": se, prefix + "_statistic": statistic,
                       prefix + "_p": p, prefix + "_bh_q": bh_adjust(p),
                       prefix + "_ci95_low": lower, prefix + "_ci95_high": upper})
        if prefix == "welch":
            result["welch_t"] = statistic
            result["welch_df"] = np.where(good, dfs, np.nan)
        else:
            result["hc3_z"] = statistic
    return result


def analyse_contrasts(
    gene_ids: list[str], values: np.ndarray, groups: list[str], universe: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if values.shape != (len(gene_ids), len(groups)) or gene_ids != [row["entrez_gene_id"] for row in universe]:
        raise ValueError("Gene/sample axes and complete family are not exactly joined")
    group_array = np.asarray(groups)
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for case, reference in CONTRASTS:
        case_mask = group_array == case
        reference_mask = group_array == reference
        result = pairwise_statistics(values[:, case_mask], values[:, reference_mask])
        contrast = case + "_vs_" + reference
        for index, gene in enumerate(gene_ids):
            rows.append({"study_accession": STUDY, "entrez_id": gene, "entrez_gene_id": gene,
                         "estimate": result["log2_intensity_difference"][index],
                         "se": result["welch_se"][index], "statistic": result["welch_t"][index],
                         "pvalue": result["welch_p"][index], "padj": result["welch_bh_q"][index],
                         "primary_test": "welch", "hc3_reference_distribution": "asymptotic_standard_normal",
                         "gene_id_namespace": "NCBI_Entrez_Gene",
                         "gene_symbols_as_published_json": universe[index]["gene_symbols_as_published_json"],
                         "eligible_probe_count": universe[index]["eligible_probe_count"],
                         "contrast": contrast, "case_group": case, "reference_group": reference,
                         "n_case": int(np.count_nonzero(case_mask)), "n_reference": int(np.count_nonzero(reference_mask)),
                         **{key: value[index] for key, value in result.items()},
                         "multiple_testing_family_size": len(gene_ids),
                         "comparison_role": "primary_replication" if reference == "Control" else "disease_control_comparison",
                         "confounding_status": "unadjusted_whole_blood_established_disease"})
        summaries.append({"contrast": contrast, "mapped_gene_family_size": len(gene_ids),
                          "n_case": int(np.count_nonzero(case_mask)), "n_reference": int(np.count_nonzero(reference_mask)),
                          "welch_estimable": int(np.count_nonzero(np.isfinite(result["welch_p"]))),
                          "welch_bh_q_le_0_05": int(np.count_nonzero(result["welch_bh_q"] <= 0.05)),
                          "hc3_estimable": int(np.count_nonzero(np.isfinite(result["hc3_p"]))),
                          "hc3_bh_q_le_0_05": int(np.count_nonzero(result["hc3_bh_q"] <= 0.05))})
    return rows, summaries


def _json_dump(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def _csv_cell(value: Any) -> Any:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g") if np.isfinite(value) else ""
    if isinstance(value, (bool, np.bool_)):
        return "true" if value else "false"
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty schema-less CSV: {path.name}")
    def write_to(stream):
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_cell(value) for key, value in row.items()})
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                    write_to(stream)
    else:
        with path.open("w", encoding="utf-8", newline="") as stream:
            write_to(stream)


def _verified_source(record: Mapping[str, Any], root: Path) -> tuple[Path, dict[str, Any]]:
    for key in ("path", "sha256", "bytes", "url", "retrieved_at_utc"):
        if key not in record:
            raise ValueError(f"Source is missing pin/provenance field {key}")
    expected_hash = record["sha256"]
    if not isinstance(expected_hash, str) or SHA_RE.fullmatch(expected_hash) is None:
        raise ValueError("Source SHA-256 must be 64 lowercase hexadecimal characters")
    if type(record["bytes"]) is not int or record["bytes"] <= 0:
        raise ValueError("Source byte count must be a positive integer")
    if not record["url"] or not record["retrieved_at_utc"]:
        raise ValueError("Source URL and retrieval date must be nonempty")
    path = (root / record["path"]).resolve()
    size, digest = path.stat().st_size, sha256_file(path)
    if digest != expected_hash or size != record["bytes"]:
        raise ValueError(f"Source SHA-256/byte pin mismatch: {record['path']}")
    return path, dict(record)


def validate_plan(plan: Mapping[str, Any], *, validate_only: bool, code_path: Path) -> None:
    if (plan.get("schema_version") != 1 or plan.get("study_accession") != STUDY
            or plan.get("platform_accession") != PLATFORM):
        raise ValueError("Unsupported array plan schema/study/platform")
    if plan.get("source_quantity") != "source_normalized_linear_intensity":
        raise ValueError("The frozen plan must explicitly establish the deposited linear intensity quantity")
    if plan.get("methods") != METHODS:
        raise ValueError("Array plan methods do not match this implementation")
    if "contrasts" in plan and plan["contrasts"] != [case + "_vs_" + reference for case, reference in CONTRASTS]:
        raise ValueError("Array contrast family differs from the prespecified four comparisons")
    if plan.get("expected_total_samples") != 370 or plan.get("expected_adult_counts") != ADULT_COUNTS or plan.get("expected_child_counts") != CHILD_COUNTS:
        raise ValueError("Array plan must retain the full prespecified 275-adult/95-child source cohort")
    if set(plan.get("sources", {})) != set(SOURCE_KEYS):
        raise ValueError("Array plan must pin all three source inputs, with no extra unexecuted source keys")
    annotation = plan.get("annotation", {})
    if set(annotation) != {"table_format", "probe_column", "entrez_column", "symbol_column"}:
        raise ValueError("Array plan must declare exact annotation columns and table format")
    if "code_sha256" in plan and plan["code_sha256"] != sha256_file(code_path):
        raise ValueError("Array implementation SHA-256 differs from the plan")
    if not validate_only:
        if plan.get("status") != "frozen" or plan.get("allow_effects") is not True or "code_sha256" not in plan:
            raise ValueError("Effects are blocked: need a frozen, explicitly authorised, code-pinned plan")
        try:
            freeze = datetime.fromisoformat(plan["frozen_at_utc"].replace("Z", "+00:00"))
        except (KeyError, ValueError, TypeError, AttributeError) as exc:
            raise ValueError("Frozen plan lacks a valid UTC freeze timestamp") from exc
        if freeze.tzinfo is None or freeze.utcoffset() != timezone.utc.utcoffset(freeze):
            raise ValueError("Freeze timestamp must explicitly be UTC")


def _prepare_destinations(output_dir: Path, work_dir: Path, root: Path) -> None:
    output = output_dir.resolve()
    work = work_dir.resolve()
    if not work.is_relative_to((root / "work").resolve()) or work == (root / "work").resolve():
        raise ValueError("Patient/sample intermediates require a dedicated ignored work/ subdirectory")
    if output == work or output.is_relative_to(work) or work.is_relative_to(output):
        raise ValueError("Public output and patient work directories must be disjoint")
    for path in (output, work):
        if path.exists() and any(path.iterdir()):
            raise ValueError(f"Refusing to mix or overwrite an existing run directory: {path}")
    output.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)


def _sample_numeric_qc(matrix: SeriesMatrix) -> tuple[list[dict[str, Any]], list[list[str]]]:
    rows = []
    signatures: dict[str, list[str]] = defaultdict(list)
    for index, sample in enumerate(matrix.sample_ids):
        values = matrix.values[:, index]
        finite_values = values[np.isfinite(values)]
        quantiles = np.quantile(finite_values, [0.01, 0.5, 0.99]) if finite_values.size else [None] * 3
        signature = hashlib.sha256(np.ascontiguousarray(values, dtype="<f8").tobytes()).hexdigest()
        signatures[signature].append(sample)
        rows.append({"sample_accession": sample, "matrix_column_number_1based": index + 2,
                     "finite_min": float(np.min(finite_values)) if finite_values.size else None,
                     "finite_max": float(np.max(finite_values)) if finite_values.size else None,
                     "quantile_01": quantiles[0], "median": quantiles[1], "quantile_99": quantiles[2],
                     "nonfinite_count": int(np.count_nonzero(~np.isfinite(values))),
                     "zero_count": int(np.count_nonzero(values == 0)),
                     "negative_count": int(np.count_nonzero(values < 0)),
                     "exact_numeric_column_sha256": signature,
                     "qc_causes_sample_exclusion": False})
    duplicates = [group for group in signatures.values() if len(group) > 1]
    return rows, duplicates


def run_analysis(
    plan_path: Path, output_dir: Path, work_dir: Path, *, validate_only: bool = False,
    expected_plan_sha256: str | None = None, root: Path = ROOT,
) -> dict[str, Any]:
    plan_hash = sha256_file(plan_path)
    if expected_plan_sha256 is not None and plan_hash != expected_plan_sha256:
        raise ValueError("Outer plan SHA-256 mismatch")
    outer = json.loads(plan_path.read_text(encoding="utf-8"))
    plan = outer.get("array_analysis", outer)
    code_path = Path(__file__).resolve()
    validate_plan(plan, validate_only=validate_only, code_path=code_path)
    sources, provenance = {}, {}
    for name in SOURCE_KEYS:
        sources[name], provenance[name] = _verified_source(plan["sources"][name], root)
    matrix = load_series_matrix(sources["series_matrix"])
    annotation, annotation_audit = load_platform_annotation(sources["platform_annotation"], plan["annotation"])
    independent_metadata = load_sample_metadata(sources["sample_metadata"])
    sample_rows, sample_audit = join_samples(matrix, independent_metadata, plan["expected_total_samples"],
                                            plan["expected_adult_counts"], plan["expected_child_counts"])
    if "expected_probe_count" in plan and len(matrix.probe_ids) != plan["expected_probe_count"]:
        raise ValueError("Matrix probe count differs from the frozen plan")
    if "expected_annotation_rows" in plan and len(annotation) != plan["expected_annotation_rows"]:
        raise ValueError("Platform annotation row count differs from the frozen plan")
    feature_rows, universe, feature_audit = build_feature_join(matrix.probe_ids, annotation, plan["annotation"])
    if not universe:
        raise ValueError("No eligible mapped gene family")
    numeric_audit = linear_value_audit(matrix.values)
    sample_qc, duplicate_vectors = _sample_numeric_qc(matrix)
    _prepare_destinations(output_dir, work_dir, root)
    # Per-sample labels, QC, and complete metadata are deliberately NOT public outputs.
    write_csv(work_dir / "array-sample-join.csv", sample_rows)
    write_csv(work_dir / "array-sample-numeric-qc.csv", sample_qc)
    _json_dump(work_dir / "array-duplicate-sample-vectors.json", duplicate_vectors)
    write_csv(output_dir / "array-feature-join.csv.gz", feature_rows)
    write_csv(output_dir / "array-gene-universe.csv.gz", universe)
    audit = {"study_accession": STUDY, "platform_accession": PLATFORM,
             "outer_plan_sha256": plan_hash, "code_sha256": sha256_file(code_path),
             "plan_frozen_at_utc": plan.get("frozen_at_utc"), "sources": provenance,
             "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__},
             "source_quantity": plan["source_quantity"], "methods": METHODS,
             "sample_audit": sample_audit, "annotation_audit": annotation_audit,
             "feature_audit": feature_audit, "numeric_audit": numeric_audit,
             "exact_duplicate_numeric_sample_vector_groups": len(duplicate_vectors),
             "duplicate_numeric_vectors_cause_exclusion": False,
             "limitations": LIMITATIONS,
             "complete_sample_joins_and_expression_location": "ignored work/ only",
             "group_effects_computed": False, "stage": "validated_without_effects"}
    if not numeric_audit["log2_without_shift_valid"]:
        audit["stage"] = "blocked_nonpositive_or_nonfinite_source_intensity"
    elif not validate_only:
        selected = [i for i, row in enumerate(sample_rows) if row["selected_adult"]]
        adult_rows = [sample_rows[i] for i in selected]
        gene_ids, gene_values = aggregate_probes(log2_no_shift(matrix.values[:, selected]), matrix.probe_ids, feature_rows)
        if gene_ids != [row["entrez_gene_id"] for row in universe]:
            raise ValueError("Aggregation changed the complete mapped gene family")
        np.save(work_dir / "array-adult-entrez-log2.npy", gene_values, allow_pickle=False)
        _json_dump(work_dir / "array-adult-entrez-log2-axes.json",
                   {"gene_ids": gene_ids, "gene_id_namespace": "NCBI_Entrez_Gene",
                    "sample_accessions": [row["sample_accession"] for row in adult_rows],
                    "analysis_groups": [row["analysis_group"] for row in adult_rows],
                    "matrix_shape": list(gene_values.shape), "quantity": "median_of_log2_source_normalized_linear_probes"})
        effects, summaries = analyse_contrasts(gene_ids, gene_values, [row["analysis_group"] for row in adult_rows], universe)
        write_csv(output_dir / "array-gene-effects.csv.gz", effects)
        _json_dump(output_dir / "array-contrast-summary.json", summaries)
        audit["group_effects_computed"] = True
        audit["stage"] = "completed_all_adult_contrasts"
        audit["effect_rows"] = len(effects)
    audit["public_artifacts"] = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                                  for path in sorted(output_dir.iterdir()) if path.is_file()]
    audit["ignored_work_artifacts"] = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                                        for path in sorted(work_dir.iterdir()) if path.is_file()]
    _json_dump(output_dir / "array-analysis-audit.json", audit)
    return audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--plan-sha256")
    parser.add_argument("--validate-only", action="store_true", help="No group effects or patient expression export")
    args = parser.parse_args(argv)
    try:
        audit = run_analysis(args.plan, args.output_dir, args.work_dir, validate_only=args.validate_only,
                             expected_plan_sha256=args.plan_sha256)
    except (ValueError, OSError, KeyError, TypeError, csv.Error) as exc:
        print(f"Array analysis blocked: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({"stage": audit["stage"], "group_effects_computed": audit["group_effects_computed"],
                      "mapped_gene_family_size": audit["feature_audit"]["mapped_gene_family_size"]}, sort_keys=True))
    return 0 if not audit["stage"].startswith("blocked_") else 2


if __name__ == "__main__":
    raise SystemExit(main())
