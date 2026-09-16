#!/usr/bin/env python3
"""Bounded liver input qualification only. No target/program/effect analysis.

Use the repository .venv/bin/python. Public summaries contain aggregate QC;
expression and sample joins remain in ignored work/psc-liver-qualification/.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/psc-liver-qualification"
WORK = ROOT / "work/psc-liver-qualification"
OUTPUT = ROOT / "data/derived/psc-liver-input-qualification.json"
REPORT = ROOT / "reports/psc-liver-input-qualification.md"
NEW_TOTAL_LIMIT = 80_000_000
MAIN_NEW_LIMIT = 58_000_000  # Other allocations: Methods 12 MB; atlas 10 MB.


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    partial.replace(path)



def copy_bounded_body(response, path, limit):
    """Never request a body read larger than the remaining authorized budget."""
    size, digest, error = 0, hashlib.sha256(), None
    with Path(path).open("wb") as stream:
        try:
            while size < limit:
                block = response.raw.read(min(65_536, limit-size), decode_content=False)
                if not block:
                    break
                if len(block) > limit-size:
                    raise ValueError("HTTP reader returned more bytes than requested")
                stream.write(block)
                digest.update(block)
                size += len(block)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
    return size, digest.hexdigest(), error


def acquire(raw=RAW):
    import requests

    raw = Path(raw)
    raw.mkdir(parents=True, exist_ok=True)
    manifest_path = raw / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "schema_version": 1, "scope": "processed/metadata qualification only", "sources": [],
        "main_new_bytes_limit": MAIN_NEW_LIMIT, "overall_new_bytes_limit": NEW_TOTAL_LIMIT,
        "other_allocations": {"methods": 12_000_000, "atlas": 10_000_000}, "raw_reads_acquired": False}
    discovery = ROOT / "data/raw/public-data-discovery-options"
    original_receipts_path = discovery / "source-receipts.json"
    original_receipts = json.loads(original_receipts_path.read_text())
    reused = ["GSE303271-series.txt", "GSE303271-family.soft.gz", "GSE159676-series.txt",
              "GSE159676-samples-quick.txt", "GSE159676-matrix-listing.html"]
    for name in reused:
        if any(row["file"] == name and row.get("usable") for row in manifest["sources"]):
            continue
        receipts = [row for row in original_receipts if row.get("file") == name and row.get("status") == 200]
        if len(receipts) != 1:
            raise ValueError("Expected one original cache receipt")
        receipt = receipts[0]
        source = discovery / name
        if source.stat().st_size != receipt["bytes"] or sha256(source) != receipt["sha256"]:
            raise ValueError("Reused source does not match original receipt")
        shutil.copyfile(source, raw / name)
        manifest["sources"].append({"file": name, "url": receipt["url"], "bytes": receipt["bytes"],
                                     "sha256": receipt["sha256"], "retrieved_at_utc": receipt["retrieved_at_utc"],
                                     "copied_at_utc": now(), "acquisition": "verified_cache_copy", "new_bytes": 0,
                                     "copied_from": str(source.relative_to(ROOT)), "usable": True,
                                     "original_receipts_sha256": sha256(original_receipts_path), "original_receipt": receipt})
        dump_json(manifest_path, manifest)
    urls = [
        ("GSE303271_raw_counts.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303271/suppl/GSE303271_raw_counts.txt.gz", 3_021_812),
        ("GSE159676_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159676/matrix/GSE159676_series_matrix.txt.gz", 2_252_130),
        ("GPL6244.annot.gz", "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz", None),
        ("GPL6244-platform-quick.txt", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244&targ=self&form=text&view=quick", None)]
    for name, url, expected_bytes in urls:
        matched = [row for row in manifest["sources"] if row["file"] == name and row.get("usable")]
        if matched:
            if len(matched) != 1 or sha256(raw/name) != matched[0]["sha256"] or (raw/name).stat().st_size != matched[0]["bytes"]:
                raise ValueError("Existing acquired source failed receipt verification")
            continue
        remaining = MAIN_NEW_LIMIT - sum(row["new_bytes"] for row in manifest["sources"])
        started = now()
        with requests.get(url, stream=True, timeout=(30, 120), headers={"Accept-Encoding": "identity"}) as response:
            receipt = {"file": name, "url": url, "response_url": response.url, "status": response.status_code,
                       "started_at_utc": started, "headers": {key: value for key, value in response.headers.items() if key.lower() != "set-cookie"},
                       "acquisition": "complete_GET", "usable": False}
            length = response.headers.get("Content-Length")
            if length is not None and int(length) > remaining:
                receipt.update({"bytes": 0, "new_bytes": 0, "retrieved_at_utc": now(), "failure": "Content-Length exceeds remaining authorized allocation"})
                manifest["sources"].append(receipt)
                dump_json(manifest_path, manifest)
                raise ValueError(receipt["failure"])
            partial = raw / (name + ".partial")
            size, digest, error = copy_bounded_body(response, partial, remaining)
            receipt.update({"bytes": size, "new_bytes": size, "sha256": digest, "retrieved_at_utc": now(),
                            "read_error": error, "body_read_strategy": "bounded raw entity reads; Accept-Encoding identity"})
            usable = (error is None and response.status_code == 200 and "Content-Range" not in response.headers
                      and response.headers.get("Content-Encoding", "identity").lower() == "identity"
                      and (length is None or int(length) == size)
                      and (size < remaining or (length is not None and int(length) == size))
                      and (expected_bytes is None or size == expected_bytes))
            if name.endswith(".gz"):
                with partial.open("rb") as stream:
                    usable = usable and stream.read(2) == b"\x1f\x8b"
            if usable:
                partial.replace(raw / name)
                receipt["usable"] = True
            else:
                failed_name = name + f".attempt-{len(manifest['sources'])+1}.failed-body"
                partial.replace(raw / failed_name)
                receipt["file"] = failed_name
                receipt["failure"] = "Noncomplete/non-gzip response or byte-size mismatch; not a qualified input"
            manifest["sources"].append(receipt)
            manifest["new_download_bytes"] = sum(row["new_bytes"] for row in manifest["sources"])
            dump_json(manifest_path, manifest)
            print(json.dumps({key: receipt[key] for key in ["file", "status", "bytes", "sha256", "usable"]}), flush=True)
    return manifest



def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def geo_records(path, include_tables=False):
    """Read SOFT metadata literally; never infer donor identity or normalize labels."""
    records, current, in_table, table_header = [], None, False, None
    with open_text(path) as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^"):
                entity, accession = line[1:].split(" = ", 1)
                current = {"entity": entity, "accession": accession, "fields": {}, "column_descriptions": {}, "table": []}
                records.append(current)
                in_table = False
            elif current is not None and re.fullmatch(r"!\w+_table_begin", line):
                in_table, table_header = True, None
            elif current is not None and re.fullmatch(r"!\w+_table_end", line):
                in_table = False
            elif in_table:
                if include_tables:
                    row = next(csv.reader([line], delimiter="\t"))
                    if table_header is None:
                        table_header = row
                    elif row:
                        if len(row) != len(table_header):
                            raise ValueError("Ragged SOFT table")
                        current["table"].append(dict(zip(table_header, row)))
            elif current is not None and " = " in line:
                key, value = line.split(" = ", 1)
                if key.startswith("!"):
                    current["fields"].setdefault(key[1:], []).append(value)
                elif key.startswith("#"):
                    current["column_descriptions"][key[1:]] = value
    if in_table:
        raise ValueError("Unclosed SOFT table")
    return records


def one_field(record, key, required=True):
    values = record["fields"].get(key, [])
    if len(values) != 1:
        if not values and not required:
            return None
        raise ValueError(f"Expected one {key}; received {len(values)}")
    return values[0]


def characteristics(record):
    result = {}
    for item in record["fields"].get("Sample_characteristics_ch1", []):
        if ": " not in item:
            raise ValueError("Characteristic without exact colon-space separator")
        key, value = item.split(": ", 1)
        if key in result:
            raise ValueError("Duplicate characteristic key")
        result[key] = value
    return result


def exact_sample_join(columns, samples, key):
    index = {}
    for sample in samples:
        value = sample["accession"] if key == "accession" else one_field(sample, key)
        if value in index:
            raise ValueError("Duplicate metadata join key")
        index[value] = sample
    if len(columns) != len(set(columns)):
        raise ValueError("Duplicate expression column identifier")
    if set(columns) != set(index):
        raise ValueError("Expression/metadata samples do not form an exact bijection")
    return [index[column] for column in columns]


def numeric_table(path, series_matrix=False):
    """Consume every row; retain source features and audit every numeric token."""
    import numpy as np
    from decimal import Decimal, InvalidOperation

    metadata, features, rows, anomalies = {}, [], [], []
    lexical = Counter()
    columns, in_table, ended = None, not series_matrix, False
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            row = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
            if not row:
                continue
            if series_matrix and row[0] == "!series_matrix_table_begin":
                if in_table or columns is not None:
                    raise ValueError("Duplicate matrix start marker")
                in_table = True
                continue
            if series_matrix and row[0] == "!series_matrix_table_end":
                if not in_table or columns is None:
                    raise ValueError("Unmatched matrix end marker")
                in_table, ended = False, True
                continue
            if not in_table:
                if ended:
                    if any(cell.strip() for cell in row):
                        raise ValueError("Unexpected content after matrix end")
                else:
                    metadata.setdefault(row[0], []).append(row[1:])
                continue
            if columns is None:
                if row[0] != ("ID_REF" if series_matrix else ""):
                    raise ValueError("Unexpected feature header")
                columns = row[1:]
                if not columns or len(columns) != len(set(columns)) or any(not x for x in columns):
                    raise ValueError("Empty or duplicate expression column names")
                continue
            if len(row) != len(columns) + 1:
                raise ValueError(f"Ragged numerical matrix at source line {line_number}")
            if not row[0]:
                raise ValueError("Empty feature label")
            features.append(row[0])
            converted = []
            for column_index, token in enumerate(row[1:]):
                lexical["tokens"] += 1
                lexical["unsigned_integer_literals"] += bool(re.fullmatch(r"[0-9]+", token))
                lexical["decimal_point_literals"] += "." in token
                lexical["exponent_literals"] += bool(re.search(r"[eE][+-]?[0-9]+$", token))
                missing = token in ("", "NA", "N/A", "null")
                invalid = False
                try:
                    value = float(token) if not missing else float("nan")
                except ValueError:
                    value, invalid = float("nan"), True
                lexical["missing_literals"] += missing
                lexical["invalid_literals"] += invalid
                if not series_matrix and not missing and not invalid:
                    try:
                        exact = Decimal(token)
                        if exact.is_finite():
                            lexical["exact_literal_noninteger_values"] += exact != exact.to_integral_value()
                            lexical["exact_literal_negative_values"] += exact < 0
                            lexical["exact_literal_abs_gt_float64_integer_limit"] += abs(exact) > 2**53
                    except InvalidOperation:
                        lexical["invalid_decimal_literals"] += 1
                converted.append(value)
                # Only a numerical integrity flag, never an exclusion or a gene effect.
                if missing or invalid or not np.isfinite(value) or value < 0 or (series_matrix and abs(value) > 100):
                    anomalies.append({"source_line_1based": line_number, "feature_index_0based": len(features)-1,
                                      "feature_label_as_submitted": row[0], "column_index_0based": column_index,
                                      "column_id_as_submitted": columns[column_index], "literal_as_submitted": token,
                                      "reason": "missing_or_invalid" if missing or invalid else "nonfinite" if not np.isfinite(value)
                                                else "negative" if value < 0 else "array_abs_gt_100_numerical_audit"})
            rows.append(converted)
    if columns is None or not rows or (series_matrix and not ended):
        raise ValueError("Incomplete or empty numerical matrix")
    array = np.asarray(rows, dtype=np.float64)
    if series_matrix:
        accessions = metadata.get("!Sample_geo_accession", [])
        if accessions != [columns]:
            raise ValueError("Series metadata accession order differs from matrix columns")
        for key, values in metadata.items():
            if key.startswith("!Sample_") and any(len(value) != len(columns) for value in values):
                raise ValueError("Sample matrix metadata width differs from expression width")
    return {"features": features, "columns": columns, "values": array, "metadata": metadata,
            "lexical_qc": dict(lexical), "anomalies": anomalies, "table_end_verified": ended if series_matrix else None}


def numeric_qc(table, quantity):
    import numpy as np

    values = table["values"]
    finite_mask = np.isfinite(values)
    finite = values[finite_mask]
    quantiles = [0, .001, .01, .25, .5, .75, .99, .999, 1]
    feature_counts = Counter(table["features"])
    fingerprints = defaultdict(list)
    for index, column in enumerate(table["columns"]):
        fingerprints[hashlib.sha256(values[:, index].astype("<f8").tobytes()).hexdigest()].append(column)
    duplicate_columns = [group for group in fingerprints.values() if len(group) > 1]
    duplicate_rows = Counter(hashlib.sha256(row.astype("<f8").tobytes()).hexdigest() for row in values)
    result = {
        "quantity": quantity, "features": values.shape[0], "columns": values.shape[1], "cells": values.size,
        "finite_cells": int(finite.size), "nonfinite_cells": int(values.size-finite.size),
        "negative_cells": int(np.count_nonzero(finite < 0)), "zero_cells": int(np.count_nonzero(finite == 0)),
        "noninteger_cells": int(np.count_nonzero(finite != np.floor(finite))),
        "absolute_values_over_float64_exact_integer_limit": int(np.count_nonzero(np.abs(finite) > 2**53)),
        "global_quantiles": {str(q): float(np.quantile(finite, q)) for q in quantiles} if finite.size else {},
        "global_quantile_method": "numpy.quantile(method=linear); finite cells only; all features/samples pooled",
        "unique_feature_labels": len(feature_counts), "duplicate_feature_label_groups": sum(n > 1 for n in feature_counts.values()),
        "duplicate_feature_label_extra_rows": sum(n-1 for n in feature_counts.values()),
        "exact_duplicate_column_groups": len(duplicate_columns), "exact_duplicate_columns_beyond_first": sum(len(g)-1 for g in duplicate_columns),
        "all_zero_feature_rows": int(np.count_nonzero(np.all(values == 0, axis=1))),
        "constant_finite_feature_rows": int(np.count_nonzero(np.all(finite_mask, axis=1) & np.all(values == values[:, :1], axis=1))),
        "exact_duplicate_numeric_row_groups": sum(n > 1 for n in duplicate_rows.values()),
        "numeric_row_duplicates_are_not_assumed_to_be_duplicate_features": True,
        "all_zero_columns": int(np.count_nonzero(np.all(values == 0, axis=0))),
        "literal_audit": table["lexical_qc"],
        "feature_order_sha256": hashlib.sha256(("\n".join(table["features"])+"\n").encode()).hexdigest(),
        "column_order_sha256": hashlib.sha256(("\n".join(table["columns"])+"\n").encode()).hexdigest(),
        "float64_numeric_sha256": hashlib.sha256(values.astype("<f8").tobytes()).hexdigest(),
    }
    per_sample = []
    for index, column in enumerate(table["columns"]):
        vector = values[:, index]
        good = vector[np.isfinite(vector)]
        record = {"column_id_as_submitted": column, "finite_features": int(good.size), "zero_features": int(np.count_nonzero(vector == 0)),
                  "positive_features": int(np.count_nonzero(vector > 0)), "min": float(good.min()) if good.size else None,
                  "max": float(good.max()) if good.size else None, "median": float(np.median(good)) if good.size else None,
                  "exact_float64_fingerprint": hashlib.sha256(vector.astype("<f8").tobytes()).hexdigest()}
        if quantity == "submitted_raw_gene_counts":
            record["total_counts"] = float(good.sum())
        else:
            record["abs_gt_100_cells"] = int(np.count_nonzero(np.abs(vector) > 100))
        per_sample.append(record)
    if quantity == "submitted_raw_gene_counts":
        totals = [row["total_counts"] for row in per_sample]
        detected = [row["positive_features"] for row in per_sample]
        result["library_total_counts_min_median_max"] = [float(np.min(totals)), float(np.median(totals)), float(np.max(totals))]
        result["positive_features_per_library_min_median_max"] = [int(min(detected)), float(np.median(detected)), int(max(detected))]
        result["count_compatible"] = (all(result[key] == 0 for key in ["nonfinite_cells", "negative_cells", "noninteger_cells", "absolute_values_over_float64_exact_integer_limit"])
                                      and all(table["lexical_qc"].get(key, 0) == 0 for key in ["invalid_literals", "missing_literals", "invalid_decimal_literals", "exact_literal_noninteger_values", "exact_literal_negative_values", "exact_literal_abs_gt_float64_integer_limit"]))
    else:
        result["absolute_value_threshold_audit"] = {
            str(threshold): {"cells": int(np.count_nonzero(np.abs(values) > threshold)),
                             "features": int(np.count_nonzero(np.any(np.abs(values) > threshold, axis=1))),
                             "columns": int(np.count_nonzero(np.any(np.abs(values) > threshold, axis=0)))}
            for threshold in [20, 100, 1000, 100_000, 1_000_000]}
        integer_mask = finite == np.floor(finite)
        result["integer_value_cells"] = int(np.count_nonzero(integer_mask))
        result["integer_value_cells_gt_100"] = int(np.count_nonzero(integer_mask & (finite > 100)))
        for label, subset in [("abs_le_100", finite[np.abs(finite) <= 100]), ("abs_gt_100", finite[np.abs(finite) > 100])]:
            result[label + "_min_max"] = [float(subset.min()), float(subset.max())] if subset.size else None
        result["abs_gt_100_is_numerical_integrity_flag_not_exclusion_rule"] = True
        result["expression_scale_qualified"] = bool(result["nonfinite_cells"] == 0 and not np.any(np.abs(values) > 100))
    return result, per_sample, duplicate_columns


def annotation_table(path):
    metadata, header, rows, in_table, ended = {}, None, [], False, False
    with open_text(path) as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
            elif line == "!platform_table_end":
                in_table, ended = False, True
            elif in_table:
                row = next(csv.reader([line], delimiter="\t"))
                if header is None:
                    header = row
                else:
                    if len(row) != len(header):
                        raise ValueError("Ragged annotation row")
                    rows.append(dict(zip(header, row)))
            elif " = " in line:
                key, value = line.split(" = ", 1)
                metadata[key] = value
    if not ended or not rows or not header or "ID" not in header:
        raise ValueError("Incomplete annotation")
    ids = [row["ID"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate annotation ID")
    return {"metadata": metadata, "columns": header, "rows": rows}


def split_annotation_ids(value):
    # The historical GEO annotation uses the literal multi-value delimiter ///.
    tokens = [token.strip() for token in value.split("///")]
    return list(dict.fromkeys(token for token in tokens if token not in ("", "---")))


def annotation_join(features, annotation):
    index = {row["ID"]: row for row in annotation["rows"]}
    joined, statuses, id_use = [], Counter(), Counter()
    for feature in features:
        row = index.get(feature)
        if row is None:
            status, ids, symbols = "absent_from_annotation", [], []
        else:
            ids, symbols = split_annotation_ids(row["Gene ID"]), split_annotation_ids(row["Gene symbol"])
            if any(not re.fullmatch(r"[0-9]+", value) for value in ids):
                raise ValueError("Unrecognized Entrez identifier syntax")
            status = "no_entrez_id" if not ids else "one_entrez_id" if len(ids) == 1 else "multiple_entrez_ids"
        statuses[status] += 1
        if len(ids) == 1:
            id_use[ids[0]] += 1
        joined.append({"feature_id_as_submitted": feature, "mapping_status": status,
                       "historical_entrez_ids": ids, "historical_gene_symbols": symbols,
                       "annotation_raw": row})
    summary = {"annotation_feature_rows": len(index), "measurement_feature_rows": len(features), "join_status_counts": dict(statuses),
               "unique_entrez_ids_among_one_id_features": len(id_use), "entrez_ids_with_multiple_measurement_features": sum(n > 1 for n in id_use.values()),
               "multiple_features_per_gene_not_collapsed": True, "annotation_date_as_submitted": annotation["metadata"].get("!Annotation_date"),
               "annotation_platform_as_submitted": annotation["metadata"].get("!Annotation_platform"), "annotation_columns": annotation["columns"],
               "current_gene_identifier_reconciliation_performed": False, "unreleased_measurements_not_imputed": True}
    return summary, joined


def quick_table_crosscheck(quick_samples, table):
    features = {feature: index for index, feature in enumerate(table["features"])}
    columns = {column: index for index, column in enumerate(table["columns"])}
    checked, disagreements, large_checked = 0, [], 0
    for sample in quick_samples:
        col = columns[sample["accession"]]
        for row in sample["table"]:
            expected = float(row["VALUE"])
            observed = float(table["values"][features[row["ID_REF"]], col])
            checked += 1
            large_checked += abs(expected) > 100
            if expected != observed:
                disagreements.append({"sample": sample["accession"], "feature": row["ID_REF"], "quick_literal": row["VALUE"], "matrix_float64": observed})
    return {"cached_quick_cells_checked": checked, "numerical_disagreements": len(disagreements),
            "quick_cells_abs_gt_100_also_in_full_matrix": large_checked, "no_quick_preview_used_as_full_feature_universe": True}, disagreements


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def matrix_qc(raw=RAW, work=WORK):
    import numpy as np

    raw, work = Path(raw), Path(work)
    work.mkdir(parents=True, exist_ok=True)
    rna = numeric_table(raw / "GSE303271_raw_counts.txt.gz")
    array = numeric_table(raw / "GSE159676_series_matrix.txt.gz", series_matrix=True)
    rna_qc, rna_sample_qc, rna_duplicate_columns = numeric_qc(rna, "submitted_raw_gene_counts")
    array_qc, array_sample_qc, array_duplicate_columns = numeric_qc(array, "submitted_processed_array_RMA_declared_by_methods")
    annotation = annotation_table(raw / "GPL6244.annot.gz")
    mapping_qc, mapping = annotation_join(array["features"], annotation)
    quick = [record for record in geo_records(raw / "GSE159676-samples-quick.txt", include_tables=True) if record["entity"] == "SAMPLE"]
    crosscheck, disagreements = quick_table_crosscheck(quick, array)
    for accession, table, sample_qc, duplicates in [
            ("GSE303271", rna, rna_sample_qc, rna_duplicate_columns),
            ("GSE159676", array, array_sample_qc, array_duplicate_columns)]:
        np.save(work / (accession + "-submitted-float64.npy"), table["values"], allow_pickle=False)
        write_tsv(work / (accession + "-features-as-submitted.tsv"),
                  [{"row_index_0based": index, "feature_label_as_submitted": feature} for index, feature in enumerate(table["features"])],
                  ["row_index_0based", "feature_label_as_submitted"])
        dump_json(work / (accession + "-column-qc.json"), sample_qc)
        dump_json(work / (accession + "-numeric-anomalies.json"), table["anomalies"])
        dump_json(work / (accession + "-exact-duplicate-column-groups.json"), duplicates)
    dump_json(work / "GSE159676-probe-annotation-join.json", mapping)
    dump_json(work / "GSE159676-quick-matrix-disagreements.json", disagreements)
    summary = {"schema_version": 1, "GSE303271": rna_qc, "GSE159676": array_qc,
               "GSE159676_annotation": mapping_qc, "GSE159676_quick_matrix_crosscheck": crosscheck,
               "gene_effects_or_program_scores_computed": False}
    dump_json(work / "matrix-qc.json", summary)
    return summary, rna, array



def metadata_structure(rna, array, raw=RAW, work=WORK):
    raw, work = Path(raw), Path(work)
    rna_records = [row for row in geo_records(raw / "GSE303271-family.soft.gz") if row["entity"] == "SAMPLE"]
    array_records = [row for row in geo_records(raw / "GSE159676-samples-quick.txt") if row["entity"] == "SAMPLE"]
    rna_join = exact_sample_join(rna["columns"], rna_records, "Sample_title")
    array_join = exact_sample_join(array["columns"], array_records, "accession")
    summaries = {}
    for accession, records, table in [("GSE303271", rna_join, rna), ("GSE159676", array_join, array)]:
        rows, groups, char_keys, duplicate_keys = [], Counter(), Counter(), {}
        for index, record in enumerate(records):
            chars = characteristics(record)
            char_keys.update(chars.keys())
            condition = chars["sample group"] if accession == "GSE303271" else chars["condition"]
            groups[condition] += 1
            title = one_field(record, "Sample_title")
            library_description = [item.removeprefix("Library name: ") for item in record["fields"].get("Sample_description", []) if item.startswith("Library name: ")]
            if accession == "GSE303271" and library_description != [table["columns"][index]]:
                raise ValueError("RNA count label/title/library-description identity does not match")
            source_number, suffix = None, None
            if accession == "GSE159676" and condition == "Primary sclerosing cholangitis":
                match = re.fullmatch(r"Primary_sclerosing_cholangitis_([0-9]+)_([ab]) PSC_([0-9]+)_([ab])", title)
                if not match or match.group(1, 2) != match.group(3, 4):
                    raise ValueError("PSC source-number title parsing is not an exact redundant match")
                source_number, suffix = match.group(1, 2)
            record_row = {
                "matrix_column_index_0based": index, "matrix_column_id_as_submitted": table["columns"][index],
                "geo_accession": record["accession"], "title_as_submitted": title, "condition_as_submitted": condition,
                "source_name_as_submitted": one_field(record, "Sample_source_name_ch1"), "characteristics_as_submitted": chars,
                "library_name_as_submitted": library_description[0] if library_description else None,
                "platform_as_submitted": one_field(record, "Sample_platform_id"),
                "relations_as_submitted": record["fields"].get("Sample_relation", []),
                "supplementary_files_as_submitted": record["fields"].get("Sample_supplementary_file", []),
                "psc_number_in_title": source_number, "psc_ab_suffix_in_title": suffix,
                "separate_published_patient_identifier": None, "region_interpretation_from_title_alone": None,
                "patient_age": None, "patient_sex": None, "fibrosis_stage": None, "treatment": None,
                "ibd_status": None, "assay_batch": None,
            }
            rows.append(record_row)
        if accession == "GSE159676":
            for key, repeated_values in array["metadata"].items():
                if not key.startswith("!Sample_"):
                    continue
                if all(key[1:] in record["fields"] for record in records):
                    expected = [record["fields"][key[1:]] for record in records]
                    actual = [[value[index] for value in repeated_values] for index in range(len(records))]
                    if actual != expected:
                        raise ValueError(f"Quick-source and full-matrix sample metadata differ at {key}")
        for key in ["geo_accession", "title_as_submitted", "library_name_as_submitted"]:
            values = [row[key] for row in rows if row[key] is not None]
            counts = Counter(values)
            duplicate_keys[key] = {"nonmissing_rows": len(values), "distinct_values": len(counts), "duplicate_groups": sum(n > 1 for n in counts.values())}
        relations = [value for row in rows for value in row["relations_as_submitted"]]
        relation_types = {label: [value for value in relations if value.startswith(label + ": ")] for label in ["BioSample", "SRA"]}
        original_processing = sorted(set(value for record in records for value in record["fields"].get("Sample_data_processing", [])))
        groups_summary = {"measurement_columns": len(rows), "exact_metadata_join_rows": len(rows), "original_condition_counts": dict(groups),
                          "join_rule": "count header = GEO title = explicit library name" if accession == "GSE303271" else "matrix GSM header = source GSM accession; shared metadata values agree exactly",
                          "identifier_uniqueness": duplicate_keys, "sample_characteristic_key_counts": dict(char_keys),
                          "source_name_counts": dict(Counter(row["source_name_as_submitted"] for row in rows)),
                          "platform_counts": dict(Counter(row["platform_as_submitted"] for row in rows)),
                          "source_processing_declarations": original_processing,
                          "value_column_declarations": sorted(set(record["column_descriptions"].get("VALUE", "") for record in records)-{""}),
                          "relation_identity_counts": {key: {"records": len(values), "distinct": len(set(values))} for key, values in relation_types.items()},
                          "covariates_linked_to_matrix_columns_from_GEO": ["condition", "tissue_or_source_name", "platform"],
                          "clinical_age_sex_stage_treatment_ibd_batch_values_in_GEO": 0,
                          "group_level_covariates_assigned_to_individuals": False}
        if accession == "GSE159676":
            blocks = defaultdict(list)
            for row in rows:
                if row["psc_number_in_title"] is not None:
                    blocks[row["psc_number_in_title"]].append(row["psc_ab_suffix_in_title"])
            if len(blocks) != 6 or any(sorted(values) != ["a", "b"] for values in blocks.values()):
                raise ValueError("Adult PSC title structure does not form six complete a/b blocks")
            groups_summary["psc_source_number_blocks"] = len(blocks)
            groups_summary["psc_complete_ab_blocks"] = sum(sorted(values) == ["a", "b"] for values in blocks.values())
            groups_summary["ab_region_meaning_provided_by_GEO"] = False
        summaries[accession] = groups_summary
        dump_json(work / (accession + "-sample-source-join.json"), rows)
    dump_json(work / "sample-structure.json", summaries)
    return summaries



def verify_pin(path, expected_bytes, expected_sha256):
    path = Path(path)
    if path.stat().st_size != expected_bytes or sha256(path) != expected_sha256:
        raise ValueError(f"Input pin mismatch: {path}")
    return {"path": str(path.relative_to(ROOT)), "bytes": expected_bytes, "sha256": expected_sha256}


def verify_sources(raw=RAW):
    raw = Path(raw)
    manifest = json.loads((raw / "source-manifest.json").read_text())
    pins = {}
    def check(path, size, digest):
        record = verify_pin(path, size, digest)
        pins[record["path"]] = record
    for receipt in manifest["sources"]:
        if receipt.get("sha256"):
            check(raw / receipt["file"], receipt["bytes"], receipt["sha256"])
        if receipt.get("copied_from"):
            check(ROOT / receipt["copied_from"], receipt["bytes"], receipt["sha256"])
    methods_receipts = json.loads((raw / "methods/new-source-receipts.json").read_text())
    for receipt in methods_receipts:
        check(raw / "methods" / receipt["file"], receipt["response_body_bytes_saved"], receipt["sha256"])
    for receipt in json.loads((raw / "methods/reused-source-receipts.json").read_text()):
        check(ROOT / receipt["source_path"], receipt["source_bytes"], receipt["source_sha256"])
        if sha256(ROOT / receipt["original_receipt_file"]) != receipt["original_receipt_file_sha256"]:
            raise ValueError("Reused Methods receipt file changed")
    import zlib
    for receipt in json.loads((raw / "methods/pediatric-xlsx-member-receipts.json").read_text()):
        path = ROOT / receipt["extracted_path"]
        check(path, receipt["extracted_bytes"], receipt["sha256"])
        if zlib.crc32(path.read_bytes()) != receipt["member_metadata"]["crc32"]:
            raise ValueError("ZIP-extracted metadata CRC mismatch")
    atlas_receipts = [json.loads(line) for line in (raw / "atlas/request-ledger.jsonl").read_text().splitlines()]
    for receipt in atlas_receipts:
        if receipt.get("body_path"):
            check(ROOT / receipt["body_path"], receipt["received_body_bytes"], receipt["body_sha256"])
    for receipt in json.loads((raw / "atlas/cache-reuse.json").read_text()):
        check(ROOT / receipt["path"], receipt["bytes"], receipt["sha256"])
    main_new = sum(row["new_bytes"] for row in manifest["sources"])
    methods_new = sum(row["response_body_bytes_saved"] for row in methods_receipts)
    atlas_new = sum(row["received_body_bytes"] for row in atlas_receipts)
    atlas_bound = sum(max(row["received_body_bytes"], row.get("declared_body_bytes") or 0)
                      if row["method"] == "GET" and row.get("body_state", "").startswith("stopped_")
                      else row["received_body_bytes"] for row in atlas_receipts)
    if main_new > MAIN_NEW_LIMIT or methods_new > 12_000_000 or atlas_bound > 10_000_000 or main_new+methods_new+atlas_bound > NEW_TOTAL_LIMIT:
        raise ValueError("Authorized source acquisition allocation exceeded")
    budget = {"overall_limit_bytes": NEW_TOTAL_LIMIT, "main_limit_bytes": MAIN_NEW_LIMIT, "main_new_body_bytes": main_new,
              "methods_limit_bytes": 12_000_000, "methods_new_body_bytes": methods_new,
              "atlas_limit_bytes": 10_000_000, "atlas_new_retained_body_bytes": atlas_new,
              "atlas_conservative_full_advertised_prefix_GET_charge": atlas_bound,
              "new_body_bytes_read_and_retained_total": main_new+methods_new+atlas_new,
              "conservative_atlas_inclusive_budget_charge": main_new+methods_new+atlas_bound,
              "within_all_allocations": True, "copied_cache_bodies_charged_as_new_bytes": False,
              "oversize_header_only_refusals": [{"url": row["url"], "declared_bytes": int(row["headers"]["Content-Length"]),
                                                 "retained_body_bytes": row["bytes"]} for row in manifest["sources"] if row.get("failure") == "Content-Length exceeds remaining authorized allocation"]}
    return {"verified_unique_input_paths": len(pins), "input_pins": sorted(pins.values(), key=lambda row: row["path"]), "acquisition_budget": budget}


def verify_closed_outputs(work=WORK):
    scope = json.loads((Path(work) / "scope.json").read_text())
    expected = scope["organoid_closed_hashes"]
    actual = {path: sha256(ROOT / path) for path in expected}
    if actual != expected:
        raise ValueError("Closed organoid outputs changed")
    material_pins = scope.get("organoid_prior_material_pins", [])
    for record in material_pins:
        verify_pin(ROOT/record["path"], record["bytes"], record["sha256"])
    return {"closed_versioned_organoid_pins_unchanged": True, "pins": expected,
            "prior_organoid_source_and_ignored_artifact_pins_checked": len(material_pins),
            "prior_organoid_material_pins_unchanged": True,
            "no_organoid_files_written_by_this_script": True}


def make_report(result):
    pediatric, adult = result["datasets"]["GSE303271"], result["datasets"]["GSE159676"]
    rna, array = pediatric["numeric_qc"], adult["numeric_qc"]
    budget, annotation = result["provenance"]["acquisition_budget"], adult["feature_annotation"]
    flag = array["absolute_value_threshold_audit"]["100"]
    atlas = result["NoPSC_atlas"]
    n = lambda value: f"{value:,}"
    return f"""# PSC liver input qualification

## Decision

This is input qualification only. No gene contrasts, target/program scores, models or biological effects were computed. The organoid handoff remains closed and unchanged. The root analysis owns any later panel and analysis freeze.

| Resource | Qualified interface | Decision |
|---|---|---|
| GSE303271 | Complete {n(rna['features'])}-label by {rna['columns']}-sample submitted count matrix; exact GEO/library joins | Count input usable with source-annotation and clinical-metadata limits |
| GSE159676 | Complete {n(array['features'])}-feature by {array['columns']}-array matrix; historical platform join; paired-biopsy structure | **Hold expression analysis:** mixed-scale source literals and unresolved cohort labels |
| NoPSC liver atlas | Bounded browser label/schema metadata only | **Park:** no qualified original-count/feature/barcode/donor export |

## GSE303271: pediatric whole-biopsy bulk RNA

- The full raw-count download is 3,021,812 bytes. Its {n(rna['cells'])} values are finite, nonnegative, exact integer literals. All {n(rna['features'])} submitted labels are unique; all {rna['columns']} columns join exactly to GEO titles and explicit library names. There are 64 distinct GSM, BioSample and SRA accessions, no repeated column identifiers, and no exactly duplicated numeric columns. These observations alone do not prove independent patients.
- Original groups are **17 PSC, 17 ASC and 30 AIH**. The paper describes 64 patients, including 34 labelled PSC, while the GEO series uses PSC n=34. The sample-level deposit splits these into PSC and ASC. Retain ASC separately: no explicit ASC recoding rule was recovered.
- Tissue is bulk RNA from cryopreserved whole-biopsy liver tissue obtained during clinically indicated diagnostic biopsies. This is not a whole organ or sorted portal fibroblasts, macrophages or cholangiocytes. No paired-region scheme is stated for this bulk deposit; the paper's separate two-patient spatial assay is not extra bulk replication.
- Keep all {n(rna['all_zero_feature_rows'])} all-zero rows in the qualified released universe. There are {n(rna['zero_cells'])} zero cells overall. Per-library count sums range from {n(int(rna['library_total_counts_min_median_max'][0]))} to {n(int(rna['library_total_counts_min_median_max'][2]))}; median {n(rna['library_total_counts_min_median_max'][1])}. These are sums of deposited assignments, not independently verified read totals. No library or feature was dropped, rescaled or normalized.
- GEO names **GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4**. The paper's DESeq2 1.26 normalization/variance stabilization is downstream analysis, not the downloaded raw-count quantity. The GTF/provider/release/checksum, assignment/strandedness/multimapping options and executed command are unavailable in inspected sources. The source read-length/depth wordings also differ; neither is silently corrected.
- The 41,096 row labels are preserved as submitted. This is the complete **released** feature universe, not a reconstructed annotation universe or a current validated gene-ID mapping. No aliases were silently updated.
- Exact sample-linked age, sex, fibrosis/stage, treatment, IBD status and assay batch were not recovered. The main paper mentions clinical measurements, but linked figure-supplement/workbook metadata did not provide a verified sample crosswalk. Uninspected figure sheets are not evidence that such information cannot exist. Cohort summaries, library-number order and submission dates were not assigned as individual covariates.

**Feasible later comparisons:** PSC versus AIH; ASC versus AIH; PSC versus ASC, with the original labels fixed. These are pediatric disease-control comparisons. There is no healthy-control group, no stage-adjusted interface without additional exact covariates, and no cell-specific causal interpretation. Sample independence remains supported by the paper's 64-patient statement rather than a separate person/visit crosswalk. No comparison was run.

## GSE159676: adult liver-tissue arrays

### Samples, donors and source conflicts

The complete 2,252,130-byte series matrix contains 33 exact GSM-matched records: **12 PSC, 7 NASH, 3 Primary biliary cirrhosis, 3 Autoimmune hepatitis, 1 Haemochromatosis, 1 Alcohol related and 6 Liver tissue healthy**. These are source labels, not a corrected diagnosis table.

1. **PSC is six patients with two biopsies per patient, not 12 independent patients.** The public HAL author proof states this explicitly (physical PDF p44, printed p42, Supplementary Methods; repeated on physical p101). GEO titles form six complete numeric a/b blocks. The title-based blocks are consistent with that design, but a separate GSM-to-person key and the anatomical meaning/date of a/b were not recovered. Do not label them as named regions or independent replication.
2. Both Methods and GEO prose describe **five** normal/control patients, but the deposit has **six** control-style records. There is no recovered explanation for the extra record, duplicated donor, or historical exclusion. All six records remain preserved; no control was dropped.
3. The HAL proof reports **PBC 2 plus sarcoidosis 1**; GEO instead labels **Primary biliary cirrhosis 3** and has no sarcoidosis label. This is a source-version conflict, not permission to relabel any sample. HAL is an author proof, not silently represented as the final publisher supplement.
4. Microarray RNA comes from frozen liver tissue obtained from explants or diagnostic biopsies. The human array Methods do not describe sorting or microdissection. They are distinct from the paper's portal-fibroblast experiments. The control tissue is tumor-free liver from patients with colorectal cancer metastases, not healthy-volunteer biopsies. Patient-specific procurement, region, fibrosis/stage, age, sex and batch are not linked to the deposited arrays. The separate French RT-qPCR cohort's fibrosis strata must not be transferred to them.
5. The Norwegian PSC biobank is named. Overlap with the later Norwegian atlas remains unresolved, not proven absent.

### Numerical integrity: analysis hold

The GEO pre-table metadata declares **Quantile normalized** for all 33 arrays. The HAL Methods name `oligo`, `rma` and `hugene10sttranscriptcluster.db` (physical PDF p46, printed p44, Supplementary Methods). Exact versions, executed options, explicit deposited log base and export/selection history are not supplied. Conventional RMA defaults are not an executed run record. Do not log-transform this release again by assumption.

All {n(array['cells'])} source values are finite and positive, but **{n(flag['cells'])} integer literals are greater than 100**, spanning **{n(flag['features'])} features and all {flag['columns']} arrays**. Their range is **{n(int(array['abs_gt_100_min_max'][0]))} to {n(int(array['abs_gt_100_min_max'][1]))}**. The other values range from **{array['abs_le_100_min_max'][0]} to {array['abs_le_100_min_max'][1]}**. The same 4,829 values exceed 20; the >100 flag is a numerical-integrity diagnostic, not a gene or sample exclusion rule. No target list was used.

All 660 available cells in the earlier cached GEO quick preview match the complete matrix numerically, including nine large literals. This traces the issue to released source values rather than this parser; it does not identify the historical cause or establish that the authors' underlying analysis was wrong. **No decimal insertion, division, imputation, extra log, sample exclusion or feature exclusion was applied.** A source-corrected processed release or a separately authorized raw-array reconstruction is needed before expression use. Gene effects cannot be rescued by simply treating these values as counts.

### Features and annotation

- All {n(array['features'])} submitted feature IDs are unique and exactly match GPL6244's complete **33,297-row** GEO annotation (dated **Aug 09 2016**). No numeric columns are exactly duplicated. Identical numeric feature vectors can occur; they were not assumed to be duplicate IDs or collapsed.
- The platform is `[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]`. Platform metadata records NetAffx build 35 in July 2016; that is an annotation release label, not a claimed genome build. No coordinates were inferred or lifted.
- Historical Entrez mapping has **{n(annotation['join_status_counts']['one_entrez_id'])} one-ID features**, **{annotation['join_status_counts']['multiple_entrez_ids']} multi-ID features**, and **{annotation['join_status_counts']['no_entrez_id']} without an Entrez ID**. The one-ID set represents **{n(annotation['unique_entrez_ids_among_one_id_features'])} unique Entrez IDs**; **{annotation['entrez_ids_with_multiple_measurement_features']} IDs have more than one measured feature**. No mapping ambiguity or multiple-feature aggregation was silently resolved.
- The historical GEO annotation is not claimed to be the exact `hugene10sttranscriptcluster.db` snapshot used by the authors. Current gene-ID reconciliation and the source rule selecting 17,046 of 33,297 platform rows remain unavailable. Missing platform measurements were not filled in.

**Conditional later comparisons:** PSC versus NASH is a disease-control interface with six paired-biopsy PSC blocks, not n=12 independent PSC. PSC versus control-style tissue additionally needs the five/six control reconciliation. Comparisons using PBC/other labels need their source-version conflict resolved. Every expression comparison remains on hold for numerical integrity and requires a frozen donor/block and feature-aggregation rule. No comparison was run.

## NoPSC atlas: one bounded linked-export check

The normal public browser serves search-label lists, annotation/display bundles and five labelled differential-expression spreadsheet links. Those spreadsheet payloads and all individual-gene endpoints were **not** requested. No original-count bulk inventory, stable feature/barcode axis or donor-to-observation join was established.

| Browser metadata | Nuclear | Spatial |
|---|---:|---:|
| Unique search-menu labels | 25,539 | 17,319 |
| Annotation-prefix observations | 72,990 | 45,998 |
| Complete categorical fields inspected | 2 | 6 |
| Declared whole annotation-bundle bytes | 1,314,400 | 1,749,037 |
| Retained prefix bytes | 146,317 | 276,568 |

The label lists include `Total RNA count` and `Total RNA features`, so their lengths are **not** verified gene-universe counts. The annotation streams were closed before the first continuous payload. Nuclear prefixes contain disease and cell-class labels; spatial prefixes contain disease, subtype, stage, location and spatial-class labels. No donor-ID or barcode field occurs in these prefixes, but uninspected tails are not evidence of absence. Observation row counts are not independent donor counts.

The paper reports **12 PSC + 4 cirrhotic-control nuclear donors**, overlapping **23 PSC + 7 cirrhotic-control spatial donors**. Modalities are not replication. The browser's measurement quantity is not established by the article's SCTransform description. Norwegian-biobank overlap with GSE159676 remains unresolved.

**Reopen only** for a normal public original-count release with assay/units, complete feature IDs, a stable barcode axis, anonymous donor mapping and file URLs/sizes. Report that inventory and obtain root authorization before bulk transfer. No usable bulk inventory was found, so no bulk request is pending. Exact endpoint URLs, sizes and partial-schema limits are in the [machine-readable result](../data/derived/psc-liver-input-qualification.json).

## Provenance, boundaries and reproduction

- Complete inputs, failed/size-refused requests, reused cache receipts and SHA-256 pins: `data/raw/psc-liver-qualification/`. Source joins, unchanged-value numeric arrays, full feature/annotation tables and integrity flags stay in ignored `work/psc-liver-qualification/`.
- New bodies read/retained total **{n(budget['new_body_bytes_read_and_retained_total'])} bytes**. Charging the full advertised bodies of the two stopped atlas-prefix GETs gives **{n(budget['conservative_atlas_inclusive_budget_charge'])} bytes**, below the 80,000,000-byte overall cap and each worker allocation. Oversize requests were stopped from headers without body reads. No FASTQ, CEL/raw-array archive, atlas matrix, researcher contact, account or environment change was used.
- Matrix gzip CRC/end markers, all numerical cells, sample joins, original labels, feature uniqueness and mapping ambiguity were checked. Synthetic tests cover truncation, malformed joins, duplicate IDs, missing/fractional/nonfinite values, float-rounding loss and bounded acquisition. An independent stdlib audit rechecks complete matrices and metadata without importing the qualifier. The source-level audit is separate from inference.
- Only lossless decompression, TSV/CSV parsing, exact metadata/annotation joins, documented a/b title parsing, and float64 copies for QC were applied. Source files/literals remain unchanged. Public outputs contain aggregate qualification only. No RNA normalization, gene aggregation, gene contrast, target/program score, model, effect or clinical recommendation was produced.
- Python/package versions, source URLs/retrieval times, Methods locations, all assay/annotation caveats, code hashes and ignored-artifact pins are in the [JSON qualification](../data/derived/psc-liver-input-qualification.json). All four closed versioned organoid pins and {result['closed_organoid_preservation']['prior_organoid_source_and_ignored_artifact_pins_checked']} preceding source/ignored-artifact pins remain unchanged; no organoid file was written.

Offline commands:

```sh
.venv/bin/python scripts/qualify_psc_liver_inputs.py --qualify
.venv/bin/python scripts/qualify_psc_liver_inputs.py --verify
.venv/bin/python -m unittest discover -s tests -p test_psc_liver_inputs.py -v
```

Acquisition is a separate, explicit `--acquire` action. Qualification and verification make no network requests. Historical unrelated repository tests may still require their own absent ETS2/GSE84161 caches or R runtime; this task does not alter those closed branches.

### Main primary sources

- [GSE303271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303271); [pediatric paper](https://doi.org/10.1172/jci.insight.199226), Methods/RNA-Seq acquisition and Data availability.
- [GSE159676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159676); [adult article](https://doi.org/10.1002/hep.32456); [public HAL author proof](https://hal.science/hal-03667547/file/HEP-21-2383.R1_Proof_hi.pdf), physical pp44 and46.
- [GPL6244 historical annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz).
- [NoPSC nuclear browser](https://riim-data.no/NoPSC_Liver_Atlas/SINGLE_NUC/), [spatial browser](https://riim-data.no/NoPSC_Liver_Atlas/SPATIAL/) and [atlas article](https://doi.org/10.1097/HEP.0000000000001432).
"""


def qualify(raw=RAW, work=WORK, output=OUTPUT, report=REPORT):
    import importlib.metadata
    import numpy as np

    raw, work, output, report = map(Path, (raw, work, output, report))
    provenance = verify_sources(raw)
    closed = verify_closed_outputs(work)
    matrices, rna, array = matrix_qc(raw, work)
    structures = metadata_structure(rna, array, raw, work)
    methods_path = raw / "methods/findings.json"
    methods = json.loads(methods_path.read_text())
    atlas_path = work / "atlas/findings.json"
    atlas = json.loads(atlas_path.read_text())
    if not matrices["GSE303271"]["count_compatible"] or matrices["GSE303271"]["duplicate_feature_label_groups"]:
        raise ValueError("Raw RNA count integrity did not qualify")
    if methods["GSE159676"]["independence"]["reported_PSC_patients"] != structures["GSE159676"]["psc_source_number_blocks"]:
        raise ValueError("Adult source design does not match anonymous source-number blocks")
    if matrices["GSE159676_quick_matrix_crosscheck"]["numerical_disagreements"]:
        raise ValueError("Adult source previews disagree with full release")
    methods_provenance = {"findings_path": str(methods_path.relative_to(ROOT)), "findings_sha256": sha256(methods_path),
                          "source_catalog": methods["sources"], "provenance": methods["provenance"],
                          "bounded_search_stopping_point": methods["bounded_search_stopping_point"]}
    platform_record = [record for record in geo_records(raw / "GPL6244-platform-quick.txt") if record["entity"] == "PLATFORM"][0]
    result = {
        "schema_version": 1,
        "qualification_type": "public_liver_input_only_no_effect_analysis",
        "datasets": {}, "NoPSC_atlas": atlas,
        "provenance": {**provenance, "main_manifest": {"path": str((raw/"source-manifest.json").relative_to(ROOT)), "sha256": sha256(raw/"source-manifest.json")},
                       "main_sources": [{key: row[key] for key in ["file", "url", "response_url", "status", "bytes", "sha256", "retrieved_at_utc", "acquisition", "new_bytes", "usable", "failure"] if key in row} for row in json.loads((raw/"source-manifest.json").read_text())["sources"]],
                       "methods": methods_provenance,
                       "atlas_findings_path": str(atlas_path.relative_to(ROOT)), "atlas_findings_sha256": sha256(atlas_path)},
        "software": {"Python": platform.python_version(), "numpy": np.__version__, "requests": importlib.metadata.version("requests"),
                     "runtime": "repository .venv", "new_dependencies_installed": False,
                     "qualifier_script_sha256": sha256(Path(__file__)), "tests_sha256": sha256(ROOT/"tests/test_psc_liver_inputs.py")},
        "scope": {"primary_gene_contrasts_computed": False, "target_or_program_scores_computed": False, "models_fit": False,
                  "source_values_repaired": False, "features_or_samples_excluded": False, "genes_aggregated": False,
                  "new_gene_alias_corrections": False, "raw_reads_or_CEL_arrays_acquired": False,
                  "source_joins_and_expression_publicly_exported": False, "clinical_recommendations": False,
                  "allele_conventions": "not applicable; no variants assessed", "analysis_freeze_owned_by_root": True},
        "transformations": ["Lossless gzip decompression and delimited-text parsing; complete gzip CRC and end-marker checks",
                            "Exact matrix-to-GEO and probe-to-historical-annotation joins; no fuzzy or case-folded mapping",
                            "RNA original numeric literals checked with Decimal before float64 QC copies",
                            "Original PSC a/b title blocks parsed without assigning anatomical regions",
                            "All released features/samples retained; no normalization, aggregation, imputation or correction"],
        "closed_organoid_preservation": closed,
    }
    result["datasets"]["GSE303271"] = {
        "status": "qualified_submitted_count_interface_with_metadata_and_annotation_limits",
        "eligible_for_automatic_expression_analysis": False, "root_analysis_freeze_still_required": True,
        "numeric_qc": matrices["GSE303271"], "metadata_structure": structures["GSE303271"],
        "source_methods": methods["GSE303271"],
        "feature_universe": {"complete_released_labels": len(rna["features"]), "labels_preserved_verbatim": True,
                             "current_gene_identifier_reconciliation_performed": False, "historical_annotation_reconstructed": False,
                             "identifier_namespace": "Submitted gene-name-style labels; no independent ID reconciliation",
                             "genome_build_declared": "GRCh38", "GTF_release_and_count_options_recovered": False},
        "feasible_future_interfaces": [{"contrast": "PSC versus AIH", "source_samples": {"PSC": 17, "AIH": 30}},
                                       {"contrast": "ASC versus AIH", "source_samples": {"ASC": 17, "AIH": 30}},
                                       {"contrast": "PSC versus ASC", "source_samples": {"PSC": 17, "ASC": 17}}],
        "interface_conditions": ["Freeze original source diagnoses without pooling ASC into PSC",
                                 "Use patient-level bulk-biopsy interpretation, not cell-specific causality or healthy-control inference",
                                 "Do not invent missing patient crosswalk, clinical covariates, annotation release or batch",
                                 "Freeze exact feature/gene reconciliation, QC and normalization rules before later scoring or contrasts"],
        "submitted_count_sum_not_equated_to_reads_or_unique_fragments": True,
    }
    result["datasets"]["GSE159676"] = {
        "status": "hold_expression_use_numeric_integrity_and_unresolved_source_identity",
        "eligible_for_automatic_expression_analysis": False,
        "numeric_qc": matrices["GSE159676"], "metadata_structure": structures["GSE159676"],
        "source_methods": methods["GSE159676"], "feature_annotation": matrices["GSE159676_annotation"],
        "quick_matrix_crosscheck": matrices["GSE159676_quick_matrix_crosscheck"],
        "platform_metadata": {"title": one_field(platform_record, "Platform_title"),
                              "data_row_count_as_submitted": int(one_field(platform_record, "Platform_data_row_count")),
                              "annotation_update_declarations": [value for value in platform_record["fields"].get("Platform_description", []) if "annotation table updated" in value],
                              "genome_build_used_for_analysis": None, "coordinates_used_or_lifted": False},
        "blocked_future_interfaces": ["PSC versus NASH: six PSC paired-biopsy source blocks, never twelve independent PSC donors",
                                      "PSC versus control-style liver: also requires control five/six reconciliation",
                                      "PSC versus other source diagnoses: also requires PBC/sarcoidosis source-version reconciliation"],
        "reopen_requirements": ["Source-corrected processed values or separately authorized raw-array reconstruction; no guessed decimal repair",
                                "Explicit anonymous donor/block and inclusion rules consistent with six PSC patients/two biopsies",
                                "Resolve or predefine the unresolved control count and diagnosis conflicts without source relabeling",
                                "Freeze platform-feature selection, gene-ID reconciliation and multiple-feature aggregation",
                                "Preserve unavailable clinical/stage/batch covariates and uncertain atlas overlap"],
    }
    audit_path = work/"independent-audit/independent-verification.json"
    audit = json.loads(audit_path.read_text())
    if set(audit["check_status_counts"]) != {"PASS"}:
        raise ValueError("Independent source verification has nonpassing checks")
    for key in ["source_manifest_pin", "primary_comparison_file_pin", "sample_structure_comparison_file_pin", "audit_script_pin"]:
        record = audit[key]
        verify_pin(ROOT/record["path"], record["bytes"], record["sha256"])
    result["validation"] = {"independent_audit_path": str(audit_path.relative_to(ROOT)),
                            "independent_audit_sha256": sha256(audit_path),
                            "independent_audit_check_status_counts": audit["check_status_counts"],
                            "independent_audit_scope": "Separate stdlib source parser, no import/inspection of qualifier script",
                            "independent_audit_script_pin": audit["audit_script_pin"],
                            "primary_QC_not_an_effect_analysis": True}
    public_artifacts = [work/"matrix-qc.json", work/"sample-structure.json"]
    for accession in ["GSE303271", "GSE159676"]:
        public_artifacts.extend(work/(accession+suffix) for suffix in ["-submitted-float64.npy", "-features-as-submitted.tsv", "-column-qc.json", "-numeric-anomalies.json", "-exact-duplicate-column-groups.json", "-sample-source-join.json"])
    public_artifacts.extend([work/"GSE159676-probe-annotation-join.json", work/"GSE159676-quick-matrix-disagreements.json"])
    result["ignored_artifact_pins"] = [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in public_artifacts]
    dump_json(output, result)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(make_report(result))
    dump_json(work/"qualification-output-pins.json", {"outputs": [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in [output, report]],
                                                    "qualifier_script_sha256": sha256(Path(__file__)), "tests_sha256": sha256(ROOT/"tests/test_psc_liver_inputs.py")})
    return result


def verify_qualification(raw=RAW, work=WORK, output=OUTPUT):
    result = json.loads(Path(output).read_text())
    current = verify_sources(raw)
    if current["input_pins"] != result["provenance"]["input_pins"] or current["acquisition_budget"] != result["provenance"]["acquisition_budget"]:
        raise ValueError("Qualification input/byte ledger changed")
    verify_closed_outputs(work)
    if sha256(Path(__file__)) != result["software"]["qualifier_script_sha256"] or sha256(ROOT/"tests/test_psc_liver_inputs.py") != result["software"]["tests_sha256"]:
        raise ValueError("Qualification code changed; rebuild required")
    for record in result["ignored_artifact_pins"]:
        verify_pin(ROOT/record["path"], record["bytes"], record["sha256"])
    if sha256(ROOT/result["validation"]["independent_audit_path"]) != result["validation"]["independent_audit_sha256"]:
        raise ValueError("Independent audit output changed")
    audit_script = result["validation"]["independent_audit_script_pin"]
    verify_pin(ROOT/audit_script["path"], audit_script["bytes"], audit_script["sha256"])
    for record in json.loads((Path(work)/"qualification-output-pins.json").read_text())["outputs"]:
        verify_pin(ROOT/record["path"], record["bytes"], record["sha256"])
    for record in [result["provenance"]["main_manifest"], {"path": result["provenance"]["methods"]["findings_path"], "sha256": result["provenance"]["methods"]["findings_sha256"]},
                   {"path": result["provenance"]["atlas_findings_path"], "sha256": result["provenance"]["atlas_findings_sha256"]}]:
        if sha256(ROOT/record["path"]) != record["sha256"]:
            raise ValueError("Qualification source summary changed")
    return {"verified": True, "input_paths": current["verified_unique_input_paths"], "ignored_artifacts": len(result["ignored_artifact_pins"]),
            "closed_organoid_pins": 4, "no_network_requests": True, "no_new_effects": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    parser.add_argument("--matrix-qc", action="store_true")
    parser.add_argument("--qualify", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if sum([args.acquire, args.matrix_qc, args.qualify, args.verify]) != 1:
        parser.error("Choose exactly one explicit action")
    if args.acquire:
        acquire()
    elif args.matrix_qc:
        summary, rna, array = matrix_qc()
        summary["sample_structure"] = metadata_structure(rna, array)
        print(json.dumps(summary, indent=2, sort_keys=True))
    elif args.qualify:
        result = qualify()
        print(json.dumps({"statuses": {key: value["status"] for key, value in result["datasets"].items()},
                          "atlas": result["NoPSC_atlas"]["decision"], "budget": result["provenance"]["acquisition_budget"]}, indent=2, sort_keys=True))
    elif args.verify:
        print(json.dumps(verify_qualification(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
