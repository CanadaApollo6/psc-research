#!/usr/bin/env python3
"""Bounded complete filtered-cohort qualification, not expression analysis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import threading
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, HTTPRedirectHandler, build_opener

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import qualify_psc_bach2_counts as count_parser
from scripts import qualify_psc_bach2_inputs as metadata_parser

RAW = ROOT / "data/raw/psc-bach2-qualification/cohort"
WORK = ROOT / "work/psc-bach2-cohort-qualification"
ARCHIVE_BUDGET = 298_421_791
EXPECTED = {"sample01_filtered_feature_bc_matrix.tar.gz": 54_449_959,
            "sample02_filtered_feature_bc_matrix.tar.gz": 44_540_169,
            "sample03_filtered_feature_bc_matrix.tar.gz": 43_425_341,
            "sample04_filtered_feature_bc_matrix.tar.gz": 36_719_700,
            "sample06_filtered_feature_bc_matrix.tar.gz": 36_674_386,
            "sample07_filtered_feature_bc_matrix.tar.gz": 29_629_283,
            "sample08_filtered_feature_bc_matrix.tar.gz": 52_982_953}
HOSTS = {"www.ebi.ac.uk", "ftp.ebi.ac.uk"}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def json_write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def check_url(url: str, filename: str):
    p = urlsplit(url)
    if (p.scheme != "https" or p.hostname not in HOSTS or p.port not in {None, 443}
            or p.username or p.password or p.query or p.fragment):
        raise ValueError("URL outside allowed ordinary public HTTPS file routes")
    allowed = {( "www.ebi.ac.uk", "/biostudies/files/E-MTAB-14013/" + filename),
               ( "ftp.ebi.ac.uk", "/biostudies/fire/E-MTAB-/013/E-MTAB-14013/Files/" + filename)}
    if filename not in EXPECTED or (p.hostname, p.path) not in allowed:
        raise ValueError("Only the exact remaining seven filtered archives and public file routes are authorized")


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def acquire_one(item: dict, raw_dir: Path, stop: threading.Event) -> dict:
    filename = item["filename"]
    if filename not in EXPECTED or item["expected_bytes"] != EXPECTED[filename]:
        raise ValueError("Filename or byte allowance differs from the frozen seven-file list")
    receipt_path = raw_dir / (filename + ".receipt.json")
    destination = raw_dir / filename
    if receipt_path.exists():
        old = json.loads(receipt_path.read_text())
        if (old["complete"] and destination.exists() and destination.stat().st_size == EXPECTED[filename]
                and count_parser.file_sha256(destination) == old["sha256"]):
            return old
        stop.set()
        raise ValueError("Prior incomplete attempt exists; no silent retry or renewed allowance")
    if destination.exists():
        stop.set()
        raise ValueError("Unverified source file exists; no overwrite")
    receipt = {"filename": filename, "expected_bytes": EXPECTED[filename], "url": item["url"],
               "started_at_utc": utc_now(), "requests": [], "complete": False, "body_bytes": 0, "error": None}
    response = None
    digest = hashlib.sha256()
    url = item["url"]
    opener = build_opener(NoRedirect())
    try:
        if stop.is_set():
            raise ValueError("Acquisition stopped after another uncertain request")
        for _ in range(6):
            if stop.is_set():
                raise ValueError("Acquisition stopped before another redirect request")
            check_url(url, filename)
            request = Request(url, method="GET", headers={"User-Agent": "PSC-public-cohort-qualification/1.0", "Accept-Encoding": "identity"})
            try:
                response = opener.open(request, timeout=30)
            except HTTPError as exc:
                response = exc
            receipt["requests"].append({"url": url, "method": "GET", "status": response.status, "headers": dict(response.headers.items())})
            if response.status in {301, 302, 303, 307, 308}:
                location = response.headers.get("Location")
                response.close(); response = None
                if not location:
                    raise ValueError("Redirect lacks Location")
                url = urljoin(url, location)
                continue
            if response.status != 200:
                raise ValueError("Non-200 response; no archive body read")
            length = response.headers.get("Content-Length")
            if length is None or int(length) != EXPECTED[filename]:
                raise ValueError("Content-Length differs from exact allowance; no body read")
            if response.headers.get("Content-Encoding", "identity") != "identity":
                raise ValueError("Unexpected transport encoding; no body read")
            with destination.open("xb") as output:
                while receipt["body_bytes"] < EXPECTED[filename]:
                    if stop.is_set():
                        raise ValueError("Another request became uncertain; partial source preserved")
                    part = response.read(min(65_536, EXPECTED[filename] - receipt["body_bytes"]))
                    if not part:
                        break
                    output.write(part); digest.update(part); receipt["body_bytes"] += len(part)
            if receipt["body_bytes"] != EXPECTED[filename]:
                raise ValueError("Short archive body; no retry")
            receipt["complete"] = True
            receipt["response_url"] = response.geturl()
            break
        else:
            raise ValueError("Redirect limit reached")
    except Exception as exc:
        receipt["error"] = type(exc).__name__ + ": " + str(exc)
        stop.set()
    finally:
        if response is not None:
            response.close()
    receipt["sha256"] = digest.hexdigest()
    receipt["completed_at_utc"] = utc_now()
    json_write(receipt_path, receipt)
    print(json.dumps({"filename": filename, "complete": receipt["complete"], "body_bytes": receipt["body_bytes"], "error": receipt["error"]}), flush=True)
    return receipt


def verify_closed_files(plan):
    for rel, pin in plan["closed_eight_deliverable_pins"].items():
        p = ROOT / rel
        if p.stat().st_size != pin["bytes"] or count_parser.file_sha256(p) != pin["sha256"]:
            raise ValueError("Closed deliverable changed: " + rel)


def acquire(raw_dir: Path = RAW, workers: int = 3):
    plan = json.loads((raw_dir / "plan.json").read_text())
    verify_closed_files(plan)
    items = plan["remaining_filtered_archives"]
    if {x["filename"]: x["expected_bytes"] for x in items} != EXPECTED or len(items) != 7 or sum(EXPECTED.values()) != ARCHIVE_BUDGET:
        raise ValueError("Frozen cohort manifest differs from approved files and total")
    prior = [json.loads(p.read_text()) for p in raw_dir.glob("*.tar.gz.receipt.json")]
    if (any(x.get("filename") not in EXPECTED or not x["complete"] for x in prior)
            or sum(x["body_bytes"] for x in prior) > ARCHIVE_BUDGET):
        raise ValueError("Prior incomplete/unrecognized attempt or budget uncertainty requires parent decision")
    stop = threading.Event()
    results = []
    with ThreadPoolExecutor(max_workers=min(max(workers, 1), 3)) as pool:
        futures = [pool.submit(acquire_one, item, raw_dir, stop) for item in items]
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as exc:
                stop.set()
                results.append({"complete": False, "error": type(exc).__name__ + ": " + str(exc)})
    summary = {"approved_archive_bytes": ARCHIVE_BUDGET,
               "actual_new_body_bytes": sum(x.get("body_bytes", 0) for x in results),
               "additional_metadata_body_bytes": 0,
               "complete": len(results) == 7 and all(x["complete"] for x in results),
               "sample05_redownloaded": False,
               "results": sorted(results, key=lambda x: x.get("filename", ""))}
    json_write(raw_dir / "acquisition-summary.json", summary)
    if not summary["complete"] or summary["actual_new_body_bytes"] != ARCHIVE_BUDGET:
        raise ValueError("Incomplete acquisition; stop pending parent decision")
    return summary



def verified_metadata_sources() -> tuple[dict, dict]:
    closed = json.loads((ROOT / "data/derived/psc-bach2-input-qualification.json").read_text())
    names = {"E-MTAB-14013.sdrf.txt", "scrna_meta.tsv", "bach2-e-mtab-14013-metadata"}
    pins = {}
    bodies = {}
    for item in closed["source_pins"]:
        if item["name"] not in names:
            continue
        if item["name"] in pins:
            raise ValueError("Ambiguous pinned metadata source")
        path = ROOT / item["path"]
        expected_sha = item.get("sha256", item.get("response_body_sha256"))
        if path.stat().st_size != item["bytes"] or count_parser.file_sha256(path) != expected_sha:
            raise ValueError("Fixed metadata source changed: " + item["name"])
        pins[item["name"]] = {"path": item["path"], "bytes": item["bytes"], "sha256": expected_sha,
                              "url": item["url"], "retrieved_at_utc": item["retrieved_at_utc"]}
        bodies[item["name"]] = path.read_text()
    if set(bodies) != names:
        raise ValueError("Required source pins are absent")
    return pins, bodies


def source_interfaces(sdrf_text: str, study_text: str, expected_filenames: set[str]) -> dict:
    header, records = metadata_parser.parse_sdrf(sdrf_text)
    inventory = metadata_parser.file_inventory(json.loads(study_text))
    by_sample = metadata_parser.sample_interfaces(records, inventory)
    if {r["filtered_archive"] for r in by_sample.values()} != expected_filenames:
        raise ValueError("SDRF and approved full filtered cohort differ")
    if header.count("Derived Array Data File") < 2:
        raise ValueError("Repeated SDRF file columns were not retained")
    return by_sample


def fixed_cohort_metadata(text: str, by_sample: dict) -> tuple[dict, dict]:
    header, rows = metadata_parser.read_tsv(text)
    required = {"sample_name", "barcode", "condition", "sex"}
    if len(header) != len(set(header)) or not required.issubset(header):
        raise ValueError("Fixed metadata lacks unique exact join fields")
    lookup = {name: header.index(name) for name in required}
    retained = {sample: set() for sample in by_sample}
    bare_counts = Counter()
    for values in rows:
        row = {name: values[index] for name, index in lookup.items()}
        sample = row["sample_name"]
        if sample not in by_sample:
            raise ValueError("Retained sample lacks an explicitly linked filtered archive")
        if row["condition"] != by_sample[sample]["genotype"] or row["sex"] != by_sample[sample]["sex"]:
            raise ValueError("Retained source condition/sex disagrees with linked SDRF sample")
        prefix = sample + "_"
        if not row["barcode"].startswith(prefix):
            raise ValueError("Retained barcode lacks its exact sample prefix")
        bare = row["barcode"][len(prefix):]
        if not bare or not re.fullmatch(r"[ACGT]{16}-[0-9]+", bare):
            raise ValueError("Unexpected exact retained barcode form; no suffix repair")
        if bare in retained[sample]:
            raise ValueError("Duplicate retained barcode within a source sample")
        retained[sample].add(bare)
        bare_counts[bare] += 1
    if any(not values for values in retained.values()):
        raise ValueError("An explicitly linked cohort sample has no fixed retained metadata")
    return retained, {"rows": len(rows), "samples": len(retained),
                      "source_metadata_columns_as_published": header,
                      "unique_sample_barcode_pairs": sum(map(len, retained.values())),
                      "bare_barcode_strings_repeated_across_samples": sum(n > 1 for n in bare_counts.values()),
                      "RNA_SCT_margin_columns_present": [name for name in header if re.search(r"nCount|nFeature|RNA|SCT|UMI", name)],
                      "join_rule": "Published SDRF filtered-archive basename to sample; remove only exact sample_name + '_' from that sample's retained barcode; retain sample in every cell key",
                      "no_clinical_join_or_state_selection": True}


def cell_membership(barcodes: list[str], retained: set[str]) -> dict:
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("Duplicate source matrix barcode")
    source = set(barcodes)
    return {"source_cells": len(barcodes), "fixed_retained_cells": len(retained),
            "retained_cells_matched": len(retained & source),
            "retained_cells_missing": len(retained - source),
            "source_cells_outside_fixed_retained_metadata": len(source - retained)}


def technical_totals(column_sums: dict, barcodes: list[str], retained: set[str]) -> dict:
    selected = [barcode in retained for barcode in barcodes]
    result = {}
    for kind, values in sorted(column_sums.items()):
        if len(values) != len(barcodes):
            raise ValueError("Technical margin and barcode axes disagree")
        result[kind] = {"all_source_cells_total_units_exact": str(sum(values)),
                        "fixed_retained_cells_total_units_exact": str(sum(value for value, keep in zip(values, selected) if keep)),
                        "all_source_cells_zero_total_columns": sum(value == 0 for value in values),
                        "fixed_retained_cells_zero_total_columns": sum(value == 0 and keep for value, keep in zip(values, selected)),
                        "unit_label": "released source count-format units; not an original-UMI or absolute-protein certificate"}
    return result


def axis_sha256(rows: list[list[str]]) -> str:
    canonical = json.dumps(rows, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def compare_rna_axes(axes: dict[str, list[list[str]]]) -> dict:
    if not axes or any(not rows for rows in axes.values()):
        raise ValueError("Every cohort source needs a nonempty declared Gene Expression axis")
    baseline = sorted(axes)[0]
    reference = axes[baseline]
    reference_ids = {row[0] for row in reference}
    reference_by_id = {row[0]: row for row in reference}
    comparisons = {}
    for sample, rows in sorted(axes.items()):
        ids = {row[0] for row in rows}
        by_id = {row[0]: row for row in rows}
        if len(ids) != len(rows):
            raise ValueError("Duplicate source RNA IDs; no implicit merge")
        comparisons[sample] = {
            "ordered_ID_name_type_axis_equal_to_reference": rows == reference,
            "row_count": len(rows), "canonical_ordered_axis_sha256": axis_sha256(rows),
            "IDs_only_in_this_sample": len(ids - reference_ids),
            "reference_IDs_absent_from_this_sample": len(reference_ids - ids),
            "shared_IDs_with_name_or_type_difference": sum(by_id[x] != reference_by_id[x] for x in ids & reference_ids),
            "same_ID_set_as_reference": ids == reference_ids,
        }
    equal = all(item["ordered_ID_name_type_axis_equal_to_reference"] for item in comparisons.values())
    return {"reference_sample_selected_by_filename_order_only": baseline,
            "full_original_ID_name_type_order_checked": True, "all_eight_ordered_RNA_axes_identical": equal,
            "RNA_feature_universe_rows_if_identical": len(reference) if equal else None,
            "comparisons": comparisons, "no_intersection_union_renaming_collapsing_or_imputation": True,
            "assembly_performed": False,
            "assembly_status": "not performed; qualification only" if equal else "STOP: RNA axes differ; parent choice required before assembly"}


def antibody_annotations(axes: dict[str, list[list[str]]]) -> dict:
    targets = {}
    for sample, rows in sorted(axes.items()):
        for row in rows:
            if row[2] != "Antibody Capture":
                continue
            target = tuple(row)
            targets.setdefault(target, []).append(sample)
    enumerated = [{"source_feature_id": key[0], "source_feature_name": key[1],
                   "source_feature_type": key[2], "present_in_samples": samples,
                   "explicit_BACH2_token_in_source_ID_or_name": any(re.search(r"(?<![A-Za-z0-9])BACH2(?![A-Za-z0-9])", value, re.I) for value in key[:2])}
                  for key, samples in sorted(targets.items())]
    explicit = [row for row in enumerated if row["explicit_BACH2_token_in_source_ID_or_name"]]
    return {"distinct_original_ID_name_type_triples": len(enumerated), "targets_as_published": enumerated,
            "all_target_names_equal_their_source_IDs": all(r["source_feature_name"] == r["source_feature_id"] for r in enumerated),
            "all_source_names_use_ADT_prefix": bool(enumerated) and all(r["source_feature_name"].startswith("ADT_") for r in enumerated),
            "explicit_BACH2_named_Antibody_Capture_feature_present": bool(explicit),
            "explicit_BACH2_target_rows": explicit,
            "interpretation": "Source annotation only: a named BACH2 feature is not a specificity validation" if explicit else "No deposited Antibody Capture ID or name explicitly labels BACH2; RNA BACH2 is not an ADT measurement",
            "antibody_reagent_identity_specificity_independently_reconciled": False,
            "tag_units_are_absolute_protein_molecules": False,
            "RNA_to_protein_translation_effect_tested": False}


def qualify_one(sample: str, archive_path: Path, receipt: dict, retained: set[str], work_dir: Path) -> tuple[dict, list[list[str]], list[str]]:
    if not receipt["complete"] or archive_path.stat().st_size != receipt["body_bytes"] or count_parser.file_sha256(archive_path) != receipt["sha256"]:
        raise ValueError("Archive does not match its complete acquisition receipt")
    tar_path, outer = count_parser.decompress_outer(archive_path, work_dir)
    members = count_parser.inventory_tar(tar_path)
    expected_directory = sample + "_filtered_feature_bc_matrix"
    if any(item["canonical_path"].split("/")[0] != expected_directory for item in members):
        raise ValueError("Internal source member sample prefix differs from the explicitly linked archive")
    decoded = count_parser.prepare_members(tar_path, members, work_dir)
    features, feature_summary = count_parser.parse_features(ROOT / decoded["features"]["decoded_path"])
    barcodes = count_parser.parse_barcodes(ROOT / decoded["barcodes"]["decoded_path"])
    matrix_summary, column_sums = count_parser.parse_matrix(ROOT / decoded["matrix"]["decoded_path"], features, barcodes)
    membership = cell_membership(barcodes, retained)
    totals = technical_totals(column_sums, barcodes, retained)
    kinds = sorted({row[2] for row in features})
    axis_types = {kind: {"rows": sum(row[2] == kind for row in features),
                         "canonical_ordered_ID_name_type_axis_sha256": axis_sha256([row for row in features if row[2] == kind])}
                  for kind in kinds}
    json_write(work_dir / "full-original-feature-axis.json", features)
    with (work_dir / "technical-column-margins.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_sample", "source_column_1based", "source_barcode", "in_fixed_retained_metadata", "source_feature_type", "total_source_units_exact"])
        for col, barcode in enumerate(barcodes):
            for kind, values in sorted(column_sums.items()):
                writer.writerow([sample, col + 1, barcode, int(barcode in retained), kind, str(values[col])])
    exporter = matrix_summary["export_metadata_as_published"]
    count_format = ("Gene Expression" in kinds and matrix_summary["fractional_entries"] == 0
                    and isinstance(exporter, dict) and str(exporter.get("software_version", "")).startswith("cellranger-"))
    result = {"source_archive": archive_path.name, "source_sample": sample,
              "source": {key: receipt[key] for key in ("url", "response_url", "body_bytes", "sha256", "started_at_utc", "completed_at_utc")},
              "source_path": str(archive_path.relative_to(ROOT)) if archive_path.is_relative_to(ROOT) else str(archive_path),
              "source_body_reused_from_closed_first_archive_phase": sample == "sample05",
              "archive_integrity": {"outer_gzip": outer, "members": members,
                                     "tar_headers_paths_types_and_zero_trailer_passed": True,
                                     "source_member_paths_never_used_as_output_paths": True},
              "nested_member_integrity": decoded, "feature_axis": feature_summary,
              "feature_axes_by_type": axis_types, "matrix": matrix_summary,
              "cell_membership": membership, "technical_totals_by_feature_type": totals,
              "released_RNA_source_count_format_supported": count_format,
              "original_UMI_count_likelihood_gate_certified": False}
    json_write(work_dir / "qualification.json", result)
    return result, features, barcodes


def qualify(raw_dir: Path = RAW, work_dir: Path = WORK,
            output_path: Path = ROOT / "data/derived/psc-bach2-cohort-qualification.json") -> dict:
    plan = json.loads((raw_dir / "plan.json").read_text())
    verify_closed_files(plan)
    acquisition = json.loads((raw_dir / "acquisition-summary.json").read_text())
    if not acquisition["complete"] or acquisition["actual_new_body_bytes"] != ARCHIVE_BUDGET or acquisition["additional_metadata_body_bytes"] != 0:
        raise ValueError("Cohort acquisition is incomplete or differs from its closed scope")
    receipts = {item["filename"]: item for item in acquisition["results"]}
    if set(receipts) != set(EXPECTED) or len(acquisition["results"]) != 7:
        raise ValueError("Exactly seven new filtered-archive receipts required")
    for filename, expected in EXPECTED.items():
        if receipts[filename]["body_bytes"] != expected:
            raise ValueError("New source body length differs from its authorized file allowance")
        if json.loads((raw_dir / (filename + ".receipt.json")).read_text()) != receipts[filename]:
            raise ValueError("Individual acquisition receipt differs from the acquisition summary")
    reuse = plan["reuse_archive"]
    reused_path = ROOT / reuse["path"]
    if reused_path.stat().st_size != reuse["bytes"] or count_parser.file_sha256(reused_path) != reuse["sha256"]:
        raise ValueError("Closed sample05 source body changed")
    first_receipt = json.loads((reused_path.parent / "acquisition-receipt.json").read_text())
    closed_first = json.loads((ROOT / "data/derived/psc-bach2-count-qualification.json").read_text())
    if first_receipt != closed_first["source"]:
        raise ValueError("First-archive source receipt differs from its closed qualified source pin")
    receipts[reused_path.name] = first_receipt
    pins, bodies = verified_metadata_sources()
    interfaces = source_interfaces(bodies["E-MTAB-14013.sdrf.txt"], bodies["bach2-e-mtab-14013-metadata"], set(receipts))
    retained, retained_summary = fixed_cohort_metadata(bodies["scrna_meta.tsv"], interfaces)
    json_write(work_dir / "source-sample-donor-genotype-joins.json", interfaces)
    results = {}
    all_features = {}
    source_bare_barcodes = Counter()
    for sample, interface in sorted(interfaces.items()):
        filename = interface["filtered_archive"]
        archive_path = reused_path if filename == reused_path.name else raw_dir / filename
        result, features, barcodes = qualify_one(sample, archive_path, receipts[filename], retained[sample], work_dir / sample)
        results[sample] = result
        all_features[sample] = features
        source_bare_barcodes.update(barcodes)
        print(json.dumps({"sample": sample, "cells": result["cell_membership"], "feature_types": result["feature_axis"]["feature_types_as_published"]}), flush=True)
    rna_axes = {sample: [row for row in rows if row[2] == "Gene Expression"] for sample, rows in all_features.items()}
    axes = compare_rna_axes(rna_axes)
    missing = sum(row["cell_membership"]["retained_cells_missing"] for row in results.values())
    count_format_passed = all(row["released_RNA_source_count_format_supported"] for row in results.values())
    full_interface = axes["all_eight_ordered_RNA_axes_identical"] and not missing and count_format_passed
    types = sorted({kind for row in results.values() for kind in row["technical_totals_by_feature_type"]})
    from fractions import Fraction
    totals = {}
    for kind in types:
        totals[kind] = {key: str(sum(Fraction(row["technical_totals_by_feature_type"].get(kind, {}).get(key, "0")) for row in results.values()))
                        for key in ("all_source_cells_total_units_exact", "fixed_retained_cells_total_units_exact")}
    all_identical = len({axis_sha256(rows) for rows in all_features.values()}) == 1
    result = {
        "scope": "Complete released filtered-cohort RNA-feature source-count-format and retained-cell qualification; no expression analysis",
        "qualified_complete_released_RNA_feature_interface_for_fixed_retained_cells": full_interface,
        "ready_for_expression_analysis": False, "analysis_method_frozen_or_authorized": False,
        "new_archive_body_bytes": acquisition["actual_new_body_bytes"],
        "reused_sample05_archive_bytes": reuse["bytes"], "all_eight_archive_bytes": acquisition["actual_new_body_bytes"] + reuse["bytes"],
        "additional_metadata_body_bytes": 0, "additional_metadata_budget_bytes": 1_000_000,
        "source_pins": pins, "sample_count": len(results),
        "source_genotype_donor_counts_as_published": dict(Counter(row["genotype"] for row in interfaces.values())),
        "unique_accession_local_donors": len({row["individual"] for row in interfaces.values()}),
        "fixed_retained_metadata": retained_summary,
        "source_cells": sum(row["cell_membership"]["source_cells"] for row in results.values()),
        "fixed_retained_cells": retained_summary["rows"],
        "retained_cells_matched": sum(row["cell_membership"]["retained_cells_matched"] for row in results.values()),
        "retained_cells_missing": missing,
        "source_cells_outside_fixed_retained_metadata": sum(row["cell_membership"]["source_cells_outside_fixed_retained_metadata"] for row in results.values()),
        "bare_source_barcode_strings_repeated_across_samples": sum(n > 1 for n in source_bare_barcodes.values()),
        "RNA_axis_comparison": axes, "all_complete_feature_axes_ordered_identical": all_identical,
        "technical_totals_by_feature_type": totals, "antibody_target_annotation": antibody_annotations(all_features),
        "samples": results,
        "source_quantity_limits": {
            "source_format": "Direct Cell Ranger exporter metadata plus original Gene Expression type and exact sparse source values",
            "Methods_IDF_software_as_published": "Cellranger 3.0.2",
            "Methods_IDF_reference_as_published": "GRCh38 human reference genome (release 93)",
            "actual_reference_bundle_or_export_command_recovered": False,
            "historical_RNA_SCT_margin_comparison_available": False,
            "original_UMI_count_likelihood_gate_certified": False,
            "repeated_schema_is_not_historical_export_provenance": True,
            "complete_deposited_RNA_axis_is_not_a_reconstructed_prefilter_reference_universe": True,
            "no_additional_search_for_unavailable_history": True},
        "closed_eight_deliverables_sha256_preserved": {rel: pin["sha256"] for rel, pin in plan["closed_eight_deliverable_pins"].items()},
        "parser_reused_without_edits": {"path": plan["parser_reference"], "sha256": plan["parser_sha256"]},
        "private_outputs": "Full axes, per-cell type-specific margins and explicit sample/donor/genotype joins stay under ignored work/",
        "expression_effects_normalization_scores_models_or_state_gene_selection": False,
        "variant_inference_or_clinical_join": False, "assembly_performed": False,
    }
    verify_closed_files(plan)
    json_write(output_path, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["acquire", "qualify"])
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--output", type=Path, default=ROOT / "data/derived/psc-bach2-cohort-qualification.json")
    args = parser.parse_args()
    if args.command == "acquire":
        print(json.dumps(acquire(workers=args.workers), sort_keys=True))
    elif args.command == "qualify":
        result = qualify(output_path=args.output)
        print(json.dumps({key: result[key] for key in ("sample_count", "source_cells", "fixed_retained_cells", "retained_cells_missing", "qualified_complete_released_RNA_feature_interface_for_fixed_retained_cells")}, sort_keys=True))


if __name__ == "__main__":
    main()
