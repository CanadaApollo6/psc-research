#!/usr/bin/env python3
"""GSE303271 identifier-only program preparation. No expression analysis.

Use the native repository .venv. All source labels and externally fixed source
memberships are retained. This module cannot normalize, match, score or fit.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/psc-liver-program-inputs"
OUTPUT = ROOT / "data/derived/psc-liver-program-inputs.json"
REPORT = ROOT / "reports/psc-liver-program-inputs.md"
MAP_PATH = "data/derived/gse84161-gene-map.csv.gz"
ETS2 = "ENSG00000157557"
COMPARATORS = ("inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling")
INTERFACES = ("exact_symbol_ensembl_bijection", "plus_reciprocal_entrez")
SOURCE_PINS = {'config/ets2-program-plan.json': {'bytes': 10598,
                                   'sha256': '234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9'},
 'config/ets2-validation-amended-rules.json': {'bytes': 11396,
                                               'sha256': 'cbbaa6e626dad9f494367cf1a8c45c0db270f4c330f9bcd61a72280fe1b8ef6b'},
 'config/gse84161-gene-map-sources.json': {'bytes': 7879,
                                           'sha256': 'e1518e126fcd8a61e262ad68afe09e43025c66f2c21d49b26820128ad23c7154'},
 'data/derived/ets2-program-membership.csv': {'bytes': 1669240,
                                              'sha256': 'b6e468f000802e6f1ec6baa53722e8a7d3e0f98d6df4aa2bb49ef886efc199a3'},
 'data/derived/ets2-source-gene-sets.csv': {'bytes': 965247,
                                            'sha256': 'd700811a139cb0dc012f5af8fe22ced452b4d1f76fe5dabb3afebf7ea61e2524'},
 'data/derived/gse84161-gene-map.csv.gz': {'bytes': 1328327,
                                           'sha256': '9db59e5c3ccf1c7f70cbc98aa59e08c0b1168e15ce16ee91c8158e3816cf1aab'},
 'data/derived/psc-liver-input-qualification.json': {'bytes': 78519,
                                                     'sha256': '3c10ec1199725c19a3de37a0e532a66dac15b9e106fc11c9973a209760bfa4de'},
 'data/raw/psc-liver-qualification/GSE303271_raw_counts.txt.gz': {'bytes': 3021812,
                                                                  'sha256': '966be15cf8cd893d0da3ec62d73ee0a56a34932a9d1c990c3d0d58b903a24046'},
 'reports/ets2-program-audit.json': {'bytes': 37753,
                                     'sha256': 'c2589587f3f217d44623c52db9787d8c8d71834fbf1f070799235d42a993d69a'},
 'reports/gse84161-mapping-audit.json': {'bytes': 35244,
                                         'sha256': 'dcd792e85681d49fddbac46a02f29c4a9fe617ae4f7657ad77607ac267b8c83c'},
 'reports/gse84161-normalization-verification.json': {'bytes': 15132,
                                                      'sha256': 'd84273d454da66a9e4c47c517eacd601503e474c3805d2b33d5e3275368edd88'},
 'work/psc-liver-qualification/GSE303271-features-as-submitted.tsv': {'bytes': 554657,
                                                                      'sha256': '19b297c4c860c1b3fd956fba7c0351159aa17c1d866d51426a91908a7fd623ab'}}


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def fingerprint(values):
    return hashlib.sha256(("\n".join(values)+"\n").encode()).hexdigest()


def dump_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name+".partial")
    partial.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False)+"\n")
    partial.replace(path)


def read_rows(path, delimiter=","):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Missing/duplicate table header")
        result = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Ragged table")
            result.append(row)
    return result


def write_rows(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_sources(root=ROOT):
    for relative, pin in SOURCE_PINS.items():
        path = Path(root)/relative
        if path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
            raise ValueError(f"Source pin mismatch: {relative}")
    old_verification = json.loads((Path(root)/"reports/gse84161-normalization-verification.json").read_text())
    # The fixed output checksum is present both in its replay and output inventory.
    if old_verification["artifact_sha256"].get(MAP_PATH) != SOURCE_PINS[MAP_PATH]["sha256"]:
        raise ValueError("Whole-human reference not pinned by historical verification")
    old_audit = json.loads((Path(root)/"reports/ets2-program-audit.json").read_text())
    if old_audit["outputs_sha256"]["data/derived/ets2-program-membership.csv"] != SOURCE_PINS["data/derived/ets2-program-membership.csv"]["sha256"]:
        raise ValueError("Historical complete source-membership checksum mismatch")
    return {"sources": SOURCE_PINS, "network_requests": 0, "new_download_bytes": 0}


def stable_ensembl(value):
    if not re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", value):
        raise ValueError(f"Not a literal Ensembl gene ID: {value!r}")
    return re.sub(r"\.[0-9]+$", "", value)


def build_reference(relationships):
    """Build all directions from the complete reference BEFORE any restriction."""
    valid_ens, ens_symbols, symbol_ens = set(), defaultdict(set), defaultdict(set)
    ens_entrez, entrez_ens = defaultdict(set), defaultdict(set)
    seen, normalized = set(), []
    for row in relationships:
        ens = stable_ensembl(row["ensembl_gene_id"])
        if ens != row["ensembl_gene_id"]:
            raise ValueError("Reference exported normalized Ensembl ID has a version")
        entrez, symbol = row["entrez_gene_id"], row["hgnc_symbol"]
        if entrez and not re.fullmatch(r"[0-9]+", entrez):
            raise ValueError("Nonnumeric reference Entrez ID")
        key = (ens, entrez, symbol)
        if key in seen:
            raise ValueError("Duplicated complete reference relationship")
        seen.add(key)
        valid_ens.add(ens)
        if symbol:
            symbol_ens[symbol].add(ens)
            ens_symbols[ens].add(symbol)
        if entrez:
            ens_entrez[ens].add(entrez)
            entrez_ens[entrez].add(ens)
        normalized.append(key)
    return {"valid_ensembl": valid_ens, "ens_symbols": ens_symbols, "symbol_ens": symbol_ens,
            "ens_entrez": ens_entrez, "entrez_ens": entrez_ens, "relationships": normalized}


def entrez_status(ens, reference):
    ids = reference["ens_entrez"].get(ens, set())
    if not ids:
        return "no_entrez_relation"
    if len(ids) > 1:
        return "multiple_entrez_for_ensembl"
    if len(reference["entrez_ens"][next(iter(ids))]) > 1:
        return "entrez_has_multiple_ensembl"
    return "reciprocal_one_to_one"


def reference_summary(reference):
    ens, symbols = reference["valid_ensembl"], reference["symbol_ens"]
    return {"relationship_rows": len(reference["relationships"]), "unique_ensembl_ids": len(ens),
            "distinct_nonblank_hgnc_symbols": len(symbols),
            "symbols_with_one_ensembl": sum(len(values) == 1 for values in symbols.values()),
            "symbols_with_multiple_ensembl": sum(len(values) > 1 for values in symbols.values()),
            "ensembl_symbol_cardinality": dict(Counter(str(len(reference["ens_symbols"].get(value, set()))) for value in ens)),
            "ensembl_entrez_status_counts": dict(Counter(entrez_status(value, reference) for value in ens)),
            "whole_reference_built_before_program_or_source_feature_restriction": True,
            "reference_release": "Ensembl 116, June 2026, pinned full-human BioMart relationships",
            "biological_locus_or_primary_assembly_allocation_known_from_crosswalk": False}


def load_reference(path, root=ROOT):
    relationships = read_rows(path)
    for index, row in enumerate(relationships, 1):
        if int(row["crosswalk_row_number"]) != index:
            raise ValueError("Reference source row-order gap")
        if stable_ensembl(row["ensembl_gene_id_original"]) != row["ensembl_gene_id"]:
            raise ValueError("Reference Ensembl normalized field does not follow literal numeric-version rule")
        if row["entrez_gene_id_original"].strip() != row["entrez_gene_id"] or row["hgnc_symbol_original"].strip() != row["hgnc_symbol"]:
            raise ValueError("Reference source annotation transformation differs from preserved original fields")
    reference = build_reference(relationships)
    manifest = json.loads((Path(root)/"config/gse84161-gene-map-sources.json").read_text())
    expected = manifest["completeness_check"]
    if len(relationships) != expected["observed_relationship_rows"] or len(reference["valid_ensembl"]) != expected["observed_unique_ensembl_gene_ids"]:
        raise ValueError("Full-human reference export is incomplete")
    for row in relationships:
        ens, entrez, symbol = row["ensembl_gene_id"], row["entrez_gene_id"], row["hgnc_symbol"]
        if int(row["ensembl_entrez_distinct_count"]) != len(reference["ens_entrez"].get(ens, set())):
            raise ValueError("Reference Ensembl cardinality audit mismatch")
        if entrez and int(row["entrez_ensembl_distinct_count"]) != len(reference["entrez_ens"][entrez]):
            raise ValueError("Reference reciprocal cardinality audit mismatch")
        if symbol and int(row["hgnc_symbol_distinct_ensembl_count"]) != len(reference["symbol_ens"][symbol]):
            raise ValueError("Reference symbol cardinality audit mismatch")
    return reference


def load_submitted_labels(root=ROOT):
    """Read header and first field only. Numerical count fields are not parsed."""
    labels = []
    count_path = Path(root)/"data/raw/psc-liver-qualification/GSE303271_raw_counts.txt.gz"
    with gzip.open(count_path, "rt", encoding="utf-8") as stream:
        header = next(stream).rstrip("\n").split("\t")
        if header[0] != "":
            raise ValueError("Unexpected count header")
        for line in stream:
            label, separator, _unparsed_count_fields = line.partition("\t")
            if not separator or not label or label.strip() != label:
                raise ValueError("Missing or whitespace-altered source feature label")
            labels.append(label)
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate submitted labels; no collapsing permitted")
    retained = read_rows(Path(root)/"work/psc-liver-qualification/GSE303271-features-as-submitted.tsv", delimiter="\t")
    if [int(row["row_index_0based"]) for row in retained] != list(range(len(retained))):
        raise ValueError("Retained source feature axis has gaps")
    if labels != [row["feature_label_as_submitted"] for row in retained]:
        raise ValueError("Submitted count labels differ from closed full-file qualification")
    qualification = json.loads((Path(root)/"data/derived/psc-liver-input-qualification.json").read_text())
    qc = qualification["datasets"]["GSE303271"]["numeric_qc"]
    if fingerprint(labels) != qc["feature_order_sha256"] or fingerprint(header[1:]) != qc["column_order_sha256"]:
        raise ValueError("Closed count axes differ")
    return labels


def map_features(labels, reference):
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate labels cannot be collapsed")
    rows, candidates_by_ens = [], defaultdict(list)
    for index, label in enumerate(labels):
        candidate_ids = reference["symbol_ens"].get(label, set())
        unique_ens = next(iter(candidate_ids)) if len(candidate_ids) == 1 else ""
        if not candidate_ids:
            core = "no_exact_reference_symbol"
        elif len(candidate_ids) > 1:
            core = "symbol_has_multiple_ensembl"
        elif len(reference["ens_symbols"][unique_ens]) != 1:
            core = "ensembl_has_multiple_symbols"
        else:
            core = "eligible"
        if unique_ens:
            candidates_by_ens[unique_ens].append(index)
        status_entrez = entrez_status(unique_ens, reference) if unique_ens else "not_resolved_to_one_ensembl"
        strict = core if core != "eligible" else "eligible" if status_entrez == "reciprocal_one_to_one" else status_entrez
        rows.append({"source_row_index_0based": index, "label_as_submitted": label,
                     "candidate_ensembl_ids": sorted(candidate_ids),
                     "candidate_entrez_ids": sorted(set().union(*(reference["ens_entrez"].get(value, set()) for value in candidate_ids))) if candidate_ids else [],
                     "reference_symbols_for_unique_ensembl": sorted(reference["ens_symbols"].get(unique_ens, set())),
                     "ensembl_id_if_unique_symbol_resolution": unique_ens,
                     "entrez_relationship_status": status_entrez,
                     INTERFACES[0]: core, INTERFACES[1]: strict})
    for ens, indices in candidates_by_ens.items():
        if len(indices) > 1:
            for index in indices:
                rows[index]["source_feature_collision_for_ensembl"] = True
                for interface in INTERFACES:
                    if rows[index][interface] == "eligible":
                        rows[index][interface] = "multiple_submitted_features_for_ensembl"
    for row in rows:
        row.setdefault("source_feature_collision_for_ensembl", False)
    return rows


def annotation_only(root=ROOT, output_directory=WORK):
    source_check = verify_sources(root)
    labels = load_submitted_labels(root)
    reference = load_reference(Path(root)/MAP_PATH, root)
    feature_rows = map_features(labels, reference)
    summary = {"source_check": source_check, "reference": reference_summary(reference),
               "released_source_features": len(labels), "feature_label_order_sha256": fingerprint(labels),
               "literal_Ensembl_ID_rows": sum(bool(re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", label)) for label in labels),
               "gene_name_style_rows": sum(not re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", label) for label in labels),
               "feature_status_counts": {interface: dict(Counter(row[interface] for row in feature_rows)) for interface in INTERFACES},
               "no_numeric_counts_parsed": True, "no_alias_or_case_folding": True}
    dump_json(Path(output_directory)/"annotation-summary.json", summary)
    dump_json(Path(output_directory)/"feature-identifier-audit.json", feature_rows)
    return summary, reference, feature_rows



def selected_identifier_rows(path, fields):
    """Retain only declared identifier/provenance fields from an existing audit."""
    result = []
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if not set(fields).issubset(reader.fieldnames or []):
            raise ValueError("Required identifier columns absent")
        for row in reader:
            result.append({field: row[field] for field in fields})
    return result


def load_original_programs(root=ROOT):
    root = Path(root)
    plan = json.loads((root/"config/ets2-program-plan.json").read_text())
    amended = json.loads((root/"config/ets2-validation-amended-rules.json").read_text())
    definitions = [row for row in plan["programs"] if "source_member" in row]
    amended_definitions = [row for row in amended["programs"] if "source_member" in row]
    if definitions != amended_definitions or len(definitions) != 8:
        raise ValueError("External fixed source-program definitions changed")
    fields = ["dataset_key", "program_id", "source_index", "source_value"]
    membership = selected_identifier_rows(root/"data/derived/ets2-program-membership.csv", fields)
    # The other columns of this legacy source-gene-set table include published
    # effect annotations. They are never retained, interpreted, or used here.
    ensembl_rows = selected_identifier_rows(root/"data/derived/ets2-source-gene-sets.csv",
                                           ["set_name", "gene_id_original", "source_member", "source_member_bytes", "source_member_sha256"])
    old_audit = json.loads((root/"reports/gse84161-mapping-audit.json").read_text())
    old_definitions = {row["program_id"]: row for row in old_audit["source_programs"]["definitions"]}
    programs = {}
    for definition in definitions:
        program = definition["id"]
        copies = []
        for modality in ["sc", "sn"]:
            rows = sorted((row for row in membership if row["dataset_key"] == modality and row["program_id"] == program),
                          key=lambda row: int(row["source_index"]))
            if [int(row["source_index"]) for row in rows] != list(range(len(rows))):
                raise ValueError("Missing/repeated source-membership index")
            copies.append([row["source_value"] for row in rows])
        if not copies[0] or copies[0] != copies[1] or len(copies[0]) != len(set(copies[0])):
            raise ValueError("Complete original source lists disagree or have duplicate members")
        namespace = old_audit["source_programs"]["member_fingerprints"][program]["source_namespace"]
        if namespace == "Ensembl_gene_id":
            alternate = [row for row in ensembl_rows if row["set_name"] == definition["column"]]
            if [row["gene_id_original"] for row in alternate] != copies[0]:
                raise ValueError("Complete original Ensembl source-set exports disagree")
            if any(row["source_member"] != definition["source_member"] for row in alternate):
                raise ValueError("Source-member archive path mismatch")
            normalized = [stable_ensembl(value) for value in copies[0]]
        elif namespace == "HGNC_symbol":
            normalized = list(copies[0])
        else:
            raise ValueError("Unsupported original source namespace")
        if len(normalized) != len(set(normalized)):
            raise ValueError("Source members collide after numeric version removal; no silent collapse")
        intended = sum(value not in [ETS2, "ETS2"] for value in normalized)
        expected = old_definitions[program]
        if len(copies[0]) != expected["source_member_count_raw"] or intended != expected["intended_member_count_after_ets2_exclusion"]:
            raise ValueError("Original program denominator changed")
        programs[program] = {"definition": definition, "source_namespace": namespace,
                             "source_values": copies[0], "normalized_values": normalized,
                             "intended_members": intended,
                             "source_member_fingerprint": old_audit["source_programs"]["member_fingerprints"][program]}
    if plan["minimum_mapped_fraction"] != amended["minimum_mapped_fraction"] or plan["minimum_mapped_genes"] != amended["minimum_mapped_genes"]:
        raise ValueError("External mapping gates disagree")
    return programs, {"minimum_mapped_fraction": amended["minimum_mapped_fraction"], "minimum_mapped_genes": amended["minimum_mapped_genes"]}


def resolve_source_member(value, namespace, reference):
    normalized = stable_ensembl(value) if namespace == "Ensembl_gene_id" else value
    if normalized in (ETS2, "ETS2"):
        return normalized, [], "excluded_ETS2"
    if namespace == "Ensembl_gene_id":
        candidates = [normalized] if normalized in reference["valid_ensembl"] else []
        absent = "ensembl_id_absent_from_reference_snapshot"
    elif namespace == "HGNC_symbol":
        candidates = sorted(reference["symbol_ens"].get(normalized, set()))
        absent = "no_exact_reference_symbol"
    else:
        raise ValueError("Unsupported source namespace")
    if not candidates:
        status = absent
    elif len(candidates) > 1:
        status = "symbol_has_multiple_ensembl"
    else:
        ens = candidates[0]
        symbols = reference["ens_symbols"].get(ens, set())
        if not symbols:
            status = "ensembl_present_without_hgnc_symbol"
        elif len(symbols) > 1:
            status = "ensembl_has_multiple_symbols"
        elif len(reference["symbol_ens"][next(iter(symbols))]) != 1:
            status = "symbol_has_multiple_ensembl"
        else:
            status = "resolved_exact_symbol_ensembl_bijection"
    return normalized, candidates, status


def map_program_members(programs, gates, reference, features):
    features_by_ens = defaultdict(list)
    labels = {row["label_as_submitted"]: row["source_row_index_0based"] for row in features}
    for row in features:
        if row["ensembl_id_if_unique_symbol_resolution"]:
            features_by_ens[row["ensembl_id_if_unique_symbol_resolution"]].append(row)
    mapping, inventory = [], []
    mapped_sets = {interface: {} for interface in INTERFACES}
    for program, program_data in programs.items():
        rows = []
        for source_index, source_value in enumerate(program_data["source_values"]):
            normalized, candidates, resolution = resolve_source_member(source_value, program_data["source_namespace"], reference)
            candidate_symbols = sorted(set().union(*(reference["ens_symbols"].get(value, set()) for value in candidates))) if candidates else []
            candidate_entrez = sorted(set().union(*(reference["ens_entrez"].get(value, set()) for value in candidates))) if candidates else []
            audit = {"program_id": program, "source_index_0based": source_index, "source_value": source_value,
                     "normalized_source_value": normalized, "source_namespace": program_data["source_namespace"],
                     "source_reference_status": resolution, "candidate_ensembl_ids": candidates,
                     "candidate_current_reference_symbols": candidate_symbols, "candidate_entrez_ids": candidate_entrez,
                     "reference_symbols_present_as_submitted": [value for value in candidate_symbols if value in labels],
                     "literal_source_symbol_present_as_submitted": source_value in labels if program_data["source_namespace"] == "HGNC_symbol" else None,
                     "source_locus_assignment_reconstructed": False}
            for interface in INTERFACES:
                status, feature_index, ens_id = resolution, None, None
                if resolution == "resolved_exact_symbol_ensembl_bijection":
                    ens = candidates[0]
                    possible = features_by_ens.get(ens, [])
                    if not possible:
                        status = "exact_reference_symbol_absent_from_submitted_labels"
                    elif len(possible) > 1:
                        status = "multiple_submitted_features_for_ensembl"
                    elif possible[0][interface] != "eligible":
                        status = possible[0][interface]
                    else:
                        status = "mapped"
                        feature_index = possible[0]["source_row_index_0based"]
                        ens_id = ens
                audit[interface] = {"status": status, "source_feature_index_0based": feature_index, "ensembl_gene_id": ens_id}
            rows.append(audit)
        mapping.extend(rows)
        for interface in INTERFACES:
            mapped = [row[interface]["ensembl_gene_id"] for row in rows if row[interface]["status"] == "mapped"]
            if len(mapped) != len(set(mapped)):
                raise ValueError("Original distinct source members collapse to one mapped gene; no automatic merge")
            mapped_sets[interface][program] = set(mapped)
            denominator = program_data["intended_members"]
            mapped_count = len(mapped)
            inventory.append({"interface": interface, "program_id": program, "source_namespace": program_data["source_namespace"],
                              "raw_original_source_members": len(rows), "intended_source_members_after_ETS2_exclusion": denominator,
                              "excluded_ETS2_members": sum(row[interface]["status"] == "excluded_ETS2" for row in rows),
                              "mapped_source_members": mapped_count, "mapped_unique_ensembl_ids": len(set(mapped)),
                              "not_mapped_intended_members": denominator-mapped_count,
                              "mapped_fraction_of_original_intended_members": mapped_count/denominator,
                              "mapping_gate_pass": mapped_count >= gates["minimum_mapped_genes"] and mapped_count/denominator >= gates["minimum_mapped_fraction"],
                              "gate_basis": "original intended source-member denominator; annotation mapping only",
                              "status_counts": dict(Counter(row[interface]["status"] for row in rows)),
                              "mapped_ensembl_ids_sha256": fingerprint(sorted(mapped)), "effects_authorized": False})
    for interface in INTERFACES:
        parent = mapped_sets[interface]["ets2_g1_dn"]
        comparators = set().union(*(mapped_sets[interface][program] for program in COMPARATORS))
        retained = parent-comparators
        mapped_sets[interface]["ets2_g1_dn_without_comparators"] = retained
        parent_row = next(row for row in inventory if row["interface"] == interface and row["program_id"] == "ets2_g1_dn")
        inventory.append({"interface": interface, "program_id": "ets2_g1_dn_without_comparators", "derived_from": "ets2_g1_dn",
                          "intended_source_members_after_ETS2_exclusion": programs["ets2_g1_dn"]["intended_members"],
                          "mapped_unique_ensembl_ids": len(retained), "parent_mapped_unique_ensembl_ids": len(parent),
                          "mapped_comparator_union_ensembl_ids": len(comparators), "overlap_removed": len(parent & comparators),
                          "retained_fraction_of_original_intended_parent": len(retained)/programs["ets2_g1_dn"]["intended_members"],
                          "parent_mapping_fraction_used_for_gate": parent_row["mapped_fraction_of_original_intended_members"],
                          "mapping_gate_pass": parent_row["mapping_gate_pass"] and len(retained) >= gates["minimum_mapped_genes"],
                          "gate_basis": "unchanged parent coverage gate plus separate ten-gene minimum; no new coverage denominator",
                          "mapped_ensembl_ids_sha256": fingerprint(sorted(retained)), "effects_authorized": False})
    return mapping, inventory, mapped_sets


def candidate_pools(features, mapped_sets):
    summaries, rows = {}, []
    for interface in INTERFACES:
        union = set().union(*(genes for program, genes in mapped_sets[interface].items() if program != "ets2_g1_dn_without_comparators"))
        eligible = [row for row in features if row[interface] == "eligible"]
        counts = Counter()
        background_genes, background_labels = [], []
        for row in features:
            ens = row["ensembl_id_if_unique_symbol_resolution"]
            if row[interface] != "eligible":
                status = "identity_not_qualified"
            elif ens == ETS2 or row["label_as_submitted"] == "ETS2":
                status = "excluded_ETS2"
            elif ens in union:
                status = "excluded_original_program_union"
            else:
                status = "eligible_background_candidate"
                background_genes.append(ens)
                background_labels.append(row["label_as_submitted"])
            counts[status] += 1
            rows.append({"interface": interface, "source_row_index_0based": row["source_row_index_0based"],
                         "label_as_submitted": row["label_as_submitted"], "ensembl_gene_id_if_resolved": ens,
                         "identity_status": row[interface], "pool_status": status,
                         "source_feature_row_retained_for_original_count_universe": True})
        if len(background_genes) != len(set(background_genes)):
            raise ValueError("Candidate pool repeats an Ensembl gene; no hidden collapse")
        summaries[interface] = {"released_feature_axis_rows": len(features), "ID_qualified_feature_rows": len(eligible),
                                "pool_status_counts": dict(counts), "eight_original_mapped_program_union_genes": len(union),
                                "background_candidate_genes": len(background_genes),
                                "background_ensembl_ids_source_order_sha256": fingerprint(background_genes),
                                "background_source_labels_source_order_sha256": fingerprint(background_labels),
                                "includes_ID_eligible_rows_regardless_of_count_or_detection": True,
                                "abundance_bins_or_draws_created": False, "matching_pool_sufficiency_observed": False,
                                "counts_normalization_denominator_not_filtered_by_this_pool": True}
    return summaries, rows


def prepare(root=ROOT, output_directory=WORK):
    root, output_directory = Path(root), Path(output_directory)
    annotation, reference, features = annotation_only(root, output_directory)
    programs, gates = load_original_programs(root)
    mapping, inventory, mapped_sets = map_program_members(programs, gates, reference, features)
    pool_summary, pool_rows = candidate_pools(features, mapped_sets)
    dump_json(output_directory/"complete-original-source-programs.json", programs)
    dump_json(output_directory/"program-member-identifier-audit.json", mapping)
    dump_json(output_directory/"program-coverage.json", inventory)
    write_rows(output_directory/"complete-reference-pool-audit.csv", pool_rows, list(pool_rows[0]))
    program_feature_rows = []
    by_ens = {row["ensembl_id_if_unique_symbol_resolution"]: row for row in features if row[INTERFACES[0]] == "eligible"}
    for interface in INTERFACES:
        for program, genes in mapped_sets[interface].items():
            for ens in sorted(genes):
                feature = by_ens[ens]
                program_feature_rows.append({"interface": interface, "program_id": program, "ensembl_gene_id": ens,
                                             "source_feature_index_0based": feature["source_row_index_0based"], "label_as_submitted": feature["label_as_submitted"]})
    write_rows(output_directory/"qualified-program-feature-indices.csv", program_feature_rows, list(program_feature_rows[0]))
    summary = {"schema_version": 1, "status": "identifier_only_preparation_no_analysis_authorization",
               "annotation": annotation, "mapping_gates": gates, "program_coverage": inventory,
               "candidate_matching_reference_pools": pool_summary,
               "scope": {"no_numeric_counts_parsed": True, "normalization_performed": False, "expression_matching_performed": False,
                         "scores_or_effects_computed": False, "network_requests": 0, "aliases_guessed": False,
                         "ambiguous_IDs_collapsed": False, "missing_genes_zero_filled": False,
                         "atlas_mapped_program_subsets_reused": False, "organoid_derived_genes_used": False,
                         "historical_counting_GTF_reconstructed": False, "root_method_and_analysis_freeze_required": True}}
    dump_json(output_directory/"preparation-summary.json", summary)
    return summary



PREPARATION_ID = "gse303271-external-fixed-programs-identifier-preparation-v1"
SOURCE_AUDIT_FINDINGS_SHA256 = "a11fae6df4a3fb42ce547d95605be4782ba5dab3d17317b8e4762d29f5eaf1fb"
PREPARATION_ARTIFACTS = (
    "annotation-summary.json", "feature-identifier-audit.json", "complete-original-source-programs.json",
    "program-member-identifier-audit.json", "program-coverage.json", "complete-reference-pool-audit.csv",
    "qualified-program-feature-indices.csv", "preparation-summary.json",
)


def check_pin(root, pin):
    path = Path(root)/pin["path"]
    if path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
        raise ValueError(f"Pinned material changed: {pin['path']}")


def pin_record(path, root=ROOT):
    path = Path(path)
    return {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def verify_closed_materials(root=ROOT):
    root = Path(root)
    scope = json.loads((root/"work/psc-liver-program-inputs/scope.json").read_text())
    closed_liver = json.loads((root/"data/derived/psc-liver-input-qualification.json").read_text())
    old_scope = json.loads((root/"work/psc-liver-qualification/scope.json").read_text())
    versioned = {**scope["closed_liver_versioned_pins"], **scope["closed_organoid_versioned_pins"]}
    for relative, digest in versioned.items():
        if sha256(root/relative) != digest:
            raise ValueError(f"Closed versioned result changed: {relative}")
    groups = {"closed_liver_source_pins": closed_liver["provenance"]["input_pins"],
              "closed_liver_ignored_artifact_pins": closed_liver["ignored_artifact_pins"],
              "closed_organoid_source_and_artifact_pins": old_scope["organoid_prior_material_pins"]}
    for group in groups.values():
        for pin in group:
            check_pin(root, pin)
    return {"all_checked_pins_unchanged": True, "closed_liver_versioned_files": len(scope["closed_liver_versioned_pins"]),
            "closed_organoid_versioned_files": len(scope["closed_organoid_versioned_pins"]),
            **{key: len(pins) for key, pins in groups.items()}, "closed_versioned_pins": versioned,
            "no_closed_tree_writes_by_this_script": True}


def source_provenance(root=ROOT):
    root = Path(root)
    path = root/"work/psc-liver-program-inputs/source-audit/findings.json"
    if sha256(path) != SOURCE_AUDIT_FINDINGS_SHA256:
        raise ValueError("Independent original-source audit changed")
    findings = json.loads(path.read_text())
    if not findings["checks"]["all_passed"] or findings["checks"]["failed"]:
        raise ValueError("Original-source audit did not pass")
    for pin in findings["input_receipts"].values():
        # Parent-owned context documents/readiness were background reading, not
        # identifier inputs or a lock on unrelated parent progress. Keep their
        # audited snapshots, but do not freeze live contextual reports here.
        if pin["path"] not in {"docs/evidence-log.md", "reports/macromap-program-design.md",
                               "work/macromap-design/preparation/source-only-readiness.json"}:
            check_pin(root, pin)
    complete = findings["full_human_relation_export"]
    for key in ["source_response_reconstruction", "count_response_reconstruction"]:
        reconstruction = complete[key]
        if not reconstruction["exact_pin_match"]:
            raise ValueError("Retained original reference fields do not reconstruct pinned historical bytes")
        check_pin(root, reconstruction["reconstructed_output"])
    check_pin(root, findings["outputs"]["source_membership_audit"])
    source_programs = findings["source_membership_recovery"]
    fields = ["program_id", "role", "source_namespace", "source_member", "source_column",
              "source_unique_before_exclusion", "intended_source_denominator", "excluded_target_gene_count",
              "ordered_original_values_sha256_LF_terminated", "ordered_intended_values_sha256_LF_terminated"]
    program_inventory = [{field: row[field] for field in fields} for row in source_programs["per_program"]]
    count_qualification = json.loads((root/"data/derived/psc-liver-input-qualification.json").read_text())
    count_source = next(row for row in count_qualification["provenance"]["main_sources"] if row.get("file") == "GSE303271_raw_counts.txt.gz")
    verification_path = root/"work/psc-liver-program-inputs/source-audit/verification.json"
    verification = json.loads(verification_path.read_text())
    if verification["native_audit_and_independent_tests_exit_code"] != 0 or not verification["replay_byte_identical"]:
        raise ValueError("Source audit validation or offline replay failed")
    return {"fixed_source_programs": program_inventory, "source_audit": pin_record(path, root),
            "source_audit_verification": pin_record(verification_path, root),
            "source_audit_checks_passed": findings["checks"]["total_individual_checks"],
            "source_audit_tests_passed": verification["independent_tests_passed"],
            "complete_original_source_memberships_including_ETS2": source_programs["full_original_membership_rows"],
            "complete_intended_memberships_after_ETS2_exclusion": source_programs["intended_membership_rows_after_ETS2_exclusion"],
            "original_ETS2_archive": source_programs["original_archive"],
            "original_archive_member_pins": source_programs["original_source_member_pins"],
            "complete_original_memberships_recovered_from_verified_exports_not_old_mapped_fields": True,
            "SC_SN_copies_agree_but_are_not_independent_biological_sources": True,
            "reference_historical_raw_response": complete["historical_raw_response"],
            "reference_historical_count_response": complete["historical_count_response"],
            "reference_raw_response_reconstruction": complete["source_response_reconstruction"],
            "reference_count_response_reconstruction": complete["count_response_reconstruction"],
            "reference_query": complete["query"], "reference_completeness_limit": complete["completeness_limit"],
            "reference_genome_build_status": complete["source_genome_build"],
            "count_file_previous_acquisition": count_source,
            "no_fresh_network_requests": True,
            "context_only_audit_snapshots_not_identifier_inputs": ["docs/evidence-log.md", "reports/macromap-program-design.md", "work/macromap-design/preparation/source-only-readiness.json"],
            "context_snapshot_boundary": "Historical read receipts remain in source-audit/findings.json; live contextual documents are not inputs to this mapping and do not freeze unrelated parent updates.",
            "legacy_source_statistic_preview_disclosure": "Initial schema previews in this task and its source audit also displayed existing legacy source-statistic strings. They were not converted, summarized, selected on or used for new conclusions. All preparation/audit code uses identifier and provenance fields only. No new liver expression quantities were inspected."}


def public_payload(summary, root=ROOT, output_directory=WORK):
    root, output_directory = Path(root), Path(output_directory)
    scope = json.loads((root/"work/psc-liver-program-inputs/scope.json").read_text())
    provenance = source_provenance(root)
    validation_path = root/"work/psc-liver-program-inputs/validation-run.json"
    validation = json.loads(validation_path.read_text()) if validation_path.exists() else {"status": "final_validation_not_yet_recorded"}
    return {"schema_version": 1, "preparation_id": PREPARATION_ID,
            "status": "identifier_only_preparation_root_analysis_freeze_required", "source_study": "GSE303271",
            "scope": summary["scope"], "prospective_identifier_interfaces": scope["prospective_annotation_interfaces"],
            "policy_record": pin_record(root/"work/psc-liver-program-inputs/scope.json", root),
            "source_feature_axis": {"released_source_rows": summary["annotation"]["released_source_features"],
                                    "feature_label_order_sha256": summary["annotation"]["feature_label_order_sha256"],
                                    "gene_name_style_rows": summary["annotation"]["gene_name_style_rows"],
                                    "literal_Ensembl_ID_rows": summary["annotation"]["literal_Ensembl_ID_rows"],
                                    "all_original_source_rows_preserved": True,
                                    "ID_eligibility_not_a_count_normalization_universe": True,
                                    "complete_raw_count_values_not_parsed_in_this_task": True},
            "source_provenance": provenance, "source_input_pins": SOURCE_PINS,
            "whole_human_reference": summary["annotation"]["reference"],
            "feature_identifier_status_counts": summary["annotation"]["feature_status_counts"],
            "mapping_gates": {**summary["mapping_gates"], "source_denominators_not_redefined": True,
                              "nonoverlap_requires_parent_coverage_gate_and_separate_minimum_ten_retained_genes": True,
                              "mapping_gate_pass_is_not_expression_or_matching_readiness": True},
            "program_coverage": summary["program_coverage"],
            "candidate_matching_reference_pools": summary["candidate_matching_reference_pools"],
            "transformations": ["Verify fixed source byte pins before any preparation.",
                                "Read only raw gzip header and first field; verify exact complete feature axis against the closed qualification.",
                                "Build both directions of every identifier relationship from all 91,748 human reference rows before any program or measured-feature restriction.",
                                "Preserve all original program source_value strings and SC/SN source copies; strip only terminal numeric versions from valid source Ensembl IDs, if any.",
                                "Exclude the original fixed ETS2 Ensembl ID or literal ETS2 symbol, preserving excluded membership records.",
                                "Require exact case-sensitive symbol/Ensembl bijection and one source row; report the additional reciprocal Entrez interface separately.",
                                "Retain missing/ambiguous mappings with null mapped feature indexes, never zero-fill or choose one ambiguous locus.",
                                "Subtract the four mapped generic-comparator unions from the mapped primary under each interface, preserving the parent denominator and coverage gate.",
                                "Enumerate every source row in each ID-only pool audit; exclude all eight original mapped programs and ETS2 from candidate backgrounds without abundance/detection/protein-coding filters."],
            "inference_boundaries": {
                "root_identity_policy_choice_required_before_matching_or_effects": True,
                "root_normalization_matching_and_analysis_freeze_not_written_here": True,
                "historical_counting_annotation": "Paper declares GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4; exact GTF/provider/release/checksum and count-assignment settings remain unavailable.",
                "current_reference_not_historical_GTF_or_biological_locus_reconstruction": True,
                "reference_has_no_coordinates_biotypes_or_assembly_region_fields": True,
                "multiple_Ensembl_relations_do_not_prove_independent_biological_loci": True,
                "sample_design_retained_from_closed_qualification": {"PSC": 17, "ASC": 17, "AIH": 30, "tissue": "whole-biopsy bulk liver, not sorted macrophages or cholangiocytes"},
                "donor_and_covariate_limits": "64 patients are asserted in the article. Independent person/visit crosswalk and sample-linked age/sex/treatment/stage/batch covariates remain unavailable. ASC remains separate from PSC.",
                "GSE159676": "HOLD expression use; unreconciled mixed-scale integer literals and six paired PSC patients, not twelve independent patients. No adult expression use or repair here.",
                "NoPSC_atlas": "PARK; no qualified full-count/full-feature/barcode/donor export. No atlas-derived gene subset used here.",
                "no_specificity_causality_transfer_replication_or_clinical_claim": True},
            "closed_material_preservation": verify_closed_materials(root),
            "ignored_artifact_pins": [pin_record(output_directory/name, root) for name in PREPARATION_ARTIFACTS],
            "software": {"python": platform.python_version(), "native_command": ".venv/bin/python",
                         "script": "scripts/prepare_psc_liver_program_inputs.py", "script_sha256": sha256(root/"scripts/prepare_psc_liver_program_inputs.py"),
                         "tests": "tests/test_psc_liver_program_inputs.py", "tests_sha256": sha256(root/"tests/test_psc_liver_program_inputs.py"),
                         "new_dependencies_installed": False, "inference_model_version": "not applicable; no model used"},
            "validation": validation}


def report_text(data):
    core, strict = INTERFACES
    features = data["feature_identifier_status_counts"]
    pools = data["candidate_matching_reference_pools"]
    source = data["source_provenance"]
    inventory = {(row["interface"], row["program_id"]): row for row in data["program_coverage"]}
    lines = ["# GSE303271 fixed-program identifier preparation", "", "## Status and boundary", "",
             "Completed offline **identifier-only preparation**, not an analysis freeze. Root must choose the identity, normalization and matching rules before effects. No normalization, expression matching, scoring, effects or models were run. No network requests, aliases, ambiguity collapse, zero-filling or organoid-derived genes were used.", "",
             "All **41,096 released RNA rows** remain in original order. Candidate identifier pools are **not normalization universes**. No abundance, detection or protein-coding filter was applied. Detailed labels, source memberships, ambiguous candidates and indexes remain under `work/psc-liver-program-inputs/`; this report and the public JSON contain aggregates.", "",
             "## Source lineage", "",
             "The complete original `source_value` memberships were recovered from the pinned `ets2-program-membership.csv`, including previously unmapped and ambiguous atlas members. SC/SN copies agree in order; they are duplicate exports, not independent biological sources. The four original Ensembl lists also agree with the complete source-gene-set export. Old mapped fields, MacroMap/atlas subsets and legacy source-statistic columns were not used to define programs.", "",
             f"There are **{source['complete_original_source_memberships_including_ETS2']:,} original memberships** across eight programs and **{source['complete_intended_memberships_after_ETS2_exclusion']:,} intended memberships** after the single fixed ETS2 exclusion. The original ETS2 ZIP remains absent; this is verified-export reuse, not a new archive read.", "",
             "The full Ensembl 116 / June 2026 reference has **91,748 relationships and 86,411 Ensembl IDs**. All relationships are used before any feature/program restriction, including 54,919 Ensembl IDs without Entrez. The reference gzip SHA-256 is `9db59e5c3ccf1c7f70cbc98aa59e08c0b1168e15ce16ee91c8158e3816cf1aab`; the historical normalization verification independently pins it.", "",
             "The absent historical BioMart response was reconstructed exactly from preserved original fields/order, its pinned header, tabs/LF and terminal `[success]`. The owned reconstructed copy is **2,236,970 bytes**, SHA-256 `a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167`. This is not a fresh server retrieval or a reread/restoration of the absent original cache. The count response also reconstructs to its historical pin. Success/count agreement supports the returned export but cannot rule out an unobserved remote backend omission.", "",
             "Historical source URLs and retrieval times:", "",
             f"- RNA count file: {source['count_file_previous_acquisition']['url']} — {source['count_file_previous_acquisition']['retrieved_at_utc']}.",
             f"- BioMart: {source['reference_historical_raw_response']['url']} — {source['reference_historical_raw_response']['retrieved_at_utc']}; independent count response {source['reference_historical_count_response']['retrieved_at_utc']}.",
             f"- ETS2 archive: {source['original_ETS2_archive']['historical_source_record']['url']} — {source['original_ETS2_archive']['historical_source_record']['retrieved_at_utc']} (historical manifest only; ZIP absent).", "",
             "Full source/member/export pins, query XML and reconstruction receipts are in the JSON and `source-audit/findings.json`. Source previews incidentally displayed existing legacy statistic strings; these were not converted, summarized, selected on or used for conclusions. Preparation/audit code retains only identifier/provenance fields. No new liver expression quantities were inspected.", "",
             "## Two prospective identity interfaces", "",
             "Both policies were recorded in `work/psc-liver-program-inputs/scope.json` before mapping. Neither was selected by these coverage results:", "",
             "1. **Core** (`exact_symbol_ensembl_bijection`): exact case-sensitive submitted symbol ↔ one current-reference Ensembl ID, plus one submitted row for that ID. Entrez is reported, not required as an RNA bridge.",
             "2. **Strict** (`plus_reciprocal_entrez`): additionally require globally reciprocal one-to-one Ensembl ↔ Entrez relations. This is a separate optional cross-namespace interface, not the chosen analysis policy.", "",
             "All released labels are gene-name-style; none are literal Ensembl IDs. Ambiguity checks use the complete human reference, not just measured genes or programs. The reference has no coordinates, biotypes or assembly-region fields. Multiple Ensembl links do not by themselves prove separate biological loci; no primary-locus choice was made.", "",
             "The table shows mutually exclusive eligibility statuses, not all annotation relationships. Entrez conditions do not exclude core rows; the core zeros below do not mean those relations are absent.", "",
             "| Feature status | Core | Strict |", "|---|---:|---:|"]
    for status in ["eligible", "no_exact_reference_symbol", "symbol_has_multiple_ensembl", "ensembl_has_multiple_symbols", "no_entrez_relation", "multiple_entrez_for_ensembl", "entrez_has_multiple_ensembl"]:
        lines.append(f"| `{status}` | {features[core].get(status, 0):,} | {features[strict].get(status, 0):,} |")
    lines += ["", "Missing Entrez is not absence of an RNA row or biological absence. Missing current symbols are not automatically retired genes. Every unresolved relationship remains explicit; no alias, zero substitute or arbitrary ambiguous ID was used.", "",
              "## Original denominators and coverage", "",
              "The unchanged mapping gate is at least **80% of the original intended members and ten mapped genes**. All eight programs pass that annotation gate under both interfaces. This does not establish matching-pool sufficiency or authorize scoring.", "",
              "| Fixed program | Original intended | Core mapped | Strict mapped |", "|---|---:|---:|---:|"]
    for program in source["fixed_source_programs"]:
        name = program["program_id"]
        first, second = inventory[(core, name)], inventory[(strict, name)]
        lines.append(f"| `{name}` | {program['intended_source_denominator']} | {first['mapped_unique_ensembl_ids']} ({first['mapped_fraction_of_original_intended_members']:.2%}) | {second['mapped_unique_ensembl_ids']} ({second['mapped_fraction_of_original_intended_members']:.2%}) |")
    lines += ["", "The strict opposite-direction program has 536/668 members (80.24%); its original denominator is not reduced. Per-program missing-symbol, missing-current-Ensembl, absent-row, multiple-symbol and Entrez failure counts are separate in the JSON and complete ignored audit.", "",
              "**Nonoverlap sensitivity:** subtract the mapped union of inflammation, interferon-gamma, oxidative-stress and apoptotic-signaling programs. Core retains **695** of 848 mapped primary genes after removing **153**; strict retains **673** of 825 after removing **152**. The gate remains the parent's 927-member coverage gate plus the separate ten-retained-gene minimum. The retained fractions 695/927 and 673/927 do not create a new denominator or a new 80% test. Both sensitivities pass only this annotation gate.", "",
              "## Complete candidate matching-reference pools", "",
              "| Interface | ID-qualified RNA rows | Eight-program union excluded | ETS2 excluded | Background candidates |", "|---|---:|---:|---:|---:|"]
    for interface, label in [(core, "Core"), (strict, "Strict")]:
        pool = pools[interface]
        lines.append(f"| {label} | {pool['ID_qualified_feature_rows']:,} | {pool['eight_original_mapped_program_union_genes']:,} | {pool['pool_status_counts']['excluded_ETS2']} | {pool['background_candidate_genes']:,} |")
    lines += ["", "Every source row appears in each pool audit, including identity-ineligible rows. All eight mapped original programs are excluded, regardless of their coverage gates. These are complete ID-only candidate lists, not expression-matched sets. No abundance bins/draws, pool-sufficiency result, normalization, scores or effects exist. Any later normalization retains the separately qualified full 41,096-row count universe; these lists must not replace it.", "",
              "## Remaining limits", "",
              "- The article declares GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4. Exact historical GTF/provider/release/checksum and counting options remain unavailable. Current-reference symbol identity does not reconstruct that GTF or establish historical locus equivalence. Allele conventions and inference-model versions are not applicable here.",
              "- Preserve **17 PSC, 17 ASC and 30 AIH** as released; ASC is not folded into PSC. This is whole-biopsy bulk liver, not sorted macrophages or cholangiocytes. The article asserts 64 patients, but an independent person/visit crosswalk and sample-linked clinical/stage/batch covariates remain unavailable.",
              "- GSE159676 remains **HOLD** for expression: unreconciled mixed-scale literals and **six paired PSC patients**, not twelve independent patients. No repair or adult expression use occurred.",
              "- NoPSC atlas remains **PARK**: no qualified complete count/feature/barcode/donor export. No atlas subset was used.",
              "- External program coverage is not evidence for ETS2 specificity, causal mechanism, transfer, replication or clinical utility. Root owns the next method choice and prospective analysis freeze.", "",
              "## Reproduce and verify", "", "```text", ".venv/bin/python scripts/prepare_psc_liver_program_inputs.py --prepare", ".venv/bin/python -m unittest discover -s tests -p test_psc_liver_program_inputs.py -v", ".venv/bin/python scripts/prepare_psc_liver_program_inputs.py --verify", "```", "",
              "`--prepare` writes only the four owned outputs/ignored work area; it does not acquire data. `--verify` checks pins read-only. The JSON schema is version 1, with source lineage, both interfaces, full feature-status counts, original program coverage/status denominators, nonoverlap gates, candidate pools, transformation/scope boundaries, preservation pins and validation receipts.", "",
              f"Independent source audit: **{source['source_audit_checks_passed']:,} checks and {source['source_audit_tests_passed']} tests passed**, with byte-identical reconstruction replay. Closed liver/organoid versioned files and prior source/artifact pins are checked unchanged.", "",
              "Final validation receipt:", "", "```json", json.dumps(data["validation"], sort_keys=True, indent=2), "```", ""]
    return "\n".join(lines)


def publish(summary, root=ROOT, output_directory=WORK):
    data = public_payload(summary, root, output_directory)
    dump_json(Path(root)/"data/derived/psc-liver-program-inputs.json", data)
    (Path(root)/"reports/psc-liver-program-inputs.md").write_text(report_text(data))
    return data


def verify_preparation(root=ROOT):
    root = Path(root)
    verify_sources(root)
    preservation = verify_closed_materials(root)
    data = json.loads((root/"data/derived/psc-liver-program-inputs.json").read_text())
    if data["schema_version"] != 1 or data["preparation_id"] != PREPARATION_ID:
        raise ValueError("Unexpected preparation schema")
    for pin in data["ignored_artifact_pins"]+[data["policy_record"], data["source_provenance"]["source_audit"], data["source_provenance"]["source_audit_verification"]]:
        check_pin(root, pin)
    for pin in data["validation"].get("checked_artifacts", []):
        check_pin(root, pin)
    for name in ["script", "tests"]:
        if sha256(root/data["software"][name]) != data["software"][name+"_sha256"]:
            raise ValueError("Prepared software bytes changed")
    summary = json.loads((root/"work/psc-liver-program-inputs/preparation-summary.json").read_text())
    expected = public_payload(summary, root, root/"work/psc-liver-program-inputs")
    if data != expected:
        raise ValueError("Public JSON differs from retained preparation and pinned provenance")
    if (root/"reports/psc-liver-program-inputs.md").read_text() != report_text(expected):
        raise ValueError("Aggregate report differs from prepared JSON")
    return {"schema_version": 1, "status": "all_preparation_and_closed_material_pins_verified",
            "all_closed_material_pins_unchanged": preservation["all_checked_pins_unchanged"],
            "ignored_preparation_artifacts_checked": len(data["ignored_artifact_pins"]),
            "network_requests": 0, "expression_analysis_performed": False}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-only", action="store_true")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if sum([args.annotation_only, args.prepare, args.verify]) != 1:
        parser.error("Choose exactly one explicit offline action")
    if args.annotation_only:
        summary, _, _ = annotation_only()
        print(json.dumps(summary, indent=2, sort_keys=True))
    elif args.prepare:
        summary = prepare()
        data = publish(summary)
        print(json.dumps({"status": data["status"], "program_coverage": summary["program_coverage"], "candidate_pools": summary["candidate_matching_reference_pools"]}, indent=2, sort_keys=True))
    elif args.verify:
        print(json.dumps(verify_preparation(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
