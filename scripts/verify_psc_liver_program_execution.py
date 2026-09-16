#!/usr/bin/env python3
"""Independent native-output audit. Default: metadata only. No author imports.

The real CLI requires explicit root acknowledgement and exact external hashes.
Tests inject artificial contracts into private helpers, never into the real CLI.
The frozen synthetic bundle entrypoint is never called or relabeled.
"""
from __future__ import annotations
import argparse
import csv
import decimal
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
import zipfile
from collections import Counter

import numpy as np
import scipy
from scipy import special

ROOT = Path(__file__).resolve().parents[1]
SELF_REL = "scripts/verify_psc_liver_program_execution.py"
LOADED_ADAPTER_SHA = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
MATH_REL = "scripts/verify_psc_liver_programs.py"
MATH_SHA = "82403f2cb63ae1eeb7e4f4a1232b21ae6772515d0e1b3ead1ca7397737a62391"
CONTRACT_REL = "work/psc-liver-program-independent-review/execution-adapter/production-contract.json"
CONTRACT_SHA = "af87402ee46bee56751f8acbc54bdf932d92db0924b434bf588666868bb53785"
OUTPUT_REL = "work/psc-liver-program-independent-review/execution-adapter/runs"
RUNS_REL = "work/psc-liver-program-plan/runs"
REAL_ACK = "ROOT_AUTHORIZES_COMPLETE_GSE303271_NATIVE_OUTPUT_AUDIT"
CORE = "exact_symbol_ensembl_bijection"
STRICT = "plus_reciprocal_entrez"
RTOL, ATOL = 1e-10, 1e-12
RELATIVE_ONLY_FIELDS = frozenset(("p_raw", "p_for_correction", "q_bh", "q_by", "variance_a", "variance_b",
    "standard_error", "degrees_of_freedom", "p_vector_for_correction", "q_bh_numeric_including_placeholders",
    "q_by_numeric_including_placeholders"))
ARTIFACTS = frozenset((
    "aggregate-public-summary.json", "private-frozen-plan-original.json",
    "private-source-units.json", "private-scores.json", "private-all-unit-omissions.json",
    "private-library-totals.json", "private-correction-vector.json", "private-compiled-programs.json",
    "private-external-identifier-audit.json", "private-bootstrap-draws.npz",
    "private-bootstrap-differences.npz"))
COMPLETE = "completion-receipt.json"

class AuditError(ValueError):
    """Messages are fixed error codes, never source values or identifiers."""

def require(condition, code):
    if not condition:
        raise AuditError(code)

def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def checked_path(root, value, *, file=True, exists=True):
    root = Path(root).absolute()
    value = Path(value)
    require(".." not in value.parts, "path_parent_traversal")
    path = value if value.is_absolute() else root / value
    require(path.is_relative_to(root), "path_outside_root")
    for candidate in (path, *path.parents):
        require(not candidate.is_symlink(), "path_symlink")
        if candidate == root:
            break
    if exists:
        require(path.exists(), "missing_path")
        mode = path.stat()
        require(stat.S_ISREG(mode.st_mode) if file else stat.S_ISDIR(mode.st_mode), "wrong_path_type")
        if file:
            require(mode.st_nlink == 1, "hardlink_input_refused")
    return path

def pin_file(root, pin):
    require(isinstance(pin, dict) and re.fullmatch(r"[0-9a-f]{64}", str(pin.get("sha256", ""))), "invalid_pin")
    path = checked_path(root, pin["path"])
    if "bytes" in pin:
        require(type(pin["bytes"]) is int and path.stat().st_size == pin["bytes"], "pinned_size_changed")
    require(digest(path) == pin["sha256"], "pinned_hash_changed")
    return path

def json_read(path):
    def pairs(entries):
        result = {}
        for key, value in entries:
            require(key not in result, "duplicate_json_key")
            result[key] = value
        return result
    def invalid(_):
        raise AuditError("nonfinite_json")
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs, parse_constant=invalid)
    except (UnicodeError, json.JSONDecodeError):
        raise AuditError("invalid_json") from None

class Compare:
    def __init__(self):
        self.scalars = 0
    def equal(self, actual, expected, field=None):
        self.scalars += 1
        if isinstance(expected, dict):
            require(type(actual) is dict and set(actual) == set(expected), "object_schema_mismatch")
            for key in expected:
                self.equal(actual[key], expected[key], key)
        elif isinstance(expected, list):
            require(type(actual) is list and len(actual) == len(expected), "array_inventory_mismatch")
            for a, e in zip(actual, expected):
                self.equal(a, e, field)
        elif type(expected) is float:
            require(type(actual) in (int, float) and math.isfinite(actual), "invalid_numeric_field")
            absolute = 0.0 if field in RELATIVE_ONLY_FIELDS else ATOL
            require(math.isclose(actual, expected, rel_tol=RTOL, abs_tol=absolute), "numeric_agreement_failure")
        else:
            require(type(actual) is type(expected) and actual == expected, "exact_field_mismatch")

def load_math(root):
    path = pin_file(root, {"path": MATH_REL, "sha256": MATH_SHA})
    spec = importlib.util.spec_from_file_location("_frozen_psc_liver_independent_math", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def runtime():
    return {"python": ".".join(map(str, sys.version_info[:3])), "numpy": np.__version__, "scipy": scipy.__version__}

def production_contract(root):
    contract = json_read(pin_file(root, {"path": CONTRACT_REL, "sha256": CONTRACT_SHA}))
    require(contract["synthetic_fixture"] is False, "artificial_contract_in_real_entrypoint")
    require(contract["runtime"] == runtime(), "runtime_changed")
    return contract

def check_source_metadata(root, contract):
    for relative, pin in contract["immutable_pins"].items():
        pin_file(root, {"path": relative, **pin})
    for pin in contract["source_metadata_pins"].values():
        pin_file(root, pin)
    review = json_read(pin_file(root, contract["accepted_review_pin"]))
    require(review["status"] == "pass_synthetic_and_guard_review" and review["real_execution_reviewed"] is False,
            "independent_acceptance_not_valid")
    require(review["compiled_programs_sha256"] == contract["compiled_sha256"], "accepted_mapping_changed")
    require(review["implementation_sha256"] == {p["path"]: p["sha256"] for p in contract["fixed_plan_fields"]["implementation_pins"]},
            "accepted_implementation_changed")

def authenticate_plan(root, path, expected_sha, contract):
    path = checked_path(root, path)
    require(str(path.relative_to(root)) in contract["plan_paths"], "not_designated_plan")
    pin_file(root, {"path": str(path.relative_to(root)), "sha256": expected_sha})
    plan = json_read(path)
    Compare().equal(sorted(plan), contract["required_plan_fields"])
    for key, value in contract["fixed_plan_fields"].items():
        require(canonical(plan.get(key)) == canonical(value), "fixed_plan_field_changed")
    require(plan["state"] == "frozen", "plan_not_frozen")
    auth = plan["root_authorization"]
    require(type(auth) is dict and set(auth) == {"approver_role", "real_execution_authorized", "reference", "approved_at"}, "root_authorization_schema")
    require(auth["approver_role"] == "root" and auth["real_execution_authorized"] is True and
            type(auth["reference"]) is str and bool(auth["reference"]) and
            type(auth["approved_at"]) is str and bool(auth["approved_at"]), "root_execution_not_authorized")
    review = plan["independent_review"]
    require(set(review) == {"status", "reference", "artifact_pin"} and
            review["status"] == "accepted_synthetic_and_guard_review" and type(review["reference"]) is str and review["reference"],
            "plan_review_not_accepted")
    require(canonical(review["artifact_pin"]) == canonical(contract["accepted_review_pin"]), "review_pin_changed")
    require(plan["expected_compiled_programs_sha256"] == contract["compiled_sha256"], "plan_mapping_not_frozen")
    return plan, path

def inventory(root, completion_path, expected_sha, plan_sha, plan, contract, *, hash_payloads):
    completion_path = checked_path(root, completion_path)
    require(completion_path.name == COMPLETE and completion_path.parent.is_relative_to(root / RUNS_REL), "not_native_run_receipt")
    run = checked_path(root, completion_path.parent, file=False)
    require(run != root / RUNS_REL, "native_run_collection_root_refused")
    require(run.stat().st_mode & 0o077 == 0 and run.stat().st_uid == os.getuid(), "run_directory_not_owner_private")
    actual = list(run.iterdir())
    require({p.name for p in actual} == ARTIFACTS | {COMPLETE} and len(actual) == 12, "incomplete_or_extra_run_inventory")
    for path in actual:
        checked_path(root, path)
    pin_file(root, {"path": str(completion_path.relative_to(root)), "sha256": expected_sha})
    receipt = json_read(completion_path)
    expected_keys = {"schema_version", "mode", "status", "plan_sha256", "root_authorization_reference", "access_audit",
                     "real_analysis_executed", "root_freeze_written", "network_used", "artifacts", "endpoint_rows", "no_clinical_or_causal_conclusion"}
    require(set(receipt) == expected_keys, "completion_schema")
    fixed = {"schema_version": 1, "mode": "real_execution", "status": "complete", "plan_sha256": plan_sha,
             "root_authorization_reference": plan["root_authorization"]["reference"], "real_analysis_executed": True,
             "root_freeze_written": False, "network_used": False, "endpoint_rows": 60, "no_clinical_or_causal_conclusion": True,
             "access_audit": {k: True for k in ("expression_payload_numeric_loading_attempted", "real_expression_values_loaded",
                 "real_scores_computed", "real_computation_attempted", "all_score_statistics_completed", "correction_completed")}}
    for key, value in fixed.items():
        Compare().equal(receipt[key], value)
    entries = receipt["artifacts"]
    require(type(entries) is list and len(entries) == 11 and [x.get("path") for x in entries] == sorted(ARTIFACTS), "artifact_manifest_inventory")
    for entry in entries:
        require(set(entry) == {"path", "sha256", "bytes", "visibility"}, "artifact_pin_schema")
        require(entry["visibility"] == ("aggregate_public_candidate" if entry["path"].startswith("aggregate-") else "private"), "artifact_visibility_changed")
        require(type(entry["bytes"]) is int and entry["bytes"] >= 0 and re.fullmatch(r"[a-f0-9]{64}", str(entry["sha256"])), "artifact_pin_invalid")
        target = checked_path(root, run / entry["path"])
        require(target.stat().st_size == entry["bytes"], "artifact_size_changed")
        if hash_payloads:
            require(digest(target) == entry["sha256"], "artifact_hash_changed")
    return run, receipt

def axes(root, contract):
    features = json_read(checked_path(root, contract["feature_audit_path"]))
    units_raw = json_read(checked_path(root, contract["unit_source_path"]))
    nr, nc = contract["shape"]
    require(len(features) == nr and [r["source_row_index_0based"] for r in features] == list(range(nr)), "feature_axis_inventory")
    require(len(units_raw) == nc and [r["matrix_column_index_0based"] for r in units_raw] == list(range(nc)), "unit_axis_inventory")
    units = [{"unit_id": r["matrix_column_id_as_submitted"], "group": r["condition_as_submitted"],
              "column_index": r["matrix_column_index_0based"], "geo_accession": r["geo_accession"],
              "unit_type": "source_patient_sample", "verified_person_id": r["separate_published_patient_identifier"]} for r in units_raw]
    require(dict(Counter(u["group"] for u in units)) == contract["groups"], "diagnosis_axis_changed")
    for key in ("unit_id", "geo_accession"):
        require(len({u[key] for u in units}) == nc, "repeated_source_unit")
    persons = [u["verified_person_id"] for u in units if u["verified_person_id"] is not None]
    require(len(persons) == len(set(persons)), "known_repeated_person")
    return [r["label_as_submitted"] for r in features], units

def count_integer(token):
    require(type(token) is str and re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", token), "count_lexical_form")
    try:
        number = decimal.Decimal(token)
        require(number.is_finite() and 0 <= number <= 2**53 and number == number.to_integral_value(), "count_not_exact_safe_integer")
        return int(number)
    except decimal.InvalidOperation:
        raise AuditError("count_not_exact_safe_integer") from None

def parse_original_gzip(path, feature_labels, units):
    nr, nc = len(feature_labels), len(units)
    data = np.empty((nr, nc), dtype=np.uint64)
    seen, header = 0, False
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.reader(stream, delimiter="\t", strict=True):
            if not row:
                continue
            if not header:
                require(row == [""] + [u["unit_id"] for u in units], "raw_header_axis_changed")
                header = True
                continue
            require(seen < nr and len(row) == nc + 1, "raw_row_inventory_changed")
            require(row[0] == feature_labels[seen], "raw_feature_axis_changed")
            data[seen] = [count_integer(value) for value in row[1:]]
            seen += 1
    require(header and seen == nr, "raw_incomplete_source")
    return data

def array_read(path, shape, dtype):
    # Header validation precedes allocation; object arrays and giant headers fail.
    with Path(path).open("rb") as stream:
        version = np.lib.format.read_magic(stream)
        require(version == (1, 0), "unexpected_npy_version")
        declared_shape, fortran, declared_dtype = np.lib.format.read_array_header_1_0(stream, max_header_size=4096)
        require(tuple(declared_shape) == tuple(shape) and not fortran and declared_dtype == np.dtype(dtype), "npy_header_changed")
        require(Path(path).stat().st_size == stream.tell() + math.prod(shape) * np.dtype(dtype).itemsize, "npy_length_changed")
    return np.load(path, allow_pickle=False)

def npz_read(path, shapes, dtype):
    with zipfile.ZipFile(path) as archive:
        infos = archive.infolist()
        require(len(infos) == len(shapes) and {x.filename for x in infos} == {key + ".npy" for key in shapes}, "npz_member_inventory")
        for info in infos:
            require(info.filename.endswith(".npy") and not info.is_dir(), "npz_unsafe_member")
            shape = shapes[info.filename[:-4]]
            require(info.file_size <= math.prod(shape) * np.dtype(dtype).itemsize + 4096, "npz_size_limit")
            with archive.open(info) as stream:
                require(np.lib.format.read_magic(stream) == (1, 0), "npz_npy_version")
                sh, fortran, dt = np.lib.format.read_array_header_1_0(stream, max_header_size=4096)
                require(tuple(sh) == tuple(shape) and not fortran and dt == np.dtype(dtype), "npz_header_changed")
                require(info.file_size == stream.tell() + math.prod(shape) * np.dtype(dtype).itemsize, "npz_member_length")
    with np.load(path, allow_pickle=False) as archive:
        result = {key: archive[key] for key in shapes}
    require(all(np.isfinite(value).all() for value in result.values()), "npz_nonfinite")
    return result

def unavailable(n_a, n_b):
    return {"n_a": n_a, "n_b": n_b, "mean_a": None, "mean_b": None, "variance_a": None, "variance_b": None,
        "effect": None, "effect_status": "unavailable_score", "welch_status": "unavailable_score",
        "variance_a_status": "unavailable_score", "variance_b_status": "unavailable_score", "numeric_inference_failure": None,
        "standard_error": None, "degrees_of_freedom": None, "t_statistic": None, "p_raw": None, "log_p_raw": None,
        "p_underflow_to_zero": False, "ci_low": None, "ci_high": None, "confidence": 0.95,
        "interval_type": "pointwise_not_multiplicity_adjusted", "zero_variance_groups": []}

def moment_status(v, values):
    require(bool(values) and all(type(x) in (int, float) and math.isfinite(x) for x in values), "invalid_score_vector")
    constant = min(values) == max(values)
    mean = float(values[0]) if constant else math.fsum(x / len(values) for x in values)
    if constant:
        return mean, 0.0, 0.0, "exact_zero"
    try:
        mean, sd = v.moments(values)
    except ValueError:
        # Only label a nonrepresentable SD/variance. Never replace it with an
        # epsilon or calculate different valid inference from altered scores.
        centered = [x - mean for x in values]
        scale = max(map(abs, centered))
        if not math.isfinite(scale):
            return mean, math.inf, None, "unavailable_numeric_variance"
        sd = scale * math.sqrt(math.fsum((x / scale)**2 for x in centered) / (len(values) - 1))
        require(sd == 0 or not math.isfinite(sd), "unclassified_independent_moment_failure")
    variance = sd * sd
    status = "finite_positive" if 0 < variance < math.inf else "unavailable_numeric_variance_underflow" if variance == 0 else "unavailable_numeric_variance"
    return mean, sd, variance if status == "finite_positive" else None, status


def expected_welch(v, a, b):
    result = unavailable(len(a), len(b))
    result.update(effect_status="available" if a and b else "unavailable_no_group_observations",
                  welch_status="unavailable_n_below_2", variance_a_status="not_estimated", variance_b_status="not_estimated")
    for label, values in (("a", a), ("b", b)):
        result["mean_" + label] = moment_status(v, values)[0] if values else None
    if a and b:
        result["effect"] = result["mean_a"] - result["mean_b"]
    if min(len(a), len(b)) < 2:
        return result
    moments = [moment_status(v, values) for values in (a, b)]
    means_sd = [(mean, sd) for mean, sd, _, _ in moments]
    for label, (_, sd, variance, status) in zip(("a", "b"), moments):
        result["variance_" + label + "_status"] = status
        result["variance_" + label] = variance
        if status == "exact_zero":
            result["zero_variance_groups"].append(label)
    for label in ("a", "b"):
        if result["variance_" + label] is None:
            status = result["variance_" + label + "_status"]
            result.update(welch_status=status, numeric_inference_failure=status)
            return result
    ref = v.welch_reference(a, b)
    if ref["status"] == "zero_standard_error":
        result["welch_status"] = "unavailable_zero_variance_policy"
    elif ref["status"] == "ok":
        result.update(welch_status="available", standard_error=ref["se"], degrees_of_freedom=ref["df"],
                      t_statistic=ref["t"], p_raw=ref["p_raw"], log_p_raw=math.log(ref["p_raw"]), ci_low=ref["ci_low"], ci_high=ref["ci_high"])
    else:
        se = math.hypot(means_sd[0][1] / math.sqrt(len(a)), means_sd[1][1] / math.sqrt(len(b)))
        terms = [sd / math.sqrt(n) / se for (_, sd), n in zip(means_sd, (len(a), len(b)))]
        df = 1 / sum(w**4 / (n - 1) for w, n in zip(terms, (len(a), len(b))))
        t = result["effect"] / se
        tail = math.isfinite(t) and 2 * special.stdtr(df, -abs(t)) == 0
        result.update(welch_status="unavailable_numeric_tail_underflow" if tail else "unavailable_numeric_inference",
                      numeric_inference_failure="positive_tail_probability_not_representable" if tail else "nonfinite_or_invalid_inference_arithmetic",
                      p_underflow_to_zero=bool(tail))
    return result

def expected_bootstrap(v, values, replicates):
    low = None if values is None else v.percentile_linear(values, 0.025)
    high = None if values is None else v.percentile_linear(values, 0.975)
    return {"status": "unavailable_score" if values is None else "available_degenerate" if low == high else "available",
            "requested_replicates": replicates, "valid_replicates": 0 if values is None else replicates,
            "ci_low": low, "ci_high": high, "method": "percentile_linear_95", "p_value_computed": False}

def omissions(v, scores, units, contrast, base):
    contrast_id, left, right = contrast
    output = []
    for index, unit in enumerate(units):
        involved = unit["group"] in (left, right)
        stats = None if scores is None else base if not involved else expected_welch(v,
            [score for i, score in enumerate(scores) if i != index and units[i]["group"] == left],
            [score for i, score in enumerate(scores) if i != index and units[i]["group"] == right])
        effect = None if stats is None else stats["effect"]
        baseline = base["effect"]
        known = effect is not None and baseline is not None
        output.append({"unit_id": unit["unit_id"], "unit_group": unit["group"], "omitted_index": index,
            "in_contrast": involved, "contrast_id": contrast_id, "effect_delta": effect - baseline if known else None,
            "strict_sign_reversal": ((effect < 0 < baseline) or (baseline < 0 < effect)) if known else None,
            "sensitivity_only": True, "discovery_q_values": None,
            "status": "unavailable_score" if scores is None else "available" if involved else "unchanged_outside_contrast", "statistics": stats})
    relevant = [r for r in output if r["in_contrast"]]
    defined = [r for r in relevant if r["statistics"] is not None and r["statistics"]["effect"] is not None]
    effects = [r["statistics"]["effect"] for r in defined]
    changes = [r["effect_delta"] for r in defined if r["effect_delta"] is not None]
    summary = {"status": "available" if effects else "unavailable_effects", "relevant_unit_denominator": len(relevant),
        "in_contrast_rows": len(relevant), "denominator_kind": "in_contrast_units_only", "defined_effect_rows": len(defined),
        "effect_min": min(effects) if effects else None, "effect_max": max(effects) if effects else None,
        "max_absolute_effect_change": max(map(abs, changes)) if changes else None,
        "strict_sign_reversal_rows": sum(r["strict_sign_reversal"] is True for r in defined) if effects else None, "sensitivity_only": True}
    return output, summary

def numeric_domain(value, *, positive=False, probability=False):
    require(type(value) in (float, int) and math.isfinite(value), "nonfinite_native_number")
    require(not positive or value > 0, "native_positive_number_not_positive")
    require(not probability or 0 < value <= 1, "native_probability_out_of_domain")


def validate_native_stats(row):
    require(type(row) is dict, "native_statistics_schema")
    for arm in ("a", "b"):
        value, status = row["variance_" + arm], row["variance_" + arm + "_status"]
        if status == "finite_positive":
            numeric_domain(value, positive=True)
        elif status == "exact_zero":
            require(type(value) in (int, float) and value == 0, "exact_zero_variance_changed")
        else:
            require(value is None, "unavailable_variance_is_not_null")
    if row["welch_status"] == "available":
        numeric_domain(row["p_raw"], probability=True)
        numeric_domain(row["standard_error"], positive=True)
        numeric_domain(row["degrees_of_freedom"], positive=True)
        for key in ("t_statistic", "effect", "ci_low", "ci_high"):
            numeric_domain(row[key])
        require(row["ci_low"] <= row["effect"] <= row["ci_high"], "welch_interval_order")
        require(row["p_underflow_to_zero"] is False, "available_p_flagged_underflow")
        if row["log_p_raw"] is not None:
            numeric_domain(row["log_p_raw"])
            require(row["log_p_raw"] <= 0, "log_probability_out_of_domain")
    else:
        require(all(row[key] is None for key in ("p_raw", "log_p_raw", "standard_error", "degrees_of_freedom", "t_statistic", "ci_low", "ci_high")),
                "unavailable_inference_not_null")


def native_domains_and_exact_copies(public, private_scores, private_omissions, correction):
    require(type(public) is dict and type(public["endpoint_rows"]) is list, "native_public_schema")
    require(len(public["endpoint_rows"]) == 60, "native_endpoint_row_count")
    endpoints = {}
    stat_keys = set(unavailable(0, 0))
    for row in public["endpoint_rows"]:
        validate_native_stats(row)
        key = (row["interface"], row["program_id"], row["contrast_id"])
        require(key not in endpoints, "duplicate_native_endpoint")
        endpoints[key] = row
        for name in ("q_bh", "q_by", "p_for_correction"):
            if row[name] is not None:
                numeric_domain(row[name], probability=True)
        if row["interface"] == CORE and row["p_raw"] is None:
            require(row["p_for_correction"] == 1 and row["q_bh"] is None and row["q_by"] is None, "missing_core_p_policy")
        boot = row["bootstrap"]
        if boot["status"] != "unavailable_score":
            numeric_domain(boot["ci_low"])
            numeric_domain(boot["ci_high"])
            require(boot["ci_low"] <= boot["ci_high"], "bootstrap_interval_order")
            require((boot["ci_low"] == boot["ci_high"]) == (boot["status"] == "available_degenerate"), "bootstrap_degenerate_status")
    for row in private_scores:
        if row["values"] is not None:
            for value in row["values"]:
                numeric_domain(value)
                require(value >= 0, "unsigned_score_negative")
    for row in private_omissions:
        if row["statistics"] is not None:
            validate_native_stats(row["statistics"])
        if not row["in_contrast"] and row["statistics"] is not None:
            baseline = endpoints[(row["interface"], row["program_id"], row["contrast_id"])]
            require(row["statistics"] == {key: baseline[key] for key in stat_keys}, "outside_contrast_statistics_not_exact_copy")
            require(row["effect_delta"] == 0 and row["strict_sign_reversal"] is False, "outside_contrast_change_not_exact_zero")
    for key in ("p_vector_for_correction", "q_bh_numeric_including_placeholders", "q_by_numeric_including_placeholders"):
        for value in correction[key]:
            numeric_domain(value, probability=True)
    for index, key in enumerate(correction["ordered_keys"]):
        row = endpoints.get((CORE, *key))
        require(row is not None, "unknown_correction_endpoint")
        if row["p_raw"] is None:
            require(all(correction[name][index] == 1 for name in ("p_vector_for_correction", "q_bh_numeric_including_placeholders", "q_by_numeric_including_placeholders")), "private_missing_placeholders_not_exact_one")


def numerical_audit(root, run, plan, contract, v):
    """Private helper. Caller has authenticated all sources/artifacts and permission."""
    compare = Compare()
    native_public = json_read(run / "aggregate-public-summary.json")
    native_scores = json_read(run / "private-scores.json")
    native_omissions = json_read(run / "private-all-unit-omissions.json")
    native_correction = json_read(run / "private-correction-vector.json")
    native_domains_and_exact_copies(native_public, native_scores, native_omissions, native_correction)
    labels, units = axes(root, contract)
    compare.equal(json_read(run / "private-source-units.json"), units)
    specs = json_read(run / "private-compiled-programs.json")
    require(canonical(specs) == contract["compiled_sha256"], "independently_accepted_mapping_mismatch")
    require(digest(run / "private-external-identifier-audit.json") == contract["external_audit_sha256"], "external_identity_audit_changed")
    programs = [s["program_id"] for s in specs if s["interface"] == CORE]
    keys = [(s["interface"], s["program_id"]) for s in specs]
    require(len(programs) == 10 and len(set(programs)) == 10 and keys == [(tier, p) for tier in (CORE, STRICT) for p in programs], "compiled_inventory_changed")
    raw = parse_original_gzip(root / contract["raw_pin"]["path"], labels, units)
    cache = array_read(root / contract["cache_pin"]["path"], contract["shape"], "float64")
    require(np.isfinite(cache).all() and np.array_equal(raw, cache), "complete_source_cache_mismatch")
    normalized = v.logcpm_reference(raw.tolist())
    native_totals = json_read(run / "private-library-totals.json")
    require(type(native_totals) is list and len(native_totals) == len(units), "library_total_inventory")
    require(all(type(a) in (float, int) and math.isfinite(a) and a == e
                for a, e in zip(native_totals, normalized["library_totals"])), "library_totals_not_exact")
    groups = [u["group"] for u in units]
    replicates = contract["method"]["bootstrap"]["replicates"]
    draws = v.stratified_indices(groups, replicates, contract["method"]["bootstrap"]["seed"])
    native_draws = npz_read(run / "private-bootstrap-draws.npz", {g: (replicates, n) for g, n in contract["groups"].items()}, "int64")
    for group in draws:
        require(np.array_equal(native_draws[group], np.array(draws[group], dtype=np.int64)), "shared_rng_indices_changed")
    available = [s for s in specs if s["coverage"]["score_status"] == "available"]
    vectors = npz_read(run / "private-bootstrap-differences.npz", {
        "__".join((s["interface"], s["program_id"], c[0])): (replicates,) for s in available for c in contract["contrasts"]}, "float64")
    score_rows, endpoint_rows, omission_rows = [], [], []
    bootstrap_values = 0
    for spec in specs:
        require(spec["score_kind"] == "unsigned_equal_gene" and len(spec["components"]) == 1, "new_or_signed_endpoint_refused")
        indices = spec["components"][0]["indices"]
        scores = None if spec not in available else v.score_reference(normalized["logcpm"], indices)
        score_rows.append({"interface": spec["interface"], "program_id": spec["program_id"], "values": scores})
        for contrast in contract["contrasts"]:
            cid, left, right = contrast
            stats = unavailable(groups.count(left), groups.count(right)) if scores is None else expected_welch(v,
                [x for x, group in zip(scores, groups) if group == left], [x for x, group in zip(scores, groups) if group == right])
            bootstrap = None if scores is None else v.bootstrap_reference(scores, groups, draws, left, right)["draws"]
            if bootstrap is not None:
                compare.equal(vectors["__".join((spec["interface"], spec["program_id"], cid))].tolist(), bootstrap)
                bootstrap_values += len(bootstrap)
            omitted, summary = omissions(v, scores, units, contrast, stats)
            for entry in omitted:
                entry.update(interface=spec["interface"], program_id=spec["program_id"])
            omission_rows.extend(omitted)
            endpoint_rows.append({"program_id": spec["program_id"], "interface": spec["interface"], "contrast_id": cid,
                "group_a": left, "group_b": right, "contrast_role": "primary" if cid == contract["contrasts"][0][0] else "secondary",
                "inferential_role": "sole_discovery_family" if spec["interface"] == CORE else "identity_sensitivity_only",
                "score_kind": spec["score_kind"], "coverage": spec["coverage"], "target_gene_membership": spec.get("target_gene_membership"),
                **stats, "bootstrap": expected_bootstrap(v, bootstrap, replicates), "omission": summary})
    compare.equal(native_scores, score_rows)
    compare.equal(native_omissions, omission_rows)
    core_rows = [r for r in endpoint_rows if r["interface"] == CORE]
    raw_p = [r["p_raw"] for r in core_rows]
    adjusted = v.adjust_reference(raw_p)
    correction = {"family_size": 30, "missing_raw_p": sum(p is None for p in raw_p),
        "p_vector_for_correction": [1.0 if p is None else p for p in raw_p],
        "q_bh_numeric_including_placeholders": adjusted["q_bh"], "q_by_numeric_including_placeholders": adjusted["q_by"],
        "ordered_keys": [[r["program_id"], r["contrast_id"]] for r in core_rows]}
    for i, row in enumerate(core_rows):
        missing = row["p_raw"] is None
        row.update(p_for_correction=1.0 if missing else row["p_raw"], p_was_substituted=missing,
            q_bh=None if missing else adjusted["q_bh"][i], q_by=None if missing else adjusted["q_by"][i],
            family_id="Core_full_panel_three_contrasts", family_size=30, discovery_eligible=not missing)
    for row in endpoint_rows[30:]:
        row.update(p_for_correction=None, p_was_substituted=False, q_bh=None, q_by=None, family_id=None, family_size=None, discovery_eligible=False)
    compare.equal(native_correction, correction)
    public = {"schema_version": 1, "status": "complete_for_declared_analysis", "analysis_units": "source_patient_samples_not_verified_unique_donors",
        "independence_status": "article_claim_only_crosswalk_unavailable", "group_sizes": contract["groups"],
        "released_normalization_rows": contract["shape"][0], "normalization": "all_released_rows_CPM_then_log2_plus_1",
        "endpoint_order": programs, "endpoint_rows": endpoint_rows, "expected_endpoint_rows": 60,
        "expected_all_unit_omission_rows": 60 * len(units), "core_family_size": 30, "unavailable_core_p_values": correction["missing_raw_p"],
        "public_omission_denominators": {cid: groups.count(a) + groups.count(b) for cid, a, b in contract["contrasts"]},
        "all_unit_omission_row_count_is_private_inventory_not_stability_denominator": True, "strict_discoveries_permitted": False,
        "does_not_fulfill_held_IL17_response_transfer": True, "not_independent_mechanistic_confirmation_of_MacroMap": True,
        "interpretation_limits": ["Whole-biopsy observational transcript-abundance composites, not pathway activity or cell-intrinsic effects.",
            "Pointwise intervals are not simultaneous. BH/BY do not repair invalid marginal tests or unverified independence.",
            "No clinical/stage/batch adjustment, causal interpretation, clinical outcome, or clinical recommendation."],
        "mode": "real_execution", "real_analysis_executed": True, "plan_sha256": digest(run / "private-frozen-plan-original.json"),
        "method": contract["method"], "software": contract["runtime"], "root_freeze_written": False,
        "source_pin_status": "all_frozen_source_and_implementation_pins_verified", "external_candidates": plan["external_candidates"]}
    compare.equal(native_public, public)
    return {"raw_source_cells_and_cache_values_compared": int(raw.size), "source_rows_checked": len(labels),
        "source_columns_checked": len(units), "full_library_totals_compared": len(units), "score_vectors_checked": 20,
        "score_values_compared": sum(len(r["values"]) for r in score_rows if r["values"] is not None),
        "unavailable_score_vectors_checked": sum(r["values"] is None for r in score_rows), "endpoint_rows_checked": 60,
        "core_correction_vector_length": 30, "unavailable_core_p_values_checked": correction["missing_raw_p"],
        "shared_rng_indices_compared": sum(x.size for x in native_draws.values()), "bootstrap_differences_compared": bootstrap_values,
        "bootstrap_interval_status_rows_checked": 60, "private_omissions_checked": len(omission_rows), "public_omission_summaries_checked": 60,
        "structured_scalar_checks": compare.scalars, "native_logcpm_matrix_comparison": "not_available_not_claimed"}

def fresh_output(root, value):
    path = checked_path(root, value, file=False, exists=False)
    require(path.is_relative_to(root / OUTPUT_REL) and path != root / OUTPUT_REL, "audit_output_not_private_work")
    require(not path.exists(), "audit_output_reuse")
    return path

def run_authenticated_audit(root, contract, plan_path, plan_sha, completion_path, completion_sha, output, ack, v):
    """Auth wrapper used with a pinned production contract, or explicit test fixtures."""
    require(ack == REAL_ACK, "explicit_root_real_audit_ack_required")
    pin_file(ROOT, {"path": SELF_REL, "sha256": LOADED_ADAPTER_SHA})
    require(contract["runtime"] == runtime(), "runtime_changed")
    check_source_metadata(root, contract)
    plan, path = authenticate_plan(root, plan_path, plan_sha, contract)
    target = fresh_output(root, output)
    run, receipt = inventory(root, completion_path, completion_sha, plan_sha, plan, contract, hash_payloads=True)
    require((run / "private-frozen-plan-original.json").read_bytes() == path.read_bytes(), "native_original_plan_bytes_changed")
    # Payload hashes are checked after every permission/metadata/inventory guard.
    pin_file(root, contract["raw_pin"])
    pin_file(root, contract["cache_pin"])
    try:
        counts = numerical_audit(root, run, plan, contract, v)
    except AuditError:
        raise
    except Exception:
        raise AuditError("native_payload_unavailable_or_invalid") from None
    # Detect input/source mutation during the audit before issuing success.
    inventory(root, completion_path, completion_sha, plan_sha, plan, contract, hash_payloads=True)
    check_source_metadata(root, contract)
    pin_file(root, contract["raw_pin"])
    pin_file(root, contract["cache_pin"])
    pin_file(root, {"path": str(path.relative_to(root)), "sha256": plan_sha})
    pin_file(ROOT, {"path": SELF_REL, "sha256": LOADED_ADAPTER_SHA})
    result = {"schema_version": 1, "status": "pass_complete_native_output_verification",
        "real_data_audited": contract["synthetic_fixture"] is False, "synthetic_fixture": contract["synthetic_fixture"],
        "root_execution_authorization_created": False, "root_audit_acknowledgment_supplied": True,
        "plan_sha256": plan_sha, "execution_receipt_sha256": completion_sha, "compiled_programs_sha256": contract["compiled_sha256"],
        "frozen_independent_math_sha256": MATH_SHA, "adapter_sha256": LOADED_ADAPTER_SHA,
        "production_contract_sha256": CONTRACT_SHA if not contract["synthetic_fixture"] else None,
        "runtime": runtime(), "checks": counts,
        "agreement": {"relative_tolerance": RTOL, "absolute_tolerance": ATOL,
                      "relative_only_fields": sorted(RELATIVE_ONLY_FIELDS), "counts_axes_statuses_inventory_indices": "exact"},
        "native_payload_hashes": {e["path"]: e["sha256"] for e in receipt["artifacts"]},
        "source_payload_hashes": {p["path"]: p["sha256"] for p in (contract["raw_pin"], contract["cache_pin"])},
        "source_values_or_unit_identifiers_exported": False, "biological_or_clinical_conclusion": "not_assessed"}
    fresh_output(root, target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.mkdir(mode=0o700, exist_ok=False)
    output_path = target / "audit-receipt.json"
    partial = output_path.with_suffix(".json.partial")
    descriptor = os.open(partial, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    partial.rename(output_path)
    return result

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-real", action="store_true", help="Future root-authorized complete native run audit; default is metadata only")
    parser.add_argument("--expected-adapter-sha256")
    parser.add_argument("--plan")
    parser.add_argument("--expected-plan-sha256")
    parser.add_argument("--execution-receipt")
    parser.add_argument("--expected-execution-receipt-sha256")
    parser.add_argument("--root-real-audit-ack")
    parser.add_argument("--output-directory")
    args = parser.parse_args(argv)
    receipt_issued = False
    try:
        contract = production_contract(ROOT)
        if not args.verify_real:
            require(not any((args.plan, args.execution_receipt, args.output_directory, args.root_real_audit_ack,
                             args.expected_adapter_sha256, args.expected_plan_sha256, args.expected_execution_receipt_sha256)), "real_arguments_without_mode")
            check_source_metadata(ROOT, contract)
            print(json.dumps({"status": "source_only_contract_authenticated", "real_expression_accessed": False,
                              "real_execution_or_audit_authorized": False, "contract_sha256": CONTRACT_SHA}, sort_keys=True))
            return 0
        require(args.root_real_audit_ack == REAL_ACK, "explicit_root_real_audit_ack_required")
        require(all((args.plan, args.expected_plan_sha256, args.execution_receipt, args.expected_execution_receipt_sha256,
                     args.output_directory, args.expected_adapter_sha256)), "missing_real_audit_arguments")
        require(args.expected_adapter_sha256 == LOADED_ADAPTER_SHA, "loaded_adapter_hash_changed")
        pin_file(ROOT, {"path": SELF_REL, "sha256": args.expected_adapter_sha256})
        result = run_authenticated_audit(ROOT, contract, args.plan, args.expected_plan_sha256, args.execution_receipt,
            args.expected_execution_receipt_sha256, args.output_directory, args.root_real_audit_ack, load_math(ROOT))
        receipt_issued = True
        print(json.dumps({"status": result["status"], "real_data_audited": True, "checks": result["checks"]}, sort_keys=True))
        return 0
    except Exception as error:
        # Never echo math values, raw tokens, identifiers, or a traceback.
        reason = str(error) if isinstance(error, AuditError) else "unavailable_or_invalid_audit_input"
        print(json.dumps({"status": "refused_or_failed_not_verified", "reason": reason,
                          "no_success_receipt_issued": not receipt_issued}, sort_keys=True))
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
