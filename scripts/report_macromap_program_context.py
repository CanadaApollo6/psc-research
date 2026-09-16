#!/usr/bin/env python3
"""Authenticated aggregate-only MacroMap reporting. Does not import analysis code.

This is post-fit presentation, not a scientific amendment or a numerical verifier.
Only fixed aggregate JSON/CSV inputs are read. All destinations must be fresh.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import io
import itertools
import json
import math
import os
from pathlib import Path
import platform
import zlib

ROOT = Path(__file__).resolve().parents[1]
SCOPES = ("mapped_primary", "all_source_lines_sensitivity")
PROTOCOLS = ("Mod_SmartSeq2", "NEB")
TIMES = (6, 24)
STIMULI = ("CIL", "IFNB", "IFNG", "IL4", "LIL10", "MBP", "P3C", "R848", "sLPS", "PIC")
PROGRAMS = ("ets2_g1_dn", "ets2_g1_up", "ets2_g2_dn", "chr21_dn", "inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling", "ets2_g1_dn_without_comparators")
VARIANTS = (PROGRAMS[0], PROGRAMS[-1])
QUANTITIES = ("raw_change", "background_change", "adjusted_change")
COUNTS = {SCOPES[0]: {PROTOCOLS[0]: 114, PROTOCOLS[1]: 71}, SCOPES[1]: {PROTOCOLS[0]: 114, PROTOCOLS[1]: 75}}
VERIFIED = "independently_verified_with_declared_source_subset"
RESPONSE_INTERVAL = "pointwise_fixed_score_conditional_only"
PREDICTION_INTERVAL = "fixed_prediction_conditional_only"
NO_INTERVAL = "not_computed_prespecified"
STATUS = {"available", "available_descriptive", "unavailable_mock_identity", "missing_response", "insufficient_paired_lines", "incomplete_predictions", "missing_protocol", "unpaired_times", "missing_evaluation_stimulus", "empty_resampled_protocol", "empty_fold", "missing_predictor", "unknown_stimulus", "failed_paired_fit", "failed_prediction", "invalid_training_data", "missing_training_stimulus", "optimizer_failed", "mapping_failed", "reference_unavailable", "background_underflow", "insufficient_training_class_lines", "insufficient_training_class_runs", "insufficient_training_controls", "not_fitted", "insufficient_class_lines_or_runs", "empty_test_fold", "target_outside_bins", "empty_target", "insufficient_control_pool"}
GATE_STATUS = {"available", "empty_test_fold", "insufficient_training_class_lines", "insufficient_training_class_runs"}
OUTER_STATUS = GATE_STATUS | {"empty_fold", "missing_predictor", "unknown_stimulus", "failed_paired_fit", "failed_prediction"}
COMPONENT_STATUS = OUTER_STATUS | {"invalid_training_data", "missing_training_stimulus", "optimizer_failed"}
ENDPOINT_STATUS = {"available", "incomplete_predictions", "missing_protocol", "unpaired_times", "missing_evaluation_stimulus", "empty_resampled_protocol"}
MATCH_STATUS = {"available", "target_outside_bins", "empty_target", "insufficient_control_pool", "reference_unavailable"}
LABELS = {"ets2_g1_dn": "ETS2 g1 down", "ets2_g1_up": "ETS2 g1 up", "ets2_g2_dn": "ETS2 g2 down", "chr21_dn": "chr21 enhancer down", "inflammation": "Inflammation", "interferon_gamma": "Interferon gamma", "oxidative_stress": "Oxidative stress", "apoptotic_signaling": "Apoptotic signaling", "ets2_g1_dn_without_comparators": "ETS2 g1 down, nonoverlap"}
SCOPE_LIMITS = ["Independent raw normalization uses eight samples, not every source count.", "Independent reference/bin/draw/score reconstruction uses fold0 only; all folds have identity/X/model/loss checks.", "Response CI values independently replayed for CIL6h only, all programs/quantities in both protocols.", "Bootstrap conditions on fixed OOF scoring/models; it does not refit the complete procedure."]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def absolute(path):
    path = Path(path)
    require(".." not in path.parts, "Literal parent traversal is not allowed")
    path = path.absolute()
    require(not any(p.is_symlink() for p in (path, *path.parents)), "Symlink paths are not allowed")
    return path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "Duplicate JSON key")
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError("Nonfinite JSON token")
    return json.loads(data, object_pairs_hook=pairs, parse_constant=nonfinite)


def read_pin(path, expected):
    require(isinstance(expected, str) and len(expected) == 64 and all(c in "0123456789abcdef" for c in expected), "Invalid external hash")
    path = absolute(path)
    require(path.is_file(), "Missing pinned input")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    with os.fdopen(fd, "rb") as stream:
        before = os.fstat(stream.fileno())
        require(before.st_size < 64 * 1024 * 1024, "Aggregate input exceeds size limit")
        data = stream.read()
        after = os.fstat(stream.fileno())
    absolute(path)
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    require(identity(before) == identity(after) == identity(path.stat()), "Input changed while reading")
    require(digest(data) == expected, "Pinned input hash differs")
    return data


def contained(base, name):
    require(isinstance(name, str) and not Path(name).is_absolute() and ".." not in Path(name).parts, "Invalid relative artifact")
    path = absolute(Path(base) / name)
    require(path.is_relative_to(absolute(base)), "Artifact escapes its bundle")
    return path


def number(value, nullable=False):
    if value is None and nullable:
        return None
    require(type(value) in (int, float) and math.isfinite(value), "Expected finite numeric aggregate")
    return value


def integer(value, nullable=False):
    if value is None and nullable:
        return None
    number(value)
    require(type(value) is int and value >= 0, "Expected nonnegative integer count")
    return value


def status(value):
    require(value in STATUS, "Unrecognized status; reporting review is required")
    return value


def interval(value, nullable=True):
    if value is None and nullable:
        return (None, None)
    require(isinstance(value, list) and len(value) == 2, "Expected two interval bounds")
    a, b = map(number, value)
    require(a <= b, "Reversed interval")
    return a, b


def close(a, b):
    return math.isclose(number(a), number(b), rel_tol=1e-11, abs_tol=1e-12)


def labels(scope, variant=None):
    require(scope in SCOPES and (variant is None or variant in VARIANTS), "Unknown fixed scope/variant")
    return {"identity_scope": scope, "cohort_role": "mapped_primary" if scope == SCOPES[0] else "identity_inclusion_sensitivity_not_independent_cohort", **({"variant": variant, "variant_role": "primary_program" if variant == VARIANTS[0] else "fixed_nonoverlap_sensitivity"} if variant is not None else {})}


def unique_rows(rows, keys, expected):
    require(isinstance(rows, list), "Expected aggregate row list")
    found = {}
    for row in rows:
        require(isinstance(row, dict), "Expected aggregate object")
        for field in ("time", "fold"):
            if field in keys: require(type(row.get(field)) is int, "Noninteger time/fold identity")
        key = tuple(row.get(k) for k in keys)
        require(key not in found, "Duplicate aggregate key")
        found[key] = row
    require(set(found) == set(expected), "Missing/extra fixed aggregate keys")
    return found


def response_table(panel, scope):
    keys = ("identity_scope", "protocol", "time", "stimulus", "program_id", "quantity")
    expected = [(scope, p, t, s, g, q) for p, t, s, g in itertools.product((*PROTOCOLS, "pooled_descriptive"), TIMES, STIMULI, PROGRAMS) for q in ((None,) if s == "PIC" else QUANTITIES)]
    index = unique_rows(panel, keys, expected)
    rows, counts = [], {}
    for key in expected:
        src = index[key]; _, p, t, s, g, q = key
        st = status(src["status"])
        require(st in ("available_descriptive", "missing_response", "insufficient_paired_lines", "unavailable_mock_identity"), "Unknown response summary state")
        identity_fields = {"identity_scope", "protocol", "time", "stimulus", "program_id", "expected_paired_lines", "expected_paired_runs", "status"}
        summary_fields = {"lines", "runs", "median", "mean", "minimum", "maximum", "positive", "negative", "zero", "leave_one_line_out_median_range", "leave_one_run_out_median_range"}
        interval_fields = {"interval_status", "median_interval", "mean_interval", "replicates", "seed"}
        pooled_fields = {"median_interval_status", "mean_interval_status", "original_contrast_protocol_line_shares"} if p == "pooled_descriptive" else set()
        allowed = identity_fields if s == "PIC" else identity_fields | {"quantity", "input_statuses"} | summary_fields | interval_fields | pooled_fields
        require(set(src) <= allowed, "Unexpected response metadata field")
        nl, nr = integer(src["expected_paired_lines"]), integer(src["expected_paired_runs"])
        ck = (p, t, s)
        require(ck not in counts or counts[ck] == (nl, nr), "Response original counts vary by program/quantity")
        counts[ck] = (nl, nr)
        row = {**labels(scope), "protocol": p, "time": t, "stimulus": s, "program_id": g, "quantity": q,
               "status": st, "expected_paired_lines": nl, "expected_paired_runs": nr}
        available = st == "available_descriptive"
        if s == "PIC":
            require(st == "unavailable_mock_identity" and nl == nr == 0, "PIC must remain unavailable")
        else:
            statuses = src.get("input_statuses", {})
            require(isinstance(statuses, dict), "Invalid response status counts")
            for sk, count in statuses.items(): status(sk); integer(count)
            require(sum(statuses.values()) == nl, "Response statuses dropped expected lines")
            row["input_status_counts"] = canonical(statuses)
        row["lines"] = integer(src.get("lines"), nullable=True)
        row["runs"] = integer(src.get("runs"), nullable=True)
        for field in ("mean", "median", "minimum", "maximum"):
            row[field] = number(src.get(field), nullable=not available)
            if not available: require(row[field] is None, "Unavailable response cannot have zero or another effect value")
        for field in ("positive", "negative", "zero"):
            row[field] = integer(src.get(field), nullable=not available)
            if not available: require(row[field] is None, "Unavailable response cannot have sign counts")
        if available:
            require(row["lines"] == nl and row["runs"] == nr and sum(row[x] for x in ("positive", "negative", "zero")) == nl, "Response counts do not match expected cohort")
            require(row["minimum"] <= row["median"] <= row["maximum"] and row["minimum"] <= row["mean"] <= row["maximum"], "Response summaries outside range")
        for field in ("leave_one_line_out_median_range", "leave_one_run_out_median_range", "median_interval", "mean_interval"):
            row[field + "_low"], row[field + "_high"] = interval(src.get(field))
            if not available: require(src.get(field) is None, "Unavailable response cannot have interval/range")
        row["interval_status"] = src.get("interval_status", "not_available")
        require(row["interval_status"] in (RESPONSE_INTERVAL, NO_INTERVAL, "insufficient_runs", "not_available"), "Unknown response interval status")
        if p == "pooled_descriptive":
            require(src.get("median_interval") is None and src.get("mean_interval") is None, "Pooled RESPONSE intervals are prohibited")
            if s != "PIC":
                require(src.get("mean_interval_status") == src.get("median_interval_status") == NO_INTERVAL and row["interval_status"] == NO_INTERVAL, "Pooled response interval status differs")
                shares = src.get("original_contrast_protocol_line_shares", {})
                require(set(shares) == set(PROTOCOLS), "Missing original response protocol shares")
                for proto in PROTOCOLS: row["original_share_" + proto] = number(shares[proto])
        elif available:
            require(row["interval_status"] in (RESPONSE_INTERVAL, "insufficient_runs"), "Missing protocol-specific response interval status")
            if row["interval_status"] == RESPONSE_INTERVAL:
                interval(src.get("median_interval"), nullable=False); interval(src.get("mean_interval"), nullable=False)
                require(src.get("replicates") == 5000 and src.get("seed") == 2026091601, "Response interval rule changed")
            else: require(src.get("median_interval") is None and src.get("mean_interval") is None, "Unavailable response uncertainty cannot be numeric")
        for field in ("mean_interval_status", "median_interval_status"):
            require(src.get(field, row["interval_status"]) == row["interval_status"], "Unrecognized per-statistic interval status")
            row[field] = row["interval_status"]
        require(src.get("replicates") in (None, 5000) and src.get("seed") in (None, 2026091601), "Unexpected response bootstrap metadata")
        row["bootstrap_replicates"] = integer(src.get("replicates"), nullable=True)
        row["bootstrap_seed"] = integer(src.get("seed"), nullable=True)
        if p == "pooled_descriptive" or not available or row["interval_status"] == "insufficient_runs":
            require(row["bootstrap_replicates"] is None and row["bootstrap_seed"] is None, "Unavailable/pooled uncertainty cannot acquire bootstrap metadata")
        if not available and p != "pooled_descriptive":
            require(row["interval_status"] == "not_available", "Failed response cannot claim computed uncertainty")
        rows.append(row)
    for t, s in itertools.product(TIMES, STIMULI):
        total = tuple(sum(counts[(p, t, s)][i] for p in PROTOCOLS) for i in (0, 1))
        require(counts[("pooled_descriptive", t, s)] == total, "Pooled response original counts differ")
        if s != "PIC":
            for g, q in itertools.product(PROGRAMS, QUANTITIES):
                pooled = index[(scope, "pooled_descriptive", t, s, g, q)]
                for p in PROTOCOLS:
                    require(close(pooled["original_contrast_protocol_line_shares"][p], counts[(p, t, s)][0] / total[0]), "Pooled response shares changed")
    return rows, counts


def prediction_table(endpoint, scope, variant):
    st = status(endpoint["status"])
    require(st in ENDPOINT_STATUS, "Wrong endpoint status layer")
    allowed = {"status", "primary_gain", "pooled", "cells", "line_counts", "protocol_weights", "interval_status", "interval", "cell_intervals", "replicates", "seed"} if st == "available" else {"status", "primary_gain"}
    require(set(endpoint) <= allowed, "Unexpected endpoint metadata field")
    ns = COUNTS[scope]; total = sum(ns.values()); available = st == "available"
    cellkeys = {str((p, t)): (p, t) for p, t in itertools.product(PROTOCOLS, TIMES)}
    if available:
        require(endpoint.get("line_counts") == ns, "Prediction original cohort changed")
        require(set(endpoint.get("protocol_weights", {})) == set(PROTOCOLS), "Missing prediction protocol shares")
        for p in PROTOCOLS: require(close(endpoint["protocol_weights"][p], ns[p] / total), "Prediction original protocol shares changed")
        require(set(endpoint.get("cells", {})) == set(cellkeys), "Missing/extra prediction cells")
        require(endpoint.get("interval_status") in (PREDICTION_INTERVAL, "insufficient_runs"), "Missing conditional prediction interval status")
        if endpoint["interval_status"] == PREDICTION_INTERVAL:
            require(set(endpoint.get("cell_intervals", {})) == set(cellkeys), "Mandatory conditional prediction CIs missing")
            require(endpoint.get("replicates") == 5000 and endpoint.get("seed") == 2026091601, "Prediction CI rule changed")
        else:
            require(endpoint.get("interval") is None and not endpoint.get("cell_intervals") and endpoint.get("replicates") is None and endpoint.get("seed") is None, "Unavailable predictive CI cannot be numeric")
        require(close(endpoint["primary_gain"], endpoint["pooled"]["gain"]), "Primary/pooled gain differs")
        require(set(endpoint["pooled"]) == {"baseline_loss", "extended_loss", "gain"} and all(set(r) == {"baseline_loss", "extended_loss", "gain"} for r in endpoint["cells"].values()), "Unexpected loss summary fields")
        for field in ("baseline_loss", "extended_loss", "gain"):
            pooled = sum(ns[p] / total * sum(endpoint["cells"][str((p, t))][field] for t in TIMES) / 2 for p in PROTOCOLS)
            require(close(endpoint["pooled"][field], pooled), "Original line/time weighting changed")
    else:
        require(endpoint.get("primary_gain") is None and not endpoint.get("cells") and not endpoint.get("pooled") and endpoint.get("interval") is None, "Unavailable endpoint contains effects")
    require(endpoint.get("interval_status", "not_available") in (PREDICTION_INTERVAL, "insufficient_runs", "not_available"), "Unknown prediction interval status")
    require(endpoint.get("replicates") in (None, 5000) and endpoint.get("seed") in (None, 2026091601), "Unexpected prediction bootstrap metadata")
    integer(endpoint.get("replicates"), nullable=True); integer(endpoint.get("seed"), nullable=True)
    if not available:
        require(endpoint.get("replicates") is None and endpoint.get("seed") is None and endpoint.get("interval_status") is None and not endpoint.get("cell_intervals"), "Failed endpoint acquired uncertainty data")
    rows = []
    for ck, (p, t) in [*cellkeys.items(), ("pooled", ("pooled", None))]:
        source = endpoint.get("pooled", {}) if p == "pooled" else endpoint.get("cells", {}).get(ck, {})
        row = {**labels(scope, variant), "protocol": p, "time": t, "status": st,
               "is_primary_endpoint": scope == SCOPES[0] and variant == VARIANTS[0] and p == "pooled",
               "original_lines": total if p == "pooled" else ns[p], "original_total_lines": total,
               "original_Mod_SmartSeq2_lines": ns[PROTOCOLS[0]], "original_NEB_lines": ns[PROTOCOLS[1]],
               "original_Mod_SmartSeq2_share": ns[PROTOCOLS[0]] / total, "original_NEB_share": ns[PROTOCOLS[1]] / total,
               "time_weight_6h": 0.5, "time_weight_24h": 0.5,
               "interval_status": endpoint.get("interval_status", "not_available"),
               "bootstrap_replicates": endpoint.get("replicates"), "bootstrap_seed": endpoint.get("seed")}
        for field in ("baseline_loss", "extended_loss", "gain"): row[field] = number(source.get(field), nullable=not available)
        if available:
            require(row["baseline_loss"] >= 0 and row["extended_loss"] >= 0 and close(row["gain"], row["baseline_loss"] - row["extended_loss"]), "Invalid prediction loss/gain")
        ci = endpoint.get("interval") if p == "pooled" else endpoint.get("cell_intervals", {}).get(ck)
        row["gain_interval_low"], row["gain_interval_high"] = interval(ci, nullable=not available or endpoint.get("interval_status") == "insufficient_runs")
        rows.append(row)
    return rows


def influence_table(records, endpoint, scope, variant, expected_runs):
    require(isinstance(records, list) and bool(records), "Missing leave-run diagnostic status")
    available = endpoint["status"] == "available"
    if available:
        ids = []
        for row in records:
            require(row.get("protocol") in PROTOCOLS and isinstance(row.get("omitted_run"), str), "Invalid private leave-run key")
            require(row.get("diagnostic") == "leave_one_run_out_fixed_predictions", "Unexpected leave-run diagnostic")
            require(row["status"] in ENDPOINT_STATUS, "Wrong leave-run aggregate status layer")
            if row["status"] == "available": number(row.get("gain"))
            else: require(row.get("gain") is None, "Unavailable leave-run effect cannot be zero")
            ids.append((row["protocol"], row["omitted_run"]))
        require(len(set(ids)) == len(ids), "Duplicate omitted-run diagnostic")
        for p in PROTOCOLS: require(sum(r["protocol"] == p for r in records) == expected_runs[p], "Dropped leave-run diagnostic")
    else:
        require(len(records) == 1 and records[0].get("status") == endpoint["status"] and records[0].get("diagnostic") == "unavailable_incomplete_predictions" and records[0].get("gain") is None, "Missing unavailable leave-run status")
    result = []
    for p in (*PROTOCOLS, "pooled"):
        chosen = [r for r in records if p == "pooled" or r.get("protocol") == p] if available else []
        values = [r["gain"] for r in chosen if r["status"] == "available"]
        expected = sum(expected_runs.values()) if p == "pooled" else expected_runs[p]
        status_counts = {key: sum(r["status"] == key for r in chosen) for key in sorted({r["status"] for r in chosen})} if available else {endpoint["status"]: 1}
        result.append({**labels(scope, variant), "omission_protocol": p, "status": "available" if available else endpoint["status"], "diagnostic_status_counts": canonical(status_counts), "diagnostic_quantity": "pooled_predictive_gain_after_single_run_omission",
            "expected_omissions": expected, "source_diagnostic_records": len(chosen) if available else 1,
            "available_omissions": len(values), "unavailable_omissions": expected - len(values),
            "minimum_gain": min(values) if values else None, "maximum_gain": max(values) if values else None,
            "positive": sum(x > 0 for x in values), "negative": sum(x < 0 for x in values), "zero": sum(x == 0 for x in values),
            "weights": "fixed_original_protocol_line_shares_and_equal_times", "interpretation": "fixed_prediction_influence_not_independent_replicates"})
    return result


class Reader:
    def __init__(self, root, plan_path, plan_sha, complete_sha, receipt_path, receipt_sha, public_path, public_sha, replay_path=None, replay_sha=None, optimizer_path=None, optimizer_sha=None):
        self.root = absolute(root); self.used = {}; self.inputs = []
        def read(path, sha):
            self.inputs.append((absolute(path), sha))
            return strict_json(read_pin(path, sha))
        self.plan = read(plan_path, plan_sha)
        require(self.plan.get("status") == "frozen" and self.plan.get("execute_real_expression") is True and self.plan.get("analysis_id") == "macromap-fixed-program-context-v1", "Unfrozen/foreign plan")
        contract = self.plan["method_contract"]
        expected = {"identity_scopes": list(SCOPES), "protocols": list(PROTOCOLS), "times": list(TIMES), "programs": list(PROGRAMS), "predictive_variants": list(VARIANTS), "folds": 5, "pooled_response_mean_interval": False, "pooled_response_median_interval": False, "bootstrap_replicates": 5000, "bootstrap_seed": 2026091601, "time_weights": [0.5, 0.5]}
        for k, v in expected.items(): require(contract.get(k) == v, "Unexpected frozen reporting contract")
        self.directory = contained(self.root, self.plan["output_directory"])
        self.prepared = contained(self.root, self.plan["prepared_directory"])
        self.complete = read(self.directory / "execution-complete.json", complete_sha)
        require(self.complete.get("stage") == "execution_complete_all_prespecified_statuses_retained" and self.complete.get("plan_sha256") == plan_sha, "Incomplete/foreign completion")
        require(not (self.directory / "execution-failed.json").exists(), "Conflicting failed/completed execution")
        for k in ("script_sha256", "runner_script_sha256", "runtime", "prepared_summary_sha256", "method_contract"):
            require(self.complete.get(k) == self.plan.get(k), "Plan/completion contract differs")
        require(set(self.complete.get("endpoints", {})) == set(SCOPES) and all(set(v) == set(VARIANTS) for v in self.complete["endpoints"].values()), "Missing/extra completion endpoint scopes/variants")
        receipt = read(receipt_path, receipt_sha)
        require(receipt.get("status") == VERIFIED and receipt.get("plan_sha256") == plan_sha and receipt.get("execution_receipt_sha256") == complete_sha and receipt.get("verification_is_read_only") is True and receipt.get("effects_read") is True, "Independent verification did not pass for these outputs")
        require(receipt.get("fit_cells") == 80 and len(receipt.get("fit_checks", [])) == 80 and len(receipt.get("endpoints", [])) == 4 and len(receipt.get("responses", [])) == 2 and len(receipt.get("source_score_contexts", [])) == 8 and receipt.get("scope_limits") == SCOPE_LIMITS, "Independent receipt coverage differs")
        approved_public = read(public_path, public_sha)
        require(approved_public == receipt, "Root-approved public/private independent receipts differ")
        code_pins = self.plan.get("code_and_protocol_artifacts", [])
        require(isinstance(code_pins, list), "Invalid frozen code inventory")
        verifier_pins = [r["sha256"] for r in code_pins if r.get("path") == "scripts/verify_macromap_program_context.py"]
        require(len(verifier_pins) == 1, "Missing unique frozen verifier pin")
        self.verifier_sha = verifier_pins[0]
        require(isinstance(self.verifier_sha, str) and len(self.verifier_sha) == 64 and all(c in "0123456789abcdef" for c in self.verifier_sha), "Invalid frozen verifier hash")
        normalization = approved_public["normalization"]
        require(normalization.get("status") == "normalized_source_subset_verified" and normalization.get("gene_count") == 58243 and normalization.get("sample_count") == 8, "Unexpected normalized subset coverage")
        normalization_error = number(normalization["max_absolute_error"])
        require(normalization_error >= 0, "Invalid normalization error")
        responses = unique_rows(approved_public["responses"], ("identity_scope",), [(s,) for s in SCOPES])
        require(all(r["all_points_statuses_checked"] is True and r["complete_cells"] == 1512 and r["pooled_response_intervals"] is False and r["protocol_CIL6h_interval_cells_checked"] == 54 for r in responses.values()), "Public response verification coverage differs")
        self.public = {"schema": "macromap-public-verification-summary-v1", "status": VERIFIED,
            "plan_sha256": plan_sha, "execution_receipt_sha256": complete_sha, "verification_receipt_sha256": receipt_sha,
            "root_approved_public_receipt_sha256": public_sha, "verifier_sha256": self.verifier_sha,
            "aggregate_checks": {"fit_cells": 80, "source_score_contexts": len(approved_public["source_score_contexts"]),
                "predictive_endpoints": 4, "response_panels": 2, "response_panel_rows": 3024,
                "normalization_samples": 8, "normalization_genes": 58243, "response_CIL6h_interval_cells": 108,
                "all_response_points_statuses_checked": True},
            "normalization_max_absolute_error": normalization_error, "scope_limits": list(SCOPE_LIMITS)}
        self.summary = read(self.prepared / "source-only-readiness.json", self.plan["prepared_summary_sha256"])
        require(self.summary["artifact_sha256"] == self.plan["preparation_artifacts"], "Prepared inventory differs from plan")
        require(self.complete["source_pins"] == self.summary["source_pins"], "Source pins differ")
        self.pins = {"plan_sha256": plan_sha, "execution_receipt_sha256": complete_sha, "independent_verification_receipt_sha256": receipt_sha, "public_verification_summary_sha256": public_sha, "prepared_summary_sha256": self.plan["prepared_summary_sha256"], "verifier_sha256": self.verifier_sha}
        require((replay_path is None) == (replay_sha is None), "Supply both replay receipt and external hash")
        if replay_path is not None:
            replay = read(replay_path, replay_sha)
            require(replay.get("status") == "all_266_non_start_artifacts_byte_identical_expected_run_metadata_differences_only" and replay.get("source_plan_sha256") == plan_sha and replay.get("primary_complete_sha256") == complete_sha, "Mechanical replay receipt differs")
            require(replay.get("all_endpoints_identical") is True and replay.get("no_frozen_scientific_change") is True and replay.get("biological_replication") is False and replay.get("allowed_plan_difference") == ["output_directory"], "Mechanical replay scope changed")
            checks = replay.get("file_checks", [])
            inventory = self.complete["artifact_sha256"]
            require(len(checks) == len(inventory) and {r.get("file") for r in checks} == set(inventory), "Replay receipt omitted/duplicated artifact")
            differences = []
            for check in checks:
                require(check.get("primary_sha256") == inventory[check["file"]] and type(check.get("identical")) is bool, "Replay primary artifact pin differs")
                require(check["identical"] == (check["primary_sha256"] == check.get("replay_sha256")), "Replay identity flag differs")
                if not check["identical"]: differences.append(check["file"])
            require(differences == ["execution-start.json"] and replay.get("files_rehashed_per_execution") == len(inventory) and replay.get("byte_identical_files") == len(inventory)-1, "Unexpected mechanical replay artifact difference")
            self.pins["root_mechanical_replay_receipt_sha256"] = replay_sha
            self.public["mechanical_replay"] = {"status": replay["status"], "files_rehashed_per_execution": len(inventory), "total_files_rehashed": 2 * len(inventory), "byte_identical_files": len(inventory)-1, "all_endpoints_identical": True, "biological_replication": False, "allowed_difference": "execution-start metadata only; no scientific change"}

        require((optimizer_path is None) == (optimizer_sha is None), "Supply optimizer diagnostic and its external hash together")
        self.optimizer_qc = read(optimizer_path, optimizer_sha) if optimizer_path is not None else None
        if optimizer_path is not None: self.pins["approved_optimizer_diagnostic_sha256"] = optimizer_sha

    def aggregate(self, name):
        # No method accepts arbitrary payload paths: raw, normalized, NPZ and
        # per-line files are never read or even hashed by this reporter.
        allowed = {"fit-manifest.json"}
        for s in SCOPES:
            allowed.update(f"{s}/{f}.json" for f in ("response-panel", "fold-diagnostics"))
            for v in VARIANTS: allowed.update(f"{s}/{v}/{f}.json" for f in ("incremental-endpoint", "run-influence"))
            for p, t, f in itertools.product(PROTOCOLS, TIMES, range(5)):
                context = f"{s}/{p}-t{t}-fold{f}"
                allowed.add(context + "/score-reference-audit.json")
                for v in VARIANTS: allowed.add(context + f"/{v}/fit-audit.json")
        require(name in allowed, "Not an approved aggregate artifact")
        sha = self.complete["artifact_sha256"].get(name)
        path = contained(self.directory, name)
        self.used[name] = sha; self.inputs.append((path, sha))
        return strict_json(read_pin(path, sha))

    def aggregate_csv(self, name):
        require(name in ("mapping-inventory.csv", "paired-design-aggregate.csv"), "Not an approved prepared aggregate")
        sha = self.summary["artifact_sha256"][name]; path = contained(self.prepared, name)
        self.used["preparation/" + name] = sha; self.inputs.append((path, sha))
        reader = csv.DictReader(io.StringIO(read_pin(path, sha).decode("utf-8")))
        require(reader.fieldnames is not None and len(reader.fieldnames) == len(set(reader.fieldnames)), "Duplicate CSV columns")
        rows = list(reader)
        require(all(None not in row and None not in row.values() for row in rows), "Ragged CSV aggregate")
        return rows

    def recheck(self):
        for path, sha in self.inputs: read_pin(path, sha)


def optimizer_projection(reader, tables):
    """Project an approved aggregate diagnostic; never reopen NPZ or rerun fits."""
    qc = reader.optimizer_qc
    require(qc.get("schema") == "macromap-independent-optimizer-flag-capture-v1" and qc.get("status") == "captured_frozen_numeric_criteria_pass" and qc.get("purpose") == "post_fit_flag_capture_not_replacement_verification" and qc.get("existing_numerical_verification_status") == "unchanged", "Unapproved optimizer diagnostic state")
    require(set(qc) == {"created_at_utc", "existing_numerical_verification_status", "limitations", "metric_definitions", "provenance", "purpose", "rows", "schema", "status", "summary", "unchanged_criteria"}, "Unknown optimizer diagnostic fields")
    criteria = {"coefficient_atol": 2e-5, "coefficient_rtol": 2e-5, "constant_training_sd_threshold": 1e-12,
        "independent_gradient_limit": 1e-7, "method": "frozen independent_softmax_fit defaults",
        "raw_solver_success_is_separate_from_numeric_acceptance": True, "scaler_atol": 1e-12, "scaler_rtol": 1e-10}
    require(qc["unchanged_criteria"] == criteria, "Optimizer diagnostic thresholds/method changed")
    provenance = qc["provenance"]
    require(provenance.get("plan_sha256") == reader.pins["plan_sha256"] and provenance.get("primary_completion_sha256") == reader.pins["execution_receipt_sha256"] and provenance.get("frozen_verifier_sha256") == reader.verifier_sha and provenance.get("fit_manifest_sha256") == reader.complete["artifact_sha256"]["fit-manifest.json"], "Optimizer diagnostic is from another execution")
    require(provenance.get("runtime") == reader.plan["runtime"], "Optimizer diagnostic numerical runtime differs")
    wrapper_sha = provenance.get("diagnostic_wrapper_sha256")
    require(isinstance(wrapper_sha, str) and len(wrapper_sha) == 64 and all(c in "0123456789abcdef" for c in wrapper_sha), "Invalid diagnostic wrapper hash")
    pair_keys = ("identity_scope", "variant", "protocol", "time", "fold")
    expected_pairs = list(itertools.product(SCOPES, VARIANTS, PROTOCOLS, TIMES, range(5)))
    npz_pins = unique_rows(provenance["fit_npz_pins"], pair_keys, expected_pairs)
    manifest = unique_rows(reader.aggregate("fit-manifest.json"), pair_keys, expected_pairs)
    # Compare authenticated commitments only. No path below is opened/hashed.
    for key in expected_pairs:
        pin, entry = npz_pins[key], manifest[key]
        require(set(pin) == set(pair_keys) | {"sha256"}, "Unexpected model commitment fields")
        require(pin["sha256"] == entry["fit_npz_sha256"] == reader.complete["artifact_sha256"].get(entry["fit_npz"]), "Optimizer diagnostic model commitment differs")
    keys = (*pair_keys, "model")
    expected = [(*key, model) for key in expected_pairs for model in ("baseline", "extended")]
    index = unique_rows(qc["rows"], keys, expected)
    original = unique_rows(tables["model-fits"], keys, expected)
    qc_fields = set(keys) | {"coefficient_agreement", "coefficient_max_absolute_error", "constant_training_feature_count", "constant_training_feature_names", "frozen_numeric_criteria_pass", "gradient_inf", "iterations", "objective", "optimizer_success", "saved_constant_flags_verified", "solver_flag_status", "training_feature_count", "training_rows"}
    output, constant_name_counts = [], {}
    for key in expected:
        src, fit = index[key], original[key]
        require(set(src) == qc_fields, "Unexpected optimizer diagnostic row fields")
        scope, variant, p, t, f, model = key
        programs = ["inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling"] + ([variant] if model == "extended" else [])
        for field in ("optimizer_success", "coefficient_agreement", "frozen_numeric_criteria_pass", "saved_constant_flags_verified"):
            require(type(src[field]) is bool, "Optimizer flags must be boolean")
        require(src["solver_flag_status"] == ("reported_success" if src["optimizer_success"] else "reported_unsuccessful"), "Raw optimizer flag/status differs")
        require(integer(src["training_rows"]) == fit["train_pairs"] and integer(src["training_feature_count"]) == len(programs), "Optimizer diagnostic training dimensions differ")
        count = integer(src["constant_training_feature_count"])
        names = src["constant_training_feature_names"]
        require(isinstance(names, list) and all(type(name) is str for name in names) and len(set(names)) == len(names) == count and set(names) <= set(programs), "Constant-feature names are not fixed predictor labels")
        require(names == [name for name in programs if name in names], "Constant-feature names/order disagree with model mask")
        for name in names: constant_name_counts[name] = constant_name_counts.get(name, 0) + 1
        gradient = number(src["gradient_inf"]); objective = number(src["objective"]); error = number(src["coefficient_max_absolute_error"])
        require(0 <= gradient <= criteria["independent_gradient_limit"] and objective >= 0 and error >= 0, "Independent numeric diagnostic outside frozen criteria")
        require(src["coefficient_agreement"] is True and src["frozen_numeric_criteria_pass"] is True and src["saved_constant_flags_verified"] is True, "Independent numeric/constant criterion did not pass")
        output.append({**labels(scope, variant), "protocol": p, "time": t, "fold": f, "model": model,
            "diagnostic_scope": "same_frozen_independent_helper_rerun_not_original_SciPy_object",
            "optimizer_success": src["optimizer_success"], "solver_flag_status": src["solver_flag_status"],
            "iterations": integer(src["iterations"]), "gradient_inf": gradient, "objective": objective,
            "coefficient_max_absolute_error": error, "coefficient_agreement": src["coefficient_agreement"],
            "frozen_numeric_criteria_pass": src["frozen_numeric_criteria_pass"], "training_rows": src["training_rows"],
            "training_feature_count": len(programs), "constant_training_feature_count": count,
            "constant_training_feature_names": canonical(names), "saved_constant_flags_verified": src["saved_constant_flags_verified"],
            "optimizer_message_or_status_code": "not_returned_by_frozen_helper_no_reason_invented"})
    by_pair = {tuple(row[k] for k in pair_keys): row for row in tables["scale-status"]}
    for key in expected_pairs:
        baseline, extended = index[(*key, "baseline")], index[(*key, "extended")]
        shared = [name for name in extended["constant_training_feature_names"] if name != key[1]]
        require(baseline["constant_training_feature_names"] == shared and baseline["training_rows"] == extended["training_rows"], "Baseline/extended constant-feature masks disagree")
        by_pair[key].update(scale_status="constant_flags_independently_verified", numeric_scaler_audit="constant_flags_only_from_approved_aggregate_diagnostic",
            evidence="approved_same_helper_rerun_aggregate_no_NPZ_read", constant_feature_flags_exported=True,
            constant_training_feature_count=extended["constant_training_feature_count"], constant_training_feature_names=canonical(extended["constant_training_feature_names"]),
            baseline_constant_training_feature_count=baseline["constant_training_feature_count"], baseline_constant_training_feature_names=canonical(baseline["constant_training_feature_names"]),
            saved_constant_flags_verified=True)
    summary = {"fit_cells": 80, "models": 160, "optimizer_success_true": sum(r["optimizer_success"] for r in output),
        "optimizer_success_false": sum(not r["optimizer_success"] for r in output), "frozen_numeric_criteria_pass": sum(r["frozen_numeric_criteria_pass"] for r in output),
        "max_gradient_inf": max(r["gradient_inf"] for r in output), "max_coefficient_absolute_error": max(r["coefficient_max_absolute_error"] for r in output),
        "models_with_constant_training_features": sum(r["constant_training_feature_count"] > 0 for r in output),
        "constant_training_feature_name_counts": constant_name_counts, "all_saved_constant_flags_verified": True}
    require(qc["summary"] == summary, "Optimizer diagnostic summary dropped/changed rows or flags")
    reader.public["optimizer_diagnostic"] = {"schema": qc["schema"], "status": qc["status"], "purpose": qc["purpose"], "existing_numerical_verification_status": "unchanged",
        "diagnostic_wrapper_sha256": wrapper_sha, "summary": summary,
        "interpretation": "raw_solver_flags_separate_from_frozen_numeric_pass; identical_helper_reruns_not_original_objects; no_failure_reasons_returned"}
    tables["independent-optimizer-qc"] = output


def assemble(reader):
    tables = {k: [] for k in ("responses", "predictions", "leave-run-summary", "model-folds", "model-fits", "matching-status", "scale-status", "mapping", "paired-design", "cohort-design")}
    panels, fold_index, endpoints = {}, {}, {}
    for scope in SCOPES:
        rows, panels[scope] = response_table(reader.aggregate(f"{scope}/response-panel.json"), scope)
        tables["responses"].extend(rows)
        folds = reader.aggregate(f"{scope}/fold-diagnostics.json")
        keys = ("identity_scope", "protocol", "time", "fold", "predictor_program")
        expected = list(itertools.product((scope,), PROTOCOLS, TIMES, range(5), VARIANTS))
        idx = unique_rows(folds, keys, expected)
        for key in expected:
            src = idx[key]; _, p, t, f, variant = key
            require(set(src) == {"identity_scope", "protocol", "time", "fold", "predictor_program", "train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs", "train_class_lines", "train_class_runs", "missing_test_stimuli", "fold_gate_status", "status", "line_weighted_gain"}, "Unexpected fold audit fields")
            require(src["status"] in OUTER_STATUS and src["fold_gate_status"] in GATE_STATUS, "Wrong fold/gate status layer")
            row = {**labels(scope, variant), "protocol": p, "time": t, "fold": f, "status": status(src["status"]), "fold_gate_status": status(src["fold_gate_status"])}
            for name in ("train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs"): row[name] = integer(src[name])
            for name in ("train_class_lines", "train_class_runs"):
                require(set(src[name]) == set(STIMULI[:-1]), "Missing fixed class counts")
                for stimulus in STIMULI[:-1]: row[name + "_" + stimulus] = integer(src[name][stimulus])
            require(isinstance(src["missing_test_stimuli"], list) and set(src["missing_test_stimuli"]) <= set(STIMULI[:-1]), "Unknown test class")
            row["missing_test_stimuli"] = canonical(src["missing_test_stimuli"])
            row["line_weighted_gain"] = number(src.get("line_weighted_gain"), nullable=src["status"] != "available")
            if src["status"] != "available": require(row["line_weighted_gain"] is None, "Failed fold has a zero/other gain")
            tables["model-folds"].append(row); fold_index[(scope, variant, p, t, f)] = row
        for variant in VARIANTS:
            endpoint = reader.aggregate(f"{scope}/{variant}/incremental-endpoint.json")
            require(endpoint == reader.complete["endpoints"][scope][variant], "Endpoint differs from completion")
            endpoints[(scope, variant)] = endpoint
            tables["predictions"].extend(prediction_table(endpoint, scope, variant))
            expected_runs = {}
            for p, t in itertools.product(PROTOCOLS, TIMES):
                subset = [fold_index[(scope, variant, p, t, f)] for f in range(5)]
                require(sum(r["test_lines"] for r in subset) == COUNTS[scope][p], "Original fold test line total differs")
                nr = sum(r["test_runs"] for r in subset)
                require(p not in expected_runs or expected_runs[p] == nr, "Run counts differ across times")
                expected_runs[p] = nr
                if endpoint["status"] == "available":
                    require(all(r["status"] == "available" for r in subset), "Available endpoint silently dropped failed fold")
                    gain = sum(r["line_weighted_gain"] * r["test_lines"] for r in subset) / COUNTS[scope][p]
                    require(close(gain, endpoint["cells"][str((p, t))]["gain"]), "Fold/cell aggregate gain differs")
            tables["leave-run-summary"].extend(influence_table(reader.aggregate(f"{scope}/{variant}/run-influence.json"), endpoint, scope, variant, expected_runs))
    manifest = reader.aggregate("fit-manifest.json")
    expected = list(itertools.product(SCOPES, VARIANTS, PROTOCOLS, TIMES, range(5)))
    idx = unique_rows(manifest, ("identity_scope", "variant", "protocol", "time", "fold"), expected)
    contexts = set()
    for key in expected:
        scope, variant, p, t, f = key; entry = idx[key]
        require(set(entry) == {"identity_scope", "variant", "protocol", "time", "fold", "status", "fit_npz", "fit_npz_sha256", "fit_audit", "stratum_npz", "stratum_audit"}, "Unexpected manifest fields")
        context = f"{scope}/{p}-t{t}-fold{f}"
        require(entry["fit_audit"] == context + f"/{variant}/fit-audit.json" and entry["stratum_audit"] == context + "/score-reference-audit.json", "Unexpected audit path/identity")
        fit = reader.aggregate(entry["fit_audit"])
        require(set(fit) == {"schema", "status", "fold_gate", "identity_scope", "variant", "protocol", "time", "fold", "feature_programs", "class_order", "penalty", "loss_normalization", "intercepts", "fits", "preprocessing_scope"}, "Unexpected fit audit fields")
        require(set(fit["fold_gate"]) == {"status", "train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs", "train_class_lines", "train_class_runs", "missing_test_stimuli"}, "Unexpected fit gate fields")
        require(tuple(fit[k] for k in ("identity_scope", "variant", "protocol", "time", "fold")) == key, "Fit audit identity differs")
        diag = fold_index[key]
        require(fit["status"] == entry["status"] == diag["status"], "Failed fit status differs across aggregates")
        require(fit.get("schema") == "macromap-private-audit-v1" and fit.get("penalty") == 0.1 and fit.get("loss_normalization") == "sum_weights_one_equal_lines" and fit.get("intercepts") == "unpenalized_symmetric_softmax" and fit.get("preprocessing_scope") == "both_train_and_test_X_use_this_outer_fold_training_control_reference", "Fit audit method differs")
        for field in ("train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs"):
            require(fit["fold_gate"][field] == diag[field], "Fit/gate original counts differ")
        require(fit["fold_gate"]["status"] == diag["fold_gate_status"], "Fit/gate status differs")
        for field in ("train_class_lines", "train_class_runs"):
            require(fit["fold_gate"][field] == {s: diag[field + "_" + s] for s in STIMULI[:-1]}, "Fit/gate class counts differ")
        require(fit["fold_gate"]["missing_test_stimuli"] == strict_json(diag["missing_test_stimuli"]), "Fit/gate missing classes differ")
        require(fit["feature_programs"] == ["inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling", variant] and fit["class_order"] == list(STIMULI[:-1]), "Fit feature/class identity differs")
        fits = fit.get("fits", {})
        require(set(fits) <= {"baseline", "extended"}, "Unexpected extra model")
        require(set(fits) == {"baseline", "extended"}, "Paired fit omitted model/status")
        if fit["status"] == "available": require(all(r.get("status") == "available" for r in fits.values()), "Available outer fit contains failed component")
        base = {**labels(scope, variant), "protocol": p, "time": t, "fold": f, "paired_status": status(fit["status"])}
        for model in ("baseline", "extended"):
            one = fits.get(model, {})
            require(set(one) in ({"status"}, {"status", "iterations", "gradient_inf", "objective", "optimizer_message"}), "Unexpected component diagnostics fields")
            require(one.get("status") in COMPONENT_STATUS, "Wrong component status layer")
            row = {**base, "model": model, "status": status(one.get("status", "not_fitted")), "predictors": 4 if model == "baseline" else 5}
            for field in ("train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs"): row[field] = diag[field]
            row["iterations"] = integer(one.get("iterations"), nullable=True)
            row["objective"] = number(one.get("objective"), nullable=True)
            row["gradient_inf"] = number(one.get("gradient_inf"), nullable=True)
            require(row["gradient_inf"] is None or row["gradient_inf"] >= 0, "Negative gradient norm")
            if row["status"] == "available": require(row["objective"] is not None and row["gradient_inf"] is not None and 0 <= row["gradient_inf"] <= 1e-6, "Invalid available optimizer diagnostics")
            tables["model-fits"].append(row)
        tables["scale-status"].append({**base, "expected_predictors": 5, "scale_status": "inferred_from_computed_fit_diagnostics" if any("iterations" in v for v in fits.values()) else "not_recorded_in_aggregate_audit", "numeric_scaler_audit": "not_available_in_allowed_aggregate_inputs", "method": "equal_line_weight_training_fold_paired_changes_only", "evidence": "paired_fit_audit_no_scale_arrays_read", "numeric_scale_values_exported": False, "constant_feature_flags_exported": False})
        if context not in contexts:
            contexts.add(context); audit = reader.aggregate(entry["stratum_audit"])
            require(set(audit) == {"scope", "protocol", "time", "heldout_fold", "program_status", "reference_error", "reference", "matching"}, "Unexpected matching audit fields")
            require(set(audit["reference"]) == {"reference_lines", "reference_runs", "control_indices"}, "Unexpected reference audit fields")
            require((audit["scope"], audit["protocol"], audit["time"], audit["heldout_fold"]) == (scope, p, t, f), "Matching context identity differs")
            require(set(audit["program_status"]) == set(PROGRAMS) and set(audit["matching"]) == set(PROGRAMS), "Missing matching program status")
            reference = audit.get("reference", {})
            for program in PROGRAMS:
                info = audit["matching"][program]
                base_fields = {"available", "reason", "seed", "target_bin_counts", "pool_bin_counts", "underflow_bins"}
                require(set(info) in ({"available", "reason"}, base_fields, base_fields | {"selection_sha256", "baseline_target", "baseline_background"}), "Unexpected matching entry fields")
                require(audit["program_status"][program] in MATCH_STATUS | {"mapping_failed"}, "Wrong program status layer")
                pools = info.get("pool_bin_counts", {}); targets = info.get("target_bin_counts", {}); under = info.get("underflow_bins", [])
                require(all(str(k).isdigit() and 0 <= int(k) < 20 for k in pools) and all(str(k).isdigit() and 0 <= int(k) < 20 for k in targets), "Unexpected bin key")
                for count in [*pools.values(), *targets.values()]: integer(count)
                require(isinstance(under, list) and all(type(x) is int and 0 <= x < 20 for x in under), "Unexpected underflow bins")
                require(type(info.get("available")) is bool, "Missing matching availability")
                require(info.get("reason") in ("available", "target_outside_bins", "empty_target", "insufficient_control_pool", "reference_unavailable"), "Unknown matching reason")
                require((info["reason"] == "available") == info["available"], "Matching reason/availability differs")
                require(audit["program_status"][program] in ("mapping_failed", info["reason"]), "Matching failure was dropped from program status")
                tables["matching-status"].append({**labels(scope), "protocol": p, "time": t, "fold": f, "program_id": program, "status": status(audit["program_status"][program]), "matching_available": info["available"], "matching_reason": info["reason"], "reference_status": "available" if audit.get("reference_error") is None else "failed", "reference_lines": integer(reference.get("reference_lines"), nullable=True), "reference_runs": integer(reference.get("reference_runs"), nullable=True), "target_bins": len(targets), "underflow_bin_count": len(under), "minimum_pool_bin_size": min(pools.values()) if pools else None})
    mapping = reader.aggregate_csv("mapping-inventory.csv")
    require(len(mapping) == len(PROGRAMS) and {r["program_id"] for r in mapping} == set(PROGRAMS), "Mapping inventory keys differ")
    integer_fields = ("source_unique_before_exclusion", "excluded_target_gene_count", "intended_unique_source_genes", "mapped_unique_genes", "missing_source_genes", "ambiguous_source_genes", "parent_mapped_unique_genes", "comparator_union_mapped_genes", "parent_comparator_overlap_genes", "derived_nonoverlap_mapped_genes")
    for g in PROGRAMS:
        src = next(r for r in mapping if r["program_id"] == g)
        require(set(src) == set(COLUMNS["mapping"]) | {"dataset_key", "role"} and src["dataset_key"] == "macromap", "Unexpected mapping aggregate schema")
        row = {"program_id": g}
        for field in integer_fields:
            try:
                value = Decimal(src[field]) if src.get(field) else None
            except InvalidOperation as exc:
                raise ValueError("Invalid mapping count token") from exc
            require(value is None or (value.is_finite() and 0 <= value < 2**53 and value == value.to_integral_value()), "Invalid mapping count")
            row[field] = int(value) if value is not None else None
        for field in ("mapping_fraction", "derived_mapping_fraction", "parent_mapping_fraction"): row[field] = number(float(src[field]), nullable=False) if src.get(field) else None
        for field in ("mapping_gate_pass", "parent_mapping_gate_pass"):
            require(src.get(field) in ("True", "False", ""), "Invalid mapping gate")
            row[field] = {"True": True, "False": False, "": None}[src[field]]
        row["parent_program_id"] = src.get("parent_program_id") or None
        require(row["parent_program_id"] in (None, VARIANTS[0]), "Unexpected parent program")
        tables["mapping"].append(row)
    design = reader.aggregate_csv("paired-design-aggregate.csv")
    require(all(set(r) == {"protocol", "time", "stimulus", "paired_lines", "status"} for r in design), "Unexpected paired-design schema")
    design_index = unique_rows([{**r, "time": int(r["time"])} for r in design], ("protocol", "time", "stimulus"), itertools.product(PROTOCOLS, TIMES, STIMULI))
    for scope, p, t, s in itertools.product(SCOPES, PROTOCOLS, TIMES, STIMULI):
        n, runs = panels[scope][(p, t, s)]
        original = design_index[(p, t, s)]
        if scope == SCOPES[0]: require(n == int(original["paired_lines"]), "Response/source design counts differ")
        require(original["status"] == ("unavailable_mock_identity" if s == "PIC" else "source_qualified"), "Design availability changed")
        tables["paired-design"].append({**labels(scope), "protocol": p, "time": t, "stimulus": s, "paired_lines": n, "paired_runs": runs, "status": original["status"], "unit": "source_line_not_treatment_row_or_gene"})
    for scope, p in itertools.product(SCOPES, PROTOCOLS):
        src = reader.summary["protocol_design"][p]
        tables["cohort-design"].append({**labels(scope), "protocol": p, "source_lines": integer(src["source_lines"]), "mapped_lines": integer(src["mapped_lines"]), "source_runs": integer(src["runs"]), "prediction_lines_both_times": COUNTS[scope][p], "prediction_total_lines": sum(COUNTS[scope].values()), "fixed_prediction_protocol_share": COUNTS[scope][p] / sum(COUNTS[scope].values())})
    if getattr(reader, "optimizer_qc", None) is not None:
        optimizer_projection(reader, tables)
    return tables


# Explicit column whitelist: dict/list source payloads are never serialized.
COLUMNS = {
 "responses": ("identity_scope", "cohort_role", "protocol", "time", "stimulus", "program_id", "quantity", "status", "expected_paired_lines", "expected_paired_runs", "lines", "runs", "input_status_counts", "median", "mean", "minimum", "maximum", "positive", "negative", "zero", "leave_one_line_out_median_range_low", "leave_one_line_out_median_range_high", "leave_one_run_out_median_range_low", "leave_one_run_out_median_range_high", "interval_status", "median_interval_status", "mean_interval_status", "median_interval_low", "median_interval_high", "mean_interval_low", "mean_interval_high", "bootstrap_replicates", "bootstrap_seed", "original_share_Mod_SmartSeq2", "original_share_NEB"),
 "predictions": ("identity_scope", "cohort_role", "variant", "variant_role", "protocol", "time", "status", "is_primary_endpoint", "original_lines", "original_total_lines", "original_Mod_SmartSeq2_lines", "original_NEB_lines", "original_Mod_SmartSeq2_share", "original_NEB_share", "time_weight_6h", "time_weight_24h", "baseline_loss", "extended_loss", "gain", "interval_status", "gain_interval_low", "gain_interval_high", "bootstrap_replicates", "bootstrap_seed"),
 "leave-run-summary": ("identity_scope", "cohort_role", "variant", "variant_role", "omission_protocol", "status", "diagnostic_status_counts", "diagnostic_quantity", "expected_omissions", "source_diagnostic_records", "available_omissions", "unavailable_omissions", "minimum_gain", "maximum_gain", "positive", "negative", "zero", "weights", "interpretation"),
 "model-folds": ("identity_scope", "cohort_role", "variant", "variant_role", "protocol", "time", "fold", "status", "fold_gate_status", "train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs", *[name + "_" + stimulus for name in ("train_class_lines", "train_class_runs") for stimulus in STIMULI[:-1]], "missing_test_stimuli", "line_weighted_gain"),
 "model-fits": ("identity_scope", "cohort_role", "variant", "variant_role", "protocol", "time", "fold", "paired_status", "model", "status", "predictors", "train_pairs", "test_pairs", "train_lines", "test_lines", "train_runs", "test_runs", "iterations", "objective", "gradient_inf"),
 "matching-status": ("identity_scope", "cohort_role", "protocol", "time", "fold", "program_id", "status", "matching_available", "matching_reason", "reference_status", "reference_lines", "reference_runs", "target_bins", "underflow_bin_count", "minimum_pool_bin_size"),
 "scale-status": ("identity_scope", "cohort_role", "variant", "variant_role", "protocol", "time", "fold", "paired_status", "expected_predictors", "scale_status", "numeric_scaler_audit", "method", "evidence", "numeric_scale_values_exported", "constant_feature_flags_exported", "constant_training_feature_count", "constant_training_feature_names", "baseline_constant_training_feature_count", "baseline_constant_training_feature_names", "saved_constant_flags_verified"),
 "independent-optimizer-qc": ("identity_scope", "cohort_role", "variant", "variant_role", "protocol", "time", "fold", "model", "diagnostic_scope", "optimizer_success", "solver_flag_status", "iterations", "gradient_inf", "objective", "coefficient_max_absolute_error", "coefficient_agreement", "frozen_numeric_criteria_pass", "training_rows", "training_feature_count", "constant_training_feature_count", "constant_training_feature_names", "saved_constant_flags_verified", "optimizer_message_or_status_code"),
 "mapping": ("program_id", "source_unique_before_exclusion", "excluded_target_gene_count", "intended_unique_source_genes", "mapped_unique_genes", "missing_source_genes", "ambiguous_source_genes", "mapping_fraction", "mapping_gate_pass", "derived_mapping_fraction", "parent_mapping_fraction", "parent_mapping_gate_pass", "parent_program_id", "parent_mapped_unique_genes", "comparator_union_mapped_genes", "parent_comparator_overlap_genes", "derived_nonoverlap_mapped_genes"),
 "paired-design": ("identity_scope", "cohort_role", "protocol", "time", "stimulus", "paired_lines", "paired_runs", "status", "unit"),
 "cohort-design": ("identity_scope", "cohort_role", "protocol", "source_lines", "mapped_lines", "source_runs", "prediction_lines_both_times", "prediction_total_lines", "fixed_prediction_protocol_share"),
}


def csv_gzip(rows, columns):
    text = io.StringIO(newline=""); writer = csv.DictWriter(text, fieldnames=columns, lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    for row in rows:
        require(set(row) <= set(columns), "Non-whitelisted public column")
        require(all(v is None or type(v) in (str, int, float, bool) for v in row.values()), "Detailed array/object in public table")
        writer.writerow(row)
    target = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=target, compresslevel=9, mtime=0) as stream: stream.write(text.getvalue().encode("utf-8"))
    return target.getvalue()


def response_matrix(rows, protocol):
    import numpy as np
    matrix = np.full((len(PROGRAMS), len(TIMES) * len(STIMULI)), np.nan)
    for i, g in enumerate(PROGRAMS):
        for j, (t, s) in enumerate(itertools.product(TIMES, STIMULI)):
            matches = [r for r in rows if r["identity_scope"] == SCOPES[0] and r["protocol"] == protocol and r["time"] == t and r["stimulus"] == s and r["program_id"] == g and r["quantity"] == (None if s == "PIC" else "adjusted_change")]
            require(len(matches) == 1, "Heatmap dropped/duplicated registered response")
            row = matches[0]
            if row["status"] == "available_descriptive": matrix[i, j] = number(row["median"])
    return matrix


def render_figures(tables):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.hashsalt": "macromap-aggregate-report-v1", "svg.fonttype": "none", "axes.unicode_minus": True})
    figures = {}
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    display_cells = [(p, t) for p, t in itertools.product(PROTOCOLS, TIMES)] + [("pooled", None)]
    for ax, scope in zip(axes, SCOPES):
        for v, variant in enumerate(VARIANTS):
            for j, (p, t) in enumerate(display_cells):
                row = next(r for r in tables["predictions"] if (r["identity_scope"], r["variant"], r["protocol"], r["time"]) == (scope, variant, p, t))
                y = j + (-0.12 if v == 0 else 0.12); color = ("#2166ac", "#b35806")[v]
                label = ("Fixed ETS2 g1-down program", "Fixed comparator-nonoverlap sensitivity")[v] if j == 0 else None
                if row["gain"] is None:
                    ax.text(0.02, y, "unavailable: " + row["status"], transform=ax.get_yaxis_transform(), va="center", fontsize=7, color=color)
                else:
                    ax.plot(row["gain"], y, "o" if v == 0 else "s", color=color, label=label, markersize=5)
                    if row["gain_interval_low"] is not None:
                        ax.hlines(y, row["gain_interval_low"], row["gain_interval_high"], color=color, linewidth=1.8)
                    else:
                        ax.annotate("CI unavailable", (row["gain"], y), xytext=(6, 0), textcoords="offset points", fontsize=7, va="center", color=color)
        ax.axvline(0, color="0.4", linewidth=0.8)
        ax.set_yticks(range(5), [f"{p} {t} h" if t is not None else "Pooled (equal time; fixed line shares)" for p, t in display_cells])
        ax.invert_yaxis(); ax.set_xlabel("Baseline − extended held-out log loss (nats per line); positive = lower loss")
        ax.set_title("Mapped primary: 185 lines (114 + 71)" if scope == SCOPES[0] else "Identity-inclusion sensitivity: 189 lines (114 + 75); not an independent cohort")
        ax.grid(axis="x", alpha=0.2)
        ax.legend(handles=[Line2D([], [], color="#2166ac", marker="o", label="Fixed ETS2 g1-down program"), Line2D([], [], color="#b35806", marker="s", label="Fixed comparator-nonoverlap sensitivity")], loc="best", fontsize=8)
    fig.suptitle("Fixed program incremental prediction\n95% conditional fixed-OOF run-bootstrap intervals; not full-procedure uncertainty")
    figures["prediction"] = fig
    matrices = [response_matrix(tables["responses"], p) for p in PROTOCOLS]
    finite = np.concatenate([a[np.isfinite(a)] for a in matrices]); vmax = float(np.max(np.abs(finite))) if len(finite) else 1.
    vmax = vmax or 1.
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), layout="constrained")
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#bdbdbd")
    for ax, p, matrix in zip(axes, PROTOCOLS, matrices):
        image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto", interpolation="nearest")
        ax.set_yticks(range(len(PROGRAMS)), [LABELS[g] for g in PROGRAMS])
        ax.set_xticks(range(20), [f"{s}\n{t}h" for t, s in itertools.product(TIMES, STIMULI)])
        for i, j in zip(*np.where(~np.isfinite(matrix))): ax.text(j, i, "N/A", ha="center", va="center", fontsize=7, color="black")
        ax.axvline(9.5, color="black", linewidth=1); ax.set_title(p + ": mapped-primary paired lines; grey = unavailable, not zero")
    fig.colorbar(image, ax=axes, shrink=0.8, label="Median adjusted paired response (arbitrary relative score units)")
    fig.suptitle("Fixed nine-program response panel: all 20 registered conditions\nPIC has no qualified mock identity. No response significance screen; complete aggregates exported.")
    figures["responses"] = fig
    result = {}
    for name, figure in figures.items():
        for fmt in ("svg", "png"):
            buf = io.BytesIO()
            metadata = {"Date": None, "Creator": "psc-research aggregate-only MacroMap reporter"} if fmt == "svg" else {"Software": "psc-research aggregate-only MacroMap reporter"}
            figure.savefig(buf, format=fmt, dpi=160, metadata=metadata)
            result[f"{name}_{fmt}"] = buf.getvalue()
        plt.close(figure)
    return result


def write_exclusive(path, data):
    path = absolute(path)
    require(path.parent.is_dir(), "Output parent does not exist")
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o644)
    with os.fdopen(fd, "wb") as stream: stream.write(data)


def publish(reader, data_output, figure_prefix):
    import matplotlib
    import numpy
    import PIL
    require(reader.optimizer_qc is not None, "Publication requires the root-approved aggregate optimizer diagnostic")
    data_output = absolute(data_output); figure_prefix = absolute(figure_prefix)
    paths = {f"{name}_{fmt}": figure_prefix.with_name(figure_prefix.name + f"-{name}.{fmt}") for name, fmt in itertools.product(("prediction", "responses"), ("png", "svg"))}
    protected = (reader.directory, reader.prepared, reader.root / "data/raw")
    for destination in [data_output, *paths.values()]:
        require(not any(destination.is_relative_to(p) or p.is_relative_to(destination) for p in protected), "Output overlaps immutable inputs")
        require(not destination.exists() and not destination.is_symlink(), "Publication destinations must be fresh")
    require(not any(p.is_relative_to(data_output) for p in paths.values()), "Separate figure and table destinations")
    tables = assemble(reader)
    files = {name + ".csv.gz": csv_gzip(rows, COLUMNS[name]) for name, rows in tables.items()}
    files["verification-summary.json"] = (json.dumps(reader.public, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    figures = render_figures(tables)
    reader.recheck()
    manifest = {"schema": "macromap-aggregate-report-v1", "stage": "post_fit_reporting_not_scientific_amendment", "input_pins": reader.pins,
        "reporter_sha256": digest(absolute(__file__).read_bytes()), "authenticated_aggregate_payload_sha256": reader.used,
        "runtime": {"python": platform.python_version(), "numpy": numpy.__version__, "matplotlib": matplotlib.__version__, "pillow": PIL.__version__, "zlib": zlib.ZLIB_RUNTIME_VERSION},
        "table_rows": {name: len(rows) for name, rows in tables.items()}, "table_columns": COLUMNS,
        "data_artifact_sha256": {name: digest(data) for name, data in files.items()}, "figure_artifact_sha256": {name: digest(data) for name, data in figures.items()},
        "privacy": "aggregate whitelist only; no source-line/sample/run identities, matching realizations, scales, coefficients or probabilities",
        "uncertainty": "fixed-OOF conditional predictive CIs; pointwise protocol response CIs; no pooled response CIs",
        "limits": ["Healthy-derived stimulus context is not PSC validation or evidence of ETS2 activity/dependence or clinical benefit.", "Response scores have arbitrary relative units, not calibrated ETS2 activity.", "Sensitivity counterparts are not independent cohorts or best-result selection.", "Numeric scaler means/SDs are not exported; constant-feature names/counts come only from the approved independent aggregate diagnostic. Raw solver flags refer to identical helper reruns, not original SciPy result objects."]}
    require(data_output.parent.is_dir() and figure_prefix.parent.is_dir(), "Create destination parent directories explicitly")
    data_output.mkdir(exist_ok=False)
    for name, data in files.items(): write_exclusive(data_output / name, data)
    for name, data in figures.items(): write_exclusive(paths[name], data)
    reader.recheck()
    write_exclusive(data_output / "reporting-manifest.json", (json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    for name in ("plan", "verification-receipt", "public-verification-summary", "replay-verification", "optimizer-diagnostics", "data-output", "figure-prefix"): parser.add_argument("--" + name, type=Path, required=True)
    for name in ("plan-sha256", "execution-complete-sha256", "verification-receipt-sha256", "public-verification-summary-sha256", "replay-verification-sha256", "optimizer-diagnostics-sha256"): parser.add_argument("--" + name, required=True)
    parser.add_argument("--publish", action="store_true", help="Explicit root-controlled publication after reporter review; otherwise validate aggregates only")
    args = parser.parse_args()
    reader = Reader(args.root, args.plan, args.plan_sha256, args.execution_complete_sha256, args.verification_receipt, args.verification_receipt_sha256, args.public_verification_summary, args.public_verification_summary_sha256, args.replay_verification, args.replay_verification_sha256, args.optimizer_diagnostics, args.optimizer_diagnostics_sha256)
    if args.publish:
        result = publish(reader, args.data_output, args.figure_prefix)
        print(json.dumps({"stage": result["stage"], "table_rows": result["table_rows"]}, sort_keys=True))
    else:
        tables = assemble(reader); reader.recheck()
        print(json.dumps({"stage": "reporting_validation_only_no_publication", "table_rows": {k: len(v) for k, v in tables.items()}}, sort_keys=True))


if __name__ == "__main__":
    main()
