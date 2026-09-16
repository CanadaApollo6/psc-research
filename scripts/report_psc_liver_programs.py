#!/usr/bin/env python3
"""Bounded, hash-specific aggregate reporting. No scientific pipeline imports.

This reporter accepts ONE authenticated, completed 60-row aggregate snapshot.
Unobserved status branches fail closed; changed data need a separate review.
Only projections are exported. No source units, scientific reruns or network.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import stat
import sys
import tempfile

ROOT = Path(__file__).absolute().parent.parent
WORK = "work/psc-liver-program-reporting"
DATA = "data/derived/psc-liver-programs"
FIGURE = "reports/figures/psc-liver-programs"
SCRIPT = "scripts/report_psc_liver_programs.py"
TESTS = "tests/test_report_psc_liver_programs.py"
DESIGN = "reports/psc-liver-program-reporting-design.md"
PUBLISH_ACK = "ROOT_APPROVES_PSC_LIVER_AGGREGATE_PUBLICATION"
COMPLETION = DATA + "/completion-receipt.json"
MANIFEST = DATA + "/report-manifest.json"
CHECKPOINT = "c1ce7f25934c3797eed4660cfed414749c8b46bf"
NATIVE_COMPLETION = "b3b4e5d35dc48c2fa52604dadc0f45dd6c10f2c7a45567e924935140d49ce524"
COMPILED = "17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698"
TIERS = ("exact_symbol_ensembl_bijection", "plus_reciprocal_entrez")
PROGRAMS = (
    "ets2_g1_dn", "ets2_g1_up", "ets2_g2_dn", "chr21_dn", "inflammation",
    "interferon_gamma", "oxidative_stress", "apoptotic_signaling",
    "ets2_g1_dn_without_comparators", "reactome_ecm_organization",
)
CONTRASTS = (("PSC_minus_AIH", "PSC", "AIH"),
             ("PSC_minus_ASC", "PSC", "ASC"),
             ("ASC_minus_AIH", "ASC", "AIH"))
GROUPS = {"PSC": 17, "ASC": 17, "AIH": 30}
DENOMINATORS = {"PSC_minus_AIH": 47, "PSC_minus_ASC": 34, "ASC_minus_AIH": 47}
INTENDED = (927, 668, 875, 141, 815, 139, 436, 588, 927, 321)
MAPPED = ((848, 612, 809, 126, 750, 121, 403, 544, 695, 301),
          (825, 536, 789, 122, 744, 118, 400, 542, 673, 299))
LABELS = ("ETS2 gRNA1-down (focal)", "ETS2 gRNA1-up (not inverted)",
          "ETS2 gRNA2-down", "CHR21-down", "Inflammation", "Interferon gamma",
          "Oxidative stress", "Apoptotic signaling",
          "ETS2 gRNA1-down\nwithout four comparators", "ECM organization (unsigned)")
ROW_COUNTS = {"endpoints.csv": 60, "group-summaries.csv": 60,
              "bootstrap-intervals.csv": 60, "omission-stability.csv": 60,
              "mapping-summary.csv": 20}
PAYLOADS = tuple(DATA + "/" + name for name in ROW_COUNTS) + (
    DATA + "/held-endpoints.json", DATA + "/verification-provenance.json",
    FIGURE + ".png", FIGURE + ".svg",
)
ALL_FILES = PAYLOADS + (MANIFEST, COMPLETION)
FONT_NAMES = ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf")

INPUT_PINS = {'aggregate': {'bytes': 161191,
               'path': 'work/psc-liver-program-plan/runs/primary-v1/aggregate-public-summary.json',
               'sha256': '350ee5d1e447b1cab4621595611de99f7a787af40102187ca0377257c9a598d0'},
 'handoff': {'bytes': 2027,
             'path': 'work/psc-liver-program-root-freeze-v1/root-reporting-handoff.json',
             'sha256': 'd148b76e4e19fc74174fb5c2693cf9caad35ab21483ae224ac251a1364b0397a'},
 'independent_verification': {'bytes': 3712,
                              'path': 'work/psc-liver-program-independent-review/execution-adapter/runs/primary-v1/audit-receipt.json',
                              'sha256': '9a5827803aa44a9124e59bb491d9b9f6d021b1c5fdd79d82cfbed73c1b856004'},
 'mechanical_replay': {'bytes': 5071,
                       'path': 'work/psc-liver-program-root-freeze-v1/root-mechanical-replay-verification.json',
                       'sha256': '00c4e01b10c01f0345933b08bf0a178ba6d94aee8581ac4242f8406c6e6b79e8'},
 'outer_seal': {'bytes': 6859,
                'path': 'config/psc-liver-program-execution-seal.json',
                'sha256': '59eb70cb6a539d84799300332c1ece54bf5384e9655d26884dacf06b62df7c53'},
 'plan': {'bytes': 9247,
          'path': 'config/psc-liver-program-freeze.json',
          'sha256': '5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d'},
 'protocol': {'bytes': 15387,
              'path': 'docs/psc-liver-program-protocol.md',
              'sha256': 'a36faad8bffee1e5f4f0789c8ce77827e1378477751c7e9e29ecd6bb1d3dcf90'},
 'verification_design': {'bytes': 9803,
                         'path': 'reports/psc-liver-program-execution-verification-design.md',
                         'sha256': '4c2550b6e92ef5b0e02db26ff90d9c29b6c646eb10034ef771ad743f25d03eef'}}

# Exact, recursively typed metadata for the accepted aggregate snapshot.
FIXED_TOP = {'all_unit_omission_row_count_is_private_inventory_not_stability_denominator': True,
 'analysis_units': 'source_patient_samples_not_verified_unique_donors',
 'core_family_size': 30,
 'does_not_fulfill_held_IL17_response_transfer': True,
 'endpoint_order': ['ets2_g1_dn',
                    'ets2_g1_up',
                    'ets2_g2_dn',
                    'chr21_dn',
                    'inflammation',
                    'interferon_gamma',
                    'oxidative_stress',
                    'apoptotic_signaling',
                    'ets2_g1_dn_without_comparators',
                    'reactome_ecm_organization'],
 'expected_all_unit_omission_rows': 3840,
 'expected_endpoint_rows': 60,
 'external_candidates': [{'name': 'epithelial_il17a_alone_response',
                          'not_an_empty_endpoint': True,
                          'original_source_members': None,
                          'status': 'HOLD_no_qualified_source'}],
 'group_sizes': {'AIH': 30, 'ASC': 17, 'PSC': 17},
 'independence_status': 'article_claim_only_crosswalk_unavailable',
 'interpretation_limits': ['Whole-biopsy observational transcript-abundance composites, not '
                           'pathway activity or cell-intrinsic effects.',
                           'Pointwise intervals are not simultaneous. BH/BY do not repair invalid '
                           'marginal tests or unverified independence.',
                           'No clinical/stage/batch adjustment, causal interpretation, clinical '
                           'outcome, or clinical recommendation.'],
 'method': {'bootstrap': {'bit_generator': 'PCG64',
                          'group_order': ['PSC', 'ASC', 'AIH'],
                          'joint_across_endpoints_tiers_contrasts': True,
                          'p_values': False,
                          'quantile_method': 'linear',
                          'replicates': 10000,
                          'retain_constant_draws': True,
                          'seed': 2026091602,
                          'within_group_order': 'original_source_column_order'},
            'companion_inversion': False,
            'confidence': 0.95,
            'coverage': {'fraction_denominator': 5,
                         'fraction_numerator': 4,
                         'minimum_mapped_total': 10,
                         'nonoverlap_parent_gate': True},
            'cpm_scale': 1000000,
            'expression_filter': False,
            'external_target_exclusion': 'none_unless_separate_root_amendment',
            'learned_weights': False,
            'multiplicity': {'family': 'approved_panel_x_three_contrasts',
                             'methods': ['BH', 'BY'],
                             'missing_p_for_correction': 1,
                             'strict_q_values': False,
                             'tier': 'exact_symbol_ensembl_bijection',
                             'unavailable_raw_and_q': None},
            'nonfinite_input_policy': 'fail_run',
            'normalization': 'all_released_rows_CPM_then_log2_plus_1',
            'numeric_inference_failure_policy': 'retain_finite_effect_p_ci_df_unavailable',
            'omissions': {'all_source_units': True,
                          'outside_contrast_unchanged_rows': 'private_only',
                          'p_values_are_diagnostics_only': True,
                          'public_stability_units': 'in_contrast_only',
                          'recompute_scores': False,
                          'recompute_welch': True,
                          'sign_reversal': 'strict_opposite_nonzero_signs'},
            'pseudocount': 1,
            'score_orientation': 'positive_equal_member_all_ten_programs',
            'variance_standardization': False,
            'zero_variance_policy': 'welch_if_one_positive'},
 'mode': 'real_execution',
 'normalization': 'all_released_rows_CPM_then_log2_plus_1',
 'not_independent_mechanistic_confirmation_of_MacroMap': True,
 'plan_sha256': '5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d',
 'public_omission_denominators': {'ASC_minus_AIH': 47, 'PSC_minus_AIH': 47, 'PSC_minus_ASC': 34},
 'real_analysis_executed': True,
 'released_normalization_rows': 41096,
 'root_freeze_written': False,
 'schema_version': 1,
 'software': {'numpy': '2.5.3', 'python': '3.12.14', 'scipy': '1.18.1'},
 'source_pin_status': 'all_frozen_source_and_implementation_pins_verified',
 'status': 'complete_for_declared_analysis',
 'strict_discoveries_permitted': False,
 'unavailable_core_p_values': 0}


AUDIT_PUBLIC_SPEC = {'agreement': {'absolute_tolerance': 1e-12,
               'counts_axes_statuses_inventory_indices': 'exact',
               'relative_only_fields': ['degrees_of_freedom',
                                        'p_for_correction',
                                        'p_raw',
                                        'p_vector_for_correction',
                                        'q_bh',
                                        'q_bh_numeric_including_placeholders',
                                        'q_by',
                                        'q_by_numeric_including_placeholders',
                                        'standard_error',
                                        'variance_a',
                                        'variance_b'],
               'relative_tolerance': 1e-10},
 'checks': {'bootstrap_differences_compared': 600000,
            'bootstrap_interval_status_rows_checked': 60,
            'core_correction_vector_length': 30,
            'endpoint_rows_checked': 60,
            'full_library_totals_compared': 64,
            'native_logcpm_matrix_comparison': 'not_available_not_claimed',
            'private_omissions_checked': 3840,
            'public_omission_summaries_checked': 60,
            'raw_source_cells_and_cache_values_compared': 2630144,
            'score_values_compared': 1280,
            'score_vectors_checked': 20,
            'shared_rng_indices_compared': 640000,
            'source_columns_checked': 64,
            'source_rows_checked': 41096,
            'structured_scalar_checks': 748519,
            'unavailable_core_p_values_checked': 0,
            'unavailable_score_vectors_checked': 0},
 'runtime': {'numpy': '2.5.3', 'python': '3.12.14', 'scipy': '1.18.1'}}


class ReportingError(ValueError):
    """Fixed-code error. Never echo input values, IDs, paths or tracebacks."""


def require(condition, code):
    if not condition:
        raise ReportingError(code)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def parse_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate_json_key")
            result[key] = value
        return result

    def bad_constant(_):
        raise ReportingError("nonfinite_json_literal")

    try:
        return json.loads(data, object_pairs_hook=pairs, parse_constant=bad_constant)
    except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise ReportingError("invalid_json") from exc


def exact(actual, expected):
    """Recursive exact keys/types/values; bool is never accepted as an int."""
    require(type(actual) is type(expected), "exact_type")
    if isinstance(expected, dict):
        require(actual.keys() == expected.keys(), "exact_keys")
        for key in expected:
            exact(actual[key], expected[key])
    elif isinstance(expected, list):
        require(len(actual) == len(expected), "exact_list_length")
        for left, right in zip(actual, expected):
            exact(left, right)
    else:
        require(actual == expected, "exact_value")


def keys(value, expected):
    require(type(value) is dict, "object_required")
    require(set(value) == set(expected), "schema_keys")


def number(value, *, low=None, high=None, positive=False):
    require(type(value) in (int, float), "numeric_scalar_required")
    try:
        finite = math.isfinite(value)
    except OverflowError as exc:
        raise ReportingError("numeric_overflow") from exc
    require(finite, "nonfinite_number")
    if low is not None:
        require(value >= low, "numeric_lower_bound")
    if high is not None:
        require(value <= high, "numeric_upper_bound")
    if positive:
        require(value > 0, "positive_number_required")


def integer(value, low, high):
    require(type(value) is int and low <= value <= high, "integer_range")


def safe_path(root, relative):
    """Canonical relative POSIX spelling; reject aliases and existing symlinks."""
    require(type(relative) is str and bool(relative), "invalid_path")
    parsed = PurePosixPath(relative)
    require(not parsed.is_absolute() and str(parsed) == relative and
            not any(part in (".", "..") for part in parsed.parts) and
            "\\" not in relative, "noncanonical_path")
    root = Path(root)
    require(root.is_absolute(), "absolute_root_required")
    # Check root ancestors too; do not silently resolve a symlinked root.
    current = Path(root.anchor)
    for part in root.parts[1:] + parsed.parts:
        current = current / part
        try:
            info = current.lstat()
        except FileNotFoundError:
            continue
        require(not stat.S_ISLNK(info.st_mode), "symlink_refused")
        if current != root / relative:
            require(stat.S_ISDIR(info.st_mode), "ancestor_not_directory")
    return root / relative


def read_regular(root, relative, maximum=4_000_000):
    path = safe_path(root, relative)
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    except OSError as exc:
        raise ReportingError("input_open_refused") from exc
    try:
        info = os.fstat(fd)
        require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1,
                "nonregular_or_hardlinked_input")
        require(info.st_size <= maximum, "input_size_limit")
        with os.fdopen(fd, "rb", closefd=False) as handle:
            data = handle.read(maximum + 1)
        require(len(data) == info.st_size, "input_size_changed")
        after = os.fstat(fd)
        stable = ("st_dev", "st_ino", "st_mode", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
        require(all(getattr(after, key) == getattr(info, key) for key in stable),
                "input_changed_while_reading")
        return data
    finally:
        os.close(fd)


def authenticated(root, pin):
    data = read_regular(root, pin["path"])
    require(len(data) == pin["bytes"] and digest(data) == pin["sha256"], "input_hash_mismatch")
    return data


def coverage_spec(tier_index, program_index):
    nonoverlap = program_index == 8
    return {
        "components": [{"coverage_gate": not nonoverlap,
                        "mapped_unique_members": MAPPED[tier_index][program_index],
                        "name": "all", "original_source_members": INTENDED[program_index]}],
        "gate_basis": "parent_gate_plus_10_retained" if nonoverlap else
                      "original_component_denominators_plus_10_total",
        "mapped_unique_total": MAPPED[tier_index][program_index],
        "parent_gate": True if nonoverlap else None,
        "reasons": [], "score_status": "available",
    }


NUMERIC_FIELDS = ("ci_high", "ci_low", "degrees_of_freedom", "effect", "log_p_raw",
                  "mean_a", "mean_b", "p_raw", "standard_error", "t_statistic",
                  "variance_a", "variance_b")


def validate_row(row, tier_index, program_index, contrast_index):
    tier = TIERS[tier_index]
    program = PROGRAMS[program_index]
    contrast, group_a, group_b = CONTRASTS[contrast_index]
    core = tier_index == 0
    target = {"excluded": False, "literal_symbol": "ETS2",
              "present_in_original_source": False} if program_index == 9 else None
    fixed = {
        "confidence": 0.95, "contrast_id": contrast,
        "contrast_role": "primary" if contrast_index == 0 else "secondary",
        "coverage": coverage_spec(tier_index, program_index),
        "discovery_eligible": core, "effect_status": "available",
        "family_id": "Core_full_panel_three_contrasts" if core else None,
        "family_size": 30 if core else None, "group_a": group_a, "group_b": group_b,
        "inferential_role": "sole_discovery_family" if core else "identity_sensitivity_only",
        "interface": tier, "interval_type": "pointwise_not_multiplicity_adjusted",
        "n_a": GROUPS[group_a], "n_b": GROUPS[group_b],
        "numeric_inference_failure": None, "p_underflow_to_zero": False,
        "p_was_substituted": False, "program_id": program,
        "score_kind": "unsigned_equal_gene", "target_gene_membership": target,
        "variance_a_status": "finite_positive", "variance_b_status": "finite_positive",
        "welch_status": "available", "zero_variance_groups": [],
    }
    keys(row, set(fixed) | set(NUMERIC_FIELDS) |
         {"p_for_correction", "q_bh", "q_by", "bootstrap", "omission"})
    for key, value in fixed.items():
        exact(row[key], value)
    for key in NUMERIC_FIELDS:
        number(row[key])
    for key in ("mean_a", "mean_b"):
        number(row[key], low=0)
    for key in ("variance_a", "variance_b", "standard_error", "degrees_of_freedom"):
        number(row[key], positive=True)
    number(row["p_raw"], positive=True, high=1)
    number(row["log_p_raw"], high=0)
    require(row["ci_low"] <= row["effect"] <= row["ci_high"], "welch_interval_order")
    if core:
        number(row["p_for_correction"], positive=True, high=1)
        require(row["p_for_correction"] == row["p_raw"], "core_correction_input")
        for key in ("q_bh", "q_by"):
            number(row[key], low=row["p_raw"], high=1, positive=True)
        require(row["q_by"] >= row["q_bh"], "correction_order")
    else:
        for key in ("p_for_correction", "q_bh", "q_by"):
            exact(row[key], None)

    boot = row["bootstrap"]
    fixed_boot = {"method": "percentile_linear_95", "p_value_computed": False,
                  "requested_replicates": 10000, "status": "available", "valid_replicates": 10000}
    keys(boot, set(fixed_boot) | {"ci_low", "ci_high"})
    for key, value in fixed_boot.items():
        exact(boot[key], value)
    number(boot["ci_low"])
    number(boot["ci_high"])
    require(boot["ci_low"] < boot["ci_high"], "bootstrap_interval_order_or_status")

    omit = row["omission"]
    denominator = DENOMINATORS[contrast]
    fixed_omit = {
        "defined_effect_rows": denominator, "denominator_kind": "in_contrast_units_only",
        "in_contrast_rows": denominator, "relevant_unit_denominator": denominator,
        "sensitivity_only": True, "status": "available",
    }
    keys(omit, set(fixed_omit) | {"effect_min", "effect_max", "max_absolute_effect_change",
                                "strict_sign_reversal_rows"})
    for key, value in fixed_omit.items():
        exact(omit[key], value)
    for key in ("effect_min", "effect_max", "max_absolute_effect_change"):
        number(omit[key])
    require(omit["effect_min"] <= omit["effect_max"], "omission_range_order")
    number(omit["max_absolute_effect_change"], low=0)
    integer(omit["strict_sign_reversal_rows"], 0, denominator)


def validate_aggregate(value):
    """Validate the entire tree before any table, plot or destination creation.

    This is deliberately snapshot-specific: all 60 Welch/bootstrap results are
    available. An unobserved unavailable state is refused, not imputed. Exact
    nulls (Strict q/family/correction, targets, parent gates) remain null.
    """
    keys(value, set(FIXED_TOP) | {"endpoint_rows"})
    for key, expected in FIXED_TOP.items():
        exact(value[key], expected)
    rows = value["endpoint_rows"]
    require(type(rows) is list and len(rows) == 60, "complete_grid_required")
    index = 0
    for tier_index in range(2):
        for program_index in range(10):
            for contrast_index in range(3):
                validate_row(rows[index], tier_index, program_index, contrast_index)
                index += 1
    # Equality, not recomputation: each reused group/mapping must be consistent.
    seen = {}
    for row in rows:
        for arm in ("a", "b"):
            key = row["interface"], row["program_id"], row["group_" + arm]
            summary = [row["n_" + arm], row["mean_" + arm],
                       row["variance_" + arm], row["variance_" + arm + "_status"]]
            if key in seen:
                exact(summary, seen[key])
            seen[key] = summary
    require(len(seen) == 60, "group_grid_required")
    return value


def validate_receipts(documents):
    handoff, plan, audit, replay = (documents[k] for k in
                                   ("handoff", "plan", "independent_verification", "mechanical_replay"))
    for key in ("aggregate", "plan", "outer_seal", "independent_verification", "mechanical_replay"):
        expected = {k: INPUT_PINS[key][k] for k in handoff[key]}
        exact(handoff[key], expected)
    exact(handoff["primary_completion_sha256"], NATIVE_COMPLETION)
    exact(handoff["actual_source_identity_token_hits_in_aggregate"], 0)
    exact(handoff["all73outerseals_unchanged_after_primary_verification_and_replay"], True)
    exact(plan["state"], "frozen")
    exact(plan["method"], FIXED_TOP["method"])
    exact(plan["contrast_order"], [list(c) for c in CONTRASTS])
    exact([p["program_id"] for p in plan["panel"]], list(PROGRAMS))
    exact(plan["identity"], {"primary": TIERS[0], "sensitivity": TIERS[1]})
    exact(plan["expected_compiled_programs_sha256"], COMPILED)
    exact(audit["plan_sha256"], INPUT_PINS["plan"]["sha256"])
    exact(audit["execution_receipt_sha256"], NATIVE_COMPLETION)
    exact(audit["compiled_programs_sha256"], COMPILED)
    exact(audit["status"], "pass_complete_native_output_verification")
    exact(audit["real_data_audited"], True)
    exact(audit["source_values_or_unit_identifiers_exported"], False)
    for key, expected in AUDIT_PUBLIC_SPEC.items():
        exact(audit[key], expected)
    exact(audit["native_payload_hashes"]["aggregate-public-summary.json"], INPUT_PINS["aggregate"]["sha256"])
    expected_checks = {
        "raw_source_cells_and_cache_values_compared": 2630144,
        "full_library_totals_compared": 64, "score_values_compared": 1280,
        "endpoint_rows_checked": 60, "core_correction_vector_length": 30,
        "shared_rng_indices_compared": 640000, "bootstrap_differences_compared": 600000,
        "private_omissions_checked": 3840, "public_omission_summaries_checked": 60,
        "native_logcpm_matrix_comparison": "not_available_not_claimed",
    }
    for key, expected in expected_checks.items():
        exact(audit["checks"][key], expected)
    exact(replay["status"], "all12_native_artifacts_byte_identical")
    exact(replay["actual_hashes_checked"], 24)
    exact(replay["plan_sha256"], INPUT_PINS["plan"]["sha256"])
    exact(replay["primary_completion_sha256"], NATIVE_COMPLETION)
    exact(replay["replay_completion_sha256"], NATIVE_COMPLETION)
    require(len(replay["comparisons"]) == 12, "native_replay_inventory")
    expected_native = dict(audit["native_payload_hashes"], **{"completion-receipt.json": NATIVE_COMPLETION})
    require(len(expected_native) == 12 and
            [x["artifact"] for x in replay["comparisons"]] == list(expected_native),
            "native_replay_inventory")
    for item in replay["comparisons"]:
        exact(item["byte_identical"], True)
        exact(item["primary"], item["replay"])
        exact(item["primary"]["sha256"], expected_native[item["artifact"]])
    exact(documents["outer_seal"]["native_plan"]["sha256"], INPUT_PINS["plan"]["sha256"])


def own_pins(root):
    return {name: {"path": path, "sha256": digest(data), "bytes": len(data)}
            for name, path in (("reporter", SCRIPT), ("tests", TESTS), ("design", DESIGN))
            for data in [read_regular(root, path)]}


def load_context(root=ROOT, expected_reporter_sha256=None):
    require(type(expected_reporter_sha256) is str and
            re.fullmatch(r"[0-9a-f]{64}", expected_reporter_sha256) is not None,
            "expected_reporter_hash_required")
    local = own_pins(root)
    require(local["reporter"]["sha256"] == expected_reporter_sha256, "reporter_hash_mismatch")
    raw = {name: authenticated(root, pin) for name, pin in INPUT_PINS.items()}
    documents = {name: parse_json(data) for name, data in raw.items()
                 if name not in ("protocol", "verification_design")}
    validate_receipts(documents)
    validate_aggregate(documents["aggregate"])
    return {"root": Path(root), "aggregate": documents["aggregate"],
            "audit": documents["independent_verification"], "pins": local,
            "expected_reporter_sha256": expected_reporter_sha256}


def reauthenticate(context):
    for pin in INPUT_PINS.values():
        authenticated(context["root"], pin)
    exact(own_pins(context["root"]), context["pins"])


def setup_matplotlib(root):
    # No user matplotlibrc/cache dependence. This is the only runtime cache.
    cache = safe_path(root, WORK + "/.matplotlib")
    cache.mkdir(mode=0o700, parents=False, exist_ok=True)
    require(cache.is_dir(), "matplotlib_cache_not_directory")
    os.environ["MPLCONFIGDIR"] = str(cache)
    os.environ["MPLBACKEND"] = "Agg"
    import matplotlib
    matplotlib.use("Agg", force=True)
    matplotlib.rcdefaults()
    matplotlib.rcParams.update({"svg.hashsalt": "psc-liver-program-report-v1",
        "svg.fonttype": "path", "font.family": "DejaVu Sans", "font.size": 10,
        "text.usetex": False, "axes.unicode_minus": True, "figure.dpi": 160,
        "savefig.dpi": 160, "path.simplify": False})
    return matplotlib


def runtime_record(root):
    matplotlib = setup_matplotlib(root)
    from matplotlib import ft2font
    fonts = Path(matplotlib.get_data_path()) / "fonts/ttf"
    return {
        "python": platform.python_version(), "python_full": sys.version,
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(), "byteorder": sys.byteorder,
        "matplotlib": matplotlib.__version__, "numpy": importlib.metadata.version("numpy"),
        "pillow": importlib.metadata.version("Pillow"), "freetype": ft2font.__freetype_version__,
        "backend": "Agg", "svg_hashsalt": "psc-liver-program-report-v1", "dpi": 160,
        "figure_inches": [17, 9], "font_family": "DejaVu Sans",
        "fonts": {name: {"sha256": digest((fonts / name).read_bytes()),
                         "bytes": (fonts / name).stat().st_size} for name in FONT_NAMES},
        "metadata_dates": "omitted", "table_float_serialization": "Python round-trip str; no rounding",
        "statistical_recomputation": False,
    }


def csv_bytes(rows):
    require(bool(rows), "empty_export_table")
    fields = list(rows[0])
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        require(list(row) == fields, "export_column_mismatch")
        converted = {}
        for key, value in row.items():
            require(type(value) in (str, int, float, bool, list, type(None)), "export_value_type")
            if value is None:
                converted[key] = ""
            elif type(value) is bool:
                converted[key] = "true" if value else "false"
            elif type(value) is list:
                converted[key] = json.dumps(value, separators=(",", ":"), allow_nan=False)
            else:
                converted[key] = value
        writer.writerow(converted)
    return output.getvalue().encode("utf-8")


ENDPOINT_FIELDS = (
    "interface", "program_id", "contrast_id", "contrast_role", "group_a", "group_b",
    "n_a", "n_b", "mean_a", "mean_b", "variance_a", "variance_b",
    "variance_a_status", "variance_b_status", "effect", "effect_status", "score_kind",
    "ci_low", "ci_high", "confidence", "interval_type", "standard_error",
    "degrees_of_freedom", "t_statistic", "p_raw", "log_p_raw", "q_bh", "q_by",
    "welch_status", "numeric_inference_failure", "p_underflow_to_zero", "zero_variance_groups",
    "p_was_substituted", "family_id", "family_size", "discovery_eligible", "inferential_role",
)


def mapping_fields(row):
    coverage = row["coverage"]
    component = coverage["components"][0]
    target = row["target_gene_membership"]
    return {
        "coverage_component": component["name"],
        "original_source_members": component["original_source_members"],
        "mapped_unique_members": component["mapped_unique_members"],
        "component_coverage_gate": component["coverage_gate"],
        "mapped_unique_total": coverage["mapped_unique_total"],
        "coverage_gate_basis": coverage["gate_basis"], "parent_gate": coverage["parent_gate"],
        "score_status": coverage["score_status"], "coverage_reasons": coverage["reasons"],
        "target_literal_symbol": None if target is None else target["literal_symbol"],
        "target_present_in_original_source": None if target is None else target["present_in_original_source"],
        "target_excluded": None if target is None else target["excluded"],
    }


def project_tables(aggregate):
    # Public function also validates; no caller can bypass validation by using it directly.
    validate_aggregate(aggregate)
    endpoints, groups, bootstrap, omission, mapping = [], [], [], [], []
    rows = aggregate["endpoint_rows"]
    group_lookup = {}
    for row in rows:
        endpoints.append({**{k: row[k] for k in ENDPOINT_FIELDS}, **mapping_fields(row)})
        identity = {k: row[k] for k in ("interface", "program_id", "contrast_id", "inferential_role")}
        bootstrap.append({**identity, "effect": row["effect"], **row["bootstrap"],
                          "interval_type": "pointwise_not_multiplicity_adjusted"})
        omission.append({**identity, "effect": row["effect"], **row["omission"]})
        if row["contrast_id"] == CONTRASTS[0][0]:
            mapping.append({"interface": row["interface"], "program_id": row["program_id"],
                            "score_kind": row["score_kind"], **mapping_fields(row)})
        for arm in ("a", "b"):
            group_lookup[row["interface"], row["program_id"], row["group_" + arm]] = {
                "interface": row["interface"], "program_id": row["program_id"],
                "group": row["group_" + arm], "n": row["n_" + arm], "mean": row["mean_" + arm],
                "sample_variance_ddof1": row["variance_" + arm],
                "variance_status": row["variance_" + arm + "_status"],
                "score_status": row["coverage"]["score_status"],
                "units": "mean_of_member_log2_1_plus_CPM",
            }
    for tier in TIERS:
        for program in PROGRAMS:
            for group in GROUPS:
                groups.append(group_lookup[tier, program, group])
    tables = dict(zip(ROW_COUNTS, (endpoints, groups, bootstrap, omission, mapping)))
    for name, table in tables.items():
        require(len(table) == ROW_COUNTS[name], "export_row_count")
    return {DATA + "/" + name: csv_bytes(table) for name, table in tables.items()}


def held_summary():
    return {
        "schema_version": 1, "source_protocol": INPUT_PINS["protocol"],
        "held_endpoint": FIXED_TOP["external_candidates"][0],
        "held_endpoint_has_no_score_or_test": True,
        "other_sources": [{"name": "GSE159676", "status": "expression_HOLD"},
                          {"name": "NoPSC_atlas", "status": "PARK"}],
        "not_zero_or_empty_endpoint": True,
    }


def verification_summary(context):
    audit = context["audit"]
    for key, expected in AUDIT_PUBLIC_SPEC.items():
        exact(audit[key], expected)
    exact(audit["status"], "pass_complete_native_output_verification")
    return {
        "schema_version": 1, "status": "authenticated_aggregate_projection_only",
        "inputs": INPUT_PINS, "reporting_artifacts": context["pins"],
        "checkpoint_commit": CHECKPOINT, "native_completion_sha256": NATIVE_COMPLETION,
        "compiled_programs_canonical_sha256": COMPILED,
        "dataset": {"accession": "GSE303271", "released_rows": 41096,
                    "source_patient_samples": 64, "group_sizes": GROUPS,
                    "independence_status": FIXED_TOP["independence_status"],
                    "healthy_group": False},
        "order": {"tiers": list(TIERS), "programs": list(PROGRAMS),
                  "contrasts": [c[0] for c in CONTRASTS], "group_summary_groups": list(GROUPS)},
        "inference": {"Core_family": 30, "Core_corrections": ["BH", "BY"],
                      "Strict_rows": 30, "Strict_role": "identity_sensitivity_only_no_discovery_q",
                      "intervals": "pointwise_95_percent_not_multiplicity_adjusted",
                      "observed_P_replaced": False, "unavailable_fields": "null_JSON_empty_CSV_never_zero_filled",
                      "omission_denominators": DENOMINATORS,
                      "omission_ranges_are_not_confidence_intervals": True},
        "independent_verification": {"status": audit["status"], "checks": audit["checks"],
                                     "agreement": audit["agreement"], "runtime": audit["runtime"],
                                     "performed_by_reporter": False},
        "native_mechanical_replay": {"files_byte_identical": 12, "actual_hashes_checked": 24,
                                     "completion_included": True, "performed_by_reporter": False,
                                     "biological_replication": False},
        "outer_seal_preservation": {"root_handoff_attests_unchanged": 73,
                                    "reporter_rehashed_all_73": False,
                                    "reporter_only_rehashed_listed_inputs": True},
        "privacy": {"bounded_to_source_aggregate_sha256": INPUT_PINS["aggregate"]["sha256"],
                    "root_handoff_source_identity_token_hits": 0,
                    "reporter_opened_source_unit_identifiers": False,
                    "recursive_snapshot_schema_validated_before_projection": True,
                    "universal_privacy_guarantee": False},
        "scientific_rerun": False, "new_statistics": False,
        "root_freeze_written_false_means": "native_writer_did_not_create_external_root_freeze",
        "limitations": [
            "Whole-biopsy observational abundance; composition and clinical/technical confounding remain.",
            "Paper-asserted source-patient independence; no independent person/visit crosswalk.",
            "No healthy group; no equivalence, generic-inflammation or PSC-cause proof.",
            "Positive equal-member unsigned means; gRNA1-up not inverted.",
            "Effects are mean log2(1+CPM) composite differences, not raw fold changes or TF activation.",
            "ECM: 321 official symbols, Core301/Strict299; unsigned organization, not fibrosis response.",
            "No cell-intrinsic regulation, induction, clinical recommendation or independent MacroMap confirmation.",
            "Held epithelial IL17 response-transfer excluded; GSE159676 HOLD and NoPSC PARK.",
        ],
    }


def figure_bytes(aggregate, root):
    validate_aggregate(aggregate)
    matplotlib = setup_matplotlib(root)
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.lines import Line2D
    fonts = Path(matplotlib.get_data_path()) / "fonts/ttf"
    font = FontProperties(fname=str(fonts / FONT_NAMES[0]))
    bold = FontProperties(fname=str(fonts / FONT_NAMES[1]))
    fig, axes = plt.subplots(1, 3, figsize=(17, 9), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.26, right=0.98, top=0.76, bottom=0.24, wspace=0.12)
    colors = ("#165B88", "#BA5913")
    lookup = {(r["interface"], r["program_id"], r["contrast_id"]): r
              for r in aggregate["endpoint_rows"]}
    for contrast_index, (contrast, group_a, group_b) in enumerate(CONTRASTS):
        ax = axes[contrast_index]
        ax.axvline(0, color="#777777", linewidth=1, linestyle="--", zorder=0)
        for program_index, program in enumerate(PROGRAMS):
            y = 9 - program_index
            if program_index % 2 == 0:
                ax.axhspan(y - .48, y + .48, color="#F1F3F4", zorder=0)
            for tier_index, tier in enumerate(TIERS):
                row = lookup[tier, program, contrast]
                offset = .13 if tier_index == 0 else -.13
                position = y + offset
                ax.plot([row["ci_low"], row["ci_high"]], [position, position],
                        color=colors[tier_index], linewidth=1.7, zorder=2)
                ax.plot([row["ci_low"], row["ci_high"]], [position, position],
                        marker="|", linestyle="none", color=colors[tier_index], markersize=5)
                ax.plot(row["effect"], position, marker="o" if tier_index == 0 else "D",
                        markersize=5, markerfacecolor=colors[tier_index] if tier_index == 0 else "white",
                        markeredgecolor=colors[tier_index], linestyle="none", zorder=3)
        ax.set_title(f"{group_a} (n={GROUPS[group_a]}) − {group_b} (n={GROUPS[group_b]})\n"
                     + ("Primary contrast" if contrast_index == 0 else "Fixed secondary contrast"),
                     fontproperties=bold, fontsize=12, pad=16)
        # Fixed full-range axis for this authenticated snapshot; clipping is refused.
        require(all(-.38 < lookup[t, p, contrast]["ci_low"] <=
                    lookup[t, p, contrast]["ci_high"] < .38 for t in TIERS for p in PROGRAMS),
                "figure_would_clip_interval")
        ax.set_xlim(-.38, .38)
        ax.set_ylim(-.65, 9.65)
        ax.set_xticks([-.3, -.15, 0, .15, .3])
        ax.set_yticks(list(range(9, -1, -1)), LABELS, fontproperties=font, fontsize=11)
        ax.tick_params(axis="y", length=0, pad=12)
        ax.grid(axis="x", alpha=.18)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        for text in ax.get_xticklabels():
            text.set_fontproperties(font)
            text.set_fontsize(10)
    fig.text(.5, .94, "Pediatric liver: all ten fixed unsigned transcript-abundance programs",
             ha="center", fontproperties=bold, fontsize=17)
    fig.text(.5, .895, "GSE303271 · whole biopsies · 64 source patient-samples · 17 PSC / 17 ASC / 30 AIH",
             ha="center", fontproperties=font, fontsize=12)
    handles = [Line2D([0], [0], color=colors[i], marker="o" if i == 0 else "D",
                      markerfacecolor=colors[i] if i == 0 else "white", linewidth=1.7)
               for i in range(2)]
    fig.legend(handles, ["Core: primary; exact symbol–Ensembl bijection",
                         "Strict: reciprocal Entrez; identity sensitivity only"],
               loc="lower center", bbox_to_anchor=(.5, .845), ncol=2, frameon=False, prop=font)
    fig.text(.62, .18, "Mean log2(1 + CPM) composite difference (first group minus second)",
             ha="center", fontproperties=bold, fontsize=12)
    footnotes = [
        "Bars: pointwise 95% Welch intervals (unadjusted, not simultaneous). No individual values shown.",
        "All 30 Core hypotheses: BH and BY reported in the complete endpoint table; Strict is not a discovery family.",
        "Positive scores do not imply TF activation or induction. ECM is unsigned organization, not a fibrosis response.",
        "Whole-biopsy composition / clinical confounding and article-assumed patient independence limit interpretation.",
    ]
    for index, text in enumerate(footnotes):
        fig.text(.02, .118 - .026 * index, text, fontproperties=font, fontsize=10)
    outputs = {}
    try:
        for extension in ("png", "svg"):
            stream = io.BytesIO()
            metadata = {"Software": "PSC liver aggregate reporter v1"} if extension == "png" else {
                "Date": None, "Creator": "PSC liver aggregate reporter v1",
                "Title": "Complete fixed PSC liver program contrasts",
                "Description": "Aggregate Core and Strict pointwise Welch intervals; no individual values.",
            }
            fig.savefig(stream, format=extension, dpi=160, metadata=metadata)
            outputs[FIGURE + "." + extension] = stream.getvalue()
    finally:
        plt.close(fig)
    return outputs


def inventory(payloads):
    return [{"path": name, "bytes": len(payloads[name]), "sha256": digest(payloads[name])}
            for name in sorted(payloads)]


def build_bundle(context):
    validate_aggregate(context["aggregate"])
    payloads = project_tables(context["aggregate"])
    payloads[DATA + "/held-endpoints.json"] = json_bytes(held_summary())
    payloads[DATA + "/verification-provenance.json"] = json_bytes(verification_summary(context))
    runtime = runtime_record(context["root"])
    payloads.update(figure_bytes(context["aggregate"], context["root"]))
    require(set(payloads) == set(PAYLOADS), "payload_inventory")
    manifest = {
        "schema_version": 1, "status": "complete_bounded_aggregate_report",
        "paths_are_logical_public_destinations": True,
        "inputs": INPUT_PINS, "reporting_artifacts": context["pins"], "runtime": runtime,
        "row_counts": ROW_COUNTS, "payloads": inventory(payloads),
        "null_policy": "JSON_null_CSV_empty; no zero-fill or observed_P_substitution",
        "native_verification_is_not_reporting_reexecution": True,
        "privacy_scope": "one_hash_authenticated_recursively_validated_aggregate_not_universal_guarantee",
    }
    payloads[MANIFEST] = json_bytes(manifest)
    completion = {
        "schema_version": 1, "status": "complete_validated_aggregate_report",
        "manifest_sha256": digest(payloads[MANIFEST]),
        "inventory": inventory(payloads), "completion_not_self_listed": True,
        "source_aggregate_sha256": INPUT_PINS["aggregate"]["sha256"],
        "reporter_sha256": context["expected_reporter_sha256"],
        "payload_file_count": 10, "whole_bundle_file_count": 11,
    }
    payloads[COMPLETION] = json_bytes(completion)
    reauthenticate(context)
    return payloads


def fresh_preview_path(root, relative):
    path = safe_path(root, relative)
    parsed = PurePosixPath(relative)
    require(parsed.parent == PurePosixPath(WORK) and
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", parsed.name) is not None,
            "preview_destination_not_dedicated_child")
    require(not path.exists(), "destination_reused")
    require(path.parent.is_dir(), "preview_parent_missing")
    return path


def actual_files(directory):
    require(directory.is_dir() and not directory.is_symlink(), "bundle_directory_required")
    found = set()
    expected_directories = {str(parent) for name in ALL_FILES
                            for parent in PurePosixPath(name).parents if str(parent) != "."}
    for path in directory.rglob("*"):
        relative = path.relative_to(directory).as_posix()
        info = path.lstat()
        require(not stat.S_ISLNK(info.st_mode), "bundle_symlink")
        if stat.S_ISDIR(info.st_mode):
            require(relative in expected_directories, "extra_bundle_directory")
        else:
            require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1, "bundle_nonregular")
            found.add(relative)
    require(found == set(ALL_FILES), "actual_bundle_inventory")


def verify_bundle(root, relative, expected_completion, context):
    directory = safe_path(root, relative)
    require(PurePosixPath(relative).parent == PurePosixPath(WORK), "preview_source_outside_work")
    actual_files(directory)
    payloads = {name: read_regular(directory, name) for name in ALL_FILES}
    require(digest(payloads[COMPLETION]) == expected_completion, "preview_completion_hash_mismatch")
    completion = parse_json(payloads[COMPLETION])
    expected = {
        "schema_version": 1, "status": "complete_validated_aggregate_report",
        "manifest_sha256": digest(payloads[MANIFEST]),
        "inventory": inventory({k: v for k, v in payloads.items() if k != COMPLETION}),
        "completion_not_self_listed": True,
        "source_aggregate_sha256": INPUT_PINS["aggregate"]["sha256"],
        "reporter_sha256": context["expected_reporter_sha256"],
        "payload_file_count": 10, "whole_bundle_file_count": 11,
    }
    exact(completion, expected)
    manifest = parse_json(payloads[MANIFEST])
    expected_manifest = {
        "schema_version": 1, "status": "complete_bounded_aggregate_report",
        "paths_are_logical_public_destinations": True,
        "inputs": INPUT_PINS, "reporting_artifacts": context["pins"],
        "runtime": runtime_record(root), "row_counts": ROW_COUNTS,
        "payloads": inventory({k: payloads[k] for k in PAYLOADS}),
        "null_policy": "JSON_null_CSV_empty; no zero-fill or observed_P_substitution",
        "native_verification_is_not_reporting_reexecution": True,
        "privacy_scope": "one_hash_authenticated_recursively_validated_aggregate_not_universal_guarantee",
    }
    exact(manifest, expected_manifest)
    # Authenticate semantic projections too, without recomputing any statistic.
    for name, data in project_tables(context["aggregate"]).items():
        require(payloads[name] == data, "table_projection_changed")
    require(payloads[DATA + "/held-endpoints.json"] == json_bytes(held_summary()), "held_projection_changed")
    require(payloads[DATA + "/verification-provenance.json"] == json_bytes(verification_summary(context)),
            "verification_projection_changed")
    reauthenticate(context)
    return payloads


def exclusive_write(path, data):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    info = None
    try:
        info = os.fstat(fd)
        with os.fdopen(fd, "wb", closefd=False) as handle:
            handle.write(data)
            handle.flush()
            os.fsync(fd)
    except BaseException:
        # Remove only this call's newly created file, never an existing O_EXCL target.
        try:
            if info is None:
                info = os.stat(fd)
            current = path.lstat()
            if (current.st_dev, current.st_ino) == (info.st_dev, info.st_ino):
                path.unlink()
        except FileNotFoundError:
            pass
        raise
    finally:
        os.close(fd)


def install_completion(temporary, destination, data):
    exclusive_write(temporary, data)
    info = temporary.lstat()
    linked = False
    try:
        # Hard-link installation is atomic and refuses replacement, unlike rename.
        os.link(temporary, destination, follow_symlinks=False)
        linked = True
        temporary.unlink()
    except BaseException:
        # A failed post-link cleanup must not leave a success marker while the
        # caller rolls back payloads. Never remove a pre-existing destination.
        if linked:
            try:
                current = destination.lstat()
                if (current.st_dev, current.st_ino) == (info.st_dev, info.st_ino):
                    destination.unlink()
            except FileNotFoundError:
                pass
        try:
            current = temporary.lstat()
            if (current.st_dev, current.st_ino) == (info.st_dev, info.st_ino):
                temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def write_preview(context, relative, reference=None):
    root = context["root"]
    destination = fresh_preview_path(root, relative)
    payloads = build_bundle(context)
    if reference is not None:
        require(payloads == reference, "report_replay_bytes_differ")
    # Rendering/validation failures create no result destination. All staging is ignored.
    stage = Path(tempfile.mkdtemp(prefix=".stage-", dir=destination.parent))
    try:
        for name in ALL_FILES:
            if name == COMPLETION:
                continue
            path = safe_path(stage, name)
            path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
            exclusive_write(path, payloads[name])
        for name in ALL_FILES[:-1]:
            require(read_regular(stage, name) == payloads[name], "actual_write_mismatch")
        reauthenticate(context)
        # Reserve destination exclusively, preventing rename-overwrite of an empty directory.
        destination.mkdir(mode=0o700)
        try:
            for child in list(stage.iterdir()):
                os.rename(child, destination / child.name)
            # Completion is atomically installed LAST, after actual final payload validation.
            for name in ALL_FILES[:-1]:
                require(read_regular(destination, name) == payloads[name], "final_write_mismatch")
            temporary = destination / ".completion.tmp"
            install_completion(temporary, destination / COMPLETION, payloads[COMPLETION])
            actual_files(destination)
        except BaseException:
            shutil.rmtree(destination)
            raise
    finally:
        shutil.rmtree(stage)
    return {"status": "preview_complete" if reference is None else "all11_reporting_files_byte_identical",
            "output_directory": relative, "files": 11,
            "completion_sha256": digest(payloads[COMPLETION]),
            "manifest_sha256": digest(payloads[MANIFEST]),
            "actual_hashes_checked": 11 if reference is None else 22}


def publish(context, preview, expected_completion, acknowledgment):
    require(acknowledgment == PUBLISH_ACK, "root_publication_ack_required")
    root = context["root"]
    payloads = verify_bundle(root, preview, expected_completion, context)
    target_dir = safe_path(root, DATA)
    figures = [safe_path(root, FIGURE + extension) for extension in (".png", ".svg")]
    require(not target_dir.exists() and all(not p.exists() for p in figures), "public_destination_reused")
    require(target_dir.parent.is_dir() and all(p.parent.is_dir() for p in figures), "public_parent_missing")
    created = []
    reserved = False
    try:
        target_dir.mkdir(mode=0o755)
        reserved = True
        for name in ALL_FILES:
            if name == COMPLETION:
                continue
            path = safe_path(root, name)
            exclusive_write(path, payloads[name])
            created.append(path)
        for name in ALL_FILES[:-1]:
            require(read_regular(root, name) == payloads[name], "public_actual_write_mismatch")
        reauthenticate(context)
        # Complete inventory of the new data directory; no partials/extras allowed.
        expected_data = {PurePosixPath(n).name for n in ALL_FILES if n.startswith(DATA + "/") and n != COMPLETION}
        require({p.name for p in target_dir.iterdir()} == expected_data, "public_data_inventory")
        temporary = target_dir / ".completion.tmp"
        install_completion(temporary, root / COMPLETION, payloads[COMPLETION])
        created.append(root / COMPLETION)
        for name in ALL_FILES:
            require(read_regular(root, name) == payloads[name], "public_final_write_mismatch")
    except BaseException:
        for path in reversed(created):
            if path.exists():
                path.unlink()
        if reserved:
            # Never remove a pre-existing directory or an unexpected foreign file.
            try:
                target_dir.rmdir()
            except OSError:
                pass
        raise
    return {"status": "published_exact_authenticated_preview", "files": 11,
            "completion_sha256": digest(payloads[COMPLETION]),
            "manifest_sha256": digest(payloads[MANIFEST])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("validate", "preview", "replay", "publish"):
        command_parser = subparsers.add_parser(command)
        command_parser.add_argument("--expected-reporter-sha256", required=True)
        if command in ("preview", "replay"):
            command_parser.add_argument("--output-directory", required=True)
        if command in ("replay", "publish"):
            command_parser.add_argument("--preview-directory", required=True)
            command_parser.add_argument("--expected-preview-completion-sha256", required=True)
        if command == "publish":
            command_parser.add_argument("--root-publication-ack", required=True)
    args = parser.parse_args(argv)
    try:
        context = load_context(ROOT, args.expected_reporter_sha256)
        if args.command == "validate":
            result = {"status": "pass_authenticated_recursive_aggregate_schema", "endpoint_rows": 60,
                      "Core_family": 30, "Strict_rows_no_q": 30, "inputs": INPUT_PINS,
                      "reporting_artifacts": context["pins"], "no_scientific_rerun": True}
        elif args.command == "preview":
            result = write_preview(context, args.output_directory)
        elif args.command == "replay":
            reference = verify_bundle(ROOT, args.preview_directory,
                                      args.expected_preview_completion_sha256, context)
            result = write_preview(context, args.output_directory, reference)
        else:
            result = publish(context, args.preview_directory,
                             args.expected_preview_completion_sha256, args.root_publication_ack)
        print(json_bytes(result).decode(), end="")
        return 0
    except (ReportingError, OSError, KeyError, TypeError, ValueError, RuntimeError) as exc:
        # Do not print arbitrary exception strings (they can contain private values).
        code = str(exc) if isinstance(exc, ReportingError) else "reporting_operation_refused"
        print(json_bytes({"status": "refused", "code": code}).decode(), end="", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
