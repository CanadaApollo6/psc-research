#!/usr/bin/env python3
"""Source-blinded GSE303271 method draft; synthetic math and guarded future execution.

No default action reads expression values. An approval JSON is a workflow control,
not a cryptographic proof of an approver's identity. Root must own the freeze.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import sys

import numpy as np
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
WORK_REL = "work/psc-liver-program-plan"
RUNS_REL = WORK_REL + "/runs"
SCRIPT_REL = "scripts/analyze_psc_liver_programs.py"
TEST_REL = "tests/test_psc_liver_programs.py"
IDENTIFIER_JSON = "data/derived/psc-liver-program-inputs.json"
FEATURE_AUDIT = "work/psc-liver-program-inputs/feature-identifier-audit.json"
LEGACY_INDICES = "work/psc-liver-program-inputs/qualified-program-feature-indices.csv"
SAMPLE_JOIN = "work/psc-liver-qualification/GSE303271-sample-source-join.json"
CORE = "exact_symbol_ensembl_bijection"
STRICT = "plus_reciprocal_entrez"
INTERFACES = (CORE, STRICT)
GROUPS = ("PSC", "ASC", "AIH")
CONTRASTS = (("PSC_minus_AIH", "PSC", "AIH"),
             ("PSC_minus_ASC", "PSC", "ASC"),
             ("ASC_minus_AIH", "ASC", "AIH"))
GROUP_COUNTS = {"PSC": 17, "ASC": 17, "AIH": 30}
SHAPE = (41096, 64)
LEGACY_DENOMINATORS = {"ets2_g1_dn": 927, "ets2_g1_up": 668,
    "ets2_g2_dn": 875, "chr21_dn": 141, "inflammation": 815,
    "interferon_gamma": 139, "oxidative_stress": 436,
    "apoptotic_signaling": 588, "ets2_g1_dn_without_comparators": 927}
NONOVERLAP = "ets2_g1_dn_without_comparators"
COMPARATORS = ("inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling")
EXECUTION_ACK = "ROOT_AUTHORIZED_REAL_EXPRESSION_EXECUTION"
ROOT_PANEL_DECISION = "agentmsg_8bf94a77-6693-40fa-b017-395f3ce33d56"
ROOT_METHOD_DECISION = "agentmsg_d84e9242-be6a-4f66-9ebb-cb9d839f0ca6"
BOOTSTRAP_SEED = 2026091602
ECM_PROGRAM = "reactome_ecm_organization"
TEN_PROGRAMS = (*LEGACY_DENOMINATORS, ECM_PROGRAM)
DECISION_NAMES = ("panel", "normalization", "score_orientation", "zero_variance",
                  "bootstrap", "omissions", "external_target_exclusion",
                  "source_patient_sample_independence")
ZERO_POLICIES = ("unavailable_if_either_zero", "welch_if_one_positive")
CLOSED_VERSIONED_PINS = {'scripts/prepare_psc_liver_program_inputs.py': 'bb7165a4e3138919b92cc06009c3212bf79097c42690e12ff273bbbc97f2b866', 'tests/test_psc_liver_program_inputs.py': '13cc0064212e7096be0ccb2cfe0f307f82375602c46da57361125bd2e1764f25', 'reports/psc-liver-program-inputs.md': '7e77ce6439da2df34fc673e7d96519623c0b6fb5fda9e515ef02f5d253c07035', 'data/derived/psc-liver-program-inputs.json': '14afd9c313138b39f039876a4b2b67e7e1cefe599340b2ea2e190f4e6914e30a', 'data/derived/psc-liver-input-qualification.json': '3c10ec1199725c19a3de37a0e532a66dac15b9e106fc11c9973a209760bfa4de', 'reports/psc-liver-input-qualification.md': '179b55b4b62f73467476d7fcb39d47da80dee46c5b3c2da056803f4aa677c3b3', 'scripts/qualify_psc_liver_inputs.py': '2e870af8031489d3bb49adf6533fb04f16cb671cf43db7be816154b4c3f88137', 'tests/test_psc_liver_inputs.py': '24059bd1eb0487b48ee2c9d2af64f77b5cdda7965ad7aa0d5fe17b8e0ef8851e', 'data/derived/psc-organoid-input-qualification.json': '1658188a59a83adcfb059161dbd0daed78fc8938d8335dce8a7fcb0bc9198e99', 'reports/psc-organoid-input-qualification.md': '169814232eda0a03673d993f231aba6851b3f1d0a48f3f2a76f72f7a431733c5', 'scripts/qualify_psc_organoid_inputs.py': '1f1bda0c8c3596e78f72c48858f52ec75acb593c9a834dcdf651ca571451489e', 'tests/test_psc_organoid_inputs.py': '49161257b8049aace73e20bebd7a4a6fb6fc42ebe463bc00cb3b2557ab146b3c'}
METADATA_PINS = {'work/psc-liver-program-inputs/feature-identifier-audit.json': {'bytes': 21875165, 'path': 'work/psc-liver-program-inputs/feature-identifier-audit.json', 'sha256': 'e469b26496ca7bee552bce9a23d20c5a9a9e2301a22a3bcb6732102efe13d556'}, 'work/psc-liver-program-inputs/program-coverage.json': {'bytes': 16315, 'path': 'work/psc-liver-program-inputs/program-coverage.json', 'sha256': '091b80c2fb6d9ad3f677232147a642fc9f9aa9b3ed95ea9a33e1e5251382f4e8'}, 'work/psc-liver-program-inputs/qualified-program-feature-indices.csv': {'bytes': 682694, 'path': 'work/psc-liver-program-inputs/qualified-program-feature-indices.csv', 'sha256': '19d80e6bcca82ebd4d782fde8986fc3aa4699aa3bddd548e0b5c8b9e0ad9de60'}, 'work/psc-liver-qualification/GSE303271-sample-source-join.json': {'bytes': 61623, 'path': 'work/psc-liver-qualification/GSE303271-sample-source-join.json', 'sha256': '30427a566b9e857a5f98dfb389eace558bb1a6c088d6da4f8bfaf423f85463f0'}}
MATRIX_PIN = {'bytes': 21041280, 'path': 'work/psc-liver-qualification/GSE303271-submitted-float64.npy', 'sha256': '4a19908237dfaeb2f6c713fd6c7f60743d9b31b556a643c9c0dfcab7159003d7'}
RAW_MATRIX_PIN = {"path": "data/raw/psc-liver-qualification/GSE303271_raw_counts.txt.gz",
                  "sha256": "966be15cf8cd893d0da3ec62d73ee0a56a34932a9d1c990c3d0d58b903a24046"}
EXTERNAL_SOURCE_REGISTRY = {'reactome_ecm_organization': {'program_id': 'reactome_ecm_organization', 'source_manifest_pin': {'path': 'config/psc-liver-external-program-sources.json', 'sha256': '2428cf428187581d992dc2b507305b259465db6259aec19b8b658614c0773305'}, 'membership_pin': {'path': 'data/derived/psc-liver-external-program-membership.csv', 'sha256': '8cf6573a97d01cc8e83558634bf62b6d48d9be015e67f50d15868071c18b7202', 'bytes': 125834}, 'schema_handoff_pin': {'path': 'work/psc-liver-external-programs/method-schema-handoff.json', 'sha256': '69d7653d5829837eed563cef2c5d22eca1394c03b62a2942909f18fef50e3652'}, 'score_kind': 'unsigned_equal_gene', 'source_identifier_namespace': 'MSigDB_geneSymbols', 'original_component_members': {'all': 321}, 'direction_literals': {'all': 'UNSIGNED'}, 'target_gene_exclusions': [], 'source_version': 'MSigDB2026.1.Hs / Reactome95, M610 / R-HSA-1474244'}}


EXTERNAL_RAW_SOURCE_PINS = {'ecm_html': {'path': 'data/raw/psc-liver-external-programs/ecm-scout/01-msigdb-ecm.html', 'sha256': '0bce8c070ad8e24961cd22fabd57da8cf26161f7ece5d3da49715812ac1895a6', 'bytes': 197895}, 'ecm_json': {'path': 'data/raw/psc-liver-external-programs/ecm-scout/02-msigdb-ecm.json', 'sha256': '2046e084f1c36c74cf107c8448d19bb4f47bb2bf37b45b2388fd3d1b837b6540', 'bytes': 3006}}


class GuardError(ValueError):
    """Approval, identity, pin or destination does not meet the contract."""


class DataIntegrityError(ValueError):
    """Invalid numeric/input data; never silently drop or substitute it."""


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode("utf-8")).hexdigest()


def file_sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    def reject(value):
        raise GuardError("Nonstandard nonfinite JSON literal: " + value)
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise GuardError("Duplicate JSON key: " + key)
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      parse_constant=reject, object_pairs_hook=unique)


def dump_json(path, data):
    path = Path(path)
    encoded = json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n"
    temp = path.with_name(path.name + ".partial")
    with temp.open("x", encoding="utf-8") as stream:
        stream.write(encoded)
    temp.replace(path)


def relative_input(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise GuardError("Input path must be a relative repository path")
    target = Path(root) / relative
    for entry in (target,*target.parents):
        if entry == Path(root):
            break
        if entry.is_symlink():
            raise GuardError("Pinned input path or ancestry cannot be a symlink")
    if not target.resolve().is_relative_to(Path(root).resolve()):
        raise GuardError("Input symlink escapes repository")
    return target


def metadata_input(root, relative):
    """Guard every ancillary input before any body read or implementation hash."""
    target = relative_input(root,relative)
    if target.exists():
        for relative_payload in (MATRIX_PIN["path"],RAW_MATRIX_PIN["path"]):
            payload = Path(root)/relative_payload
            if payload.exists() and target.samefile(payload):
                raise GuardError("Metadata input aliases a reserved expression payload")
    return target


def verify_pin(root, pin):
    if not isinstance(pin, dict) or not re.fullmatch(r"[a-f0-9]{64}", str(pin.get("sha256", ""))):
        raise GuardError("Malformed SHA-256 pin")
    payload_paths = (MATRIX_PIN["path"],RAW_MATRIX_PIN["path"])
    target = relative_input(root,pin["path"]) if pin["path"] in payload_paths else metadata_input(root,pin["path"])
    if not target.is_file():
        raise GuardError("Missing pinned input: " + pin["path"])
    if "bytes" in pin and target.stat().st_size != pin["bytes"]:
        raise GuardError("Pinned size changed: " + pin["path"])
    if file_sha(target) != pin["sha256"]:
        raise GuardError("Pinned hash changed: " + pin["path"])
    return target


def software_versions():
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__}


def draft_plan(root=ROOT):
    """Create a proposal only. Never create a freeze or an execution approval."""
    root = Path(root)
    return {
        "schema_version": 1, "analysis_id": "GSE303271-fixed-programs",
        "state": "draft_not_frozen", "panel_approved": True,
        "root_panel_decision_reference": ROOT_PANEL_DECISION,
        "root_method_decision_reference": ROOT_METHOD_DECISION,
        "research_question": "narrower_unsigned_ETS2_comparator_ECM_tissue_context",
        "does_not_fulfill_held_IL17_response_transfer": True,
        "root_authorization": {"approver_role": None, "reference": None,
            "approved_at": None, "real_execution_authorized": False},
        "decisions": {key: {"status": "approved", "reference": ROOT_METHOD_DECISION}
                      for key in DECISION_NAMES},
        "independent_review": {"status": "not_accepted", "reference": None,
                               "artifact_pin": None},
        "dataset": {"accession": "GSE303271", "shape": list(SHAPE),
                    "groups": GROUP_COUNTS, "matrix_pin": MATRIX_PIN,
                    "raw_matrix_pin": RAW_MATRIX_PIN},
        "identity": {"primary": CORE, "sensitivity": STRICT},
        "unit": {"kind": "source_patient_sample",
            "independence_status": "article_claim_only_crosswalk_unavailable",
            "assumption_approved": True,
            "known_repeat_status": "no_known_repeats_crosswalk_unavailable"},
        "panel": [{"program_id": key, "kind": "accepted_legacy", "status": "qualified",
                   "source_approval_reference": "accepted_identifier_preparation",
                   "score_kind": "unsigned_equal_gene",
                   "original_source_members": value}
                  for key, value in LEGACY_DENOMINATORS.items()] + [
                     {**EXTERNAL_SOURCE_REGISTRY[ECM_PROGRAM],"kind":"external_symbol_csv",
                      "status":"qualified","source_approval_reference":ROOT_PANEL_DECISION}],
        "external_candidates": [
            {"name": "epithelial_il17a_alone_response", "status": "HOLD_no_qualified_source",
             "original_source_members": None, "not_an_empty_endpoint": True}],
        "contrast_order": [list(x) for x in CONTRASTS],
        "method": {
            "normalization": "all_released_rows_CPM_then_log2_plus_1",
            "pseudocount": 1, "cpm_scale": 1000000,
            "expression_filter": False, "variance_standardization": False,
            "learned_weights": False, "companion_inversion": False,
            "coverage": {"fraction_numerator": 4, "fraction_denominator": 5,
                         "minimum_mapped_total": 10,
                         "nonoverlap_parent_gate": True},
            "score_orientation": "positive_equal_member_all_ten_programs",
            "external_target_exclusion": "none_unless_separate_root_amendment",
            "confidence": 0.95, "zero_variance_policy": "welch_if_one_positive",
            "nonfinite_input_policy": "fail_run",
            "numeric_inference_failure_policy": "retain_finite_effect_p_ci_df_unavailable",
            "bootstrap": {"replicates": 10000, "seed": BOOTSTRAP_SEED,
                "bit_generator": "PCG64", "quantile_method": "linear",
                "group_order": list(GROUPS), "within_group_order": "original_source_column_order",
                "joint_across_endpoints_tiers_contrasts": True,
                "retain_constant_draws": True, "p_values": False},
            "omissions": {"all_source_units": True, "recompute_scores": False,
                "recompute_welch": True, "p_values_are_diagnostics_only": True,
                "sign_reversal": "strict_opposite_nonzero_signs",
                "public_stability_units": "in_contrast_only",
                "outside_contrast_unchanged_rows": "private_only"},
            "multiplicity": {"methods": ["BH", "BY"], "tier": CORE,
                "family": "approved_panel_x_three_contrasts",
                "missing_p_for_correction": 1, "unavailable_raw_and_q": None,
                "strict_q_values": False}},
        "expected_compiled_programs_sha256": None,
        "software": software_versions(),
        "implementation_pins": [{"path": rel, "sha256": file_sha(metadata_input(root,rel))}
                                for rel in (SCRIPT_REL, TEST_REL) if (root/rel).is_file()],
    }


def coverage_status(components, score_kind, parent_gate=None):
    """Identifier-only gate, exact fractions. Missing genes never become zeros."""
    if score_kind not in ("unsigned_equal_gene", "balanced_signed"):
        raise GuardError("Unknown score kind")
    wanted = ("all",) if score_kind == "unsigned_equal_gene" else ("positive", "negative")
    if tuple(component["name"] for component in components) != wanted:
        raise GuardError("Missing, extra or reordered score components")
    seen, reasons, records = set(), [], []
    for component in components:
        denominator = component["source_members"]
        indices = component["indices"]
        if type(denominator) is not int or denominator <= 0:
            raise GuardError("Source denominator must be a fixed positive integer")
        if any(type(i) is not int or i < 0 for i in indices) or len(set(indices)) != len(indices):
            raise GuardError("Repeated/invalid measured identity, not a new membership unit")
        if seen.intersection(indices):
            raise GuardError("Opposite signed components overlap in a measured identity")
        seen.update(indices)
        if len(indices) > denominator:
            raise GuardError("Mapped members exceed original source denominator")
        passed = 5 * len(indices) >= 4 * denominator
        records.append({"name": component["name"], "original_source_members": denominator,
                        "mapped_unique_members": len(indices), "coverage_gate": passed})
        if parent_gate is None and not passed:
            reasons.append("coverage_below_80_percent:" + component["name"])
        if not indices:
            reasons.append("empty_required_component:" + component["name"])
    if parent_gate is not None:
        if score_kind != "unsigned_equal_gene" or type(parent_gate) is not bool:
            raise GuardError("Invalid nonoverlap parent gate")
        if not parent_gate:
            reasons.append("parent_coverage_gate_failed")
    if len(seen) < 10:
        reasons.append("mapped_total_below_10")
    return {"score_status": "available" if not reasons else "unavailable",
            "reasons": reasons, "components": records,
            "mapped_unique_total": len(seen), "parent_gate": parent_gate,
            "gate_basis": "parent_gate_plus_10_retained" if parent_gate is not None else "original_component_denominators_plus_10_total"}


def normalize_counts(counts, expected_shape):
    original = np.asarray(counts)
    if original.shape != tuple(expected_shape) or original.ndim != 2:
        raise DataIntegrityError("Full released matrix shape changed")
    if original.dtype.kind not in "iuf":
        raise DataIntegrityError("Counts must be numeric, not boolean, complex, string or object recodings")
    if not np.isfinite(original).all() or (original < 0).any():
        raise DataIntegrityError("Nonfinite/negative counts; no deletion or zero filling")
    limit = np.uint64(2**53)
    if original.dtype.kind in "iu":
        unsafe = (original.astype(np.uint64) > limit).any()
    else:
        unsafe = (original > np.float64(2**53)).any() or not np.equal(original,np.floor(original)).all()
    if unsafe:
        raise DataIntegrityError("Noninteger or unsafe float64 integer counts before conversion")
    values = original.astype(np.float64,copy=False)
    # Exact bounded integer totals: 1024 * 2**53 < 2**64. Reject an
    # over-limit block before accumulation; adding two allowed totals is safe.
    # A float sum would incorrectly accept [2**53,1] after rounding down.
    exact_totals = np.zeros(values.shape[1],dtype=np.uint64)
    for start in range(0,values.shape[0],1024):
        block = values[start:start+1024].astype(np.uint64).sum(axis=0,dtype=np.uint64)
        if (block > limit).any():
            raise DataIntegrityError("Unsafe exact integer library total")
        exact_totals += block
        if (exact_totals > limit).any():
            raise DataIntegrityError("Unsafe exact integer library total")
    if (exact_totals == 0).any():
        raise DataIntegrityError("Nonpositive full-universe library total")
    totals = exact_totals.astype(np.float64)
    logged = np.log1p((values / totals[None, :]) * 1000000.0) / math.log(2.0)
    if not np.isfinite(logged).all():
        raise DataIntegrityError("Nonfinite transformed expression")
    return logged, totals


def score_program(logged, specification):
    if specification["coverage"]["score_status"] != "available":
        return None
    components = specification["components"]
    scores = []
    for component in components:
        indices = component["indices"]
        if not indices or max(indices) >= logged.shape[0]:
            raise GuardError("Missing/out-of-range score members")
        scores.append(np.mean(logged[indices, :], axis=0))
    result = scores[0] if specification["score_kind"] == "unsigned_equal_gene" else 0.5*scores[0] - 0.5*scores[1]
    if not np.isfinite(result).all():
        raise DataIntegrityError("Nonfinite score: no complete-case or nanmean fallback")
    return result


def finite_vector(values):
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or not np.isfinite(result).all():
        raise DataIntegrityError("Score vectors must be finite and one-dimensional")
    return result


def welch_summary(a, b, zero_variance_policy, confidence=0.95):
    if zero_variance_policy not in ZERO_POLICIES or confidence != 0.95:
        raise GuardError("Welch policy must be explicitly fixed")
    a, b = finite_vector(a), finite_vector(b)
    n_a, n_b = len(a), len(b)
    constant_a = bool(n_a and np.all(a == a[0]))
    constant_b = bool(n_b and np.all(b == b[0]))
    mean_a = float(a[0]) if constant_a else float(np.mean(a)) if n_a else None
    mean_b = float(b[0]) if constant_b else float(np.mean(b)) if n_b else None
    result = {"n_a": n_a, "n_b": n_b, "mean_a": mean_a, "mean_b": mean_b,
        "variance_a": None, "variance_b": None, "effect": None,
        "variance_a_status": "not_estimated", "variance_b_status": "not_estimated",
        "effect_status": "unavailable_no_group_observations", "welch_status": "unavailable_n_below_2",
        "standard_error": None, "degrees_of_freedom": None, "t_statistic": None,
        "p_raw": None, "log_p_raw": None, "p_underflow_to_zero": False,
        "ci_low": None, "ci_high": None, "confidence": confidence,
        "interval_type": "pointwise_not_multiplicity_adjusted", "zero_variance_groups": [],
        "numeric_inference_failure": None}
    if n_a and n_b:
        result.update(effect=float(mean_a-mean_b), effect_status="available")
        if not math.isfinite(result["effect"]):
            raise DataIntegrityError("Nonfinite descriptive arithmetic: no finite effect to retain")
    if n_a < 2 or n_b < 2:
        return result
    with np.errstate(over="ignore",invalid="ignore",under="ignore"):
        v_a = 0.0 if constant_a else float(np.var(a, ddof=1))
        v_b = 0.0 if constant_b else float(np.var(b, ddof=1))
    for label, variance, constant in (("a",v_a,constant_a),("b",v_b,constant_b)):
        if constant:
            status = "exact_zero"
            result["zero_variance_groups"].append(label)
        elif variance == 0:
            status = "unavailable_numeric_variance_underflow"
        elif not math.isfinite(variance) or variance < 0:
            status = "unavailable_numeric_variance"
        else:
            status = "finite_positive"
        result["variance_"+label+"_status"] = status
        result["variance_"+label] = variance if status in ("exact_zero","finite_positive") else None
    unavailable_variances = [result["variance_"+label+"_status"] for label in ("a","b")
                             if result["variance_"+label] is None]
    if unavailable_variances:
        result.update(welch_status=unavailable_variances[0],numeric_inference_failure=unavailable_variances[0])
        return result
    if (v_a == 0 and v_b == 0) or (zero_variance_policy == ZERO_POLICIES[0] and (v_a == 0 or v_b == 0)):
        result["welch_status"] = "unavailable_zero_variance_policy"
        return result
    # Scale before squaring. One exact zero group is ordinary Welch with the
    # other group's actual df. Numeric inference failures do not erase effects.
    try:
        scale = max(v_a,v_b)
        sa, sb = (v_a/scale)/n_a, (v_b/scale)/n_b
        summed = sa+sb
        weight_a, weight_b = sa/summed, sb/summed
        denominator = weight_a*weight_a/(n_a-1) + weight_b*weight_b/(n_b-1)
        df = 1.0/denominator
        se = math.sqrt(scale)*math.sqrt(summed)
        if not math.isfinite(se) or se <= 0 or not math.isfinite(df) or df <= 0:
            raise FloatingPointError("Invalid Welch standard error or df")
        t_value = result["effect"]/se
        p = float(2*stats.t.sf(abs(t_value), df))
        log_p = float(math.log(2) + stats.t.logsf(abs(t_value), df))
        critical = float(stats.t.ppf((1+confidence)/2, df))
        interval = [result["effect"]-critical*se, result["effect"]+critical*se]
        if not all(math.isfinite(v) for v in (df,se,t_value,p,critical,*interval)) or not 0 <= p <= 1:
            raise FloatingPointError("Invalid Welch arithmetic")
    except (FloatingPointError, OverflowError, ZeroDivisionError):
        result.update(welch_status="unavailable_numeric_inference",
                      numeric_inference_failure="nonfinite_or_invalid_inference_arithmetic")
        return result
    if p == 0:
        result.update(welch_status="unavailable_numeric_tail_underflow",p_underflow_to_zero=True,
                      numeric_inference_failure="positive_tail_probability_not_representable")
        return result  # Raw P/CI/df remain null, never epsilon or observed zero.
    result.update(welch_status="available", standard_error=se, degrees_of_freedom=df,
                  t_statistic=t_value, p_raw=p, log_p_raw=log_p if math.isfinite(log_p) else None,
                  ci_low=interval[0], ci_high=interval[1])
    return result


def adjust_p_values(values, method):
    if method not in ("BH", "BY"):
        raise GuardError("Only the prespecified BH and BY corrections exist")
    raw = []
    for value in values:
        if value is None:
            raw.append(1.0)
        elif isinstance(value, bool) or not math.isfinite(float(value)) or not 0 <= float(value) <= 1:
            raise DataIntegrityError("Invalid P must not be disguised as planned missingness")
        else:
            raw.append(float(value))
    if not raw:
        raise GuardError("Empty multiplicity family")
    m = len(raw)
    order = np.argsort(np.asarray(raw), kind="stable")
    harmonic = math.fsum(1/i for i in range(1,m+1)) if method == "BY" else 1.0
    ranked = np.asarray(raw)[order] * (m*harmonic/np.arange(1,m+1))
    ranked = np.minimum(1.0, np.minimum.accumulate(ranked[::-1])[::-1])
    restored = np.empty(m)
    restored[order] = ranked
    return [float(x) for x in restored]


def correction_family(rows, endpoint_order):
    expected = [(program, contrast[0]) for program in endpoint_order for contrast in CONTRASTS]
    if any(row["interface"] not in INTERFACES for row in rows):
        raise GuardError("Unknown identity tier in correction inventory")
    strict_keys = [(row["program_id"],row["contrast_id"]) for row in rows if row["interface"] == STRICT]
    if strict_keys and (len(strict_keys)!=len(set(strict_keys)) or set(strict_keys)!=set(expected)):
        raise GuardError("Incomplete Strict diagnostic inventory")
    core = [row for row in rows if row["interface"] == CORE]
    keys = [(row["program_id"], row["contrast_id"]) for row in core]
    if len(keys) != len(set(keys)) or set(keys) != set(expected):
        raise GuardError("Incomplete, duplicated or extra Core hypothesis family")
    by_key = dict(zip(keys, core))
    ordered = [by_key[key] for key in expected]
    p_values = [row["p_raw"] for row in ordered]
    bh, by = adjust_p_values(p_values,"BH"), adjust_p_values(p_values,"BY")
    for row, q_bh, q_by in zip(ordered,bh,by):
        missing = row["p_raw"] is None
        row.update(p_for_correction=1.0 if missing else row["p_raw"], p_was_substituted=missing,
            q_bh=None if missing else q_bh, q_by=None if missing else q_by,
            family_id="Core_full_panel_three_contrasts", family_size=len(expected),
            discovery_eligible=not missing)
    for row in rows:
        if row["interface"] == STRICT:
            row.update(p_for_correction=None, p_was_substituted=False, q_bh=None, q_by=None,
                       family_id=None, family_size=None, discovery_eligible=False)
    return {"family_size":len(expected), "missing_raw_p":sum(value is None for value in p_values),
            "p_vector_for_correction":[1.0 if value is None else value for value in p_values],
            "q_bh_numeric_including_placeholders":bh, "q_by_numeric_including_placeholders":by,
            "ordered_keys":[list(key) for key in expected]}


def stratified_draws(units, replicates, seed):
    if type(replicates) is not int or replicates <= 0 or type(seed) is not int or seed < 0:
        raise GuardError("Bootstrap replicate count and seed must be fixed nonnegative integers")
    rng = np.random.Generator(np.random.PCG64(seed))
    draws = {}
    for group in GROUPS:
        indices = np.asarray([i for i, unit in enumerate(units) if unit["group"] == group], dtype=int)
        if not len(indices):
            raise GuardError("Bootstrap diagnostic group is empty")
        draws[group] = indices[rng.integers(0,len(indices),size=(replicates,len(indices)))]
    return draws


def bootstrap_differences(scores, draws):
    if scores is None:
        return None
    scores = finite_vector(scores)
    means = {}
    for group,indices in draws.items():
        selected = scores[indices]
        group_means = np.mean(selected,axis=1)
        constant_rows = np.all(selected == selected[:,0,None],axis=1)
        group_means[constant_rows] = selected[constant_rows,0]
        means[group] = group_means
    differences = {key:means[a]-means[b] for key,a,b in CONTRASTS}
    if any(not np.isfinite(values).all() for values in differences.values()):
        raise DataIntegrityError("Invalid bootstrap draw; no draw deletion")
    return differences


def bootstrap_summary(values, replicates):
    if values is None:
        return {"status":"unavailable_score", "requested_replicates":replicates,
                "valid_replicates":0, "ci_low":None, "ci_high":None,
                "method":"percentile_linear_95", "p_value_computed":False}
    values = finite_vector(values)
    if len(values) != replicates:
        raise DataIntegrityError("Incomplete bootstrap draws")
    lo,hi = np.quantile(values,[0.025,0.975],method="linear")
    return {"status":"available_degenerate" if lo == hi else "available",
            "requested_replicates":replicates,"valid_replicates":len(values),
            "ci_low":float(lo),"ci_high":float(hi),"method":"percentile_linear_95",
            "p_value_computed":False}


def omission_summaries(scores, units, contrast, zero_policy):
    key, a, b = contrast
    if scores is not None:
        scores = finite_vector(scores)
    base = None if scores is None else welch_summary(
        scores[[i for i,u in enumerate(units) if u["group"]==a]],
        scores[[i for i,u in enumerate(units) if u["group"]==b]], zero_policy)
    results = []
    for omitted, unit in enumerate(units):
        involved = unit["group"] in (a,b)
        result = {"unit_id":unit["unit_id"], "unit_group":unit["group"],
                  "omitted_index":omitted,"in_contrast":involved,
                  "contrast_id":key,"effect_delta":None,"strict_sign_reversal":None,
                  "sensitivity_only":True, "discovery_q_values":None}
        if scores is None:
            result.update(status="unavailable_score", statistics=None)
        else:
            stat = dict(base) if not involved else welch_summary(
                scores[[i for i,u in enumerate(units) if i != omitted and u["group"]==a]],
                scores[[i for i,u in enumerate(units) if i != omitted and u["group"]==b]],zero_policy)
            if stat["effect"] is not None and base["effect"] is not None:
                result.update(effect_delta=stat["effect"]-base["effect"],
                              strict_sign_reversal=((stat["effect"] < 0 < base["effect"]) or
                                                    (base["effect"] < 0 < stat["effect"])))
            result.update(status="available" if involved else "unchanged_outside_contrast", statistics=stat)
        results.append(result)
    return results


def load_units(root=ROOT):
    rows = read_json(verify_pin(root, METADATA_PINS[SAMPLE_JOIN]))
    if len(rows) != SHAPE[1] or [r["matrix_column_index_0based"] for r in rows] != list(range(SHAPE[1])):
        raise GuardError("Accepted source-sample axis changed")
    units = [{"unit_id":r["matrix_column_id_as_submitted"],
              "group":r["condition_as_submitted"], "column_index":r["matrix_column_index_0based"],
              "geo_accession":r["geo_accession"], "unit_type":"source_patient_sample",
              "verified_person_id":r["separate_published_patient_identifier"]} for r in rows]
    if Counter(u["group"] for u in units) != Counter(GROUP_COUNTS):
        raise GuardError("No sample-group recoding is allowed")
    for key in ("unit_id","geo_accession"):
        if len({u[key] for u in units}) != len(units):
            raise GuardError("Duplicated source unit/library")
    person_ids = [u["verified_person_id"] for u in units if u["verified_person_id"] is not None]
    if len(person_ids) != len(set(person_ids)):
        raise GuardError("Known repeated persons block independent-unit inference")
    return units


def verify_closed_metadata(root=ROOT):
    for relative, digest in CLOSED_VERSIONED_PINS.items():
        verify_pin(root,{"path":relative,"sha256":digest})
    for pin in METADATA_PINS.values():
        verify_pin(root,pin)
    summary = read_json(metadata_input(root,IDENTIFIER_JSON))
    if summary["source_feature_axis"]["released_source_rows"] != SHAPE[0]:
        raise GuardError("Released row universe changed")
    for relative,pin in summary["source_input_pins"].items():
        if relative == RAW_MATRIX_PIN["path"]:
            if pin["sha256"] != RAW_MATRIX_PIN["sha256"]:
                raise GuardError("Closed raw-count lineage pin changed")
            continue  # Real payload bytes are deferred to authorized execution.
        verify_pin(root,{"path":relative,"sha256":pin["sha256"],"bytes":pin["bytes"]})
    return summary


def load_legacy_specs(summary, root=ROOT):
    with metadata_input(root,LEGACY_INDICES).open(newline="",encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    coverage_rows = {(r["interface"],r["program_id"]):r for r in summary["program_coverage"]}
    specifications = {}
    for interface in INTERFACES:
        mapped_sets = {}
        for program, denominator in LEGACY_DENOMINATORS.items():
            selected = [r for r in rows if r["interface"]==interface and r["program_id"]==program]
            indices = [int(r["source_feature_index_0based"]) for r in selected]
            genes = [r["ensembl_gene_id"] for r in selected]
            record = coverage_rows[(interface,program)]
            if len(set(genes)) != len(genes) or len(indices) != record["mapped_unique_ensembl_ids"]:
                raise GuardError("Closed legacy identity count/uniqueness mismatch")
            if record["intended_source_members_after_ETS2_exclusion"] != denominator:
                raise GuardError("Closed original denominator changed")
            fingerprint = hashlib.sha256(("\n".join(sorted(genes))+"\n").encode()).hexdigest()
            if fingerprint != record["mapped_ensembl_ids_sha256"]:
                raise GuardError("Closed mapped gene identity set changed")
            if any(i < 0 or i >= SHAPE[0] for i in indices):
                raise GuardError("Closed feature index outside full matrix")
            mapped_sets[program] = set(indices)
            component = {"name":"all","source_members":denominator,"indices":indices}
            parent_gate = None if program != NONOVERLAP else specifications[(interface,"ets2_g1_dn")]["coverage"]["score_status"]=="available"
            gate = coverage_status([component],"unsigned_equal_gene",parent_gate)
            if (gate["score_status"]=="available") != record["mapping_gate_pass"]:
                raise GuardError("Closed eligibility gate disagrees")
            specifications[(interface,program)] = {"program_id":program,"interface":interface,
                "score_kind":"unsigned_equal_gene","source_direction":"retained_as_published_not_inverted",
                "components":[component],"coverage":gate,"source_kind":"accepted_legacy"}
        expected = mapped_sets["ets2_g1_dn"] - set().union(*(mapped_sets[p] for p in COMPARATORS))
        if mapped_sets[NONOVERLAP] != expected:
            raise GuardError("Nonoverlap set differs from parent minus all comparator sets")
    return specifications


def map_external_rows(rows, endpoint, feature_rows):
    """Exact source-symbol join to the full-reference audited RNA interface.

    Published Entrez/source IDs are provenance only; never alias bridges or weights.
    The complete feature audit already tested reciprocal identities globally.
    """
    if endpoint["score_kind"] not in ("unsigned_equal_gene","balanced_signed"):
        raise GuardError("Unresolved external score convention")
    expected = endpoint["original_component_members"]
    wanted = ["all"] if endpoint["score_kind"]=="unsigned_equal_gene" else ["positive","negative"]
    if set(expected) != set(wanted) or set(endpoint["direction_literals"]) != set(wanted):
        raise GuardError("Signed source components/polarity unresolved")
    if endpoint.get("target_gene_exclusions",[]) != []:
        raise GuardError("This implementation does not extend target exclusion to external sources")
    if endpoint["source_identifier_namespace"] not in ("MSigDB_geneSymbols","HGNC_symbol"):
        raise GuardError("Unreviewed external identifier namespace")
    labels = {r["label_as_submitted"]:r for r in feature_rows}
    if len(labels) != len(feature_rows):
        raise GuardError("Duplicate submitted labels cannot be collapsed")
    components = {name:[] for name in wanted}
    source_values = set()
    order = []
    for row in rows:
        if row["program_id"] != endpoint["program_id"]:
            raise GuardError("Unexpected source program row")
        if row["source_identifier_namespace"] != endpoint["source_identifier_namespace"]:
            raise GuardError("Source namespace changed")
        if row["membership_weight"] != "1":
            raise GuardError("Only literal unit membership weights are permitted")
        value = row["source_value"]
        if not value or value in source_values or value != value.strip():
            raise GuardError("Missing/duplicate/normalized source symbol; no automatic collapse")
        source_values.add(value)
        order.append(int(row["source_member_index"]))
        matched = [name for name in wanted if row["published_direction"]==endpoint["direction_literals"][name]]
        if len(matched) != 1:
            raise GuardError("Unknown/ambiguous source polarity")
        components[matched[0]].append(value)
    if order != list(range(1,len(rows)+1)):
        raise GuardError("Official source membership order changed")
    if {name:len(values) for name,values in components.items()} != expected:
        raise GuardError("Original unique source-member denominators changed")
    specs, audit = {}, []
    for interface in INTERFACES:
        mapped = []
        for name in wanted:
            indices = []
            for value in components[name]:
                feature = labels.get(value)
                status = "no_exact_submitted_symbol" if feature is None else feature[interface]
                index = feature["source_row_index_0based"] if status=="eligible" else None
                if index is not None:
                    indices.append(index)
                audit.append({"program_id":endpoint["program_id"],"interface":interface,
                              "component":name,"source_value":value,"mapping_status":status,
                              "feature_index":index})
            mapped.append({"name":name,"source_members":expected[name],"indices":indices})
        gate = coverage_status(mapped,endpoint["score_kind"])
        specs[(interface,endpoint["program_id"])] = {"program_id":endpoint["program_id"],
            "interface":interface,"score_kind":endpoint["score_kind"],"components":mapped,
            "coverage":gate,"source_kind":"external_symbol_csv",
            "source_direction":endpoint["direction_literals"],
            "target_gene_membership": {"literal_symbol":"ETS2",
                "present_in_original_source": "ETS2" in source_values, "excluded":False}}
    return specs,audit


def validate_panel(plan):
    panel = plan.get("panel")
    if not isinstance(panel,list) or not panel:
        raise GuardError("Missing fixed endpoint panel")
    ids = [e["program_id"] for e in panel]
    if any(not isinstance(value,str) or not re.fullmatch(r"[A-Za-z0-9_]+",value) for value in ids):
        raise GuardError("Unsafe/nonliteral program identifier")
    if tuple(ids) != TEN_PROGRAMS:
        raise GuardError("The root-fixed ten-program panel/order changed; no endpoint fallback")
    for endpoint in panel:
        if endpoint.get("status") != "qualified" or not endpoint.get("source_approval_reference"):
            raise GuardError("Unqualified/HOLD/unapproved source cannot be an executable endpoint")
        if endpoint["program_id"] in LEGACY_DENOMINATORS:
            if endpoint["kind"] != "accepted_legacy" or endpoint["score_kind"] != "unsigned_equal_gene" or endpoint["original_source_members"] != LEGACY_DENOMINATORS[endpoint["program_id"]]:
                raise GuardError("A fixed legacy endpoint was changed")
        elif endpoint.get("kind") != "external_symbol_csv":
            raise GuardError("Unknown source interface; no fallback")
    return ids




def validate_external_source(endpoint, root=ROOT):
    """Bind the approved source handoff, manifest and CSV; no generic reader.

    The source worker handoff binds the exact qualified manifest/CSV versions.
    Qualification is not root panel approval. Only ECM has a qualified adapter;
    a future signed source requires a new qualified, independently reviewed entry.
    """
    if endpoint.get("program_id") not in EXTERNAL_SOURCE_REGISTRY:
        raise GuardError("No qualified/whitelisted external source adapter; HOLD is not a program")
    accepted = EXTERNAL_SOURCE_REGISTRY[endpoint["program_id"]]
    for key in ("program_id","source_manifest_pin","membership_pin","score_kind",
                "source_identifier_namespace","original_component_members","direction_literals",
                "target_gene_exclusions","source_version"):
        if endpoint.get(key) != accepted[key]:
            raise GuardError("External source semantic/pin binding changed: " + key)
    handoff = read_json(verify_pin(root,accepted["schema_handoff_pin"]))
    expected = {"schema_version":1,"program_id":accepted["program_id"],
        "config_path":accepted["source_manifest_pin"]["path"],
        "final_config_sha256":accepted["source_manifest_pin"]["sha256"],
        "membership_path":accepted["membership_pin"]["path"],
        "membership_sha256":accepted["membership_pin"]["sha256"],
        "membership_bytes":accepted["membership_pin"]["bytes"],
        "membership_rows":321,"complete_original_denominator":321,
        "original_source_identifier_count":340,"source_direction":"UNSIGNED",
        "source_version":accepted["source_version"],"scoring_or_endpoint_authorized":False,
        "source_qualification_status":"FINAL_partial_ECM_qualified_epithelial_IL17A_response_HOLD"}
    if any(handoff.get(key)!=value for key,value in expected.items()):
        raise GuardError("Qualified source handoff does not bind the manifest and membership semantics")
    manifest = read_json(verify_pin(root,accepted["source_manifest_pin"]))
    if manifest.get("schema_version") != 1 or manifest.get("qualification_status") != "partial_ECM_qualified_epithelial_IL17A_response_HOLD" or manifest.get("source_pass_closed") is not True:
        raise GuardError("Source manifest is not the closed partial qualification")
    membership = manifest.get("membership_output",{})
    if any(membership.get(key)!=value for key,value in {**accepted["membership_pin"],"rows":321}.items()):
        raise GuardError("Source manifest points to a different membership CSV")
    programs = manifest.get("programs",[])
    if len(programs)!=1 or programs[0].get("program_id")!=accepted["program_id"]:
        raise GuardError("Qualified manifest endpoint inventory changed")
    expected_program = {"status":"qualified_complete_official_pathway_membership_source_only",
        "official_program_name":"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
        "msigdb_systematic_id":"M610","exact_source_id":"R-HSA-1474244",
        "collection":"C2:CP:REACTOME","organism":"Homo sapiens",
        "source_identifier_namespace":"Human_Ensembl_Gene_ID",
        "msigdb_version":"2026.1.Hs","reactome_release":"95","member_namespace":"MSigDB_geneSymbols",
        "expected_member_count":321,"expected_source_identifier_count":340,"direction":"UNSIGNED",
        "membership_weight":"1","published_response_directions_available":False,"scoring_authorized":False,
        "source_keys":{"html":"ecm_html","json":"ecm_json"}}
    if any(programs[0].get(key)!=value for key,value in expected_program.items()):
        raise GuardError("Qualified manifest source semantics changed")
    for key,pin in EXTERNAL_RAW_SOURCE_PINS.items():
        source = manifest.get("sources",{}).get(key,{})
        if source.get("complete") is not True or source.get("status") != 200 or any(source.get(field)!=value for field,value in pin.items()):
            raise GuardError("Complete original ECM source pin changed")
        verify_pin(root,pin)
    return verify_pin(root,accepted["membership_pin"])


def compile_programs(plan, summary, root=ROOT):
    ids = validate_panel(plan)
    specifications = load_legacy_specs(summary,root)
    external_audit = []
    feature_rows = None
    for endpoint in plan["panel"]:
        if endpoint["kind"] == "accepted_legacy":
            continue
        source_path = validate_external_source(endpoint,root)
        if feature_rows is None:
            feature_rows = read_json(metadata_input(root,FEATURE_AUDIT))
            if len(feature_rows) != SHAPE[0] or [r["source_row_index_0based"] for r in feature_rows] != list(range(SHAPE[0])):
                raise GuardError("Complete feature-audit axis changed")
        with source_path.open(newline="",encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            required = {"program_id","source_member_index","source_value","source_identifier_namespace",
                        "published_direction","membership_weight"}
            if not required.issubset(reader.fieldnames or []) or len(reader.fieldnames)!=len(set(reader.fieldnames)):
                raise GuardError("External source schema changed")
            rows = [{key:row[key] for key in required} for row in reader]
        extra, audit = map_external_rows(rows,endpoint,feature_rows)
        specifications.update(extra)
        external_audit.extend(audit)
    ordered = [specifications[(interface,program)] for interface in INTERFACES for program in ids]
    for program in ids:
        core_indices = set(i for component in specifications[(CORE,program)]["components"] for i in component["indices"])
        strict_indices = set(i for component in specifications[(STRICT,program)]["components"] for i in component["indices"])
        if not strict_indices.issubset(core_indices):
            raise GuardError("Strict identity sensitivity is not a subset of Core")
    return ordered, external_audit


def source_only_preflight(plan=None, root=ROOT):
    plan = draft_plan(root) if plan is None else plan
    summary = verify_closed_metadata(root)
    units = load_units(root)
    specifications, external_audit = compile_programs(plan,summary,root)
    approved_panel = plan.get("panel_approved") is True
    records = []
    for specification in specifications:
        for key,a,b in CONTRASTS:
            records.append({"program_id":specification["program_id"],"interface":specification["interface"],
                "contrast_id":key,"group_a":a,"group_b":b,"coverage":specification["coverage"],
                "effect_status":"not_executed_source_only", "welch_status":"not_executed_source_only",
                "bootstrap_status":"not_executed_source_only","omission_status":"not_executed_source_only",
                "effect":None,"p_raw":None,"p_for_correction":None,"q_bh":None,"q_by":None})
    receipt = {"schema_version":1,"mode":"source_only_preflight",
        "status":"complete_for_source_only_mode_not_an_execution_approval",
        "panel_approved":approved_panel,"programs_inspected":len(plan["panel"]),
        "prospectively_approved_core_family_size":3*len(plan["panel"]) if approved_panel else None,
        "supplied_plan_state":plan.get("state"),"execution_authorization_checked":False,
        "root_method_decision_reference":ROOT_METHOD_DECISION,
        "external_target_gene_membership":{spec["program_id"]:spec["target_gene_membership"]
            for spec in specifications if spec["interface"]==CORE and "target_gene_membership" in spec},
        "source_patient_samples":len(units),"groups":dict(Counter(u["group"] for u in units)),
        "released_row_universe":SHAPE[0],"independence_status":"article_claim_only_crosswalk_unavailable",
        "compiled_programs_sha256":canonical_sha(specifications),"endpoint_status_rows":records,
        "external_candidates":plan.get("external_candidates",[]),
        "research_question":"narrower_unsigned_ETS2_comparator_ECM_tissue_context",
        "does_not_fulfill_held_IL17_response_transfer":True,
        "not_independent_mechanistic_confirmation_of_MacroMap":True,
        "external_membership_rows_read":len(external_audit)//2,
        "real_expression_payload_opened":False,"real_scores_computed":False,
        "real_analysis_executed":False,"root_freeze_written":False,"network_used":False,
        "expression_payload_pin_checks":"deferred_until_authorized_execution_before_array_loading",
        "biological_conclusions":"not_assessed","software":software_versions()}
    return receipt,specifications,units,external_audit


def validate_method(plan):
    method = plan["method"]
    expected = draft_plan_method_without_runtime()
    if method != expected or type(method["bootstrap"]["seed"]) is not int:
        raise GuardError("Changed/unresolved root-approved fixed statistical contract")


def draft_plan_method_without_runtime():
    # Method proposal has no input access and cannot grant authorization.
    # Supply a nonexistent root to avoid even hashing implementation files.
    return draft_plan(Path("/__psc_liver_no_project_inputs__"))["method"]


def fresh_output_path(output_directory, root=ROOT):
    root = Path(root).resolve()
    path = Path(output_directory)
    if not path.is_absolute():
        path = root/path
    allowed = root/RUNS_REL
    if path.exists() or path.is_symlink():
        raise GuardError("Output directory already exists, including empty/failed runs")
    if not path.resolve().is_relative_to(allowed.resolve()) or path.resolve()==allowed.resolve():
        raise GuardError("Real outputs must use a fresh private run subdirectory")
    for candidate in (allowed,*allowed.parents):
        if candidate == root.parent:
            break
        if candidate.is_symlink():
            raise GuardError("Private output ancestry cannot be a symlink")
    for parent in path.parents:
        if parent == root:
            break
        if parent.is_symlink():
            raise GuardError("Output symlink traversal rejected")
    return path




def checked_plan_path(path, root=ROOT):
    """Plan JSON only, never a payload or an arbitrary source/QC JSON reader."""
    root = Path(root).resolve()
    path = Path(path)
    path = path if path.is_absolute() else root/path
    allowed = {root/"config/psc-liver-program-plan.json",
               root/"config/psc-liver-program-freeze.json",
               root/WORK_REL/"proposed-plan.json"}
    if path.resolve() not in allowed or path.is_symlink():
        raise GuardError("Plan input must be a designated liver plan JSON, not a data/source payload")
    for parent in path.parents:
        if parent == root:
            break
        if parent.is_symlink():
            raise GuardError("Plan symlink traversal rejected")
    return metadata_input(root,path.relative_to(root))


def authorize_execution(plan_path, expected_sha, acknowledgement, output_directory, root=ROOT):
    """All authorization and pin gates precede numerical matrix loading."""
    if acknowledgement != EXECUTION_ACK:
        raise GuardError("Explicit future real-execution acknowledgement required")
    if not re.fullmatch(r"[a-f0-9]{64}",str(expected_sha or "")):
        raise GuardError("Root frozen-plan SHA-256 required")
    plan_path = checked_plan_path(plan_path,root)
    if file_sha(plan_path) != expected_sha:
        raise GuardError("Frozen plan hash changed")
    plan = read_json(plan_path)
    if plan.get("schema_version") != 1 or plan.get("analysis_id") != "GSE303271-fixed-programs" or plan.get("state") != "frozen":
        raise GuardError("Not the reviewed root-frozen plan schema")
    auth = plan.get("root_authorization",{})
    if auth.get("approver_role") != "root" or auth.get("real_execution_authorized") is not True or not auth.get("reference") or not auth.get("approved_at"):
        raise GuardError("No explicit root authorization for real expression execution")
    if plan.get("root_method_decision_reference") != ROOT_METHOD_DECISION:
        raise GuardError("Root prospective method decision is not preserved")
    if plan.get("root_panel_decision_reference") != ROOT_PANEL_DECISION or plan.get("research_question") != "narrower_unsigned_ETS2_comparator_ECM_tissue_context" or plan.get("does_not_fulfill_held_IL17_response_transfer") is not True:
        raise GuardError("Root narrower research question/panel decision is not preserved")
    if plan.get("panel_approved") is not True or set(plan.get("decisions",{})) != set(DECISION_NAMES):
        raise GuardError("Panel or scientific decisions are unresolved")
    if any(d.get("status")!="approved" or not d.get("reference") for d in plan["decisions"].values()):
        raise GuardError("A statistical/source decision is not approved")
    if plan.get("identity") != {"primary":CORE,"sensitivity":STRICT} or plan.get("contrast_order") != [list(c) for c in CONTRASTS]:
        raise GuardError("Fixed identity role or ordered contrasts changed")
    unit = plan.get("unit",{})
    if unit != {"kind":"source_patient_sample","independence_status":"article_claim_only_crosswalk_unavailable",
                "assumption_approved":True,"known_repeat_status":"no_known_repeats_crosswalk_unavailable"}:
        raise GuardError("Unapproved independent-source-patient-sample assumption or known repeats")
    dataset = plan.get("dataset",{})
    if dataset != {"accession":"GSE303271","shape":list(SHAPE),"groups":GROUP_COUNTS,
                   "matrix_pin":MATRIX_PIN,"raw_matrix_pin":RAW_MATRIX_PIN}:
        raise GuardError("Released expression source/shape/groups changed")
    validate_panel(plan)
    validate_method(plan)
    if plan.get("software") != software_versions():
        raise GuardError("Native Python/NumPy/SciPy versions changed")
    pins = plan.get("implementation_pins",[])
    if {pin.get("path") for pin in pins} != {SCRIPT_REL,TEST_REL} or len(pins)!=2:
        raise GuardError("Incomplete implementation pins")
    for pin in pins:
        verify_pin(root,pin)
    review = plan.get("independent_review",{})
    if review.get("status") != "accepted_synthetic_and_guard_review" or not review.get("reference"):
        raise GuardError("Independent implementation review is not accepted")
    review_record = read_json(verify_pin(root,review["artifact_pin"]))
    if review_record.get("status") != "pass_synthetic_and_guard_review" or review_record.get("implementation_sha256") != {p["path"]:p["sha256"] for p in pins} or review_record.get("real_execution_reviewed") is not False:
        raise GuardError("Stale or incompatible independent review receipt")
    output = fresh_output_path(output_directory,root)
    preflight, specifications, units, external_audit = source_only_preflight(plan,root)
    if plan.get("expected_compiled_programs_sha256") != canonical_sha(specifications):
        raise GuardError("Unfrozen or changed compiled memberships/coverage")
    # This is the first access to expression-payload bytes. Hashing is not loading
    # numeric values. Only a fully approved real-execution action reaches it.
    verify_pin(root,RAW_MATRIX_PIN)
    matrix_path = verify_pin(root,MATRIX_PIN)
    return plan,specifications,units,external_audit,output,matrix_path,preflight


def unavailable_statistics(n_a, n_b):
    return {"n_a":n_a,"n_b":n_b,"mean_a":None,"mean_b":None,"variance_a":None,"variance_b":None,
        "effect":None,"effect_status":"unavailable_score","welch_status":"unavailable_score",
        "variance_a_status":"unavailable_score","variance_b_status":"unavailable_score",
        "numeric_inference_failure":None,
        "standard_error":None,"degrees_of_freedom":None,"t_statistic":None,
        "p_raw":None,"log_p_raw":None,"p_underflow_to_zero":False,"ci_low":None,"ci_high":None,
        "confidence":0.95,"interval_type":"pointwise_not_multiplicity_adjusted","zero_variance_groups":[]}


def summarize_omissions(rows):
    # Outside-contrast unchanged rows are private inventory, never evidence of
    # stability. All public metrics and denominators use only involved units.
    relevant = [row for row in rows if row["in_contrast"]]
    finite = [row for row in relevant if row["statistics"] is not None and row["statistics"]["effect"] is not None]
    effects = [row["statistics"]["effect"] for row in finite]
    deltas = [row["effect_delta"] for row in finite if row["effect_delta"] is not None]
    return {"status":"available" if effects else "unavailable_effects",
        "relevant_unit_denominator":len(relevant),"in_contrast_rows":len(relevant),
        "denominator_kind":"in_contrast_units_only",
        "defined_effect_rows":len(finite),"effect_min":min(effects) if effects else None,
        "effect_max":max(effects) if effects else None,
        "max_absolute_effect_change":max(map(abs,deltas)) if deltas else None,
        "strict_sign_reversal_rows":sum(row["strict_sign_reversal"] is True for row in finite) if finite else None,
        "sensitivity_only":True}


def analyze_counts(counts, specifications, units, method, expected_shape):
    """Pure computation API. Current tests pass generated synthetic values only."""
    if len(units) != expected_shape[1] or len({u["unit_id"] for u in units}) != len(units):
        raise GuardError("Source-unit axis mismatch")
    if set(u["group"] for u in units) != set(GROUPS):
        raise GuardError("Unexpected/empty diagnostic group")
    keys = [(spec["interface"],spec["program_id"]) for spec in specifications]
    program_order = [spec["program_id"] for spec in specifications if spec["interface"]==CORE]
    expected_keys = [(interface,program) for interface in INTERFACES for program in program_order]
    if len(keys) != len(set(keys)) or set(keys) != set(expected_keys) or not program_order:
        raise GuardError("Incomplete Core/Strict program inventory")
    logged, totals = normalize_counts(counts,expected_shape)
    draws = stratified_draws(units,method["bootstrap"]["replicates"],method["bootstrap"]["seed"])
    groups = {group:[i for i,u in enumerate(units) if u["group"]==group] for group in GROUPS}
    score_vectors, bootstrap_vectors, rows, omissions = {}, {}, [], []
    for spec in specifications:
        key = (spec["interface"],spec["program_id"])
        scores = score_program(logged,spec)
        score_vectors[key] = scores
        boot = bootstrap_differences(scores,draws)
        for contrast in CONTRASTS:
            contrast_id,a,b = contrast
            stat = unavailable_statistics(len(groups[a]),len(groups[b])) if scores is None else welch_summary(
                scores[groups[a]],scores[groups[b]],method["zero_variance_policy"],method["confidence"])
            values = None if boot is None else boot[contrast_id]
            bootstrap_vectors[(*key,contrast_id)] = values
            loo = omission_summaries(scores,units,contrast,method["zero_variance_policy"])
            for entry in loo:
                entry.update(interface=key[0],program_id=key[1])
            omissions.extend(loo)
            row = {"program_id":key[1],"interface":key[0],"contrast_id":contrast_id,
                "group_a":a,"group_b":b,"contrast_role":"primary" if contrast_id==CONTRASTS[0][0] else "secondary",
                "inferential_role":"sole_discovery_family" if key[0]==CORE else "identity_sensitivity_only",
                "score_kind":spec["score_kind"],"coverage":spec["coverage"],**stat,
                "target_gene_membership":spec.get("target_gene_membership"),
                "bootstrap":bootstrap_summary(values,method["bootstrap"]["replicates"]),
                "omission":summarize_omissions(loo)}
            rows.append(row)
    expected_rows = {(i,p,c[0]) for i in INTERFACES for p in program_order for c in CONTRASTS}
    if {(r["interface"],r["program_id"],r["contrast_id"]) for r in rows} != expected_rows or len(rows)!=len(expected_rows):
        raise GuardError("Incomplete endpoint/contrast status rows")
    expected_omissions = len(units)*len(expected_rows)
    if len(omissions) != expected_omissions:
        raise GuardError("Incomplete all-unit omissions")
    correction = correction_family(rows,program_order)
    public = {"schema_version":1,"status":"complete_for_declared_analysis",
        "analysis_units":"source_patient_samples_not_verified_unique_donors",
        "independence_status":"article_claim_only_crosswalk_unavailable",
        "group_sizes":{group:len(indices) for group,indices in groups.items()},
        "released_normalization_rows":expected_shape[0],"normalization":"all_released_rows_CPM_then_log2_plus_1",
        "endpoint_order":program_order,"endpoint_rows":rows,
        "expected_endpoint_rows":len(expected_rows),"expected_all_unit_omission_rows":expected_omissions,
        "core_family_size":correction["family_size"],"unavailable_core_p_values":correction["missing_raw_p"],
        "public_omission_denominators":{key:len(groups[a])+len(groups[b]) for key,a,b in CONTRASTS},
        "all_unit_omission_row_count_is_private_inventory_not_stability_denominator":True,
        "strict_discoveries_permitted":False,
        "does_not_fulfill_held_IL17_response_transfer":True,
        "not_independent_mechanistic_confirmation_of_MacroMap":True,
        "interpretation_limits":["Whole-biopsy observational transcript-abundance composites, not pathway activity or cell-intrinsic effects.",
            "Pointwise intervals are not simultaneous. BH/BY do not repair invalid marginal tests or unverified independence.",
            "No clinical/stage/batch adjustment, causal interpretation, clinical outcome, or clinical recommendation."]}
    return {"public":public,"scores":score_vectors,"library_totals":totals,"draws":draws,
            "bootstrap_vectors":bootstrap_vectors,"omission_rows":omissions,"correction":correction}


def load_count_matrix(path):
    """Real numerical loading exists in one visibly guarded CLI call site."""
    return np.load(path,allow_pickle=False)


def failed_endpoint_inventory(plan, reason):
    return [{"program_id":endpoint["program_id"],"interface":interface,"contrast_id":contrast[0],
             "effect_status":"not_available_run_failed","welch_status":"not_available_run_failed",
             "bootstrap_status":"not_available_run_failed","omission_status":"not_available_run_failed",
             "effect":None,"p_raw":None,"p_for_correction":None,"q_bh":None,"q_by":None,
             "reason":reason} for interface in INTERFACES for endpoint in plan["panel"] for contrast in CONTRASTS]


def execute(plan_path, expected_sha, acknowledgement, output_directory, root=ROOT):
    plan,specs,units,external_audit,output,matrix_path,preflight = authorize_execution(
        plan_path,expected_sha,acknowledgement,output_directory,root)
    # No output directory exists on any guard refusal. Atomic mkdir prevents reuse.
    output.parent.mkdir(parents=True,exist_ok=True)
    output.mkdir(mode=0o700,exist_ok=False)
    access = {"expression_payload_numeric_loading_attempted":False,
              "real_expression_values_loaded":False,"real_scores_computed":False,
              "real_computation_attempted":False,"all_score_statistics_completed":False,
              "correction_completed":False}
    dump_json(output/"RUNNING.json",{"status":"running_not_complete","plan_sha256":expected_sha})
    try:
        plan_bytes = checked_plan_path(plan_path,root).read_bytes()
        if hashlib.sha256(plan_bytes).hexdigest() != expected_sha:
            raise GuardError("Root plan changed after guard checks")
        (output/"private-frozen-plan-original.json").write_bytes(plan_bytes)
        access["expression_payload_numeric_loading_attempted"] = True
        counts = load_count_matrix(matrix_path)
        access["real_expression_values_loaded"] = True
        access["real_computation_attempted"] = True
        # If computation fails midway, null explicitly means partial computation
        # is unknown, not an assertion that no scores were ever computed.
        access["real_scores_computed"] = None
        access["correction_completed"] = None
        result = analyze_counts(counts,specs,units,plan["method"],SHAPE)
        access.update(real_scores_computed=True,all_score_statistics_completed=True,correction_completed=True)
        # Unit-level values and every omission remain private. Public summary has
        # group-level/endpoint aggregates only, never omitted unit identities.
        private_scores = [{"interface":interface,"program_id":program,
                           "values":None if values is None else values.tolist()}
                          for (interface,program),values in result["scores"].items()]
        dump_json(output/"private-source-units.json",units)
        dump_json(output/"private-scores.json",private_scores)
        dump_json(output/"private-all-unit-omissions.json",result["omission_rows"])
        dump_json(output/"private-library-totals.json",result["library_totals"].tolist())
        dump_json(output/"private-correction-vector.json",result["correction"])
        dump_json(output/"private-compiled-programs.json",specs)
        dump_json(output/"private-external-identifier-audit.json",external_audit)
        np.savez_compressed(output/"private-bootstrap-draws.npz",**result["draws"])
        np.savez_compressed(output/"private-bootstrap-differences.npz",
            **{"__".join(key):value for key,value in result["bootstrap_vectors"].items() if value is not None})
        public = result["public"]
        public.update(mode="real_execution",real_analysis_executed=True,plan_sha256=expected_sha,
                      method=plan["method"],software=software_versions(),root_freeze_written=False,
                      source_pin_status="all_frozen_source_and_implementation_pins_verified",
                      external_candidates=plan.get("external_candidates",[]))
        dump_json(output/"aggregate-public-summary.json",public)
        (output/"RUNNING.json").unlink()
        artifacts = [{"path":p.name,"sha256":file_sha(p),"bytes":p.stat().st_size,
                      "visibility":"aggregate_public_candidate" if p.name.startswith("aggregate-") else "private"}
                     for p in sorted(output.iterdir()) if p.is_file()]
        receipt = {"schema_version":1,"mode":"real_execution","status":"complete",
            "plan_sha256":expected_sha,"root_authorization_reference":plan["root_authorization"]["reference"],
            "access_audit":access,"real_analysis_executed":True,"root_freeze_written":False,
            "network_used":False,"artifacts":artifacts,"endpoint_rows":len(public["endpoint_rows"]),
            "no_clinical_or_causal_conclusion":True}
        dump_json(output/"completion-receipt.json",receipt)
        return receipt
    except Exception as exc:
        dump_json(output/"FAILED.json",{"schema_version":1,"mode":"real_execution","status":"failed_not_complete",
            "plan_sha256":expected_sha,"error_type":type(exc).__name__,"reason":str(exc),
            "access_audit":access,"endpoint_status_rows":failed_endpoint_inventory(plan,type(exc).__name__),
            "correction_performed":access["correction_completed"],"biological_conclusions":"not_assessed"})
        raise


def main(argv=None, root=ROOT):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight",action="store_true",help="Identifier/sample metadata only; no expression payload access")
    action.add_argument("--write-proposal",action="store_true",help="Record the prospective root panel in an unfrozen method draft; never authorize execution")
    action.add_argument("--execute",action="store_true",help="Future separately root-authorized real execution only")
    parser.add_argument("--plan")
    parser.add_argument("--expected-plan-sha256")
    parser.add_argument("--acknowledge-real-execution")
    parser.add_argument("--output-directory")
    args = parser.parse_args(argv)
    try:
        if args.execute:
            if not all((args.plan,args.expected_plan_sha256,args.acknowledge_real_execution,args.output_directory)):
                raise GuardError("Execution requires root plan, its frozen hash, explicit acknowledgement and fresh private output path")
            result = execute(args.plan,args.expected_plan_sha256,args.acknowledge_real_execution,args.output_directory,root)
        else:
            if any((args.expected_plan_sha256,args.acknowledge_real_execution,args.output_directory)):
                raise GuardError("Execution flags are not preflight/proposal permissions")
            work = Path(root)/WORK_REL
            if work.is_symlink() or not work.resolve().is_relative_to(Path(root).resolve()):
                raise GuardError("Planning output path escapes its owned work directory")
            for parent in work.parents:
                if parent == Path(root):
                    break
                if parent.is_symlink():
                    raise GuardError("Planning output ancestry cannot be a symlink")
            work.mkdir(parents=True,exist_ok=True)
            if args.write_proposal:
                if args.plan:
                    raise GuardError("Proposal creation does not read or modify a root plan")
                result = draft_plan(root)
                dump_json(work/"proposed-plan.json",result)
            else:
                plan = read_json(checked_plan_path(args.plan,root)) if args.plan else None
                result,specs,units,audit = source_only_preflight(plan,root)
                dump_json(work/"source-only-preflight.json",result)
                dump_json(work/"source-only-compiled-programs.json",specs)
                dump_json(work/"source-only-external-identifier-audit.json",audit)
        print(json.dumps({"status":result.get("status",result.get("state")),
                          "mode":result.get("mode","proposal"),
                          "real_analysis_executed":result.get("real_analysis_executed",False)},sort_keys=True))
        return 0
    except (GuardError,DataIntegrityError,OSError,KeyError,TypeError,ValueError) as exc:
        print(json.dumps({"status":"refused_or_failed","error_type":type(exc).__name__,"reason":str(exc),
                          "no_execution_authorization_created":True},sort_keys=True),file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
