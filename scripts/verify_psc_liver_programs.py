#!/usr/bin/env python3
"""Independent read-only verifier for the prospective liver fixed-program method.

The CLI has ONLY source-only and built-in synthetic modes. It cannot open a real
count matrix or a real score table. Mathematical functions are independent of the
analysis author. Later real numerical use needs a separate root freeze/receipt.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Sequence

from scipy import special

CORE = "exact_symbol_ensembl_bijection"
STRICT = "plus_reciprocal_entrez"
CONTRASTS = (("PSC", "AIH"), ("PSC", "ASC"), ("ASC", "AIH"))
TOLERANCE = {"rtol": 1e-10, "atol": 1e-12}


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def number(value: Any) -> float:
    demand(not isinstance(value, bool), "boolean is not a numerical measurement")
    demand(isinstance(value, (int, float)), "measurement must be a number")
    result = float(value)
    demand(math.isfinite(result), "measurement must be finite")
    return result


def close(actual: Any, expected: Any, label: str = "value") -> None:
    if expected is None:
        demand(actual is None, f"{label}: missing status was not retained")
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        value = number(actual)
        demand(math.isclose(value, expected, rel_tol=TOLERANCE["rtol"], abs_tol=TOLERANCE["atol"]),
               f"{label}: expected {expected!r}, observed {actual!r}")
    else:
        demand(actual == expected, f"{label}: exact mismatch")


def logcpm_reference(counts: Sequence[Sequence[float]]) -> dict:
    """Full synthetic feature universe, genes by patient-samples; stable log1p."""
    demand(len(counts) > 0, "empty feature universe")
    width = len(counts[0])
    demand(width > 0, "empty sample universe")
    matrix = []
    for row in counts:
        demand(len(row) == width, "ragged count matrix")
        values = []
        for original in row:
            # Validate Python integer literals before float conversion; 2**53+1
            # would otherwise round to the apparently allowed value 2**53.
            if isinstance(original, int) and not isinstance(original, bool):
                demand(0 <= original <= 2**53, "count exceeds exact float64 integer range")
                value = original
            else:
                value_float = number(original)
                demand(value_float >= 0 and value_float.is_integer(),
                       "counts must be nonnegative integers")
                demand(value_float <= 2**53, "count exceeds exact float64 integer range")
                value = int(value_float)
            values.append(value)
        matrix.append(values)
    # Sum exact integers before any floating conversion. fsum([2**53,1]) is
    # rounded back to 2**53 and cannot enforce this safety boundary.
    integer_totals = [sum(row[j] for row in matrix) for j in range(width)]
    demand(all(0 < x <= 2**53 for x in integer_totals), "invalid/zero/unsafe library total")
    totals = [float(x) for x in integer_totals]
    values = [[math.log1p(row[j] / totals[j] * 1e6) / math.log(2.0)
               for j in range(width)] for row in matrix]
    return {"library_totals": totals, "logcpm": values, "feature_rows": len(matrix)}


def mapping_gate(original_members: int, mapped_members: int,
                 minimum_genes: int = 10) -> bool:
    demand(type(original_members) is int and original_members > 0, "invalid source denominator")
    demand(type(mapped_members) is int and 0 <= mapped_members <= original_members,
           "invalid mapped membership")
    demand(type(minimum_genes) is int and minimum_genes > 0, "invalid minimum")
    # Integer form avoids boundary rounding and never changes the denominator.
    return mapped_members >= minimum_genes and 5 * mapped_members >= 4 * original_members


def score_reference(logcpm: Sequence[Sequence[float]], positive: Sequence[int],
                    negative: Sequence[int] | None = None) -> list[float]:
    """Unsigned equal-gene mean; explicit optional 1/2-up minus 1/2-down.

    The signed option tests math only. It does not authorize a missing source or
    choose a signed endpoint. Both components must exist and must be disjoint.
    """
    demand(len(logcpm) > 0 and len(logcpm[0]) > 0, "empty expression universe")
    width = len(logcpm[0])
    demand(all(len(row) == width for row in logcpm), "ragged score matrix")
    for row in logcpm:
        for x in row:
            demand(number(x) >= 0, "logCPM must be nonnegative")
    components = [list(positive)] if negative is None else [list(positive), list(negative)]
    all_indices = [index for component in components for index in component]
    demand(all(components), "required component absent; never zero-fill")
    demand(all(type(x) is int and 0 <= x < len(logcpm) for x in all_indices),
           "program index is absent/out of range")
    demand(len(all_indices) == len(set(all_indices)), "repeated or conflicting member")
    means = [[math.fsum(logcpm[i][j] for i in component) / len(component)
              for j in range(width)] for component in components]
    if negative is None:
        return means[0]
    return [(means[0][j] - means[1][j]) / 2.0 for j in range(width)]


def moments(values: Sequence[float]) -> tuple[float, float]:
    demand(len(values) > 0, "empty group")
    data = [number(x) for x in values]
    mean = math.fsum(x / len(data) for x in data)
    if len(data) < 2:
        return mean, 0.0
    centered = [x - mean for x in data]
    scale = max(abs(x) for x in centered)
    # Test exact constancy using original values; summing a constant must not
    # invent tiny nonzero variance from arithmetic rounding.
    if min(data) == max(data):
        return data[0], 0.0
    demand(scale > 0 and math.isfinite(scale), "invalid centered scale")
    sd = scale * math.sqrt(math.fsum((x / scale)**2 for x in centered) / (len(data) - 1))
    demand(sd > 0 and math.isfinite(sd), "nonconstant numerical variance underflow/overflow")
    return mean, sd


def welch_reference(left: Sequence[float], right: Sequence[float]) -> dict:
    """Left minus right, two sided; no arbitrary positive-variance threshold."""
    demand(len(left) > 0 and len(right) > 0, "empty group")
    ml, sl = moments(left)
    mr, sr = moments(right)
    result = {"n_left": len(left), "n_right": len(right),
              "mean_left": ml, "mean_right": mr, "sd_left": sl, "sd_right": sr,
              "effect": ml - mr, "se": None, "df": None, "t": None,
              "p_raw": None, "ci_low": None, "ci_high": None, "status": "insufficient_units"}
    if min(len(left), len(right)) < 2:
        return result
    a, b = sl / math.sqrt(len(left)), sr / math.sqrt(len(right))
    se = math.hypot(a, b)
    if se == 0:
        result["status"] = "zero_standard_error"
        return result
    demand(math.isfinite(result["effect"]), "nonfinite descriptive difference")
    # This scaled df avoids squaring tiny variances into numerical zero.
    df = 1.0 / ((a / se)**4 / (len(left) - 1) + (b / se)**4 / (len(right) - 1))
    t = result["effect"] / se
    if not all(math.isfinite(x) for x in (se, df, t)):
        result["status"] = "unavailable_numerical_inference"
        return result
    p = float(2 * special.stdtr(df, -abs(t)))
    halfwidth = float(special.stdtrit(df, 0.975)) * se
    low, high = result["effect"] - halfwidth, result["effect"] + halfwidth
    if p == 0 or not all(math.isfinite(x) for x in (p, halfwidth, low, high)):
        result["status"] = "unavailable_numerical_inference"
        return result
    result.update(status="ok", se=se, df=df, t=t, p_raw=p,
                  ci_low=low, ci_high=high)
    return result


def adjust_reference(raw_p: Sequence[float | None]) -> dict:
    """Independent full-family step-up BH/BY; missing is1 only here."""
    demand(len(raw_p) > 0, "empty correction family")
    p = [1.0 if value is None else number(value) for value in raw_p]
    demand(all(0 <= value <= 1 for value in p), "P outside [0,1]")
    order = sorted(range(len(p)), key=lambda i: (p[i], i))
    harmonic = math.fsum(1.0 / rank for rank in range(1, len(p) + 1))
    q = {}
    for method, factor in (("q_bh", 1.0), ("q_by", harmonic)):
        adjusted = [1.0] * len(p)
        previous = 1.0
        for position in range(len(order) - 1, -1, -1):
            index = order[position]
            candidate = p[index] * len(p) / (position + 1) * factor
            previous = min(previous, candidate, 1.0)
            adjusted[index] = previous
        q[method] = adjusted
    return {"p_for_correction": p, "family_size": len(p), **q}


def assert_family(rows: Sequence[dict], programs: Sequence[str], interface: str = CORE) -> int:
    """Exact ordered full grid, including unavailable rows; no Strict discoveries.

    Minimal independent wire schema: interface/program_id/left/right/status,
    p_raw/q_bh/q_by. Source failures must have missing P; valid rows must have P.
    """
    demand(interface in (CORE, STRICT), "unknown identity interface")
    demand(len(programs) > 0 and len(programs) == len(set(programs)), "duplicate/empty panel")
    expected = [(interface, p, a, b) for p in programs for a, b in CONTRASTS]
    observed = [(r["interface"], r["program_id"], r["left"], r["right"]) for r in rows]
    demand(observed == expected, "missing/extra/duplicate/reordered family row")
    for row in rows:
        demand(isinstance(row["status"], str) and bool(row["status"]), "missing status")
        demand((row["p_raw"] is not None) == (row["status"] == "ok"),
               "inferential status and raw P disagree")
    if interface == STRICT:
        demand(all(r["q_bh"] is None and r["q_by"] is None for r in rows),
               "Strict cannot supply discovery q values")
    else:
        adjusted = adjust_reference([r["p_raw"] for r in rows])
        for i, row in enumerate(rows):
            # Reported q is missing when the original test is unavailable.
            # The full numeric adjustment vector still contains its placeholder.
            for method in ("q_bh", "q_by"):
                expected_q = None if row["p_raw"] is None else adjusted[method][i]
                close(row[method], expected_q, "full-family " + method)
    return len(expected)


def stratified_indices(groups: Sequence[str], draws: int, seed: int) -> dict[str, list[list[int]]]:
    """Joint unit resampling in fixed PSC, ASC, AIH order, NumPy PCG64.

    Separate from author helpers. Replay should additionally check the frozen
    array hash: the seed alone does not pin later NumPy implementation changes.
    """
    import numpy as np
    demand(type(draws) is int and draws > 0, "invalid bootstrap draw count")
    demand(type(seed) is int and seed >= 0, "invalid bootstrap seed")
    demand(set(groups) == {"PSC", "ASC", "AIH"}, "bootstrap group universe changed")
    rng = np.random.Generator(np.random.PCG64(seed))
    output = {}
    for group in ("PSC", "ASC", "AIH"):
        original = [i for i, value in enumerate(groups) if value == group]
        local = rng.integers(0, len(original), size=(draws, len(original)))
        output[group] = [[original[int(k)] for k in row] for row in local]
    return output


def percentile_linear(values: Sequence[float], probability: float) -> float:
    demand(0 <= probability <= 1 and len(values) > 0, "invalid percentile")
    ordered = sorted(number(x) for x in values)
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    fraction = position - low
    return ordered[low] + fraction * (ordered[min(low + 1, len(ordered) - 1)] - ordered[low])


def bootstrap_reference(scores: Sequence[float], groups: Sequence[str], indices: dict,
                        left: str, right: str) -> dict:
    demand(len(scores) == len(groups), "sample axis mismatch")
    values = [number(x) for x in scores]
    demand(left != right and left in indices and right in indices, "invalid contrast")
    n_draws = len(indices[left])
    demand(n_draws > 0 and len(indices[right]) == n_draws, "draw count mismatch")
    for group, draws in indices.items():
        expected_n = groups.count(group)
        demand(expected_n > 0 and len(draws) == n_draws, "invalid stratum")
        for sample in draws:
            demand(len(sample) == expected_n, "bootstrap unit count changed")
            demand(all(type(i) is int and 0 <= i < len(groups) and groups[i] == group
                       for i in sample), "bootstrap left its stratum")
    effects = []
    for a, b in zip(indices[left], indices[right]):
        effects.append(math.fsum(values[i] / len(a) for i in a) -
                       math.fsum(values[i] / len(b) for i in b))
    return {"draws": effects, "ci_low": percentile_linear(effects, 0.025),
            "ci_high": percentile_linear(effects, 0.975)}


def omissions_reference(scores: Sequence[float], groups: Sequence[str],
                        unit_ids: Sequence[str], left: str, right: str) -> list[dict]:
    demand(len(scores) == len(groups) == len(unit_ids), "omission sample axis mismatch")
    demand(len(unit_ids) == len(set(unit_ids)), "unit IDs are repeated")
    demand(left != right and left in groups and right in groups, "invalid omission contrast")
    values = [number(x) for x in scores]
    rows = []
    for omitted, unit in enumerate(unit_ids):
        a = [values[i] for i, g in enumerate(groups) if g == left and i != omitted]
        b = [values[i] for i, g in enumerate(groups) if g == right and i != omitted]
        effect = None if not a or not b else math.fsum(a) / len(a) - math.fsum(b) / len(b)
        status = "not_in_contrast" if groups[omitted] not in (left, right) else "omitted"
        if effect is None:
            status = "insufficient_units"
        rows.append({"unit_id": unit, "omitted_group": groups[omitted], "n_left": len(a),
                     "n_right": len(b), "effect": effect, "status": status})
    return rows


def omission_summary_reference(rows: Sequence[dict], baseline_effect: float | None) -> dict:
    """Public stability counts/ranges use only units involved in the contrast."""
    relevant = [r for r in rows if r["status"] != "not_in_contrast"]
    defined = [r for r in relevant if r["effect"] is not None]
    effects = [number(r["effect"]) for r in defined]
    baseline = None if baseline_effect is None else number(baseline_effect)
    changes = [] if baseline is None else [effect - baseline for effect in effects]
    reversals = None if baseline is None or not effects else sum((effect < 0 < baseline) or
                                                                  (effect > 0 > baseline) for effect in effects)
    return {"in_contrast_rows":len(relevant),
            "stability_denominator":len(relevant),"defined_in_contrast_effect_rows":len(defined),
            "effect_min":min(effects) if effects else None,
            "effect_max":max(effects) if effects else None,
            "max_absolute_effect_change":max(map(abs,changes)) if changes else None,
            "strict_sign_reversal_rows":reversals}


def no_symlink(path: Path) -> None:
    demand(not any(p.is_symlink() for p in (path, *path.parents)), "symlink path prohibited")


def relative_file(root: Path, relative: str) -> Path:
    demand(isinstance(relative, str) and relative != "", "missing path")
    rel = Path(relative)
    demand(not rel.is_absolute() and ".." not in rel.parts, "path escape prohibited")
    path = root / rel
    no_symlink(path)
    demand(path.is_file(), "pinned file missing")
    demand(path.resolve().is_relative_to(root.resolve()), "resolved source escaped root")
    return path


def check_pin(root: Path, relative: str, pin: dict) -> None:
    path = relative_file(root, relative)
    demand(path.stat().st_size == pin["bytes"], f"size changed: {relative}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    demand(digest == pin["sha256"], f"hash changed: {relative}")


def assert_fresh_private(root: Path, output: Path, sources: Sequence[Path]) -> None:
    """Read-only guard oracle, not a directory creator or privacy attestation."""
    no_symlink(output)
    target = output.resolve()
    private_base = (root / "work/psc-liver-program-plan").resolve()
    demand(target.is_relative_to(private_base) and target != private_base,
           "private output must be a new child under the owned work area")
    demand(not output.exists(), "output must be fresh, including empty directories")
    for source in sources:
        no_symlink(source)
        resolved = source.resolve()
        demand(target != resolved and not resolved.is_relative_to(target)
               and not target.is_relative_to(resolved), "source/output path overlap")


def assert_approval(plan: dict, approval: dict, plan_bytes: bytes) -> None:
    """Independent fail-closed schema oracle; not cryptographic authorization."""
    demand(plan.get("schema_version") == 1, "unsupported plan schema")
    demand(plan.get("frozen") is True and approval.get("root_approved") is True,
           "root-approved separate freeze required")
    demand(plan.get("real_expression_authorized") is True and
           approval.get("real_expression_authorized") is True, "effects not authorized")
    demand(plan == json.loads(plan_bytes), "parsed plan does not match supplied bytes")
    demand(approval.get("plan_sha256") == hashlib.sha256(plan_bytes).hexdigest(),
           "approval does not bind exact plan bytes")
    demand(plan.get("normalization_feature_rows") == 41096, "wrong full normalization universe")
    demand(plan.get("primary_interface") == CORE and plan.get("sensitivity_interface") == STRICT,
           "identity policy changed")
    demand(plan.get("contrasts") == [list(x) for x in CONTRASTS], "contrast family changed")
    programs = plan.get("programs")
    demand(isinstance(programs, list) and bool(programs), "missing fixed program panel")
    ids = [p.get("program_id") for p in programs]
    demand(all(isinstance(x, str) and x for x in ids) and len(ids) == len(set(ids)), "invalid panel")
    for p in programs:
        if p.get("source_status") == "HOLD":
            demand(p.get("members") is None and p.get("effect_status") == "unavailable_source",
                   "unavailable source cannot become an empty or fallback score")
        else:
            demand(p.get("source_status") == "qualified" and p.get("source_sha256"),
                   "program source not qualified/pinned")


def source_only(root: Path, scope: Path) -> dict:
    """Check named metadata/code/report pins only; never open an RNA count file."""
    no_symlink(scope)
    record = json.loads(scope.read_text())
    pins = record["source_summary_pins"]
    allowed = {"data/derived/psc-liver-program-inputs.json", "reports/psc-liver-program-inputs.md",
               "config/psc-liver-external-program-sources.json",
               "data/derived/psc-liver-external-program-membership.csv"}
    demand(set(pins) == allowed, "source-only whitelist changed")
    for relative, pin in pins.items():
        check_pin(root, relative, pin)
    summary = json.loads(relative_file(root, "data/derived/psc-liver-program-inputs.json").read_text())
    groups = summary["inference_boundaries"]["sample_design_retained_from_closed_qualification"]
    demand({g: groups[g] for g in ("PSC", "ASC", "AIH")} == {"PSC":17, "ASC":17, "AIH":30},
           "sample-group source labels changed")
    for policy in (CORE, STRICT):
        pool = summary["candidate_matching_reference_pools"][policy]
        demand(pool["released_feature_axis_rows"] == 41096 and
               pool["counts_normalization_denominator_not_filtered_by_this_pool"] is True,
               "identity subset substituted for normalization universe")
        coverage = [r for r in summary["program_coverage"] if r["interface"] == policy]
        expected_ids = {"ets2_g1_dn", "ets2_g1_up", "ets2_g2_dn", "chr21_dn",
                        "inflammation", "interferon_gamma", "oxidative_stress",
                        "apoptotic_signaling", "ets2_g1_dn_without_comparators"}
        demand(len(coverage) == 9 and {r["program_id"] for r in coverage} == expected_ids,
               "eight original source programs plus nonoverlap not retained")
        for row in coverage:
            demand(row["effects_authorized"] is False, "identifier stage does not authorize effects")
            denominator = row["intended_source_members_after_ETS2_exclusion"]
            mapped = row["mapped_unique_ensembl_ids"]
            if row["program_id"] == "ets2_g1_dn_without_comparators":
                parent = next(r for r in coverage if r["program_id"] == "ets2_g1_dn")
                demand(denominator == parent["intended_source_members_after_ETS2_exclusion"],
                       "nonoverlap denominator changed")
                expected_gate = parent["mapping_gate_pass"] and mapped >= 10
                demand(mapped + row["overlap_removed"] == parent["mapped_unique_ensembl_ids"],
                       "nonoverlap member accounting changed")
            else:
                expected_gate = mapping_gate(denominator, mapped)
            demand(expected_gate == row["mapping_gate_pass"], "source-only annotation gate mismatch")
    with relative_file(root, "data/derived/psc-liver-external-program-membership.csv").open(newline="") as h:
        external_rows = list(csv.DictReader(h))
    demand(len(external_rows) == 321 and len({r["source_value"] for r in external_rows}) == 321,
           "ECM official membership changed")
    demand(all(r["program_id"] == "reactome_ecm_organization" and
               r["published_direction"] == "UNSIGNED" and r["membership_weight"] == "1"
               for r in external_rows), "external source substituted or invented direction")
    demand(sum(int(r["original_source_id_count"]) for r in external_rows) == 340,
           "official source-ID multiplicity changed")
    external = json.loads(relative_file(root, "config/psc-liver-external-program-sources.json").read_text())
    hold = external["unavailable_requested_programs"]
    demand(len(hold) == 1 and hold[0]["program_id"] == "epithelial_il17a_alone_response"
           and hold[0]["member_count"] is None and hold[0]["published_directions"] is None
           and hold[0]["scoring_authorized"] is False
           and hold[0]["membership_is_not_empty_or_zero"] is True,
           "IL17 source HOLD was lost or converted to a zero program")
    return {"status":"source_only_pass", "source_metadata_pins_checked":len(pins),
            "original_and_nonoverlap_gates_checked":18, "ECM_members_checked":321,
            "no_real_count_or_score_rows_opened":True, "real_effect_verification":False,
            "effects_authorized":False,
            "limits":["reads already qualified metadata, not independent raw-source reconstruction",
                      "IL17A response stays HOLD outside the ten-program analysis; root execution freeze is still required",
                      "not an independent-person or causal-inference verification"]}


def ecm_mapping_reference(feature_rows: Sequence[dict], source_rows: Sequence[dict],
                          expected_members: int = 321) -> dict:
    """Independent exact-symbol reuse of the already qualified full feature audit.

    No count values, aliases, published Entrez or many-to-one source Ensembl IDs
    are used as a bridge. This does not redo the closed global-reference audit.
    """
    demand(len(source_rows) == expected_members, "official ECM denominator changed")
    labels = [r["label_as_submitted"] for r in feature_rows]
    demand(len(labels) == len(set(labels)), "feature labels are not unique")
    demand([r["source_row_index_0based"] for r in feature_rows] == list(range(len(feature_rows))),
           "full released feature-index order changed")
    feature_by_label = dict(zip(labels, feature_rows))
    symbols = [r["source_value"] for r in source_rows]
    demand(len(symbols) == len(set(symbols)) and all(s and s == s.strip() for s in symbols),
           "repeated/missing/nonliteral official symbols")
    demand([int(r["source_member_index"]) for r in source_rows] == list(range(1, expected_members + 1)),
           "official source membership order changed")
    demand(all(r["program_id"] == "reactome_ecm_organization" and
               r["source_identifier_namespace"] == "MSigDB_geneSymbols" and
               r["membership_weight"] == "1" and r["published_direction"] == "UNSIGNED"
               for r in source_rows), "ECM source namespace/weight/orientation changed")
    result = {}
    for tier in (CORE, STRICT):
        indices, audit = [], []
        for symbol in symbols:
            feature = feature_by_label.get(symbol)
            status = "no_exact_submitted_symbol" if feature is None else feature[tier]
            index = feature["source_row_index_0based"] if status == "eligible" else None
            if index is not None:
                indices.append(index)
            audit.append({"source_value":symbol,"mapping_status":status,"feature_index":index})
        demand(len(indices) == len(set(indices)), "multiple members collapse onto a feature")
        result[tier] = {"original_source_members":expected_members,"mapped_members":len(indices),
                        "indices":indices,"coverage_available":mapping_gate(expected_members,len(indices)),
                        "audit":audit}
    demand(set(result[STRICT]["indices"]).issubset(result[CORE]["indices"]),
           "Strict must be an identity subset of Core")
    return result


def verify_numeric_bundle(bundle: dict) -> dict:
    """Compare a complete synthetic run export, with no file access.

    This in-memory wire format is an independent audit boundary, not the author's
    native output schema. A separately pinned adapter must preserve all fields.
    Real liver rows are not accepted or authorized by this prospective API.
    """
    demand(bundle.get("schema_version") == 1 and bundle.get("synthetic_only") is True,
           "only explicitly synthetic method fixtures are authorized")
    counts, groups, unit_ids = bundle["counts"], bundle["groups"], bundle["unit_ids"]
    demand(len(groups) == len(unit_ids) == len(counts[0]), "sample axis mismatch")
    demand(len(unit_ids) == len(set(unit_ids)), "duplicate unit IDs")
    demand(set(groups) == {"PSC", "ASC", "AIH"}, "group universe changed")
    normalized = logcpm_reference(counts)
    observed = bundle["observed"]
    checks = 0

    def match(actual, expected, label):
        nonlocal checks
        if isinstance(expected, dict):
            demand(set(actual) == set(expected), label + ": keys changed")
            for key in expected:
                match(actual[key], expected[key], label + "." + key)
        elif isinstance(expected, (list, tuple)):
            demand(len(actual) == len(expected), label + ": size changed")
            for i, value in enumerate(expected):
                match(actual[i], value, label + f"[{i}]")
        else:
            close(actual, expected, label)
            checks += 1

    match(observed["normalization"], normalized, "normalization")
    program_order = bundle["program_order"]
    demand(len(program_order) == len(set(program_order)), "duplicate fixed endpoint")
    expected_keys = [(tier, p) for tier in (CORE, STRICT) for p in program_order]
    specs = bundle["programs"]
    demand([(s["interface"], s["program_id"]) for s in specs] == expected_keys,
           "complete Core/Strict fixed program grid required")
    demand(len(observed["programs"]) == len(specs), "score output grid changed")
    indices = stratified_indices(groups, bundle["replicates"], bundle["seed"])
    match(observed["bootstrap_indices"], indices, "bootstrap_indices")
    inference = []
    for spec, actual in zip(specs, observed["programs"]):
        match([actual["interface"], actual["program_id"]],
              [spec["interface"], spec["program_id"]], "score_key")
        if spec["coverage_available"]:
            scores = score_reference(normalized["logcpm"], spec["positive"], spec.get("negative"))
        else:
            scores = None
        match(actual["scores"], scores, "scores")
        demand(len(actual["contrasts"]) == len(CONTRASTS), "contrast grid changed")
        for (left, right), row in zip(CONTRASTS, actual["contrasts"]):
            match([row["left"], row["right"]], [left, right], "contrast_order")
            expected_welch = (welch_reference([scores[i] for i,g in enumerate(groups) if g==left],
                                             [scores[i] for i,g in enumerate(groups) if g==right])
                              if scores is not None else None)
            match(row["welch"], expected_welch, "Welch")
            expected_boot = (bootstrap_reference(scores, groups, indices, left, right)
                             if scores is not None else None)
            match(row["bootstrap"], expected_boot, "bootstrap")
            expected_omissions = (omissions_reference(scores, groups, unit_ids, left, right)
                                  if scores is not None else None)
            match(row["omissions"], expected_omissions, "omissions")
            if expected_omissions is None:
                n_in = groups.count(left) + groups.count(right)
                expected_summary = {"in_contrast_rows":n_in,
                    "stability_denominator":n_in,"defined_in_contrast_effect_rows":0,
                    "effect_min":None,"effect_max":None,"max_absolute_effect_change":None,
                    "strict_sign_reversal_rows":None}
            else:
                expected_summary = omission_summary_reference(expected_omissions, expected_welch["effect"])
            match(row["omission_summary"], expected_summary, "public_omission_summary")
            inference.append({"interface":spec["interface"], "program_id":spec["program_id"],
                              "left":left,"right":right,
                              "status":"unavailable_score" if expected_welch is None else expected_welch["status"],
                              "p_raw":None if expected_welch is None else expected_welch["p_raw"]})
    adjustment = adjust_reference([r["p_raw"] for r in inference if r["interface"] == CORE])
    match(observed["full_core_correction"], adjustment, "full_Core_correction")
    for i, expected in enumerate(inference):
        actual = observed["inference"][i]
        q_missing = expected["interface"] == STRICT or expected["p_raw"] is None
        expected.update(q_bh=None if q_missing else adjustment["q_bh"][i],
                        q_by=None if q_missing else adjustment["q_by"][i])
        match(actual, expected, "full_inference_row")
    demand(len(observed["inference"]) == len(inference), "extra inference rows")
    for tier in (CORE, STRICT):
        assert_family([r for r in observed["inference"] if r["interface"] == tier], program_order, tier)
    return {"status":"pass_synthetic_numerical_bundle", "scalar_checks":checks,
            "normalization_rows":len(counts), "source_patient_sample_units":len(groups),
            "program_tier_rows":len(specs), "full_Core_family_size":adjustment["family_size"],
            "full_status_rows":len(inference), "real_expression_reviewed":False}


def builtin_synthetic() -> dict:
    counts = [[1, 0, 4, 2, 5, 3], [9, 10, 6, 8, 5, 7], [90, 90, 90, 90, 90, 90]]
    normalized = logcpm_reference(counts)
    scores = score_reference(normalized["logcpm"], [0, 1])
    groups = ["PSC", "PSC", "ASC", "ASC", "AIH", "AIH"]
    indices = stratified_indices(groups, 100, 20260916)
    contrasts = [welch_reference([scores[i] for i,g in enumerate(groups) if g==a],
                                 [scores[i] for i,g in enumerate(groups) if g==b])
                 for a,b in CONTRASTS]
    adjusted = adjust_reference([r["p_raw"] for r in contrasts] + [None]*3)
    bootstrap = [bootstrap_reference(scores, groups, indices, a,b) for a,b in CONTRASTS]
    omissions = [omissions_reference(scores, groups, [f"S{i}" for i in range(6)], a,b)
                 for a,b in CONTRASTS]
    return {"status":"builtin_synthetic_only", "full_universe_totals":normalized["library_totals"],
            "score_count":len(scores), "contrast_count":len(contrasts),
            "full_family_size_including_unavailable":adjusted["family_size"],
            "bootstrap_contrast_count":len(bootstrap), "omission_rows":sum(map(len,omissions)),
            "no_real_count_or_score_rows_opened":True}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--source-only", action="store_true")
    mode.add_argument("--synthetic-only", action="store_true")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    if args.source_only:
        output = source_only(root, root/"work/psc-liver-program-independent-review/scope.json")
    else:
        output = builtin_synthetic()
    print(json.dumps(output, sort_keys=True, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError) as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
