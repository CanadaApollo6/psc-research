#!/usr/bin/env python3
"""Bounded pinned-file audit for ONE BACH2 endpoint; default is metadata only.

Real streams and result JSON require a separate externally hashed root seal,
explicit acknowledgment, author authorization, and complete last-written receipt.
This adapter does not authorize execution. Its tests use synthetic sources only.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from collections import Counter
from decimal import Context, Decimal, InvalidOperation, localcontext
from fractions import Fraction
from functools import lru_cache
from pathlib import Path, PurePosixPath
import re
import stat
import sys

ROOT = Path(__file__).resolve().parents[1]
OWN = "work/psc-bach2-naive-rna-execution-review"
AUTHOR_BASE = "work/psc-bach2-naive-rna-plan"
VERIFIER = "scripts/verify_psc_bach2_naive_rna.py"
VERIFIER_SHA = "15fcd8d2926a0bb83c076aa707c5fa50e175763f905039b90ed2e2802e73a18d"
ADAPTER = "scripts/verify_psc_bach2_naive_rna_execution.py"
TESTS = "tests/test_verify_psc_bach2_naive_rna_execution.py"
DESIGN = "reports/psc-bach2-naive-rna-execution-verification-design.md"
AUTHOR_CODE = "scripts/analyze_psc_bach2_naive_rna.py"
AUTHOR_TESTS = "tests/test_psc_bach2_naive_rna.py"
CANDIDATE = "work/psc-bach2-followup-review/candidate-specification.json"
CANDIDATE_SHA = "b54c4fc47a93afc48b6463410902e467e459049cef759dc3ba50d8fb6149d627"
CANDIDATE_CANONICAL_SHA = "57d59202ab934fda7f7713f700a556b015acf49b6861f8a391b670081f54a6cb"
SAMPLES = tuple(f"sample{i:02d}" for i in range(1, 9))
GROUPS = ("SNP", "noSNP")
TARGET = ("ENSG00000112182", "BACH2", "Gene Expression")
TARGET_ROW = 12314  # Full source-feature position, ONE-based; not genomic.
HEADER = "%%MatrixMarket matrix coordinate integer general"
EXPORTER = {"software_version": "cellranger-6.1.1", "format_version": 2}
MAX_NNZ = 50_000_000
MAX_R = MAX_NNZ * 10**100
PRECISION = 4000
PROTOCOL = {
    "schema_version": 1, "reference_decimal_digits": PRECISION,
    "source_integer_digits": 100, "maximum_matrix_coordinates": MAX_NNZ,
    "logarithms": "direct Decimal.ln of exact rational Q/products; no author log series",
    "points_relative": "1e-10", "moments_SD_SE_df_relative": "1e-10",
    "CI_standardized_critical_relative": "1e-9", "P_relative": "1e-10",
    "zero_and_sign": "exact; no absolute epsilon, no zero replacement",
    "P_reference": "positive Decimal incomplete-beta series, 100 digits, float log-gamma normalization",
    "CI_reference": "scipy.special.betainc inverted with scipy.optimize.brentq",
    "native_availability_only": "scipy.special.stdtr/stdtrit; no scipy.stats import",
    "scope": "one scientific primary, independent verification, mechanical replay; no extra endpoint",
}
PAYLOADS = {"measurements.private.json", "public-aggregate.json", "authorization.used.json"}
PLAN_FILES = {"plan.json", "selection.private.json", "public-schema.json"}
AUTH_STATUS = "ROOT_AUTHORIZED_ONE_FIXED_BACH2_NAIVE_RNA_EXTRACTION"
ACK = "ROOT_AUTHORIZED_ONE_FIXED_BACH2_REAL_OUTPUT_AUDIT"

PUBLIC_SCHEMA = {'schema_version': 'integer:1',
 'status': 'literal:fixed_endpoint_complete',
 'endpoint': 'literal:ENSG00000112182/BACH2; fixed pooled C0+C1 naive CD4; eight source donors',
 'quantity': 'literal:mean donor log2(1+1000000*released_BACH2_RNA/all36601_released_RNA_units)',
 'groups': {'SNP': {'n': 'integer:4', 'mean': 'decimal-string', 'sample_sd': 'decimal-string'},
            'noSNP': {'n': 'integer:4', 'mean': 'decimal-string', 'sample_sd': 'decimal-string'}},
 'contrast': {'SNP_minus_noSNP': 'decimal-string',
              'SE': 'decimal-string',
              'df': 'decimal-string-or-null',
              'CI95': 'two-decimal-strings-or-null',
              'P_two_sided': 'decimal-string-or-null',
              'inference_status': 'enum:unavailable_zero_SE|unavailable_nonfinite_inference|CI_available_P_unavailable_numeric_tail|available_nominal_Welch'},
 'leave_one_donor_out': {'count': 'integer:8',
                         'minimum': 'decimal-string',
                         'maximum': 'decimal-string',
                         'positive_count': 'integer',
                         'negative_count': 'integer',
                         'zero_count': 'integer',
                         'strict_sign_reversal_count': 'integer',
                         'full_point_sign': 'integer:-1/0/1'},
 'QC': {'donors': 'integer:8',
        'selected_cells': 'integer:16659',
        'retained_cells': 'integer:55460',
        'RNA_rows': 'integer:36601',
        'zero_BACH2_donors_retained': 'integer',
        'BACH2_positive_selected_cells': 'integer',
        'source_and_membership_gates': 'literal:passed'},
 'provenance': {'plan_sha256': 'sha256',
                'selection_sha256': 'sha256',
                'source_manifest_sha256': 'sha256',
                'implementation_sha256': 'sha256',
                'runtime_sha256': 'sha256',
                'public_schema_sha256': 'sha256',
                'authorization_sha256': 'sha256',
                'candidate_sha256': 'sha256'}}

class AuditError(ValueError):
    """Fail the whole audit, without disclosing source values in CLI errors."""


def require(condition, code):
    if not condition:
        raise AuditError(code)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def digest(value):
    return sha(canonical(value))


def keys(value, expected):
    require(type(value) is dict and set(value) == set(expected), "schema_keys")


def hash_string(value):
    require(type(value) is str and re.fullmatch("[a-f0-9]{64}", value) is not None, "sha256_required")
    return value


def path_at(root, relative, *, regular=False):
    """Canonical repo-relative paths, no link component or hard-linked input."""
    require(type(relative) is str and "\\" not in relative, "unsafe_path")
    pure = PurePosixPath(relative)
    require(relative == pure.as_posix() and pure.parts and not pure.is_absolute()
            and all(p not in (".", "..") for p in pure.parts), "unsafe_path")
    root = Path(root).absolute()
    require(not any(p.is_symlink() for p in (root, *root.parents)), "root_link")
    path = root
    for part in pure.parts:
        path /= part
        require(not path.is_symlink(), "path_link")
    if regular:
        st = path.stat()
        require(stat.S_ISREG(st.st_mode) and st.st_nlink == 1, "regular_nonalias_input_required")
    return path


def read_bytes(root, relative, expected=None, maximum=32_000_000):
    path = path_at(root, relative, regular=True)
    require(path.stat().st_size <= maximum, "bounded_metadata_file")
    data = path.read_bytes()
    require(len(data) <= maximum, "bounded_metadata_file")
    if expected is not None:
        require(sha(data) == hash_string(expected), "file_pin_changed")
    return data


def decode_json(data):
    def object_pairs(pairs):
        require(len({key for key, _ in pairs}) == len(pairs), "duplicate_JSON_key")
        return dict(pairs)
    value = json.loads(data, object_pairs_hook=object_pairs,
                       parse_constant=lambda _: (_ for _ in ()).throw(AuditError("nonfinite_JSON")))
    require(data == canonical(value), "noncanonical_JSON")
    return value


def json_pin(root, relative, expected=None):
    return decode_json(read_bytes(root, relative, expected))


def inventory(root, relative, names):
    directory = path_at(root, relative)
    require(directory.is_dir() and {p.name for p in directory.iterdir()} == set(names), "complete_inventory_required")
    identities = set()
    for name in names:
        p = path_at(root, relative + "/" + name, regular=True)
        identity = (p.stat().st_dev, p.stat().st_ino)
        require(identity not in identities, "aliased_inventory")
        identities.add(identity)
    return directory


def ignored_policy(root):
    patterns = [line.strip() for line in read_bytes(root, ".gitignore").decode().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
    require("work/" in patterns and all(not p.startswith("!") or p == "!.env.example" for p in patterns),
            "ignored_private_policy")


def fresh_output(root, relative, create=False):
    require(type(relative) is str and re.fullmatch(re.escape(OWN) + r"/[A-Za-z0-9][A-Za-z0-9_-]{0,79}", relative),
            "owned_fresh_output_required")
    ignored_policy(root)
    p = path_at(root, relative)
    require(not p.exists(), "output_reuse_refused")
    if create:
        p.parent.mkdir(parents=True, exist_ok=True)
        path_at(root, relative)
        p.mkdir(mode=0o700)
    return p


def load_verifier(root=ROOT):
    # Authentication precedes import; no production implementation is imported.
    read_bytes(root, VERIFIER, VERIFIER_SHA)
    spec = importlib.util.spec_from_file_location("bach2_pinned_independent_metadata", path_at(root, VERIFIER))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_plan_structure(plan):
    keys(plan, {"RNA_rows_1based", "candidate_sha256", "closed_limits", "exporter_comment", "fixed_candidate",
                "implementation", "matrix_header", "nnz_gate", "numerics", "permitted_future_run_root", "public_schema",
                "public_schema_sha256", "reference_identity", "runtime", "schema_version", "selection_sha256",
                "source_manifest", "source_manifest_sha256", "status", "future_run_inventory", "completion_status"})
    require(type(plan["schema_version"]) is int and plan["schema_version"] == 1 and
            plan["status"] == "METHOD_PREPARATION_ONLY_NOT_EXECUTION_AUTHORIZATION", "preparation_plan_required")
    require(plan["candidate_sha256"] == CANDIDATE_SHA and digest(plan["fixed_candidate"]) == CANDIDATE_CANONICAL_SHA,
            "fixed_candidate")
    require(plan["public_schema"] == PUBLIC_SCHEMA and plan["public_schema_sha256"] == digest(PUBLIC_SCHEMA), "public_schema")
    require(plan["matrix_header"] == HEADER and plan["exporter_comment"] == EXPORTER, "source_format")
    require(type(plan["RNA_rows_1based"]) is list and len(plan["RNA_rows_1based"]) == 36601
            and all(type(row) is int for row in plan["RNA_rows_1based"]), "complete_integer_RNA_axis")
    require(plan["permitted_future_run_root"] == AUTHOR_BASE + "/runs/", "author_run_root")
    require(plan["future_run_inventory"] == sorted(PAYLOADS | {"completion.json"}) and
            plan["completion_status"] == "completed_fixed_endpoint_bound_payload_inventory", "completion_contract")
    require(plan["numerics"] == {
        "Decimal_precision": 320, "Student_t": "native SciPy actual Welch df; P tail underflow is unavailable, not zero",
        "maximum_coordinates_per_matrix": MAX_NNZ, "maximum_value_digits": 100,
        "numeric_JSON": "finite Decimal strings; working precision is not biological precision",
        "point_and_LOO": "algebraically identical exact rational products of 1+1e6B/R; stable logarithm preserves true zero/sign",
        "variance": "exact Fraction moments of fixed-precision Decimal Y; no constancy tolerance"}, "native_numerics_changed")
    require(type(plan["implementation"]) is list and [p["path"] for p in plan["implementation"]] == [AUTHOR_CODE, AUTHOR_TESTS],
            "implementation_inventory")
    manifest = plan["source_manifest"]
    keys(manifest, {"annotation_and_receipt_pins", "matrix_expected_pins_NOT_read", "archive_expected_pins_NOT_read"})
    require(digest(manifest) == plan["source_manifest_sha256"], "source_manifest_hash")
    for p in manifest["annotation_and_receipt_pins"] + plan["implementation"]:
        keys(p, {"path", "sha256", "bytes"})
        hash_string(p["sha256"])
        require(type(p["bytes"]) is int and p["bytes"] > 0, "positive_file_bytes")
    for collection in ("matrix_expected_pins_NOT_read", "archive_expected_pins_NOT_read"):
        pins = manifest[collection]
        require(type(pins) is list and len(pins) == 8 and {p["sample"] for p in pins} == set(SAMPLES), "eight_source_pins")
        for p in pins:
            keys(p, {"sample", "path", "bytes", "sha256", "rows", "columns"} if collection.startswith("matrix")
                 else {"sample", "path", "bytes", "sha256", "source_URL", "retrieved_at_UTC"})
            hash_string(p["sha256"])
            require(type(p["bytes"]) is int and p["bytes"] > 0, "source_bytes")


def annotation_review(plan, selection, verifier, root=ROOT):
    """Independent source parsing; never opens archives, matrices or margins."""
    v = verifier
    expected = set(v.FIXED) | set(v.RECEIPT_PINS) | {".gitignore", v.METADATA_PIN[0], v.SDRF_PIN[0]}
    for sample in SAMPLES:
        expected |= {f"work/psc-bach2-cohort-qualification/{sample}/members/{kind}.txt" for kind in ("features", "barcodes")}
    declared = plan["source_manifest"]["annotation_and_receipt_pins"]
    require(len(declared) == len(expected) and {p["path"] for p in declared} == expected, "metadata_read_allowlist")
    # Validate every path before even the independent metadata module opens it.
    for relative in expected:
        path_at(root, relative, regular=True)
    for relative, expected_sha in {**v.FIXED, **v.RECEIPT_PINS}.items():
        read_bytes(root, relative, expected_sha)
    for relative, expected_sha, size in (v.METADATA_PIN, v.SDRF_PIN):
        require(len(read_bytes(root, relative, expected_sha)) == size, "source_metadata_size")
    for pin in declared + plan["implementation"]:
        require(len(read_bytes(root, pin["path"], pin["sha256"])) == pin["bytes"], "source_pin_bytes")
    members, archives, _ = v.receipt_contract(root)
    mapping = v.source_mapping(path_at(root, v.SDRF_PIN[0]))
    axes, barcodes = [], {}
    for sample in SAMPLES:
        receipt = members[sample]
        for kind in ("features", "barcodes"):
            p = receipt[kind]
            require(len(read_bytes(root, p["decoded_path"], p["decoded_sha256"])) == p["decoded_bytes"], "axis_pin")
        text = read_bytes(root, receipt["features"]["decoded_path"]).decode()
        axis = [tuple(row) for row in csv.reader(io.StringIO(text), delimiter="\t")]
        require(len(axis) == 36639 and all(len(r) == 3 and all(r) for r in axis), "full_feature_axis")
        require(len({r[0] for r in axis}) == len(axis), "duplicate_feature")
        require(Counter(r[2] for r in axis) == {"Gene Expression": 36601, "Antibody Capture": 38}, "modality_axis")
        require([i for i, row in enumerate(axis, 1) if row == TARGET] == [TARGET_ROW]
                and sum(r[0] == TARGET[0] or r[1] == TARGET[1] for r in axis) == 1, "unique_target_feature_row")
        axes.append(axis)
        barcodes[sample] = read_bytes(root, receipt["barcodes"]["decoded_path"]).decode().splitlines()
        require(len(set(barcodes[sample])) == len(barcodes[sample]) and
                all(re.fullmatch(r"[ACGT]+-[0-9]+", b) for b in barcodes[sample]), "barcode_axis")
    require(all(axis == axes[0] for axis in axes), "ordered_feature_axes")
    require(plan["RNA_rows_1based"] == [i for i, r in enumerate(axes[0], 1) if r[2] == "Gene Expression"], "complete_RNA_axis")
    with path_at(root, v.METADATA_PIN[0]).open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        require(all(r["barcode"].startswith(r["sample_name"] + "_") for r in rows), "exact_sample_prefix")
    retained, selected, counts, states = v.metadata_selection(path_at(root, v.METADATA_PIN[0]), mapping, barcodes)
    require(len(retained) == 55460 and len(selected) == 16659 and
            states == {k: x[2] for k, x in v.STATES.items()}, "complete_fixed_selection")
    require(set(r["sample"] for r in selected) == set(SAMPLES) and
            all({r["paper_clusters"] for r in selected if r["sample"] == s} == set(v.STATES) for s in SAMPLES), "all_eight_both_states")
    require(sum(map(len, barcodes.values())) == 59315, "complete_source_cells")
    v.reference_identity(path_at(root, v.REFERENCE_PATH))
    keys(selection, {"retained_keys", "donors"})
    require(type(selection["donors"]) is list and len(selection["donors"]) == 8
            and {d["sample"] for d in selection["donors"]} == set(SAMPLES), "selection_eight_donors")
    for d in selection["donors"]:
        keys(d, {"sample", "donor", "group", "retained_count", "selected", "state_counts"})
        keys(d["state_counts"], {"C0", "C1"})
        for cell in d["selected"]:
            keys(cell, {"barcode", "column_1based", "state"})
            require(type(cell["column_1based"]) is int, "source_column_integer")
    private = {"mapping": mapping, "retained_cell_counts": dict(counts),
               "source_cell_counts": {s: len(b) for s, b in barcodes.items()}, "selected": selected}
    return private


def metadata_only(plan_directory, plan_sha, root=ROOT):
    hash_string(plan_sha)
    require(re.fullmatch(re.escape(AUTHOR_BASE) + r"/[A-Za-z0-9][A-Za-z0-9_-]{0,79}", plan_directory)
            and plan_directory.rsplit("/", 1)[-1] not in ("runs", "independent-review"), "plan_directory")
    inventory(root, plan_directory, PLAN_FILES)
    plan = json_pin(root, plan_directory + "/plan.json", plan_sha)
    check_plan_structure(plan)
    selection = json_pin(root, plan_directory + "/selection.private.json", plan["selection_sha256"])
    require(read_bytes(root, plan_directory + "/public-schema.json", plan["public_schema_sha256"]) == canonical(PUBLIC_SCHEMA), "schema_payload")
    ignored_policy(root)
    # Verify that subsequent calculations import the installed pinned packages,
    # rather than a project-local module with the same import name.
    import importlib.metadata
    import numpy
    import scipy
    for name, module in (("numpy", numpy), ("scipy", scipy)):
        dist = importlib.metadata.distribution(name)
        require(Path(module.__file__).resolve() == Path(dist.locate_file(name+"/__init__.py")).resolve()
                and module.__version__ == dist.version, "runtime_import_shadowed")
    v = load_verifier(root)
    private = annotation_review(plan, selection, v, root)
    # This authenticated independent primitive verifies source pin/selection/runtime
    # identities again, including installed wheel bodies. It never reads matrices.
    summary = v.audit_compiled_plan(plan_directory, {}, private, root)
    return plan, selection, {"status": "METADATA_ONLY_AUTHENTICATED_NO_REAL_ACCESS", "plan_sha256": plan_sha,
        "selection_sha256": plan["selection_sha256"], "protocol_sha256": digest(PROTOCOL),
        "verifier_sha256": VERIFIER_SHA, "real_matrices_opened": 0, "numeric_payloads_parsed": 0,
        "expression_authorized": False, "compiled_source_runtime_verified": summary["status"]}


def integer_token(token):
    # Independent lexical exponent/integrality parser; no production Decimal parser.
    require(type(token) is str and len(token) <= 100, "value_token_bound")
    match = re.fullmatch(r"\+?([0-9]+)(?:\.([0-9]+))?(?:[eE]([+-]?[0-9]+))?", token)
    require(match is not None, "positive_integer_syntax")
    whole, fractional, exponent = match.groups()
    fractional = fractional or ""
    exponent = exponent or "0"
    require(len(exponent.lstrip("+-").lstrip("0")) <= 4, "value_exponent_bound")
    power = int(exponent) - len(fractional)
    coefficient = (whole + fractional).lstrip("0")
    require(coefficient, "explicit_zero_refused")
    if power < 0:
        require(-power <= len(coefficient) and coefficient.endswith("0" * -power), "fractional_source_unit")
        coefficient = coefficient[:power]
        power = 0
    require(0 < len(coefficient) + power <= 100, "source_integer_digits")
    return int(coefficient) * 10**power


def stream_matrix(handle, pin, columns, rna_rows):
    """Full decoded MatrixMarket pass, only B/R/detection accumulated.

    All coordinates, including unselected/ADT entries, undergo source-integrity
    checks. No ADT or unselected-gene endpoint or marginal total is accumulated.
    """
    require(type(columns) is list and columns and len(set(columns)) == len(columns)
            and all(type(c) is int and 1 <= c <= pin["columns"] for c in columns), "selected_columns")
    require(type(rna_rows) is list and len(rna_rows) == 36601 and len(set(rna_rows)) == 36601
            and all(type(r) is int and 1 <= r <= 36639 for r in rna_rows) and TARGET_ROW in rna_rows, "RNA_rows")
    require(type(pin["rows"]) is int and pin["rows"] == 36639 and type(pin["columns"]) is int and pin["columns"] > 0,
            "matrix_axis")
    selected, rna = set(columns), set(rna_rows)
    total_bytes, checksum = 0, hashlib.sha256()
    def get_line():
        nonlocal total_bytes
        raw = handle.readline(513)
        require(len(raw) <= 512, "matrix_line_limit")
        total_bytes += len(raw)
        require(total_bytes <= pin["bytes"], "matrix_byte_limit")
        checksum.update(raw)
        try:
            return raw.decode("ascii").rstrip("\r\n") if raw else None
        except UnicodeDecodeError as exc:
            raise AuditError("matrix_not_ASCII") from exc
    require(get_line() == HEADER, "matrix_header")
    metadata = get_line()
    require(metadata is not None and metadata.startswith("%metadata_json: "), "exporter_comment")
    def unique(pairs):
        require(len(dict(pairs)) == len(pairs), "exporter_duplicate_key")
        return dict(pairs)
    require(json.loads(metadata[16:], object_pairs_hook=unique) == EXPORTER, "exporter_identity")
    dimension = get_line()
    require(dimension is not None and re.fullmatch(r"[0-9]+ [0-9]+ [0-9]+", dimension), "dimension_syntax")
    nr, nc, declared = map(int, dimension.split())
    require((nr, nc) == (pin["rows"], pin["columns"]) and 0 <= declared <= min(MAX_NNZ, nr*nc), "matrix_dimensions")
    b = r = positive = observed = 0
    previous_index = -1
    while (line := get_line()) is not None:
        fields = line.split()
        require(len(fields) == 3 and all(re.fullmatch(r"[0-9]{1,9}", s) for s in fields[:2]), "coordinate_syntax")
        row, col = int(fields[0]), int(fields[1])
        require(1 <= row <= nr and 1 <= col <= nc, "coordinate_bounds")
        index = (col-1)*nr + row-1
        require(index > previous_index, "duplicate_or_unsorted_coordinate")
        previous_index = index
        value = integer_token(fields[2])
        observed += 1
        require(observed <= declared, "too_many_coordinates")
        if col in selected and row in rna:
            r += value
            if row == TARGET_ROW:
                b += value
                positive += 1
    require(observed == declared and total_bytes == pin["bytes"] and checksum.hexdigest() == pin["sha256"], "matrix_integrity")
    require(0 <= b <= r < MAX_R and r > 0, "whole_endpoint_bad_R")
    return {"B": b, "R": r, "B_positive_cells": positive, "selected_cells": len(selected)}


def dec(q):
    if isinstance(q, Fraction):
        return Decimal(q.numerator) / Decimal(q.denominator)
    return Decimal(q)


def sign(value):
    return int(value > 0) - int(value < 0)


@lru_cache(maxsize=256)
def direct_log2(q):
    require(isinstance(q, Fraction) and q > 0, "positive_rational_log")
    if q == 1:
        return Decimal(0)
    with localcontext(Context(prec=PRECISION)):
        result = dec(q).ln() / Decimal(2).ln()
        require(result.is_finite() and sign(result) == sign(q-1), "reference_log_resolution")
        return result


def product_point(arguments):
    require(set(arguments) == set(GROUPS) and all(len(arguments[g]) in (3, 4) for g in GROUPS)
            and sum(map(len, arguments.values())) in (7, 8), "fixed_point_sizes")
    products = {g: math.prod(arguments[g], start=Fraction(1)) for g in GROUPS}
    n1, n0 = len(arguments["SNP"]), len(arguments["noSNP"])
    # For 4/4 reduce the identity before logging. 3/4 uses powers 4 and 3.
    q, divisor = ((products["SNP"] / products["noSNP"], 4) if n1 == n0 == 4 else
                  (products["SNP"]**n0 / products["noSNP"]**n1, n1*n0))
    with localcontext(Context(prec=PRECISION)):
        result = direct_log2(q) / divisor
    require(sign(result) == sign(q-1), "exact_point_sign")
    return result


def beta_tail(statistic, df):
    """Independent nonzero P, even below float64 subnormal range.

    Regularized incomplete beta via its convergent hypergeometric series on
    x<=1/2. log-gamma normalization is float64 (disclosed, 1e-10 P criterion).
    Never substitute a float zero for the positive reference probability.
    """
    require(statistic >= 0 and 3 <= df <= 6, "Student_domain")
    if statistic == 0:
        return Decimal(1)
    with localcontext(Context(prec=100)):
        a, b = df/2, Decimal("0.5")
        x = df / (df + statistic*statistic)
        complement = x > Decimal("0.5")
        if complement:
            # Form 1-x without cancellation for tiny nonzero statistic.
            x, a, b = statistic*statistic / (df + statistic*statistic), b, a
        require(0 < x <= Decimal("0.5"), "positive_beta_argument")
        term = series = Decimal(1)
        for n in range(512):
            term *= (a+n)*(1-b+n)*x / ((a+1+n)*(n+1))
            series += term
            if not term or abs(term) <= abs(series)*Decimal("1e-80"):
                break
        else:
            raise AuditError("reference_beta_convergence")
        log_beta = Decimal.from_float(math.lgamma(float(a)) + math.lgamma(float(b)) - math.lgamma(float(a+b)))
        tail = (a*x.ln() - a.ln() - log_beta).exp() * series
        probability = 1-tail if complement else tail
        require(probability.is_finite() and 0 < probability <= 1, "positive_reference_probability")
        return probability


def student_reference(delta, se, df):
    from scipy.special import betainc, stdtr, stdtrit
    from scipy.optimize import brentq
    with localcontext(Context(prec=PRECISION)):
        statistic = abs(delta/se)
        expected_p = beta_tail(statistic, df)
    degree = float(df)
    critical = brentq(lambda t: float(betainc(degree/2, 0.5, degree/(degree+t*t))) - 0.05,
                      0.0, 10.0, xtol=5e-15)
    # These shared native functions verify availability only, not numerical P/CI.
    native_critical = float(stdtrit(degree, 0.975))
    native_statistic = float(statistic)
    native_p = 2*float(stdtr(degree, -native_statistic)) if math.isfinite(native_statistic) else float("nan")
    if not math.isfinite(native_critical) or native_critical <= 0:
        status = "unavailable_nonfinite_inference"
    elif not math.isfinite(native_p) or not 0 < native_p <= 1:
        status = "CI_available_P_unavailable_numeric_tail"
    else:
        status = "available_nominal_Welch"
    require(math.isfinite(critical) and critical > 0, "reference_critical")
    return {"critical": Decimal.from_float(critical), "P": expected_p, "status": status}


def reference_endpoint(measurements, selection):
    require(set(measurements) == set(SAMPLES) and len(selection["donors"]) == 8
            and {d["sample"] for d in selection["donors"]} == set(SAMPLES)
            and len({d["donor"] for d in selection["donors"]}) == 8
            and Counter(d["group"] for d in selection["donors"]) == {"SNP": 4, "noSNP": 4}, "whole_eight_gate")
    arguments, scores = {}, {}
    for donor in selection["donors"]:
        sample, m = donor["sample"], measurements[donor["sample"]]
        keys(m, {"B", "R", "B_positive_cells", "selected_cells"})
        require(all(type(x) is int for x in m.values()) and 0 <= m["B"] <= m["R"] < MAX_R and m["R"] > 0,
                "whole_endpoint_units")
        require(m["selected_cells"] == len(donor["selected"]) and
                0 <= m["B_positive_cells"] <= m["selected_cells"] and
                (m["B"] == 0) == (m["B_positive_cells"] == 0), "whole_endpoint_detection")
        arguments[sample] = Fraction(m["R"]+1_000_000*m["B"], m["R"])
        scores[sample] = direct_log2(arguments[sample])
    grouped = {g: [arguments[d["sample"]] for d in selection["donors"] if d["group"] == g] for g in GROUPS}
    with localcontext(Context(prec=PRECISION)):
        groups, variances = {}, {}
        for group in GROUPS:
            values = [Fraction(scores[d["sample"]]) for d in selection["donors"] if d["group"] == group]
            centre = sum(values, Fraction()) / 4
            variance = sum(((value-centre)**2 for value in values), Fraction()) / 3
            require((variance == 0) == (len(set(grouped[group])) == 1), "exact_constancy")
            variances[group] = variance
            groups[group] = {"n": 4, "mean": direct_log2(math.prod(grouped[group]))/4,
                             "sample_sd": dec(variance).sqrt()}
        delta = product_point(grouped)
        v1, v0 = variances["SNP"]/4, variances["noSNP"]/4
        se2 = v1+v0
        se = dec(se2).sqrt()
        df = dec(se2**2/(v1*v1/3+v0*v0/3)) if se2 else None
        student = student_reference(delta, se, df) if se2 else {"critical": None, "P": None, "status": "unavailable_zero_SE"}
        loo = {}
        for d in selection["donors"]:
            args = {g: [arguments[k["sample"]] for k in selection["donors"] if k["group"] == g and k["sample"] != d["sample"]]
                    for g in GROUPS}
            loo[d["sample"]] = product_point(args)
        summary = {"count": 8, "minimum": min(loo.values()), "maximum": max(loo.values()),
                   "positive_count": sum(x > 0 for x in loo.values()), "negative_count": sum(x < 0 for x in loo.values()),
                   "zero_count": sum(x == 0 for x in loo.values()), "full_point_sign": sign(delta),
                   "strict_sign_reversal_count": sum(sign(x)*sign(delta) < 0 for x in loo.values())}
    return {"scores": scores, "groups": groups, "delta": delta, "SE": se, "df": df,
            "student": student, "LOO": loo, "LOO_summary": summary}


def number(value):
    require(type(value) is str and len(value) <= 352 and
            re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", value), "finite_decimal_string")
    result = Decimal(value)
    require(result.is_finite() and (not result or -10000 <= result.adjusted() <= 10000), "bounded_numeric_output")
    return result


def agree(value, reference, tolerance="1e-10"):
    actual = number(value)
    require(sign(actual) == sign(reference), "numeric_zero_or_sign_disagreement")
    with localcontext(Context(prec=PRECISION)):
        require(actual == reference if reference == 0 else abs((actual-reference)/reference) <= Decimal(tolerance),
                "numeric_relative_disagreement")
    return actual


def public_shape(value, specification=PUBLIC_SCHEMA):
    if type(specification) is dict:
        keys(value, specification)
        for key in specification:
            public_shape(value[key], specification[key])
    elif specification.startswith("literal:"):
        require(type(value) is str and value == specification[8:], "public_literal")
    elif specification.startswith("enum:"):
        require(type(value) is str and value in specification[5:].split("|"), "public_enum")
    elif specification.startswith("integer"):
        require(type(value) is int, "public_integer")
        require(value in [int(n) for n in specification.split(":")[1].split("/")] if ":" in specification else value >= 0,
                "public_integer_range")
    elif specification == "sha256":
        hash_string(value)
    elif specification in ("decimal-string", "decimal-string-or-null"):
        if value is not None or specification == "decimal-string":
            number(value)
    elif specification == "two-decimal-strings-or-null":
        if value is not None:
            require(type(value) is list and len(value) == 2, "public_CI_shape")
            for x in value:
                number(x)
    else:
        raise AuditError("unknown_public_specification")


def inference_compare(inference, reference):
    keys(inference, {"delta", "SE", "df", "CI95", "P", "inference_status", "groups"})
    keys(inference["groups"], GROUPS)
    for group in GROUPS:
        value, expected = inference["groups"][group], reference["groups"][group]
        keys(value, {"n", "mean", "sample_sd"})
        require(type(value["n"]) is int and value["n"] == 4, "group_N")
        agree(value["mean"], expected["mean"])
        agree(value["sample_sd"], expected["sample_sd"], "1e-10")
    agree(inference["delta"], reference["delta"])
    agree(inference["SE"], reference["SE"], "1e-10")
    expected = reference["student"]
    status = expected["status"]
    require(inference["inference_status"] == status, "inference_availability")
    if reference["df"] is None:
        require(inference["df"] is None, "zero_SE_df")
    else:
        degree = agree(inference["df"], reference["df"], "1e-10")
        require(3 <= degree <= 6, "df_domain")
    if status in ("unavailable_zero_SE", "unavailable_nonfinite_inference"):
        require(inference["CI95"] is None and inference["P"] is None, "point_only_missing_inference")
        return
    ci = inference["CI95"]
    require(type(ci) is list and len(ci) == 2, "CI_shape")
    bounds = [number(x) for x in ci]
    require(bounds[0] < reference["delta"] < bounds[1], "CI_positive_width_and_centre")
    with localcontext(Context(prec=PRECISION)):
        for bound, direction in zip(bounds, (-1, 1)):
            standardized = (bound-reference["delta"])/reference["SE"]
            require(abs(standardized/(direction*expected["critical"])-1) <= Decimal("1e-9"), "CI_standardized_disagreement")
    if status == "CI_available_P_unavailable_numeric_tail":
        require(inference["P"] is None and expected["P"] > 0, "positive_tail_unavailable_not_zero")
    else:
        probability = agree(inference["P"], expected["P"], "1e-10")
        require(0 < probability <= 1 and (reference["delta"] != 0 or probability == 1), "probability_domain")


def provenance(plan, authorization_sha):
    return {"plan_sha256": digest(plan), "selection_sha256": plan["selection_sha256"],
            "source_manifest_sha256": plan["source_manifest_sha256"],
            "implementation_sha256": plan["implementation"][0]["sha256"], "runtime_sha256": digest(plan["runtime"]),
            "public_schema_sha256": plan["public_schema_sha256"], "authorization_sha256": authorization_sha,
            "candidate_sha256": CANDIDATE_SHA}


def compare_payloads(private, public, measurements, selection, plan, authorization_sha):
    """All eight private records and LOO points are checked, never returned."""
    reference = reference_endpoint(measurements, selection)
    keys(private, {"donors", "inference", "LOO_private", "LOO_summary", "zero_B_donors", "B_positive_selected_cells"})
    require(type(private["donors"]) is list and len(private["donors"]) == 8
            and [d["sample"] for d in private["donors"]] == [d["sample"] for d in selection["donors"]], "private_donor_inventory")
    for value, donor in zip(private["donors"], selection["donors"]):
        keys(value, {"sample", "donor", "group", "selected_cells", "state_cell_counts", "state_cell_fractions",
                     "B_source_units", "RNA_source_units", "B_positive_cells", "B_positive_cell_fraction", "Y"})
        m = measurements[donor["sample"]]
        for key in ("sample", "donor", "group"):
            require(value[key] == donor[key], "private_donor_metadata")
        for key in ("selected_cells", "B_positive_cells"):
            require(type(value[key]) is int and value[key] == m[key], "private_selected_quantities")
        require(value["B_source_units"] == str(m["B"]) and value["RNA_source_units"] == str(m["R"]), "independent_B_R")
        keys(value["state_cell_counts"], {"C0", "C1"})
        require(all(type(x) is int for x in value["state_cell_counts"].values()) and
                value["state_cell_counts"] == donor["state_counts"], "private_state_counts")
        keys(value["state_cell_fractions"], {"C0", "C1"})
        with localcontext(Context(prec=PRECISION)):
            for state in ("C0", "C1"):
                agree(value["state_cell_fractions"][state], dec(Fraction(donor["state_counts"][state], m["selected_cells"])))
            agree(value["B_positive_cell_fraction"], dec(Fraction(m["B_positive_cells"], m["selected_cells"])))
        agree(value["Y"], reference["scores"][donor["sample"]])
    inference_compare(private["inference"], reference)
    require(type(private["LOO_private"]) is list and len(private["LOO_private"]) == 8 and
            [x["omitted_sample"] for x in private["LOO_private"]] == [d["sample"] for d in selection["donors"]], "all_eight_LOO")
    for row in private["LOO_private"]:
        keys(row, {"omitted_sample", "delta"})
        agree(row["delta"], reference["LOO"][row["omitted_sample"]])
    keys(private["LOO_summary"], reference["LOO_summary"])
    for key, expected in reference["LOO_summary"].items():
        if key in ("minimum", "maximum"):
            agree(private["LOO_summary"][key], expected)
        else:
            require(type(private["LOO_summary"][key]) is int and private["LOO_summary"][key] == expected, "LOO_summary_counts")
    zero = sum(m["B"] == 0 for m in measurements.values())
    positive = sum(m["B_positive_cells"] for m in measurements.values())
    require(type(private["zero_B_donors"]) is int and private["zero_B_donors"] == zero and
            type(private["B_positive_selected_cells"]) is int and private["B_positive_selected_cells"] == positive, "private_QC")
    public_shape(public)
    # Exact correspondence to the already independently checked private outputs.
    inf = private["inference"]
    expected_contrast = {"SNP_minus_noSNP": inf["delta"], "SE": inf["SE"], "df": inf["df"], "CI95": inf["CI95"],
                         "P_two_sided": inf["P"], "inference_status": inf["inference_status"]}
    require(public["groups"] == inf["groups"] and public["contrast"] == expected_contrast and
            public["leave_one_donor_out"] == private["LOO_summary"], "public_private_agreement")
    require(public["QC"] == {"donors": 8, "selected_cells": 16659, "retained_cells": 55460, "RNA_rows": 36601,
            "zero_BACH2_donors_retained": zero, "BACH2_positive_selected_cells": positive,
            "source_and_membership_gates": "passed"}, "public_QC")
    require(public["provenance"] == provenance(plan, authorization_sha), "public_provenance")
    return {"donor_records_checked": 8, "whole_donor_LOO_points_checked": 8, "complete_matrix_streams_checked": 8,
            "all_private_metadata_and_units_agree": True, "all_aggregate_numbers_status_QC_provenance_agree": True,
            "public_nested_allowlist_passed": True}


def root_json(root, relative, expected, maximum):
    require(type(relative) is str and re.fullmatch(r"work/psc-bach2-naive-rna-freeze/[A-Za-z0-9][A-Za-z0-9_-]{0,79}\.json", relative),
            "dedicated_root_freeze_JSON_role")
    value = decode_json(read_bytes(root, relative, hash_string(expected), maximum))
    require(type(value) is dict, "root_JSON_object_required")
    return value


def author_authorization(plan, run_directory):
    return {"schema_version": 1, "status": AUTH_STATUS, "expression_authorized": True, "authorizer": "root",
            "one_fixed_extraction_only": True, "output_directory": run_directory,
            **{k: val for k, val in provenance(plan, "0"*64).items() if k != "authorization_sha256"}}


def seal_bindings(plan, plan_directory, run_directory, output, authorization_path, authorization_sha,
                  completion_sha, root=ROOT):
    """Expected schema only; this does not create or grant root authorization."""
    own_pins = {"adapter_sha256": sha(read_bytes(root, ADAPTER)), "adapter_tests_sha256": sha(read_bytes(root, TESTS)),
                "design_sha256": sha(read_bytes(root, DESIGN)), "verifier_sha256": VERIFIER_SHA,
                "verifier_tests_sha256": sha(read_bytes(root, "tests/test_verify_psc_bach2_naive_rna.py"))}
    return {"schema_version": 1, "status": ACK, "authorizer": "root", "real_matrix_audit_authorized": True,
            "one_fixed_endpoint_only": True, "plan_directory": plan_directory, "run_directory": run_directory,
            "audit_output_directory": output, "authorization_path": authorization_path,
            "completion_sha256": completion_sha, "protocol_sha256": digest(PROTOCOL),
            **provenance(plan, authorization_sha), **own_pins}


def authenticate_completed(plan, run_directory, completion_sha, authorization_sha, root=ROOT):
    require(type(run_directory) is str and re.fullmatch(re.escape(AUTHOR_BASE) + r"/runs/[A-Za-z0-9][A-Za-z0-9_-]{0,79}", run_directory),
            "completed_run_role")
    directory = inventory(root, run_directory, PAYLOADS | {"completion.json"})
    completion = decode_json(read_bytes(root, run_directory+"/completion.json", hash_string(completion_sha), 16384))
    keys(completion, {"schema_version", "status", "output_directory", "allowed_inventory", "provenance", "payloads"})
    keys(completion["payloads"], PAYLOADS)
    payloads = {}
    for name in sorted(PAYLOADS):
        pin = completion["payloads"][name]
        keys(pin, {"bytes", "sha256"})
        require(type(pin["bytes"]) is int and pin["bytes"] > 0, "payload_byte_type")
        maximum = 8192 if name == "authorization.used.json" else 2_000_000
        # All three opaque hashes pass before EITHER result JSON is parsed.
        data = read_bytes(root, run_directory+"/"+name, hash_string(pin["sha256"]), maximum)
        require(len(data) == pin["bytes"], "payload_byte_count")
        require((directory/"completion.json").stat().st_mtime_ns >= (directory/name).stat().st_mtime_ns, "completion_not_last_written")
        payloads[name] = data
    expected = {"schema_version": 1, "status": "completed_fixed_endpoint_bound_payload_inventory", "output_directory": run_directory,
                "allowed_inventory": sorted(PAYLOADS | {"completion.json"}), "provenance": provenance(plan, authorization_sha),
                "payloads": {name: {"bytes": len(body), "sha256": sha(body)} for name, body in payloads.items()}}
    require(canonical(completion) == canonical(expected), "completion_binding")
    require(sha(payloads["authorization.used.json"]) == authorization_sha, "used_authorization_hash")
    return payloads


def audit_real(*, acknowledge=False, plan_directory=None, plan_sha=None, run_directory=None,
               authorization_path=None, authorization_sha=None, seal_path=None, seal_sha=None,
               completion_sha=None, output=None, root=ROOT):
    """Future opt-in only. No caller can infer authorization from preparation."""
    require(acknowledge is True, "explicit_root_acknowledgment_required")
    require(all(type(x) is str and x for x in (plan_directory, plan_sha, run_directory, authorization_path,
                authorization_sha, seal_path, seal_sha, completion_sha, output)), "all_external_bindings_required")
    for h in (plan_sha, authorization_sha, seal_sha, completion_sha):
        hash_string(h)
    require(seal_path != authorization_path, "distinct_root_seal_and_authorization")
    # Role checks precede opens. A count file is never treated as a root JSON.
    seal = root_json(root, seal_path, seal_sha, 16384)
    require(seal.get("status") == ACK and seal.get("authorizer") == "root" and
            seal.get("real_matrix_audit_authorized") is True, "root_real_audit_authorization")
    auth = root_json(root, authorization_path, authorization_sha, 8192)
    fresh_output(root, output)
    plan, selection, _ = metadata_only(plan_directory, plan_sha, root)
    expected_seal = seal_bindings(plan, plan_directory, run_directory, output, authorization_path,
                                  authorization_sha, completion_sha, root)
    require(canonical(seal) == canonical(expected_seal), "exact_root_audit_seal")
    require(canonical(auth) == canonical(author_authorization(plan, run_directory)), "exact_native_author_authorization")
    payloads = authenticate_completed(plan, run_directory, completion_sha, authorization_sha, root)
    require(payloads["authorization.used.json"] == canonical(auth), "exact_used_authorization")
    # Exclusive creation consumes this audit output before any numerical decode.
    destination = fresh_output(root, output, create=True)
    try:
        measurements = {}
        for pin in plan["source_manifest"]["matrix_expected_pins_NOT_read"]:
            donor = next(d for d in selection["donors"] if d["sample"] == pin["sample"])
            path = path_at(root, pin["path"], regular=True)
            with path.open("rb") as handle:
                before = os.fstat(handle.fileno())
                measurements[pin["sample"]] = stream_matrix(handle, pin,
                    [c["column_1based"] for c in donor["selected"]], plan["RNA_rows_1based"])
                after = os.fstat(handle.fileno())
            require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
                    (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), "matrix_changed_during_stream")
        private = decode_json(payloads["measurements.private.json"])
        public = decode_json(payloads["public-aggregate.json"])
        comparisons = compare_payloads(private, public, measurements, selection, plan, authorization_sha)
        # Catch stale/replaced metadata, code, seals, inventories or outputs before PASS.
        metadata_only(plan_directory, plan_sha, root)
        require(root_json(root, seal_path, seal_sha, 16384) == expected_seal and
                root_json(root, authorization_path, authorization_sha, 8192) == auth, "root_seal_changed")
        require(seal_bindings(plan, plan_directory, run_directory, output, authorization_path,
                             authorization_sha, completion_sha, root) == expected_seal, "audit_implementation_changed")
        require(authenticate_completed(plan, run_directory, completion_sha, authorization_sha, root) == payloads, "completed_payloads_changed")
        receipt = {"schema_version": 1, "status": "INDEPENDENT_REAL_OUTPUT_AUDIT_PASS",
                   "scope": "one fixed BACH2 pooled C0+C1 endpoint; independent verification, not biological replication",
                   "comparisons": comparisons, "real_matrix_streams_checked": 8, "numeric_payloads_parsed": 2,
                   "root_seal_sha256": seal_sha, "completion_sha256": completion_sha,
                   "protocol_sha256": digest(PROTOCOL), "adapter_sha256": expected_seal["adapter_sha256"],
                   "verifier_sha256": VERIFIER_SHA, "provenance": provenance(plan, authorization_sha),
                   "donor_quantities_or_IDs_or_LOO_values_disclosed": False, "additional_analyses": 0}
        with (destination/"aggregate-verification.json").open("xb") as handle:
            handle.write(canonical(receipt))
        return receipt
    except Exception:
        (destination/"failure.json").write_bytes(canonical({"status": "whole_endpoint_audit_unavailable", "no_rescue": True}))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-directory", default=AUTHOR_BASE+"/candidate-v4")
    parser.add_argument("--plan-sha256")
    parser.add_argument("--real-audit", action="store_true")
    parser.add_argument("--acknowledge-root-authorized-real-audit", action="store_true")
    parser.add_argument("--run-directory")
    parser.add_argument("--authorization")
    parser.add_argument("--authorization-sha256")
    parser.add_argument("--root-seal")
    parser.add_argument("--root-seal-sha256")
    parser.add_argument("--completion-sha256")
    parser.add_argument("--output-directory")
    args = parser.parse_args(argv)
    try:
        if args.real_audit:
            result = audit_real(acknowledge=args.acknowledge_root_authorized_real_audit,
                plan_directory=args.plan_directory, plan_sha=args.plan_sha256, run_directory=args.run_directory,
                authorization_path=args.authorization, authorization_sha=args.authorization_sha256,
                seal_path=args.root_seal, seal_sha=args.root_seal_sha256, completion_sha=args.completion_sha256,
                output=args.output_directory)
        else:
            require(not any((args.acknowledge_root_authorized_real_audit, args.run_directory, args.authorization,
                    args.authorization_sha256, args.root_seal, args.root_seal_sha256, args.completion_sha256, args.output_directory)),
                    "real_options_require_real_audit")
            _, _, result = metadata_only(args.plan_directory, args.plan_sha256)
        print(json.dumps(result, sort_keys=True, allow_nan=False))
        return 0
    except Exception:
        # A path, donor label, quantity or raw exception message must not leak.
        print(json.dumps({"status": "whole_endpoint_audit_unavailable", "no_rescue": True}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
