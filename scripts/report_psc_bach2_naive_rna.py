#!/usr/bin/env python3
"""Presentation only for one hash-approved public BACH2 aggregate. No scientific rerun."""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import sys

ROOT = Path(__file__).resolve().parents[1]
WORK = "work/psc-bach2-naive-rna-reporting"
DATA = "data/derived/psc-bach2-naive-rna"
FIGURE = "reports/figures/psc-bach2-naive-rna"
HANDOFF_PATH = "work/psc-bach2-naive-rna-freeze/root-reporting-handoff.json"
HANDOFF_SHA = "53f636477b1f8c33c44308b3b48ab3f3a9b448f53422dc17f1e07fbdd935fab6"
AGGREGATE_SHA = "3f9383678afc4a576c12f89225e3148ca07ff65486c8187e8579f4dfd4ac9837"
MPL_VERSION = "3.11.2"
FONT_SHA = "3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22"
SVG_HASHSALT = "psc-bach2-naive-rna-public-report-v1"
SOURCE_KEYS = ("outer_freeze", "primary_completion", "public_aggregate", "public_schema",
               "independent_verification", "replay_verification")
COMMON = ("schema_version", "status", "endpoint", "quantity")
PROVENANCE = ("authorization_sha256", "candidate_sha256", "implementation_sha256",
              "plan_sha256", "public_schema_sha256", "runtime_sha256",
              "selection_sha256", "source_manifest_sha256")
PROVENANCE_COLUMNS = tuple("provenance." + key for key in PROVENANCE)
ENDPOINT_COLUMNS = COMMON + ("SNP_minus_noSNP", "SE", "df", "CI95_lower", "CI95_upper",
                             "P_two_sided", "inference_status") + PROVENANCE_COLUMNS
GROUP_COLUMNS = COMMON + ("group", "n", "mean", "sample_sd") + PROVENANCE_COLUMNS
PAYLOAD_PATHS = tuple(sorted([DATA + "/" + name for name in (
    "public-aggregate.json", "public-schema.json", "endpoint.csv", "groups.csv",
    "influence-qc.json", "aggregate-verification.json", "mechanical-replay-verification.json",
    "reporting-info.json")] + [FIGURE + ".png", FIGURE + ".svg"]))
MANIFEST = DATA + "/manifest.json"
COMPLETION = DATA + "/completion.json"
INVENTORY = tuple(sorted(PAYLOAD_PATHS + (MANIFEST, COMPLETION)))


class ReportingError(ValueError):
    """A closed reporting gate failed."""


def require(condition, message):
    if not condition:
        raise ReportingError(message)


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def json_bytes(obj):
    return (json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def parse_json(blob):
    def pairs(items):
        obj = {}
        for key, value in items:
            require(key not in obj, "duplicate JSON key")
            obj[key] = value
        return obj

    def invalid(value):
        raise ReportingError("nonfinite JSON constant: " + value)

    try:
        return json.loads(blob, object_pairs_hook=pairs, parse_constant=invalid)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ReportingError("invalid JSON") from exc


def safe_path(relative):
    """Canonical repository-relative paths; no symlinks, hardlink aliases or special files."""
    require(type(relative) is str and relative and "\\" not in relative, "invalid path")
    parts = relative.split("/")
    require(not PurePosixPath(relative).is_absolute() and all(
        p not in ("", ".", "..") and re.fullmatch(r"[A-Za-z0-9._-]+", p) for p in parts),
        "path alias or traversal")
    path = ROOT
    require(ROOT.absolute() == ROOT.resolve(), "repository root is an alias")
    for index, part in enumerate(parts):
        path = path / part
        try:
            info = path.lstat()
        except FileNotFoundError:
            continue
        require(not stat.S_ISLNK(info.st_mode), "symlink rejected")
        if stat.S_ISDIR(info.st_mode):
            continue
        require(index == len(parts) - 1 and stat.S_ISREG(info.st_mode) and info.st_nlink == 1,
                "nonregular file, hardlink alias or nondirectory parent")
    return path


def read_file(relative, maximum=1_000_000):
    path = safe_path(relative)
    require(path.is_file() and path.stat().st_size <= maximum, "missing or oversized file")
    with path.open("rb") as stream:
        blob = stream.read(maximum + 1)
    require(len(blob) <= maximum, "oversized file")
    return blob


def read_pin(pin):
    require(set(pin) == {"path", "bytes", "sha256"}, "invalid input pin")
    blob = read_file(pin["path"], pin["bytes"])
    require(len(blob) == pin["bytes"] and digest(blob) == pin["sha256"], "input pin drift")
    return blob


def decimal_string(value):
    require(type(value) is str and len(value) <= 2048 and re.fullmatch(
        r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", value),
        "invalid decimal string")
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ReportingError("invalid decimal") from exc
    require(number.is_finite() and -10000 <= number.adjusted() <= 10000,
            "nonfinite or unbounded decimal")
    return number


def validate_tree(value, schema):
    """Interpret only the pinned schema's small descriptor vocabulary, recursively."""
    if type(schema) is dict:
        require(type(value) is dict and set(value) == set(schema), "nested key/type mismatch")
        for key in schema:
            validate_tree(value[key], schema[key])
    elif schema.startswith("integer"):
        require(type(value) is int, "integer required (not boolean)")
        if ":" in schema:
            require(value in [int(v) for v in schema.split(":", 1)[1].split("/")], "fixed integer drift")
    elif schema.startswith("literal:"):
        require(type(value) is str and value == schema[8:], "literal drift")
    elif schema.startswith("enum:"):
        require(type(value) is str and value in schema[5:].split("|"), "unknown status")
    elif schema == "sha256":
        require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value), "invalid SHA256")
    elif schema in ("decimal-string", "decimal-string-or-null"):
        if value is not None or schema == "decimal-string":
            decimal_string(value)
    elif schema == "two-decimal-strings-or-null":
        if value is not None:
            require(type(value) is list and len(value) == 2, "invalid CI grid")
            for endpoint in value:
                decimal_string(endpoint)
    else:
        raise ReportingError("unknown pinned schema descriptor")


def validate_public(public, schema, expected_provenance):
    validate_tree(public, schema)
    require(public["provenance"] == expected_provenance, "provenance hash drift")
    contrast = public["contrast"]
    # This reporter accepts only the actual, hash-approved availability state.
    require(contrast["inference_status"] == "available_nominal_Welch", "unseen status: fail closed")
    require(contrast["CI95"] is not None and contrast["P_two_sided"] is not None
            and contrast["df"] is not None, "available status requires complete CI/P/df")
    point, se, df, p = [decimal_string(contrast[k]) for k in
                        ("SNP_minus_noSNP", "SE", "df", "P_two_sided")]
    lower, upper = map(decimal_string, contrast["CI95"])
    require(se > 0 and 3 <= df <= 6 and 0 < p <= 1 and lower < point < upper,
            "invalid available inference/CI")
    for group in ("SNP", "noSNP"):
        require(decimal_string(public["groups"][group]["mean"]) >= 0
                and decimal_string(public["groups"][group]["sample_sd"]) >= 0, "invalid group summary")
    loo = public["leave_one_donor_out"]
    low, high = decimal_string(loo["minimum"]), decimal_string(loo["maximum"])
    negative, positive, zero, reversal = [loo[k] for k in
        ("negative_count", "positive_count", "zero_count", "strict_sign_reversal_count")]
    sign = int(point > 0) - int(point < 0)
    require(all(0 <= c <= 8 for c in (negative, positive, zero, reversal))
            and negative + positive + zero == 8 and loo["full_point_sign"] == sign
            and reversal == (positive if sign < 0 else negative if sign > 0 else 0)
            and low <= high and (low < 0) == (negative > 0) and (high > 0) == (positive > 0)
            and (not zero or low <= 0 <= high)
            and (negative != 8 or high < 0) and (positive != 8 or low > 0)
            and (zero != 8 or low == high == 0), "invalid complete omission summary")
    qc = public["QC"]
    require(0 <= qc["BACH2_positive_selected_cells"] <= qc["selected_cells"]
            and 0 <= qc["zero_BACH2_donors_retained"] <= 8, "invalid aggregate QC")


def load_sources():
    pin = {"path": HANDOFF_PATH, "bytes": 2815, "sha256": HANDOFF_SHA}
    handoff = parse_json(read_pin(pin))
    blobs = {key: read_pin(handoff[key]) for key in SOURCE_KEYS}
    objects = {key: parse_json(blob) for key, blob in blobs.items()}
    freeze, native = objects["outer_freeze"], objects["primary_completion"]
    expected = {
        "authorization_sha256": freeze["native_per_destination_authorizations"]["primary-v1"]["sha256"],
        "candidate_sha256": freeze["candidate_source_file_sha256"],
        "implementation_sha256": next(p["sha256"] for p in freeze["public_code_protocol_review_and_schema_pins"]
                                      if p["path"] == "scripts/analyze_psc_bach2_naive_rna.py"),
        "plan_sha256": freeze["compiled_plan"]["sha256"],
        "public_schema_sha256": handoff["public_schema"]["sha256"],
        "runtime_sha256": freeze["runtime"]["runtime_sha256"],
        "selection_sha256": freeze["selection_sha256"],
        "source_manifest_sha256": freeze["source_manifest_sha256"],
    }
    validate_public(objects["public_aggregate"], objects["public_schema"], expected)
    require(digest(blobs["public_aggregate"]) == AGGREGATE_SHA, "unapproved aggregate")
    require(native["provenance"] == expected == objects["independent_verification"]["provenance"],
            "metadata provenance drift")
    require(objects["independent_verification"]["completion_sha256"] == handoff["primary_completion"]["sha256"]
            == objects["replay_verification"]["primary_completion_sha256"], "completion binding drift")
    # Only these seven approved small files are read. No following of metadata references.
    return handoff, blobs, objects


def csv_bytes(columns, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def render_figures(public):
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.font_manager import FontProperties

    require(matplotlib.__version__ == MPL_VERSION, "rendering version drift")
    font = Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans.ttf"
    require(digest(font.read_bytes()) == FONT_SHA, "rendering font drift")
    c = public["contrast"]
    numeric = [float(decimal_string(v)) for v in (c["SNP_minus_noSNP"], *c["CI95"])]
    require(all(math.isfinite(v) for v in numeric), "nonfinite display coordinate")
    point, lower, upper = numeric
    require(-1.2 < lower < point < upper < 0.4, "changed figure range: fail closed")
    creator = "PSC BACH2 aggregate reporter v1; Matplotlib " + MPL_VERSION
    rc = {"svg.hashsalt": SVG_HASHSALT, "svg.fonttype": "path", "font.family": "DejaVu Sans",
          "font.weight": "normal", "font.style": "normal", "text.usetex": False,
          "axes.unicode_minus": False}
    with matplotlib.rc_context({**matplotlib.rcParamsDefault, **rc}):
        fp = FontProperties(fname=str(font), size=10)
        fig = Figure(figsize=(6.7, 3.1), facecolor="white")
        ax = fig.add_axes((0.12, 0.34, 0.83, 0.35))
        ax.axvline(0, color="#777777", linewidth=1, linestyle="--")
        ax.hlines(0, lower, upper, color="#225577", linewidth=2)
        ax.vlines([lower, upper], -0.06, 0.06, color="#225577", linewidth=1.5)
        ax.plot([point], [0], marker="o", color="#225577", markersize=6)
        ax.set(xlim=(-1.2, 0.4), ylim=(-0.4, 0.4), yticks=[], xticks=[-1, -0.5, 0])
        ax.set_xticklabels(["-1.0", "-0.5", "0"], fontproperties=fp)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel("Difference in mean donor log2(1 + released-RNA CPM)", fontproperties=fp, labelpad=8)
        fig.text(0.12, 0.9, "BACH2 | fixed pooled C0+C1 naive CD4", fontproperties=fp, fontsize=12)
        fig.text(0.12, 0.8, "SNP - noSNP; 4 donors per group; one nominal Welch 95% CI", fontproperties=fp)
        label = (f"Delta {Decimal(c['SNP_minus_noSNP']):.3f}; 95% CI "
                 f"[{Decimal(c['CI95'][0]):.3f}, {Decimal(c['CI95'][1]):.3f}]; "
                 f"nominal two-sided P = {Decimal(c['P_two_sided']):.3f}")
        fig.text(0.12, 0.13, label, fontproperties=fp)
        fig.text(0.12, 0.04, "All 36,601 released RNA rows in the denominator; display rounding only.",
                 fontproperties=fp, fontsize=9)
        images = {}
        for extension, metadata in (("svg", {"Date": None, "Creator": creator,
                                             "Title": "Single fixed BACH2 released-RNA contrast"}),
                                    ("png", {"Software": creator})):
            buffer = io.BytesIO()
            fig.savefig(buffer, format=extension, dpi=160, metadata=metadata)
            images[FIGURE + "." + extension] = buffer.getvalue()
    return images


def prepare_payloads():
    handoff, blobs, objects = load_sources()
    public = objects["public_aggregate"]
    common = {key: public[key] for key in COMMON}
    common.update({"provenance." + k: public["provenance"][k] for k in PROVENANCE})
    contrast = public["contrast"]
    endpoint = {**common, **{k: v for k, v in contrast.items() if k != "CI95"},
                "CI95_lower": contrast["CI95"][0], "CI95_upper": contrast["CI95"][1]}
    groups = [{**common, "group": group, **public["groups"][group]} for group in ("SNP", "noSNP")]
    influence = {key: public[key] for key in COMMON + ("QC", "leave_one_donor_out", "provenance")}
    info = {
        "schema_version": 1, "status": "PRESENTATION_ONLY_APPROVED_PUBLIC_AGGREGATE",
        "new_statistical_calculations": 0, "private_or_source_quantity_payloads_read": 0,
        "source_verification_scope": "Copied pinned audit/replay receipts; reporter does not recheck matrices or private outputs.",
        "data_precision": "Original decimal strings; JSON exact source bytes. Display rounding only in figures.",
        "figure": {"contrast_count": 1, "interval": "nominal two-sided actual-df Welch 95% CI, not group SD",
                   "donor_or_cell_points": False, "matplotlib": MPL_VERSION, "font": "DejaVuSans.ttf",
                   "font_sha256": FONT_SHA, "svg_hashsalt": SVG_HASHSALT, "svg_fonts": "paths",
                   "date_metadata": None, "dpi": 160},
        "python": sys.version,
        "privacy": "Aggregate n=4 groups are not formal differential privacy or a guarantee against inferential disclosure.",
        "interpretation": "Donors, not cells, are replicates. Released RNA units are not certified untouched UMIs. No independent biological replication.",
    }
    payloads = {
        DATA + "/public-aggregate.json": blobs["public_aggregate"],
        DATA + "/public-schema.json": blobs["public_schema"],
        DATA + "/endpoint.csv": csv_bytes(ENDPOINT_COLUMNS, [endpoint]),
        DATA + "/groups.csv": csv_bytes(GROUP_COLUMNS, groups),
        DATA + "/influence-qc.json": json_bytes(influence),
        DATA + "/aggregate-verification.json": blobs["independent_verification"],
        DATA + "/mechanical-replay-verification.json": blobs["replay_verification"],
        DATA + "/reporting-info.json": json_bytes(info),
        **render_figures(public),
    }
    require(tuple(sorted(payloads)) == PAYLOAD_PATHS, "presentation inventory mismatch")
    reporter = Path(__file__).read_bytes()
    manifest = {
        "schema_version": 1, "status": "COMPLETE_FIXED_AGGREGATE_PRESENTATION_MANIFEST",
        "reporter": {"path": "scripts/report_psc_bach2_naive_rna.py", "bytes": len(reporter), "sha256": digest(reporter)},
        "authority": {"path": HANDOFF_PATH, "bytes": 2815, "sha256": HANDOFF_SHA},
        "inputs": {key: handoff[key] for key in SOURCE_KEYS},
        "schemas": {"endpoint.csv": {"rows": 1, "columns": list(ENDPOINT_COLUMNS)},
                    "groups.csv": {"rows": 2, "grid": ["SNP", "noSNP"], "columns": list(GROUP_COLUMNS)},
                    "influence-qc.json": {"keys": sorted(influence), "nested_fields": "Exact corresponding public-schema.json subtrees"},
                    "public-aggregate.json": "Exact byte copy; recursively validated against public-schema.json and frozen provenance"},
        "payloads": {key: {"bytes": len(blob), "sha256": digest(blob)} for key, blob in sorted(payloads.items())},
    }
    payloads[MANIFEST] = json_bytes(manifest)
    return payloads


def preview_path(relative, fresh=False):
    require(type(relative) is str and re.fullmatch(re.escape(WORK) + r"/[A-Za-z0-9][A-Za-z0-9_-]{0,63}", relative),
            "preview must be a simple fresh name in the ignored reporting namespace")
    path = safe_path(relative)
    require("work/" in read_file(".gitignore", 16384).decode().splitlines(), "reporting namespace is not ignored")
    require(not fresh or not path.exists(), "existing preview rejected")
    return path


def completion_bytes(relative, payloads):
    return json_bytes({"schema_version": 1, "status": "COMPLETE_APPROVED_PUBLIC_PRESENTATION",
                       "origin_preview_directory": relative, "allowed_inventory": list(INVENTORY),
                       "manifest": {"bytes": len(payloads[MANIFEST]), "sha256": digest(payloads[MANIFEST])}})


def write_new(path, blob):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(blob)
        stream.flush()
        os.fsync(stream.fileno())


def write_completion(path, blob):
    temporary = path.with_name(".completion.tmp")
    write_new(temporary, blob)
    os.replace(temporary, path)


def preview(relative):
    path = preview_path(relative, fresh=True)
    payloads = prepare_payloads()  # Validate before creating any output directory.
    safe_path(WORK).mkdir(parents=True, exist_ok=True)
    path.mkdir(mode=0o700)
    try:
        for name, blob in payloads.items():
            write_new(path / name, blob)
        complete = completion_bytes(relative, payloads)
        write_completion(path / COMPLETION, complete)  # Atomic and last.
    except BaseException:
        shutil.rmtree(path)
        raise
    return {"bundle": relative, "completion_sha256": digest(complete), "files": len(INVENTORY)}


def verified_payloads(relative, completion_sha):
    path = preview_path(relative)
    require(path.is_dir(), "missing bundle")
    files, directories = [], []
    for base, dirs, names in os.walk(path, followlinks=False):
        for name in dirs + names:
            entry = Path(base) / name
            safe_path(str(entry.relative_to(ROOT)))
            (directories if entry.is_dir() else files).append(str(entry.relative_to(path)))
    expected_dirs = {str(parent) for name in INVENTORY for parent in PurePosixPath(name).parents if str(parent) != "."}
    require(set(files) == set(INVENTORY) and set(directories) == expected_dirs, "incomplete or extra bundle inventory")
    payloads = prepare_payloads()
    complete = completion_bytes(relative, payloads)
    require(digest(complete) == completion_sha, "completion authentication failed")
    for name, expected in {**payloads, COMPLETION: complete}.items():
        require(read_file(relative + "/" + name) == expected, "presentation replay mismatch: " + name)
    return {**payloads, COMPLETION: complete}


def verify(relative, completion_sha):
    verified_payloads(relative, completion_sha)
    return {"bundle": relative, "completion_sha256": completion_sha, "status": "EXACT_PRESENTATION_REPLAY_PASS"}


def approval_template(relative, completion_sha):
    return {"schema_version": 1, "authorizer": "root", "status": "APPROVED_FIXED_PUBLIC_REPORTING_BUNDLE",
            "bundle": relative, "completion_sha256": completion_sha,
            "reporter_sha256": digest(Path(__file__).read_bytes()), "public_aggregate_sha256": AGGREGATE_SHA,
            "destinations": [DATA, FIGURE + ".png", FIGURE + ".svg"]}


def publish(relative, completion_sha, approval_path, approval_sha):
    # Root-only later action. This function is not used to create a preview or run a test replay.
    payloads = verified_payloads(relative, completion_sha)
    require(re.fullmatch(r"work/psc-bach2-naive-rna-freeze/[A-Za-z0-9][A-Za-z0-9_-]*\.json", approval_path),
            "root approval must be in the dedicated root namespace")
    approval = read_file(approval_path, 8192)
    require(digest(approval) == approval_sha and json_bytes(parse_json(approval))
            == json_bytes(approval_template(relative, completion_sha)),
            "root approval mismatch")
    destinations = [safe_path(p) for p in (DATA, FIGURE + ".png", FIGURE + ".svg")]
    require(all(not p.exists() for p in destinations), "existing publication destination rejected")
    created = []
    try:
        destinations[0].mkdir()
        created.append(destinations[0])
        for name in INVENTORY:
            if name == COMPLETION:
                continue
            destination = safe_path(name)
            if name.startswith(FIGURE + "."):
                created.append(destination)
            write_new(destination, payloads[name])
        write_completion(safe_path(COMPLETION), payloads[COMPLETION])
    except BaseException:
        for path in reversed(created):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink(missing_ok=True)
        raise
    return {"status": "APPROVED_PRESENTATION_PUBLISHED", "approval_sha256": approval_sha,
            "completion_sha256": completion_sha}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("validate", help="read approved public aggregate/metadata only; write nothing")
    p = sub.add_parser("preview", help="fresh ignored presentation only")
    p.add_argument("--output", required=True)
    for command in ("verify", "approval-template", "publish"):
        p = sub.add_parser(command)
        p.add_argument("--bundle", required=True)
        p.add_argument("--completion-sha256", required=True)
        if command == "publish":
            p.add_argument("--approval", required=True)
            p.add_argument("--approval-sha256", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command in (None, "validate"):
            load_sources()
            result = {"status": "APPROVED_AGGREGATE_SCHEMA_AND_METADATA_PINS_PASS", "outputs_written": 0}
        elif args.command == "preview":
            result = preview(args.output)
        elif args.command == "verify":
            result = verify(args.bundle, args.completion_sha256)
        elif args.command == "approval-template":
            verify(args.bundle, args.completion_sha256)
            result = approval_template(args.bundle, args.completion_sha256)
        else:
            result = publish(args.bundle, args.completion_sha256, args.approval, args.approval_sha256)
        print(json_bytes(result).decode(), end="")
        return 0
    except (ReportingError, OSError, KeyError, TypeError) as exc:
        print("REPORTING_FAILED: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
