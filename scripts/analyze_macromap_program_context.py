"""MacroMap source-only preparation and tested numerical building blocks.

The CLI does not read numerical count fields or run real program scores.
The parent must freeze a separate execution plan before calling the numerical
helpers on MacroMap. No existing ETS2 protocol is changed by this module.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

try:
    from . import analyze_ets2_program as original
except ImportError:
    import analyze_ets2_program as original

ROOT = Path(__file__).resolve().parents[1]
RAW = "data/raw/macromap-qualification"
PROGRAMS = original.PROGRAM_IDS
COMPARATORS = original.COMPARATOR_PROGRAM_IDS
STIMULI = ("CIL", "IFNB", "IFNG", "IL4", "LIL10", "MBP", "P3C", "R848", "sLPS")
REGISTRY_STIMULI = STIMULI + ("PIC",)
TIMES = (6, 24)
PROTOCOLS = ("Mod_SmartSeq2", "NEB")
FOLD_SEED = 20260916
BOOTSTRAP_SEED = 2026091601
N_FOLDS = 5
RIDGE_LAMBDA = 0.1
MATCH_SEED = 20260912
MATCH_REPLICATES = 1000
MATCH_BINS = 20
MIN_LINES = 3
MIN_RUNS = 3
MAX_ITER = 2000
GRADIENT_TOL = 1e-6
FEATURE_FIELDS = ("Chr", "Geneid", "Start", "End", "Strand", "Length")

# Exact local reuse pins. No network fallback, automatic alias repair or fresh
# program membership selection is permitted by the preparation command.
SOURCE_PINS = {
    "config/ets2-program-plan.json": {
        "bytes": 10598,
        "sha256": "234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9"
    },
    "config/ets2-validation-amended-rules.json": {
        "bytes": 11396,
        "sha256": "cbbaa6e626dad9f494367cf1a8c45c0db270f4c330f9bcd61a72280fe1b8ef6b"
    },
    "data/derived/ets2-program-membership.csv": {
        "bytes": 1669240,
        "sha256": "b6e468f000802e6f1ec6baa53722e8a7d3e0f98d6df4aa2bb49ef886efc199a3"
    },
    "data/derived/ets2-source-gene-sets.csv": {
        "bytes": 965247,
        "sha256": "d700811a139cb0dc012f5af8fe22ced452b4d1f76fe5dabb3afebf7ea61e2524"
    },
    "reports/ets2-program-audit.json": {
        "bytes": 37753,
        "sha256": "c2589587f3f217d44623c52db9787d8c8d71834fbf1f070799235d42a993d69a"
    },
    "scripts/analyze_ets2_program.py": {
        "bytes": 87516,
        "sha256": "7e7c73fab9f4bef31f535f4a3fb8bfbff2d9642b473542533a42aa15c2be6092"
    },
    "config/ets2-benchmark-sources.json": {
        "bytes": 1422,
        "sha256": "60ca96b13c02931ac14156f76df0caa0632ffadb53386759e8e425018628796a"
    },
    "data/raw/macromap-qualification/MacroMap_raw_expression.txt.gz": {
        "bytes": 174122135,
        "sha256": "eb9dd32b6f25ae9912e47dcc5e1f5170a98af919a42d48a4bca3a680fc5353d9"
    },
    "data/raw/macromap-qualification/MacroMap_raw_expression.header.txt": {
        "bytes": 109238,
        "sha256": "59e360723555df694ff50115dfb796465fbfc30347c7fc833fde83126e0408f0"
    },
    "data/raw/macromap-qualification/raw-gene-universe.tsv": {
        "bytes": 2652422,
        "sha256": "da50bdcaa04f0c6122aa77282673733fefa0f578e19ca802278607238ac1b35c"
    },
    "data/raw/macromap-qualification/matrix-download-receipt.json": {
        "bytes": 863,
        "sha256": "82bd72138271f83c922988bf69d25362efdf4e48b568e799ba818eac7a487ba1"
    },
    "data/raw/macromap-qualification/matrix-structural-audit.json": {
        "bytes": 1155,
        "sha256": "a4e8d73e7cd8628e6218246b7199d2c2bf4a380789770bae0deec239f0bed515"
    },
    "data/raw/macromap-qualification/gene-universe-audit.json": {
        "bytes": 2327,
        "sha256": "4aa0ca7f632f9ddfce9a50492a96efc5d4f81eddfdae63bb116cb2294f448f86"
    },
    "data/raw/macromap-qualification/worker-methods/sample-metadata.csv": {
        "bytes": 798498,
        "sha256": "b6f20190d43c7918f0dc31f19a26fec8020058c744fac66444991fabc1845987"
    },
    "data/raw/macromap-qualification/worker-methods/expression-header-sample-map.csv": {
        "bytes": 736731,
        "sha256": "0fa1b2cf5c0cd61c2bfb03bec722aafacd1fa7f7df2f1c3fd303dee86f37165b"
    },
    "data/raw/macromap-qualification/worker-methods/github_gencode_v27_annotation.body": {
        "bytes": 3079419,
        "sha256": "db11e855d959656bc3a25de65460ae12368ca13b1ae1a5493973f6ea433bdf41"
    },
    "data/raw/macromap-qualification/worker-methods/publication_eqtl_supp11_lines.body": {
        "bytes": 526439,
        "sha256": "97805429aaf777d8d13f625375a0dee124f130a2499fd6237bec279295bd89d8"
    },
    "data/raw/macromap-qualification/worker-methods/publication_eqtl_fulltext_xml.body": {
        "bytes": 137743,
        "sha256": "d3f69f885db807ff964f9158c7b7f581bea649225d6c76a37e4d3ef34f38cb8f"
    },
    "data/raw/macromap-qualification/worker-methods/source-ledger.json": {
        "bytes": 33437,
        "sha256": "7e6a7c818e938ad889515710b1d417dba19598d896447ab47c7b1ed584ecfd79"
    }
}


def sha256(path: Path) -> str:
    reject_path_aliases(path)
    with path.open("rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def verify_pins(root: Path, pins: Mapping[str, Mapping[str, Any]] = SOURCE_PINS) -> dict:
    checked = {}
    for relative, pin in sorted(pins.items()):
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
            raise ValueError(f"Source pin mismatch: {relative}")
        checked[relative] = dict(pin)
    return checked


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    reject_path_aliases(path)
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing/duplicated columns: {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Ragged table: {path}")
    return rows


def _absolute(path: Path) -> Path:
    # Never erase a lexical parent step before a caller opens its original
    # path: link/../file can mean a different inode to the OS than abspath.
    if ".." in Path(path).parts:
        raise ValueError("Parent-traversal components are forbidden in source/output paths")
    return Path(os.path.abspath(path))


def reject_path_aliases(path: Path) -> None:
    """Reject symlinks, including dangling links, in every path component."""
    absolute = _absolute(path)
    for component in (*reversed(absolute.parents), absolute):
        if component.is_symlink():
            raise ValueError(f"Output path contains a symlink: {component}")


def fresh_private_run(root: Path, output: Path) -> Path:
    """Create a new nonaliased directory under the literal ignored anchor.

    Resolve neither the anchor nor output through a symlink. Existing run
    directories are immutable snapshots, even if currently empty. Exclusive
    artifact writes additionally reject aliases inserted as filenames.
    """
    root, output = _absolute(root), _absolute(output)
    anchor = root / "work/macromap-design"
    if output == anchor or not output.is_relative_to(anchor):
        raise ValueError("Preparation output must be a fresh subdirectory of ignored work/macromap-design")
    reject_path_aliases(anchor)
    reject_path_aliases(output)
    anchor.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    reject_path_aliases(output)
    output.mkdir(exist_ok=False)
    return output


def exclusive_text(path: Path, text: str) -> None:
    reject_path_aliases(path)
    with path.open("x", encoding="utf-8") as fh:
        fh.write(text)


def publish_without_replacement(source: Path, destination: Path) -> None:
    """Atomic publication of a finished same-filesystem file, never overwrite.

    Unlike rename/replace, link raises if the destination exists, including a
    dangling symlink. The temporary link is removed only after publication.
    """
    reject_path_aliases(source)
    reject_path_aliases(destination)
    os.link(source, destination, follow_symlinks=False)
    source.unlink()


def write_csv(path: Path, rows: Sequence[Mapping]) -> None:
    if not rows:
        raise ValueError("Cannot write a schema-free empty table")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    reject_path_aliases(path)
    with path.open("x", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def xlsx_sheet(path: Path, sheet_name: str) -> list[dict[str, str]]:
    """Read the cached RNA-seq worksheet without an extra XLSX dependency.

    Text cells remain exact. Numeric values use the same Python float/int text
    conversion as the preceding openpyxl-to-CSV extraction. No expression
    workbook is read by this function in the source-only command.
    """
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    relns = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    with zipfile.ZipFile(path) as zf:
        strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            tree = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            strings = ["".join(item.itertext()) for item in tree.findall("s:si", ns)]
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        sheets = [item for item in workbook.findall("s:sheets/s:sheet", ns)
                  if item.attrib["name"] == sheet_name]
        if len(sheets) != 1:
            raise ValueError("Expected one exact worksheet")
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = [r.attrib["Target"] for r in rels
                  if r.attrib["Id"] == sheets[0].attrib[relns]]
        if len(target) != 1:
            raise ValueError("Worksheet relationship not unique")
        member = target[0].lstrip("/") if target[0].startswith("/") else "xl/" + target[0]
        tree = ET.fromstring(zf.read(member))
    cells = []
    for row in tree.findall("s:sheetData/s:row", ns):
        values = {}
        for cell in row.findall("s:c", ns):
            coltext = re.match(r"[A-Z]+", cell.attrib["r"]).group()
            col = 0
            for letter in coltext:
                col = col * 26 + ord(letter) - ord("A") + 1
            node = cell.find("s:v", ns)
            text = node.text if node is not None and node.text is not None else ""
            kind = cell.attrib.get("t", "n")
            if kind == "s":
                text = strings[int(text)]
            elif kind == "inlineStr":
                text = "".join(cell.find("s:is", ns).itertext())
            elif kind == "n" and text:
                value = float(text)
                text = str(int(value)) if value.is_integer() else str(value)
            values[col - 1] = text
        cells.append(values)
    header = [cells[0][i] for i in range(len(cells[0]))]
    if len(header) != len(set(header)):
        raise ValueError("Duplicated worksheet header")
    return [{key: row.get(i, "") for i, key in enumerate(header)} for row in cells[1:]]


def retained_source_programs(root: Path) -> dict[str, list[str]]:
    """Recover ALL original source values, not the atlas's mapped subset.

    The old ZIP is not required to exist. This route uses the exact hash-pinned
    complete source-value audit, whose sc/sn copies must agree, and independently
    checks each of the four Ensembl lists against the old source-set table.
    Missing or ambiguous atlas mappings are deliberately NOT filtered out.
    """
    rows = read_csv(root / "data/derived/ets2-program-membership.csv")
    old_audit = json.loads((root / "reports/ets2-program-audit.json").read_text())
    relative = "data/derived/ets2-program-membership.csv"
    if sha256(root / relative) != old_audit["outputs_sha256"][relative]:
        raise ValueError("Historical membership pin mismatch")
    source_rows = read_csv(root / "data/derived/ets2-source-gene-sets.csv")
    plan = json.loads((root / "config/ets2-program-plan.json").read_text())
    definitions = {item["id"]: item for item in plan["programs"]}
    result = {}
    for program in original.ORIGINAL_PROGRAM_IDS:
        copies = []
        for dataset in ("sc", "sn"):
            selected = sorted((r for r in rows if r["dataset_key"] == dataset
                               and r["program_id"] == program),
                              key=lambda r: int(r["source_index"]))
            if [int(r["source_index"]) for r in selected] != list(range(len(selected))):
                raise ValueError("Historical source indices not complete")
            copies.append([r["source_value"] for r in selected])
        if not copies[0] or copies[0] != copies[1] or len(copies[0]) != len(set(copies[0])):
            raise ValueError(f"Historical source membership conflict: {program}")
        result[program] = copies[0]
        if program not in COMPARATORS:
            source_set = {r["gene_id_original"] for r in source_rows
                          if r["set_name"] == definitions[program]["column"]}
            if set(result[program]) != source_set:
                raise ValueError(f"Original Ensembl list mismatch: {program}")
    if len(set(result["ets2_g1_dn"]) - {original.ETS2_STABLE_ID}) != 927:
        raise ValueError("Original primary membership count changed")
    return result


def reference_features(universe: Sequence[Mapping], annotation_lines: Sequence[str]):
    """Exact versioned joins; preserve PAR_Y and reject stable-ID collisions."""
    import pandas as pd
    annotation = {}
    for line in annotation_lines:
        fields = line.split()
        if len(fields) != 6 or fields[4] in annotation:
            raise ValueError("Reference rows not unique six-field records")
        annotation[fields[4]] = fields[5]
    ids = [r["Geneid"] for r in universe]
    if len(ids) != len(set(ids)) or not set(ids) <= annotation.keys():
        raise ValueError("Released feature IDs duplicate or not in reference")
    stable = [re.sub(r"\.\d+$", "", value) for value in ids]
    if len(stable) != len(set(stable)):
        raise ValueError("Collision after numeric-version removal")
    # original matcher recognizes ENSG[digits][.digits] only, never PAR_Y as
    # another copy of the X gene. All 45 absent PAR_Y rows stay absent.
    var = pd.DataFrame({"feature_name": [annotation[value] for value in ids]}, index=ids)
    features = original.build_feature_map(var, len(ids))
    return features, sorted(set(annotation) - set(ids))


def validate_samples(metadata: Sequence[Mapping[str, str]], header: Sequence[str]) -> list[dict]:
    if len(header) != len(set(header)) or not header:
        raise ValueError("Count sample header missing or duplicated")
    required = {"SampleID", "Stimulus_Hours", "RunID", "Library_prep", "HipsciID"}
    expected_conditions = {f"{s}_{t}" for s in REGISTRY_STIMULI + ("Ctrl",) for t in TIMES}
    expected_conditions |= {"Prec_D0", "Prec_D2"}
    rows = []
    for item in metadata:
        if not required <= item.keys() or any(not item[key] for key in required):
            raise ValueError("Missing required sample metadata")
        condition = item["Stimulus_Hours"]
        suffix = "_" + condition
        if condition not in expected_conditions or not item["SampleID"].endswith(suffix):
            raise ValueError("SampleID/condition suffix mismatch")
        if item["Library_prep"] not in PROTOCOLS:
            raise ValueError("Unknown protocol")
        line = item["SampleID"][:-len(suffix)]
        if not line:
            raise ValueError("Empty source line")
        rows.append({**item, "sample_id": item["SampleID"] + "_" + item["RunID"],
                     "line": line, "run": item["RunID"], "protocol": item["Library_prep"],
                     "condition": condition, "identity_mapped": item["HipsciID"] != "NOMATCH"})
    by_sample = {r["sample_id"]: r for r in rows}
    if len(by_sample) != len(rows) or set(by_sample) != set(header):
        raise ValueError("Count/metadata samples differ or duplicate")
    keys = [(r["line"], r["condition"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Multiple released cultures for one line/condition; no automatic collapse")
    groups = defaultdict(list)
    run_protocols = defaultdict(set)
    hipsci_lines = defaultdict(set)
    for row in rows:
        groups[row["line"]].append(row)
        run_protocols[row["run"]].add(row["protocol"])
        if row["identity_mapped"]:
            hipsci_lines[row["HipsciID"]].add(row["line"])
    # All metadata other than sample/condition is declared static by the source.
    static_fields = set(metadata[0]) - {"SampleID", "Stimulus_Hours"}
    for line, group in groups.items():
        for field in static_fields:
            if len({r[field] for r in group}) != 1:
                raise ValueError(f"Within-line metadata conflict: {line}, {field}")
    if any(len(values) != 1 for values in run_protocols.values()):
        raise ValueError("A sequencing run crosses protocols")
    if any(len(values) != 1 for values in hipsci_lines.values()):
        raise ValueError("One mapped HipSci identity has multiple source lines")
    return [{**by_sample[sample], "sample_index": i} for i, sample in enumerate(header)]


def pair_registry(samples: Sequence[Mapping], *, include_nomatch: bool = False) -> list[dict]:
    lookup = {(r["line"], r["condition"]): r for r in samples}
    if len(lookup) != len(samples):
        raise ValueError("Duplicate line/condition")
    lines = {r["line"]: r for r in samples}
    pairs = []
    for line, info in sorted(lines.items()):
        for stimulus in REGISTRY_STIMULI:
            for time in TIMES:
                treated = lookup.get((line, f"{stimulus}_{time}"))
                control = lookup.get((line, f"Ctrl_{time}"))
                if not info["identity_mapped"] and not include_nomatch:
                    status = "excluded_NOMATCH_identity"
                elif stimulus == "PIC":
                    status = "unavailable_mock_identity"
                elif treated is None and control is None:
                    status = "missing_stimulus_and_control"
                elif treated is None:
                    status = "missing_stimulus"
                elif control is None:
                    status = "missing_time_matched_control"
                else:
                    if any(treated[key] != control[key] for key in ("line", "run", "protocol")):
                        raise ValueError("Unpaired run/protocol/line")
                    status = "available"
                pairs.append({"line": line, "run": info["run"], "protocol": info["protocol"],
                              "identity_mapped": info["identity_mapped"], "stimulus": stimulus,
                              "time": time, "status": status,
                              "treatment_id": treated["sample_id"] if treated else "",
                              "control_id": control["sample_id"] if control else "",
                              "treatment_index": treated["sample_index"] if treated else "",
                              "control_index": control["sample_index"] if control else ""})
    return pairs


def build_folds(samples: Sequence[Mapping], *, folds: int = N_FOLDS,
                seed: int = FOLD_SEED) -> dict[tuple[str, str], int]:
    if folds < 2:
        raise ValueError("At least two folds required")
    line_groups = defaultdict(set)
    runs = defaultdict(set)
    for row in samples:
        line_groups[row["line"]].add((row["protocol"], row["run"]))
        runs[row["protocol"]].add(row["run"])
    if any(len(value) != 1 for value in line_groups.values()):
        raise ValueError("A line crosses run/protocol groups")
    assignment = {}
    for protocol, group in sorted(runs.items()):
        if len(group) < folds:
            raise ValueError("Too few runs for fixed folds; do not reduce K")
        ordered = sorted(group, key=lambda run: (hashlib.sha256(
            f"{seed}:{protocol}:{run}".encode()).hexdigest(), run))
        assignment.update({(protocol, run): i % folds for i, run in enumerate(ordered)})
    return assignment


def prediction_cohort(pairs: Sequence[Mapping]) -> list[dict]:
    available = [dict(r) for r in pairs if r["status"] == "available"]
    times = defaultdict(set)
    for row in available:
        times[row["line"]].add(row["time"])
    keep = {line for line, values in times.items() if values == set(TIMES)}
    return [r for r in available if r["line"] in keep]


def split_pairs(pairs: Sequence[Mapping], assignment: Mapping, protocol: str, time: int,
                heldout_fold: int) -> tuple[list[dict], list[dict], dict]:
    selected = [dict(r) for r in pairs if r["protocol"] == protocol and r["time"] == time
                and r["status"] == "available"]
    train = [r for r in selected if assignment[(r["protocol"], r["run"])] != heldout_fold]
    test = [r for r in selected if assignment[(r["protocol"], r["run"])] == heldout_fold]
    if {r["line"] for r in train} & {r["line"] for r in test}:
        raise ValueError("Training/test line leakage")
    by_class = {s: len({r["line"] for r in train if r["stimulus"] == s}) for s in STIMULI}
    by_class_runs = {s: len({r["run"] for r in train if r["stimulus"] == s}) for s in STIMULI}
    missing_test = [s for s in STIMULI if not any(r["stimulus"] == s for r in test)]
    reason = "available"
    if not test:
        reason = "empty_test_fold"
    elif min(by_class.values(), default=0) < MIN_LINES:
        reason = "insufficient_training_class_lines"
    elif min(by_class_runs.values(), default=0) < MIN_RUNS:
        reason = "insufficient_training_class_runs"
    return train, test, {"status": reason, "train_pairs": len(train), "test_pairs": len(test),
                         "train_lines": len({r["line"] for r in train}),
                         "test_lines": len({r["line"] for r in test}),
                         "train_runs": len({r["run"] for r in train}),
                         "test_runs": len({r["run"] for r in test}),
                         "train_class_lines": by_class, "train_class_runs": by_class_runs,
                         "missing_test_stimuli": missing_test}


def line_weights(lines: Sequence[str]) -> np.ndarray:
    """Each line sums to 1/N, irrespective of its available condition count."""
    counts = Counter(lines)
    if not counts:
        raise ValueError("Empty line set")
    return np.array([1 / (len(counts) * counts[line]) for line in lines], dtype=float)


def training_control_reference(log_values: np.ndarray, pairs: Sequence[Mapping],
                               assignment: Mapping, protocol: str, time: int,
                               heldout_fold: int) -> tuple[np.ndarray, dict]:
    """Features x samples input; train controls once per line, no test controls."""
    selected = [r for r in pairs if r["status"] == "available" and r["protocol"] == protocol
                and r["time"] == time and assignment[(protocol, r["run"])] != heldout_fold]
    by_line = defaultdict(set)
    runs = set()
    for row in selected:
        by_line[row["line"]].add(int(row["control_index"]))
        runs.add(row["run"])
    if len(by_line) < MIN_LINES or len(runs) < MIN_RUNS:
        raise ValueError("Insufficient training reference lines/runs")
    if any(len(indices) != 1 for indices in by_line.values()):
        raise ValueError("Multiple control cultures require an explicit procedure")
    indices = [next(iter(by_line[line])) for line in sorted(by_line)]
    if len(indices) != len(set(indices)):
        raise ValueError("A control sample is assigned to different lines")
    if min(indices) < 0 or max(indices) >= log_values.shape[1]:
        raise ValueError("Control index outside matrix")
    # Slice only training controls. Neither held-out controls nor any treatment
    # contributes even to a numerical validity check of the matching reference.
    values = np.asarray(log_values[:, indices], dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("Invalid training control expression")
    return values.mean(axis=1), {"reference_lines": len(by_line), "reference_runs": len(runs),
                                "control_indices": indices}


def matched_weights(reference: np.ndarray, feature_ids: Sequence[str],
                    programs: Mapping[str, Sequence[int]], excluded: Sequence[int],
                    protocol: str, time: int, fold: int, *, replicates: int = MATCH_REPLICATES,
                    bins: int = MATCH_BINS) -> dict[str, dict]:
    """Same random-set mean as old scoring, compressed to gene weights.

    Weighting is an algebraic reduction, not a new score or random null. It
    avoids the samples x 1000 x target-genes array. Selection hashes preserve
    each realized set. Mapping gates must be checked by the caller first.
    """
    if bins < 1 or replicates < 1:
        raise ValueError("Matching bins and replicate count must be positive")
    excluded_set = set(map(int, excluded))
    if any(not set(map(int, indices)) <= excluded_set for indices in programs.values()):
        raise ValueError("All target program genes must be excluded from background pools")
    _, bin_arrays = original.expression_bins(reference, feature_ids, bins=bins)
    result = {}
    for program, indices in programs.items():
        seedtext = f"{MATCH_SEED}:macromap:protocol={protocol}:time={time}:fold={fold}:{program}"
        seed = int(hashlib.sha256(seedtext.encode()).hexdigest()[:16], 16)
        match = original.sample_expression_matched_sets(
            bin_arrays, indices, excluded, replicates=replicates, seed=seed)
        info = {"available": match.available, "reason": match.reason, "seed": seed,
                "target_bin_counts": match.target_bin_counts,
                "pool_bin_counts": match.pool_bin_counts,
                "underflow_bins": list(match.underflow_bins)}
        if match.available:
            weights = np.bincount(match.selections.ravel(), minlength=len(reference)).astype(float)
            weights /= match.selections.size
            target = np.zeros(len(reference))
            target[match.target_indices] = 1 / len(match.target_indices)
            info.update(target_weights=target, background_weights=weights,
                        score_weights=target - weights,
                        selection_sha256=hashlib.sha256(match.selections.astype("<i8").tobytes()).hexdigest(),
                        baseline_target=float(reference @ target),
                        baseline_background=float(reference @ weights))
        result[program] = info
    return result


def score_with_weights(log_values: np.ndarray, weights: Mapping, sample_indices: Sequence[int]):
    if not weights["available"]:
        return {"status": weights["reason"], "raw": None, "background": None, "adjusted": None}
    values = np.asarray(log_values[:, sample_indices], dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite scoring values")
    raw = weights["target_weights"] @ values
    background = weights["background_weights"] @ values
    return {"status": "available", "raw": raw, "background": background,
            "adjusted": raw - background}


@dataclass(frozen=True)
class Scaler:
    mean: np.ndarray
    scale: np.ndarray
    constant: np.ndarray

    def transform(self, x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        if values.ndim != 2 or values.shape[1] != len(self.mean) or not np.isfinite(values).all():
            raise ValueError("Invalid held-out predictor matrix")
        return (values - self.mean) / self.scale


def fit_scaler(training_x: np.ndarray, training_lines: Sequence[str]) -> Scaler:
    """Fit only to training-run paired changes, with equal line weights."""
    x = np.asarray(training_x, dtype=float)
    if x.ndim != 2 or len(x) != len(training_lines) or not np.isfinite(x).all():
        raise ValueError("Invalid training predictors")
    weights = line_weights(training_lines)
    mean = weights @ x
    variance = weights @ np.square(x - mean)
    sd = np.sqrt(variance)
    constant = sd <= 1e-12
    # Constant training features remain at training mean with a unit scale.
    # The ridge optimum has zero slopes there; no held-out-based fallback.
    return Scaler(mean, np.where(constant, 1.0, sd), constant)


def ridge_objective(theta: np.ndarray, x: np.ndarray, y: np.ndarray, weights: np.ndarray,
                    n_classes: int, penalty: float = RIDGE_LAMBDA) -> tuple[float, np.ndarray]:
    """Weighted mean multinomial loss + lambda/2 * ||all class slopes||^2.

    Use a symmetric K-class softmax, not an arbitrarily unpenalized reference
    slope. Intercepts are unpenalized. The intercept shift nonidentifiability
    is fixed by zero initialization and explicitly recentering the solution.
    """
    parameters = np.asarray(theta, dtype=float).reshape(x.shape[1] + 1, n_classes)
    intercept = parameters[0]
    beta = parameters[1:]
    logits = x @ beta + intercept
    logp = logits - logsumexp(logits, axis=1, keepdims=True)
    loss = -np.sum(weights * logp[np.arange(len(y)), y]) + penalty / 2 * np.square(beta).sum()
    residual = np.exp(logp)
    residual[np.arange(len(y)), y] -= 1
    residual *= weights[:, None]
    gradient = np.vstack([residual.sum(axis=0), x.T @ residual + penalty * beta])
    return float(loss), gradient.ravel()


@dataclass
class RidgeFit:
    status: str
    coefficients: np.ndarray | None
    n_classes: int
    iterations: int = 0
    gradient_inf: float | None = None
    objective: float | None = None
    optimizer_message: str = ""

    def log_proba(self, x: np.ndarray) -> np.ndarray:
        if self.status != "available" or self.coefficients is None:
            raise ValueError(f"Prediction unavailable: {self.status}")
        x = np.asarray(x, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.coefficients.shape[0] - 1 or not np.isfinite(x).all():
            raise ValueError("Invalid prediction features")
        with np.errstate(over="ignore", invalid="ignore"):
            logits = x @ self.coefficients[1:] + self.coefficients[0]
            result = logits - logsumexp(logits, axis=1, keepdims=True)
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite held-out log probabilities; no clipping fallback")
        return result


def fit_multinomial_ridge(x: np.ndarray, y: Sequence[int], lines: Sequence[str], *,
                          n_classes: int = 9, penalty: float = RIDGE_LAMBDA,
                          max_iter: int = MAX_ITER) -> RidgeFit:
    x = np.asarray(x, dtype=float)
    yraw = np.asarray(y)
    if (x.ndim != 2 or len(x) == 0 or len(x) != len(lines) or yraw.shape != (len(x),)
            or not np.isfinite(x).all() or not np.issubdtype(yraw.dtype, np.integer)
            or penalty <= 0 or n_classes < 2):
        return RidgeFit("invalid_training_data", None, n_classes)
    y = yraw.astype(int)
    if set(y) != set(range(n_classes)):
        return RidgeFit("missing_training_stimulus", None, n_classes)
    weights = line_weights(lines)
    initial = np.zeros((x.shape[1] + 1) * n_classes)
    result = minimize(ridge_objective, initial, args=(x, y, weights, n_classes, penalty),
                      jac=True, method="L-BFGS-B",
                      options={"maxiter": max_iter, "ftol": 1e-13, "gtol": 1e-8, "maxls": 50})
    value, gradient = ridge_objective(result.x, x, y, weights, n_classes, penalty)
    norm = float(np.max(np.abs(gradient)))
    ok = bool(result.success and np.isfinite(value) and np.isfinite(result.x).all()
              and np.isfinite(gradient).all() and norm <= GRADIENT_TOL)
    coeff = result.x.reshape(x.shape[1] + 1, n_classes).copy() if ok else None
    if coeff is not None:
        coeff -= coeff.mean(axis=1, keepdims=True)
    return RidgeFit("available" if ok else "optimizer_failed", coeff, n_classes,
                    int(result.nit), norm, value, str(result.message))


def heldout_increment(train_x: np.ndarray, test_x: np.ndarray, train_pairs: Sequence[Mapping],
                      test_pairs: Sequence[Mapping], *, max_iter: int = MAX_ITER) -> dict:
    """Four fixed comparator columns then ONE prespecified ETS2 column.

    This pure helper never reads files. Callers must supply one protocol/time
    and gate the mapping, matching and fold design before fitting. No missing
    rows or failed models are dropped to rescue the endpoint.
    """
    if not train_pairs or not test_pairs:
        return {"status": "empty_fold", "predictions": []}
    training_lines = [r["line"] for r in train_pairs]
    testing_lines = [r["line"] for r in test_pairs]
    train_groups = {(r["protocol"], r["run"]) for r in train_pairs}
    test_groups = {(r["protocol"], r["run"]) for r in test_pairs}
    strata = {(r["protocol"], r["time"]) for r in list(train_pairs) + list(test_pairs)}
    if set(training_lines) & set(testing_lines) or train_groups & test_groups or len(strata) != 1:
        raise ValueError("Leakage or mixed protocol/time")
    train_x = np.asarray(train_x, dtype=float)
    test_x = np.asarray(test_x, dtype=float)
    if train_x.shape != (len(train_pairs), 5) or test_x.shape != (len(test_pairs), 5):
        raise ValueError("Exactly four comparator and one ETS2 columns required")
    if not np.isfinite(train_x).all() or not np.isfinite(test_x).all():
        return {"status": "missing_predictor", "predictions": []}
    if any(r["stimulus"] not in STIMULI for r in list(train_pairs) + list(test_pairs)):
        return {"status": "unknown_stimulus", "predictions": []}
    y = [STIMULI.index(r["stimulus"]) for r in train_pairs]
    test_y = np.array([STIMULI.index(r["stimulus"]) for r in test_pairs])
    scaler = fit_scaler(train_x, training_lines)
    tr = scaler.transform(train_x)
    te = scaler.transform(test_x)
    baseline = fit_multinomial_ridge(tr[:, :4], y, training_lines, max_iter=max_iter)
    extended = fit_multinomial_ridge(tr, y, training_lines, max_iter=max_iter)
    fits = {"baseline": baseline, "extended": extended, "scaler": scaler}
    if baseline.status != "available" or extended.status != "available":
        return {"status": "failed_paired_fit", "predictions": [], **fits}
    try:
        logs = [fit.log_proba(values) for fit, values in [(baseline, te[:, :4]), (extended, te)]]
    except ValueError:
        return {"status": "failed_prediction", "predictions": [], **fits}
    rows = []
    for i, pair in enumerate(test_pairs):
        lb, le = (-float(logp[i, test_y[i]]) for logp in logs)
        rows.append({**pair, "baseline_loss": lb, "extended_loss": le, "gain": lb - le,
                     # np.argmax picks the first in the frozen STIMULI order on exact ties.
                     "baseline_class": STIMULI[int(np.argmax(logs[0][i]))],
                     "extended_class": STIMULI[int(np.argmax(logs[1][i]))],
                     "baseline_log_probabilities": logs[0][i].tolist(),
                     "extended_log_probabilities": logs[1][i].tolist()})
    return {"status": "available", "predictions": rows, **fits}


def line_time_losses(predictions: Sequence[Mapping]) -> list[dict]:
    groups = defaultdict(list)
    seen = set()
    line_runs = defaultdict(set)
    for row in predictions:
        key = (row["line"], row["time"], row["stimulus"])
        if key in seen or row["stimulus"] not in STIMULI or row["time"] not in TIMES:
            raise ValueError("Duplicated/unknown held-out label")
        seen.add(key)
        if any(not math.isfinite(row[field]) for field in ("baseline_loss", "extended_loss")):
            raise ValueError("Missing/failed prediction cannot be dropped")
        if any(row[field] < 0 for field in ("baseline_loss", "extended_loss")):
            raise ValueError("Negative log loss")
        line_runs[row["line"]].add((row["protocol"], row["run"]))
        groups[(row["protocol"], row["run"], row["line"], row["time"])].append(row)
    if any(len(values) != 1 for values in line_runs.values()):
        raise ValueError("A predicted line crosses protocol/run")
    result = []
    for (protocol, run, line, time), rows in sorted(groups.items()):
        base = float(np.mean([r["baseline_loss"] for r in rows]))
        extended = float(np.mean([r["extended_loss"] for r in rows]))
        result.append({"protocol": protocol, "run": run, "line": line, "time": time,
                       "stimuli_observed": len(rows), "baseline_loss": base,
                       "extended_loss": extended, "gain": base - extended})
    return result


def aggregate_prediction(predictions: Sequence[Mapping], *, expected_pairs: Sequence[Mapping],
                         multiplicity: Mapping[tuple[str, str], int] | None = None,
                         protocol_weights: Mapping[str, float] | None = None) -> dict:
    """Line-weighted risk, 6h/24h equal once, then line-share protocol pooling.

    expected_pairs prevents silent dropping of failed folds or missing feature
    rows. Missing study measurements belong in that pre-frozen registry; an
    absent expected prediction is an unavailable endpoint, never a zero loss.
    """
    key = lambda r: (r["protocol"], r["run"], r["line"], r["time"], r["stimulus"])
    if Counter(map(key, predictions)) != Counter(map(key, expected_pairs)):
        return {"status": "incomplete_predictions", "primary_gain": None}
    losses = line_time_losses(predictions)
    by_protocol = defaultdict(list)
    for row in losses:
        by_protocol[row["protocol"]].append(row)
    if set(by_protocol) != set(PROTOCOLS):
        return {"status": "missing_protocol", "primary_gain": None}
    cells = {}
    line_counts = {}
    for protocol, group in by_protocol.items():
        lines_by_time = [{r["line"] for r in group if r["time"] == time} for time in TIMES]
        if not lines_by_time[0] or lines_by_time[0] != lines_by_time[1]:
            return {"status": "unpaired_times", "primary_gain": None}
        line_counts[protocol] = len(lines_by_time[0])
        for time in TIMES:
            observed = {r["stimulus"] for r in predictions if r["protocol"] == protocol and r["time"] == time}
            if observed != set(STIMULI):
                return {"status": "missing_evaluation_stimulus", "primary_gain": None}
            subset = [r for r in group if r["time"] == time]
            w = np.array([multiplicity.get((protocol, r["run"]), 0) if multiplicity is not None else 1
                          for r in subset], dtype=float)
            if not w.sum() > 0:
                return {"status": "empty_resampled_protocol", "primary_gain": None}
            w /= w.sum()
            cells[(protocol, time)] = {field: float(w @ np.array([r[field] for r in subset]))
                                       for field in ("baseline_loss", "extended_loss", "gain")}
    pw = dict(protocol_weights) if protocol_weights is not None else {
        p: line_counts[p] / sum(line_counts.values()) for p in PROTOCOLS}
    if set(pw) != set(PROTOCOLS) or any(v <= 0 for v in pw.values()) or not np.isclose(sum(pw.values()), 1):
        raise ValueError("Invalid protocol target weights")
    pooled = {field: sum(pw[p] * sum(cells[(p, t)][field] for t in TIMES) / 2 for p in PROTOCOLS)
              for field in ("baseline_loss", "extended_loss", "gain")}
    return {"status": "available", "primary_gain": pooled["gain"], "pooled": pooled,
            "cells": cells, "line_counts": line_counts, "protocol_weights": pw}


def cluster_bootstrap_gain(predictions: Sequence[Mapping], *, expected_pairs: Sequence[Mapping],
                           replicates: int = 5000, seed: int = BOOTSTRAP_SEED) -> dict:
    """Conditional uncertainty of fixed OOF predictions; no model refitting.

    Runs are sampled within protocol. All lines/times/conditions in a run move
    together. These intervals do NOT include variability in learned references
    or overlapping training fits and are NOT unconditional algorithm CIs.
    """
    point = aggregate_prediction(predictions, expected_pairs=expected_pairs)
    if point["status"] != "available":
        return point
    runs = {p: sorted({r["run"] for r in predictions if r["protocol"] == p}) for p in PROTOCOLS}
    if min(map(len, runs.values())) < MIN_RUNS:
        return {**point, "interval_status": "insufficient_runs", "interval": None}
    if replicates < 1:
        raise ValueError("Positive bootstrap replicates required")
    # Pre-aggregate the line losses once. This is equivalent to copying all
    # rows for each sampled run, but does not revisit thousands of rows 5000x.
    losses = line_time_losses(predictions)
    sums = {p: np.array([[sum(r["gain"] for r in losses if r["protocol"] == p
                             and r["run"] == run and r["time"] == time) for time in TIMES]
                         for run in runs[p]]) for p in PROTOCOLS}
    counts = {p: np.array([len({r["line"] for r in losses if r["protocol"] == p and r["run"] == run})
                          for run in runs[p]]) for p in PROTOCOLS}
    rng = np.random.Generator(np.random.PCG64(seed))
    draws = np.empty(replicates)
    cell_draws = {key: np.empty(replicates) for key in point["cells"]}
    for b in range(replicates):
        total = 0.0
        for p in PROTOCOLS:
            multiplicity = np.bincount(rng.integers(len(runs[p]), size=len(runs[p])), minlength=len(runs[p]))
            by_time = (multiplicity @ sums[p]) / (multiplicity @ counts[p])
            for j, t in enumerate(TIMES):
                cell_draws[(p, t)][b] = by_time[j]
            total += point["protocol_weights"][p] * float(by_time.mean())
        draws[b] = total
    interval = lambda values: np.quantile(values, [0.025, 0.975], method="linear").tolist()
    return {**point, "interval_status": "fixed_prediction_conditional_only", "interval": interval(draws),
            "cell_intervals": {key: interval(values) for key, values in cell_draws.items()},
            "replicates": replicates, "seed": seed}


def clustered_response_summary(values: Sequence[float], lines: Sequence[str], runs: Sequence[str], *,
                               replicates: int = 5000, seed: int = BOOTSTRAP_SEED,
                               compute_interval: bool = True) -> dict:
    """One protocol/time/program cell, one paired response per line.

    Point summaries retain the earlier median/range/sign and leave-line-out
    conventions. Run bootstrap intervals add conditional, pointwise context;
    the same explicit seed/run order lets callers reuse cluster draws across
    the panel. No genes or treatment rows become independent replicates.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) != len(lines) or len(x) != len(runs) or len(set(lines)) != len(lines):
        raise ValueError("Response cell must contain one paired value per line")
    if not np.isfinite(x).all():
        return {"status": "missing_response", "median": None}
    if len(x) < MIN_LINES:
        return {"status": "insufficient_paired_lines", "median": None, "lines": len(x)}
    unique_runs = sorted(set(runs))
    result = {"status": "available_descriptive", "lines": len(x), "runs": len(unique_runs),
              "median": float(np.median(x)), "mean": float(np.mean(x)),
              "minimum": float(np.min(x)), "maximum": float(np.max(x)),
              "positive": int(np.sum(x > 0)), "negative": int(np.sum(x < 0)), "zero": int(np.sum(x == 0))}
    loo = [float(np.median(np.delete(x, i))) for i in range(len(x))] if len(x) - 1 >= MIN_LINES else []
    result["leave_one_line_out_median_range"] = [min(loo), max(loo)] if loo else None
    lor = [float(np.median(x[np.array(runs) != run])) for run in unique_runs
           if int(np.sum(np.array(runs) != run)) >= MIN_LINES]
    result["leave_one_run_out_median_range"] = [min(lor), max(lor)] if lor else None
    if not compute_interval:
        return {**result, "interval_status": "descriptive_only_no_pooled_interval",
                "median_interval": None, "mean_interval": None}
    if len(unique_runs) < MIN_RUNS:
        return {**result, "interval_status": "insufficient_runs", "median_interval": None, "mean_interval": None}
    if replicates < 1:
        raise ValueError("Positive bootstrap replicates required")
    members = [np.flatnonzero(np.array(runs) == run) for run in unique_runs]
    rng = np.random.Generator(np.random.PCG64(seed))
    medians, means = np.empty(replicates), np.empty(replicates)
    for b in range(replicates):
        selected = rng.integers(len(unique_runs), size=len(unique_runs))
        sample = x[np.concatenate([members[i] for i in selected])]
        medians[b], means[b] = np.median(sample), np.mean(sample)
    interval = lambda v: np.quantile(v, [0.025, 0.975], method="linear").tolist()
    return {**result, "interval_status": "pointwise_fixed_score_conditional_only",
            "median_interval": interval(medians), "mean_interval": interval(means),
            "replicates": replicates, "seed": seed}


# Fast syntax proves nonnegative integrality BEFORE binary-float parsing.
# Rare exponent/signed forms use exact Decimal below, never a rounded float.
_INTEGRAL_TOKEN = r"[+]?[0-9]+(?:\.0*)?"
_INTEGRAL_ROW = re.compile(_INTEGRAL_TOKEN + r"(?:\t" + _INTEGRAL_TOKEN + r")*")
_DECIMAL_TOKEN = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def exact_count_row(text: str, expected_columns: int) -> np.ndarray:
    """Require exact finite nonnegative integer source quantities below 2^53.

    Float-first validation would accept 9.9999999999999999 as 10 and -1e-400
    as zero. A lexical integer fast path handles usual count exports without
    millions of Decimal objects. Other allowed number syntax is checked with
    exact Decimal before conversion; no rounded/underflowed token is accepted.
    """
    if expected_columns < 1 or text.count("\t") + 1 != expected_columns:
        raise ValueError("Wrong numerical field count")
    if _INTEGRAL_ROW.fullmatch(text):
        values = np.fromstring(text, dtype=np.float64, sep="\t")
        if (len(values) != expected_columns or not np.isfinite(values).all()
                or np.any(values >= 2**53)):
            raise ValueError("Raw integer count outside exact float64 range")
        return values
    integers = []
    for source_token in text.split("\t"):
        token = source_token.strip()
        if not _DECIMAL_TOKEN.fullmatch(token):
            raise ValueError("Invalid raw-count lexical token")
        try:
            value = Decimal(token)
            if (not value.is_finite() or value < 0 or value >= 2**53
                    or value != value.to_integral_value()):
                raise ValueError("Source count is not an exact nonnegative integer below 2^53")
            integers.append(int(value))
        except InvalidOperation as exc:
            raise ValueError("Invalid exact raw-count quantity") from exc
    return np.asarray(integers, dtype=np.float64)


def normalize_count_cache(source: Path, output: Path, feature_ids: Sequence[str],
                          sample_ids: Sequence[str], *, chunk_rows: int = 256) -> dict:
    """Numerical helper for AFTER execution freeze; never called by prepare().

    One gzip pass writes float64 gene x sample counts and sums every released
    gene. A bounded second memmap pass replaces counts with log2(CPM+1).
    Files are created exclusively through a new file handle, not a pathname
    that might be a dangling symlink. Finished cache/sidecar files are published
    without replacing an existing name. Source counts remain in the raw gzip.
    """
    if chunk_rows < 1 or not feature_ids or not sample_ids:
        raise ValueError("Invalid cache dimensions/chunks")
    if len(feature_ids) != len(set(feature_ids)) or len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Duplicate cache identities")
    source, output = _absolute(source), _absolute(output)
    reject_path_aliases(source)
    source_before = source.stat()
    source_signature = lambda stat: (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    partial = output.with_name(output.name + ".partial")
    sidecar = output.with_name(output.name + ".json")
    partial_sidecar = output.with_name(output.name + ".json.partial")
    reservation = output.with_name(output.name + ".lock")
    paths = (output, partial, sidecar, partial_sidecar, reservation)
    for path in paths:
        reject_path_aliases(path)
        if path.exists():
            raise FileExistsError(f"Do not overwrite cache artifacts: {path}")
    output.parent.mkdir(parents=True, exist_ok=True)
    reject_path_aliases(output.parent)
    totals = np.zeros(len(sample_ids), dtype=np.float64)
    nrows = 0
    created = set()
    matrix = None
    reservation_handle = None
    try:
        # Reserve the complete namespace before checking it a second time.
        # Keep this exclusive reservation through cache AND sidecar publication;
        # on success it remains as an immutable no-reuse marker. A crash leaves
        # a visible stale reservation, never permission to overwrite a cache.
        reservation_handle = reservation.open("x", encoding="utf-8")
        created.add(reservation)
        reservation_handle.write("Reserved immutable numerical-cache namespace.\n")
        reservation_handle.flush()
        for path in (output, partial, sidecar, partial_sidecar):
            reject_path_aliases(path)
            if path.exists():
                raise FileExistsError(f"Cache name occupied after reservation: {path}")
        # x+b atomically rejects every existing leaf, including a dangling
        # symlink. Mapping the open handle avoids a second pathname open.
        with partial.open("x+b") as cache_file:
            created.add(partial)
            cache_file.truncate(len(feature_ids) * len(sample_ids) * 8)
            matrix = np.memmap(cache_file, dtype="<f8", mode="r+", shape=(len(feature_ids), len(sample_ids)))
            with gzip.open(source, "rt") as fh:
                header = fh.readline().rstrip("\r\n").split("\t")
                if header != list(FEATURE_FIELDS) + list(sample_ids):
                    raise ValueError("Count header mismatch")
                for i, line in enumerate(fh):
                    fields = line.rstrip("\r\n").split("\t", 6)
                    if i >= len(feature_ids) or len(fields) != 7 or fields[1] != feature_ids[i]:
                        raise ValueError("Count feature order/row mismatch")
                    values = exact_count_row(fields[6], len(sample_ids))
                    matrix[i] = values
                    totals += values
                    nrows += 1
            if nrows != len(feature_ids) or np.any(totals <= 0) or np.any(totals >= 2**53):
                raise ValueError("Missing rows or invalid released-universe totals")
            for start in range(0, len(feature_ids), chunk_rows):
                stop = min(start + chunk_rows, len(feature_ids))
                values = np.asarray(matrix[start:stop]).copy()
                matrix[start:stop] = np.log2(1 + values * (1e6 / totals)[None, :])
            matrix.flush()
            del matrix
            matrix = None
        source_hash = sha256(source)
        if source_signature(source.stat()) != source_signature(source_before):
            raise ValueError("Raw source changed during normalization; cache not published")
        publish_without_replacement(partial, output)
        created.discard(partial)
        created.add(output)
        audit = {"complete": True, "quantity": "log2(released-universe CPM + 1)",
                 "shape": [len(feature_ids), len(sample_ids)], "dtype": "<f8", "source_sha256": source_hash,
                 "cache_sha256": sha256(output), "all_library_totals_positive": True,
                 "feature_id_sha256": hashlib.sha256("\n".join(feature_ids).encode()).hexdigest(),
                 "sample_id_sha256": hashlib.sha256("\n".join(sample_ids).encode()).hexdigest()}
        # Track only paths created by this invocation. Never unlink a prior
        # user artifact when rejecting an occupied destination.
        reject_path_aliases(partial_sidecar)
        with partial_sidecar.open("x", encoding="utf-8") as fh:
            created.add(partial_sidecar)
            fh.write(json.dumps(audit, indent=2) + "\n")
        publish_without_replacement(partial_sidecar, sidecar)
        created.discard(partial_sidecar)
        created.add(sidecar)
        reservation_handle.close()
        reservation_handle = None
        return audit
    except BaseException:
        if matrix is not None:
            del matrix
        if reservation_handle is not None:
            reservation_handle.close()
        for path in created:
            path.unlink(missing_ok=True)
        raise


def _bool_field(value: str) -> bool:
    if value not in {"True", "False"}:
        raise ValueError("Expected an explicit True/False field")
    return value == "True"


def read_pairs(path: Path) -> list[dict]:
    """Restore CSV types; a literal string 'False' must never mean True."""
    rows = read_csv(path)
    for row in rows:
        row["time"] = int(row["time"])
        row["identity_mapped"] = _bool_field(row["identity_mapped"])
        for key in ("treatment_index", "control_index"):
            row[key] = int(row[key]) if row[key] != "" else ""
        if row["time"] not in TIMES or row["stimulus"] not in REGISTRY_STIMULI:
            raise ValueError("Unexpected prepared pair label")
        if row["status"] == "available" and any(row[key] == "" for key in ("treatment_index", "control_index")):
            raise ValueError("Prepared available pair lacks a sample index")
    return rows


def load_prepared_design(directory: Path, *, expected_summary_sha256: str) -> dict:
    """Read the root-pinned source-only design, including explicit CSV types.

    The execution freeze must supply the expected summary hash. This loader
    checks all preparation artifacts, not a mutable path chosen at inference.
    It does not authorize numerical inference or change ready_for_inference.
    """
    summary_path = directory / "source-only-readiness.json"
    if sha256(summary_path) != expected_summary_sha256:
        raise ValueError("Preparation summary does not match execution pin")
    summary = json.loads(summary_path.read_text())
    for filename, expected in summary["artifact_sha256"].items():
        if Path(filename).name != filename or sha256(directory / filename) != expected:
            raise ValueError(f"Preparation artifact changed: {filename}")
    samples = read_csv(directory / "sample-map.csv")
    for row in samples:
        row["identity_mapped"] = _bool_field(row["identity_mapped"])
        row["sample_index"] = int(row["sample_index"])
    pairs = read_pairs(directory / "pair-registry.csv")
    predictive = read_pairs(directory / "predictive-pairs.csv")
    fold_rows = read_csv(directory / "run-folds.csv")
    folds = {(row["protocol"], row["run"]): int(row["fold"]) for row in fold_rows}
    if len(folds) != len(fold_rows) or any(f not in range(N_FOLDS) for f in folds.values()):
        raise ValueError("Invalid prepared fold registry")
    inventory = read_csv(directory / "mapping-inventory.csv")
    gates = {row["program_id"]: _bool_field(row["mapping_gate_pass"]) for row in inventory}
    memberships = read_csv(directory / "program-mapping.csv")
    indices = {p: sorted({int(row["feature_index"]) for row in memberships if row["program_id"] == p
                          and row["mapping_status"] in {"mapped", "mapped_nonoverlap"}}) for p in PROGRAMS}
    features = read_csv(directory / "features.csv")
    if [int(row["feature_index"]) for row in features] != list(range(len(features))):
        raise ValueError("Prepared feature order changed")
    excluded = set().union(*(set(indices[p]) for p in original.ORIGINAL_PROGRAM_IDS))
    excluded |= {int(row["feature_index"]) for row in features
                 if row["stable_feature_id"] == original.ETS2_STABLE_ID or row["feature_name"] == original.ETS2_SYMBOL}
    return {"summary": summary, "samples": samples, "pairs": pairs, "predictive_pairs": predictive,
            "folds": folds, "mapping_gates": gates, "program_indices": indices,
            "feature_ids": [row["feature_id"] for row in features],
            "sample_ids": [row["sample_id"] for row in samples],
            "pool_exclusion_indices": sorted(excluded)}


def prepare(root: Path, output: Path) -> dict:
    """Source-identity, mapping, pairing and deterministic fold gates ONLY."""
    # Full sample joins/line identifiers remain in a fresh, nonaliased ignored
    # directory. Prior preparation snapshots are never overwritten.
    output = fresh_private_run(root, output)
    verified = verify_pins(root)
    metadata = read_csv(root / RAW / "worker-methods/sample-metadata.csv")
    xlsx_rows = xlsx_sheet(root / RAW / "worker-methods/publication_eqtl_supp11_lines.body", "RNA-seq")
    if metadata != xlsx_rows:
        raise ValueError("Cached sample metadata disagrees with original Supplementary Data 11")
    header = (root / RAW / "MacroMap_raw_expression.header.txt").read_text().strip().split("\t")
    if tuple(header[:6]) != FEATURE_FIELDS:
        raise ValueError("Wrong qualified header schema")
    # Read the source header only. No expression-valued row is read here.
    with gzip.open(root / RAW / "MacroMap_raw_expression.txt.gz", "rt") as fh:
        if fh.readline().strip().split("\t") != header:
            raise ValueError("Qualified header no longer matches raw source")
    samples = validate_samples(metadata, header[6:])
    old_join = {r["expression_sample_id"]: r for r in read_csv(root / RAW / "worker-methods/expression-header-sample-map.csv")}
    for row in samples:
        prior = old_join[row["sample_id"]]
        if any(prior[old] != row[new] for old, new in [("source_line_key", "line"), ("RunID", "run"),
                ("Library_prep", "protocol"), ("Stimulus_Hours", "condition"), ("HipsciID_as_published", "HipsciID")]):
            raise ValueError("New sample join disagrees with previous qualification")
    universe = read_csv(root / RAW / "raw-gene-universe.tsv", "\t")
    features, reference_only = reference_features(universe, (root / RAW / "worker-methods/github_gencode_v27_annotation.body").read_text().splitlines())
    if len(features) != 58243 or len(reference_only) != 45 or any(not i.endswith("_PAR_Y") for i in reference_only):
        raise ValueError("Qualified released feature universe changed")
    source_programs = retained_source_programs(root)
    plan = json.loads((root / "config/ets2-program-plan.json").read_text())
    memberships, inventory, mapped = original.build_program_mappings(plan, source_programs, {"macromap": features})
    for filename, frame in (("program-mapping.csv", memberships), ("mapping-inventory.csv", inventory), ("features.csv", features)):
        path = output / filename
        reject_path_aliases(path)
        with path.open("x", newline="") as fh:
            frame.to_csv(fh, index=False)
    pairs = pair_registry(samples)
    folds = build_folds(samples)
    predictive = prediction_cohort(pairs)
    write_csv(output / "sample-map.csv", samples)
    write_csv(output / "pair-registry.csv", pairs)
    write_csv(output / "predictive-pairs.csv", predictive)
    write_csv(output / "run-folds.csv", [{"protocol": p, "run": r, "fold": f} for (p, r), f in sorted(folds.items())])
    gates = []
    for p in PROTOCOLS:
        for t in TIMES:
            for f in range(N_FOLDS):
                _, _, gate = split_pairs(predictive, folds, p, t, f)
                gates.append({"protocol": p, "time": t, "fold": f, **gate})
    exclusive_text(output / "fold-gates.json", json.dumps(gates, indent=2) + "\n")
    pair_counts = [{"protocol": p, "time": t, "stimulus": s,
                    "paired_lines": len({r["line"] for r in pairs if r["protocol"] == p and r["time"] == t
                                         and r["stimulus"] == s and r["status"] == "available"}),
                    "status": "unavailable_mock_identity" if s == "PIC" else "source_qualified"}
                   for p in PROTOCOLS for t in TIMES for s in REGISTRY_STIMULI]
    write_csv(output / "paired-design-aggregate.csv", pair_counts)
    import platform
    import scipy
    summary = {"stage": "source_only_preparation_not_execution_freeze", "ready_for_inference": False,
               "preparation_script_sha256": sha256(Path(__file__)),
               "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
               "real_program_scores_computed": False, "real_effects_or_predictions_computed": False,
               "source_pins": verified,
               "artifact_sha256": {filename: sha256(output / filename) for filename in (
                   "program-mapping.csv", "mapping-inventory.csv", "features.csv", "sample-map.csv",
                   "pair-registry.csv", "predictive-pairs.csv", "run-folds.csv", "fold-gates.json",
                   "paired-design-aggregate.csv")},
               "source_programs_route": "complete_pinned_previous_source_value_audit",
               "original_ets2_archive_present": (root / "data/raw/ets2-benchmark/ets2-public-release.zip").is_file(),
               "released_genes": len(features), "reference_only_PAR_Y": len(reference_only), "samples": len(samples),
               "mapped_lines": len({r["line"] for r in samples if r["identity_mapped"]}),
               "unmapped_separate_lines": len({r["line"] for r in samples if not r["identity_mapped"]}),
               "protocol_design": {p: {"source_lines": len({r["line"] for r in samples if r["protocol"] == p}),
                                       "mapped_lines": len({r["line"] for r in samples if r["protocol"] == p and r["identity_mapped"]}),
                                       "runs": len({r["run"] for r in samples if r["protocol"] == p}),
                                       "predictive_lines_both_times": len({r["line"] for r in predictive if r["protocol"] == p})}
                                   for p in PROTOCOLS},
               "mapping_inventory": json.loads(inventory.to_json(orient="records")),
               "fold_gate_status": dict(Counter(g["status"] for g in gates)),
               "fold_class_absences": sum(bool(g["missing_test_stimuli"]) for g in gates),
               "matching_pools_and_numerical_fits": "not_examined_before_freeze",
               "proposal": {"folds": N_FOLDS, "fold_seed": FOLD_SEED, "ridge_lambda": RIDGE_LAMBDA,
                            "class_order": list(STIMULI), "time_order": list(TIMES),
                            "scaling": "training_run_paired_changes_only_equal_line_weight",
                            "reference": "training_run_time_matched_controls_only",
                            "uncertainty": "stratified_run_bootstrap_fixed_predictions_conditional_only"}}
    exclusive_text(output / "source-only-readiness.json", json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def execution_contract() -> dict:
    """Exact finite method contract; parent may freeze, not tune, this template."""
    return {"contract_version": 1, "programs": list(PROGRAMS), "comparators": list(COMPARATORS),
            "stimuli": list(STIMULI), "times": list(TIMES), "protocols": list(PROTOCOLS),
            "folds": N_FOLDS, "fold_seed": FOLD_SEED, "matching_bins": MATCH_BINS,
            "matching_replicates": MATCH_REPLICATES, "matching_seed": MATCH_SEED,
            "ridge_lambda": RIDGE_LAMBDA, "ridge_maxiter": MAX_ITER, "gradient_tolerance": GRADIENT_TOL,
            "minimum_training_class_lines": MIN_LINES, "minimum_training_class_runs": MIN_RUNS,
            "bootstrap_replicates": 5000, "bootstrap_seed": BOOTSTRAP_SEED,
            "identity_scopes": ["mapped_primary", "all_source_lines_sensitivity"],
            "predictive_variants": ["ets2_g1_dn", "ets2_g1_dn_without_comparators"],
            "time_weights": [0.5, 0.5], "protocol_weights": "fixed_metadata_prediction_line_shares",
            "scaling": "equal_line_weight_training_fold_paired_changes_only",
            "reference": "equal_line_weight_training_run_matched_controls_only",
            "quantity": "log2(released_58243_feature_CPM+1)", "normalization_gene_chunk": 256,
            "scoring_gene_chunk": 256, "prediction_uncertainty": "conditional_fixed_OOF_run_bootstrap",
            "pooled_response_median_interval": False,
            "pooled_response_mean_interval": False,
            "response_p_values": False,
            "missing_fit_rule": "all_expected_predictions_required_no_fallback",
            "audit_schema": "macromap-private-audit-v1",
            "independent_source_check": {"folds": [0], "protocols": list(PROTOCOLS), "times": list(TIMES),
                "training_pairs_per_context": 4, "heldout_pairs_per_context": 4, "programs": list(PROGRAMS),
                "full_universe_samples": 8, "ordering": "sha256(macromap-independent-v1:+canonicalkey)",
                "canonical_pair_key": "compact JSON [protocol,run,line,integer_time,stimulus]",
                "canonical_sample_key": "original sample_id string"}}


def runtime_versions() -> dict:
    import platform
    import scipy
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}


def execution_plan_template(prepared_directory: Path, output_directory: Path, *, root: Path = ROOT) -> dict:
    """Return an UNFROZEN template. No numerical source field is read."""
    prepared_directory = _absolute(prepared_directory)
    summary = prepared_directory / "source-only-readiness.json"
    return {"analysis_id": "macromap-fixed-program-context-v1", "status": "DRAFT_NOT_AUTHORIZED",
            "execute_real_expression": False, "frozen_at_utc": None,
            "script_sha256": sha256(Path(__file__)),
            "runner_script_sha256": sha256(Path(__file__).with_name("run_macromap_program_context.py")),
            "runtime": runtime_versions(),
            "prepared_directory": str(prepared_directory.relative_to(_absolute(root))),
            "prepared_summary_sha256": sha256(summary),
            "output_directory": str(_absolute(output_directory).relative_to(_absolute(root))),
            "method_contract": execution_contract()}


def validate_execution_plan(plan: Mapping, *, root: Path = ROOT) -> dict:
    from datetime import datetime
    if (plan.get("analysis_id") != "macromap-fixed-program-context-v1" or plan.get("status") != "frozen"
            or plan.get("execute_real_expression") is not True):
        raise ValueError("Real expression execution requires a separately frozen, explicitly authorized plan")
    try:
        frozen_at = datetime.fromisoformat(plan["frozen_at_utc"])
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError("Missing valid execution freeze timestamp") from exc
    if frozen_at.tzinfo is None:
        raise ValueError("Execution freeze timestamp must include a timezone")
    if (plan.get("script_sha256") != sha256(Path(__file__)) or
            plan.get("runner_script_sha256") != sha256(Path(__file__).with_name("run_macromap_program_context.py"))):
        raise ValueError("Execution method/runner script differs from the parent freeze")
    if plan.get("runtime") != runtime_versions() or plan.get("method_contract") != execution_contract():
        raise ValueError("Execution runtime or fixed statistical method differs from freeze")
    for key in ("prepared_directory", "output_directory"):
        relative = Path(plan[key])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Execution directories must be literal root-relative ignored paths")
        absolute = _absolute(root) / relative
        anchor = _absolute(root) / "work/macromap-design"
        if absolute == anchor or not absolute.is_relative_to(anchor):
            raise ValueError("Execution paths must remain under ignored work/macromap-design")
        reject_path_aliases(absolute)
    design = load_prepared_design(root / plan["prepared_directory"],
                                  expected_summary_sha256=plan["prepared_summary_sha256"])
    if design["summary"]["preparation_script_sha256"] != plan["script_sha256"]:
        raise ValueError("Repeat source-only preparation with the final frozen code before execution")
    verify_pins(root)
    return design


def save_npz_exclusive(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    reject_path_aliases(path)
    if any(np.asarray(value).dtype.hasobject for value in arrays.values()):
        raise ValueError("Private audit arrays must not require pickle")
    with path.open("xb") as fh:
        np.savez_compressed(fh, **arrays)


def score_program_blockwise(log_values: np.ndarray, sample_indices: Sequence[int],
                            targets: np.ndarray, backgrounds: np.ndarray, *, chunk_rows: int = 256) -> tuple:
    """All program raw/background means and detected fractions in one pass.

    targets/backgrounds are program x gene. No expression-dependent gene
    membership or normalization is learned by this linear accumulation.
    """
    targets, backgrounds = np.asarray(targets, float), np.asarray(backgrounds, float)
    if (targets.ndim != 2 or backgrounds.shape != targets.shape or targets.shape[1] != log_values.shape[0]
            or not np.isfinite(targets).all() or not np.isfinite(backgrounds).all() or chunk_rows < 1):
        raise ValueError("Invalid program weighting matrix")
    indices = np.asarray(sample_indices, dtype=int)
    if (len(indices) != len(set(indices)) or np.any(indices < 0)
            or np.any(indices >= log_values.shape[1])):
        raise ValueError("Invalid scoring sample indices")
    raw = np.zeros((len(indices), len(targets)))
    background, detected = np.zeros_like(raw), np.zeros_like(raw)
    for start in range(0, log_values.shape[0], chunk_rows):
        stop = min(start + chunk_rows, log_values.shape[0])
        values = np.asarray(log_values[start:stop, indices], dtype=float)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("Invalid normalized expression values")
        raw += values.T @ targets[:, start:stop].T
        background += values.T @ backgrounds[:, start:stop].T
        detected += (values > 0).T @ targets[:, start:stop].T
    return raw, background, detected


def _pair_matrix(pairs: Sequence[Mapping], sample_position: Mapping[int, int],
                 sample_scores: np.ndarray) -> np.ndarray:
    if not pairs:
        return np.empty((0, sample_scores.shape[1]))
    return np.array([sample_scores[sample_position[int(r["treatment_index"])]] -
                     sample_scores[sample_position[int(r["control_index"])]] for r in pairs])


def canonical_pair_key(row: Mapping) -> str:
    return json.dumps([row["protocol"], row["run"], row["line"], int(row["time"]), row["stimulus"]],
                      separators=(",", ":"), ensure_ascii=True)


def _pair_audit_arrays(prefix: str, pairs: Sequence[Mapping]) -> dict:
    result = {prefix + "_" + key: np.asarray([r[key] for r in pairs], dtype=str)
              for key in ("line", "run", "protocol", "stimulus", "treatment_id", "control_id")}
    result[prefix + "_class_index"] = np.asarray([STIMULI.index(r["stimulus"]) for r in pairs], dtype=np.int64)
    result[prefix + "_time"] = np.asarray([r["time"] for r in pairs], dtype=np.int64)
    result[prefix + "_treatment_index"] = np.asarray([r["treatment_index"] for r in pairs], dtype=np.int64)
    result[prefix + "_control_index"] = np.asarray([r["control_index"] for r in pairs], dtype=np.int64)
    result[prefix + "_line_weights"] = line_weights([r["line"] for r in pairs]) if pairs else np.array([])
    # Explicit aliases are part of the independent verifier's agreed schema.
    result[prefix + "_lines"] = result[prefix + "_line"]
    result[prefix + "_runs"] = result[prefix + "_run"]
    result[prefix + "_y"] = result[prefix + "_class_index"]
    result[prefix + "_pair_ids"] = np.asarray([canonical_pair_key(r) for r in pairs], dtype=str)
    return result


def save_fit_bundle(directory: Path, result: Mapping, train_x: np.ndarray, test_x: np.ndarray,
                    train_pairs: Sequence[Mapping], test_pairs: Sequence[Mapping], *,
                    predictor_program: str, fold_gate: Mapping, identity_scope: str = "synthetic_or_unspecified",
                    protocol: str | None = None, time: int | None = None, fold: int | None = None) -> dict:
    """Per-fit source-linkable NPZ plus numerical/failure audit; no pickle."""
    arrays = {"train_x_raw": np.asarray(train_x), "test_x_raw": np.asarray(test_x),
              "feature_programs": np.asarray(list(COMPARATORS) + [predictor_program]),
              "class_order": np.asarray(STIMULI),
              **_pair_audit_arrays("train", train_pairs), **_pair_audit_arrays("test", test_pairs)}
    arrays["train_x"], arrays["test_x"] = arrays["train_x_raw"], arrays["test_x_raw"]
    arrays["predictor_names"] = arrays["feature_programs"]
    scaler = result.get("scaler")
    if scaler is not None:
        arrays.update(scaler_mean=scaler.mean, scaler_scale=scaler.scale, scaler_constant=scaler.constant,
                      train_x_standardized=scaler.transform(train_x), test_x_standardized=scaler.transform(test_x))
    fits = {}
    for name in ("baseline", "extended"):
        fit = result.get(name)
        if fit is None:
            fits[name] = {"status": result["status"]}
        else:
            fits[name] = {"status": fit.status, "iterations": fit.iterations, "gradient_inf": fit.gradient_inf,
                          "objective": fit.objective, "optimizer_message": fit.optimizer_message}
            if fit.coefficients is not None:
                arrays[name + "_coefficients"] = fit.coefficients
    predictions = result.get("predictions", [])
    if predictions:
        arrays["baseline_log_probabilities"] = np.array([r["baseline_log_probabilities"] for r in predictions])
        arrays["extended_log_probabilities"] = np.array([r["extended_log_probabilities"] for r in predictions])
        arrays["baseline_logp"], arrays["extended_logp"] = arrays["baseline_log_probabilities"], arrays["extended_log_probabilities"]
        arrays["baseline_loss"] = np.array([r["baseline_loss"] for r in predictions])
        arrays["extended_loss"] = np.array([r["extended_loss"] for r in predictions])
        arrays["gain"] = arrays["baseline_loss"] - arrays["extended_loss"]
    save_npz_exclusive(directory / "fit-arrays.npz", arrays)
    rows = list(train_pairs) + list(test_pairs)
    audit = {"schema": "macromap-private-audit-v1", "status": result["status"], "fold_gate": fold_gate,
             "identity_scope": identity_scope, "variant": predictor_program,
             "protocol": protocol if protocol is not None else rows[0]["protocol"] if rows else None,
             "time": time if time is not None else rows[0]["time"] if rows else None, "fold": fold,
             "feature_programs": list(COMPARATORS) + [predictor_program], "class_order": STIMULI,
             "penalty": RIDGE_LAMBDA, "loss_normalization": "sum_weights_one_equal_lines",
             "intercepts": "unpenalized_symmetric_softmax", "fits": fits,
             "preprocessing_scope": "both_train_and_test_X_use_this_outer_fold_training_control_reference"}
    exclusive_text(directory / "fit-audit.json", json.dumps(original._json_safe(audit), indent=2, allow_nan=False) + "\n")
    return audit


def _response_panel(rows: Sequence[Mapping], *, scope: str) -> list[dict]:
    """Complete protocol-first 20-condition registry, then descriptive pooling."""
    panel = []
    for protocol in (*PROTOCOLS, "pooled_descriptive"):
        for time in TIMES:
            for stimulus in REGISTRY_STIMULI:
                for program in PROGRAMS:
                    selected = [r for r in rows if r["time"] == time and r["stimulus"] == stimulus
                                and r["program_id"] == program
                                and (protocol == "pooled_descriptive" or r["protocol"] == protocol)]
                    identity = {"identity_scope": scope, "protocol": protocol, "time": time,
                                "stimulus": stimulus, "program_id": program, "expected_paired_lines": len(selected),
                                "expected_paired_runs": len({(r["protocol"], r["run"]) for r in selected})}
                    if stimulus == "PIC":
                        panel.append({**identity, "status": "unavailable_mock_identity"})
                        continue
                    for quantity in ("raw_change", "background_change", "adjusted_change"):
                        values = [r[quantity] for r in selected]
                        if protocol == "pooled_descriptive":
                            result = pooled_response_summary(values, [r["line"] for r in selected],
                                [r["run"] for r in selected], [r["protocol"] for r in selected])
                        else:
                            result = clustered_response_summary(values, [r["line"] for r in selected],
                                [r["run"] for r in selected])
                        panel.append({**identity, "quantity": quantity,
                                      "input_statuses": dict(Counter(r["status"] for r in selected)), **result})
    return panel


def pooled_response_summary(values: Sequence[float], lines: Sequence[str], runs: Sequence[str],
                            protocols: Sequence[str]) -> dict:
    """Pooled mean and median are descriptive only under the final freeze rule.

    No pooled-response bootstrap or weighted-quantile convention is defined.
    This does not change protocol-specific intervals or the primary pooled
    predictive-gain interval. Missing expected values remain unavailable.
    """
    if len(protocols) != len(values) or any(p not in PROTOCOLS for p in protocols):
        raise ValueError("Unknown/misaligned pooled-response protocol")
    result = clustered_response_summary(values, lines, [p + ":" + r for p, r in zip(protocols, runs)], compute_interval=False)
    counts = Counter(protocols)
    shares = {p: counts[p] / len(protocols) for p in PROTOCOLS} if len(protocols) else {}
    return {**result, "interval_status": "not_computed_prespecified",
            "median_interval_status": "not_computed_prespecified", "mean_interval_status": "not_computed_prespecified",
            "median_interval": None, "mean_interval": None, "original_contrast_protocol_line_shares": shares}


def _prediction_diagnostics(rows: Sequence[Mapping], expected: Sequence[Mapping], point: Mapping) -> list[dict]:
    diagnostics = []
    if point["status"] != "available":
        return [{"status": point["status"], "diagnostic": "unavailable_incomplete_predictions"}]
    for protocol in PROTOCOLS:
        for run in sorted({r["run"] for r in expected if r["protocol"] == protocol}):
            keep = lambda r: not (r["protocol"] == protocol and r["run"] == run)
            subset, targets = list(filter(keep, rows)), list(filter(keep, expected))
            result = aggregate_prediction(subset, expected_pairs=targets,
                                          protocol_weights=point["protocol_weights"])
            diagnostics.append({"diagnostic": "leave_one_run_out_fixed_predictions", "protocol": protocol,
                                "omitted_run": run, "status": result["status"], "gain": result.get("primary_gain")})
    return diagnostics


def execute_fixed_context(log_values: np.ndarray, design: Mapping, output: Path) -> dict:
    """Post-freeze engine; never invoked by source-only prepare().

    Caller must validate the parent execution freeze and create a fresh private
    output directory. This pure-input engine is also available to independent
    synthetic verification without reading a real expression source.
    """
    if log_values.shape != (len(design["feature_ids"]), len(design["samples"])):
        raise ValueError("Expression/design dimensions disagree")
    summaries = {}
    fit_manifest = []
    feature_ids, sample_ids = design["feature_ids"], design["sample_ids"]
    feature_axis_hash = hashlib.sha256("\n".join(feature_ids).encode()).hexdigest()
    sample_axis_hash = hashlib.sha256("\n".join(sample_ids).encode()).hexdigest()
    target_offsets = np.cumsum([0] + [len(design["program_indices"][p]) for p in PROGRAMS])
    target_flat = np.array([i for p in PROGRAMS for i in design["program_indices"][p]], dtype=np.int64)
    n_genes = len(feature_ids)
    program_targets = np.zeros((len(PROGRAMS), n_genes))
    for j, program in enumerate(PROGRAMS):
        target = design["program_indices"][program]
        if len(target):
            program_targets[j, target] = 1 / len(target)
    for scope in execution_contract()["identity_scopes"]:
        scope_dir = output / scope
        reject_path_aliases(scope_dir)
        scope_dir.mkdir(exist_ok=False)
        all_pairs = pair_registry(design["samples"], include_nomatch=scope != "mapped_primary")
        predictive = prediction_cohort(all_pairs)
        write_csv(scope_dir / "pair-registry.csv", all_pairs)
        write_csv(scope_dir / "expected-prediction-pairs.csv", predictive)
        predictions = {variant: [] for variant in execution_contract()["predictive_variants"]}
        response_rows, fold_diagnostics = [], []
        for protocol in PROTOCOLS:
            for time in TIMES:
                for fold in range(N_FOLDS):
                    context = f"{protocol}-t{time}-fold{fold}"
                    directory = scope_dir / context
                    directory.mkdir(exist_ok=False)
                    available_pairs = [r for r in all_pairs if r["status"] == "available"
                                       and r["protocol"] == protocol and r["time"] == time]
                    selected_indices = sorted({int(r[k]) for r in available_pairs
                                               for k in ("treatment_index", "control_index")})
                    positions = {index: j for j, index in enumerate(selected_indices)}
                    train, test, fold_gate = split_pairs(predictive, design["folds"], protocol, time, fold)
                    reference_error = None
                    try:
                        reference, ref_audit = training_control_reference(log_values, all_pairs, design["folds"], protocol, time, fold)
                    except ValueError as exc:
                        reference_error = str(exc)
                        reference, ref_audit = np.full(n_genes, np.nan), {"reference_lines": 0, "reference_runs": 0, "control_indices": []}
                    if reference_error is None:
                        match = matched_weights(reference, feature_ids, design["program_indices"],
                                                design["pool_exclusion_indices"], protocol, time, fold)
                        feature_bin, _ = original.expression_bins(reference, feature_ids, bins=MATCH_BINS)
                    else:
                        match = {p: {"available": False, "reason": "reference_unavailable"} for p in PROGRAMS}
                        feature_bin = np.full(n_genes, -1)
                    background_weights = np.array([match[p].get("background_weights", np.zeros(n_genes)) for p in PROGRAMS])
                    program_status = {p: ("mapping_failed" if not design["mapping_gates"][p] else
                                          "available" if match[p]["available"] else match[p]["reason"]) for p in PROGRAMS}
                    raw, background, detected = score_program_blockwise(
                        log_values, selected_indices, program_targets, background_weights,
                        chunk_rows=execution_contract()["scoring_gene_chunk"])
                    for j, p in enumerate(PROGRAMS):
                        if not len(design["program_indices"][p]):
                            raw[:, j] = np.nan
                            detected[:, j] = np.nan
                        if program_status[p] != "available":
                            background[:, j] = np.nan
                    adjusted = raw - background
                    ref_indices = ref_audit["control_indices"]
                    saved_background_weights = background_weights.copy()
                    for j, p in enumerate(PROGRAMS):
                        if program_status[p] != "available":
                            saved_background_weights[j] = np.nan
                    save_npz_exclusive(directory / "score-reference-arrays.npz", {
                        "reference": reference, "feature_bin": feature_bin, "program_ids": np.asarray(PROGRAMS),
                        "feature_axis_sha256": np.asarray(feature_axis_hash), "sample_axis_sha256": np.asarray(sample_axis_hash),
                        "target_gene_indices": target_flat, "target_gene_offsets": target_offsets,
                        "pool_exclusion_indices": np.asarray(design["pool_exclusion_indices"], np.int64),
                        "target_weights": program_targets, "background_weights": saved_background_weights,
                        "matching_available": np.asarray([program_status[p] == "available" for p in PROGRAMS]),
                        "selected_sample_indices": np.asarray(selected_indices, np.int64),
                        "selected_sample_ids": np.asarray([sample_ids[i] for i in selected_indices]),
                        "training_control_indices": np.asarray(ref_indices, np.int64),
                        "training_control_ids": np.asarray([sample_ids[i] for i in ref_indices], dtype=str),
                        "sample_raw_scores": raw, "sample_background_scores": background,
                        "sample_adjusted_scores": adjusted, "sample_detected_fractions": detected})
                    matching_audit = {p: {k: v for k, v in info.items() if not isinstance(v, np.ndarray)} for p, info in match.items()}
                    exclusive_text(directory / "score-reference-audit.json", json.dumps(original._json_safe({
                        "scope": scope, "protocol": protocol, "time": time, "heldout_fold": fold,
                        "program_status": program_status, "reference_error": reference_error,
                        "reference": ref_audit, "matching": matching_audit}), indent=2, allow_nan=False) + "\n")
                    response_test = [r for r in available_pairs if design["folds"][(protocol, r["run"])] == fold]
                    deltas = [_pair_matrix(response_test, positions, values) for values in (raw, background, adjusted)]
                    for i, pair in enumerate(response_test):
                        for j, p in enumerate(PROGRAMS):
                            tr, ctrl = positions[pair["treatment_index"]], positions[pair["control_index"]]
                            response_rows.append({**pair, "identity_scope": scope, "heldout_fold": fold, "program_id": p,
                                "status": program_status[p], "raw_change": float(deltas[0][i, j]),
                                "background_change": float(deltas[1][i, j]), "adjusted_change": float(deltas[2][i, j]),
                                "treatment_detected_fraction": float(detected[tr, j]),
                                "control_detected_fraction": float(detected[ctrl, j])})
                    all_train_x, all_test_x = _pair_matrix(train, positions, adjusted), _pair_matrix(test, positions, adjusted)
                    for variant in predictions:
                        variant_dir = directory / variant
                        variant_dir.mkdir(exist_ok=False)
                        columns = [PROGRAMS.index(p) for p in (*COMPARATORS, variant)]
                        train_x, test_x = all_train_x[:, columns], all_test_x[:, columns]
                        result = ({"status": fold_gate["status"], "predictions": []} if fold_gate["status"] != "available" else
                                  heldout_increment(train_x, test_x, train, test))
                        fit_audit = save_fit_bundle(variant_dir, result, train_x, test_x, train, test,
                                        predictor_program=variant, fold_gate=fold_gate, identity_scope=scope,
                                        protocol=protocol, time=time, fold=fold)
                        fit_manifest.append({**{k: fit_audit[k] for k in ("identity_scope", "variant", "protocol", "time", "fold", "status")},
                            "fit_npz": str((variant_dir / "fit-arrays.npz").relative_to(output)),
                            "fit_npz_sha256": sha256(variant_dir / "fit-arrays.npz"),
                            "fit_audit": str((variant_dir / "fit-audit.json").relative_to(output)),
                            "stratum_npz": str((directory / "score-reference-arrays.npz").relative_to(output)),
                            "stratum_audit": str((directory / "score-reference-audit.json").relative_to(output))})
                        for row in result["predictions"]:
                            row.update(identity_scope=scope, predictor_program=variant, heldout_fold=fold)
                        predictions[variant].extend(result["predictions"])
                        losses = line_time_losses(result["predictions"]) if result["predictions"] else []
                        fold_diagnostics.append({"identity_scope": scope, "protocol": protocol, "time": time,
                            "fold": fold, "predictor_program": variant, **{k: v for k, v in fold_gate.items() if k != "status"},
                            "fold_gate_status": fold_gate["status"], "status": result["status"],
                            "line_weighted_gain": float(np.mean([r["gain"] for r in losses])) if losses else None})
        write_csv(scope_dir / "line-responses.csv", response_rows)
        exclusive_text(scope_dir / "response-panel.json", json.dumps(original._json_safe(_response_panel(response_rows, scope=scope)),
                                                                     indent=2, allow_nan=False) + "\n")
        exclusive_text(scope_dir / "fold-diagnostics.json", json.dumps(original._json_safe(fold_diagnostics), indent=2, allow_nan=False) + "\n")
        summaries[scope] = {}
        for variant, rows in predictions.items():
            variant_dir = scope_dir / variant
            variant_dir.mkdir(exist_ok=False)
            if rows:
                write_csv(variant_dir / "heldout-loss-records.csv", rows)
            point = cluster_bootstrap_gain(rows, expected_pairs=predictive)
            summaries[scope][variant] = point
            diagnostics = _prediction_diagnostics(rows, predictive, point)
            exclusive_text(variant_dir / "incremental-endpoint.json", json.dumps(original._json_safe(point), indent=2, allow_nan=False) + "\n")
            exclusive_text(variant_dir / "run-influence.json", json.dumps(original._json_safe(diagnostics), indent=2, allow_nan=False) + "\n")
    exclusive_text(output / "fit-manifest.json", json.dumps(fit_manifest, indent=2, allow_nan=False) + "\n")
    return summaries


def run_frozen_analysis(plan_path: Path, *, expected_plan_sha256: str, root: Path = ROOT) -> dict:
    """Only real-data execution entry: separate freeze, fresh ignored outputs."""
    from datetime import datetime, timezone
    if sha256(plan_path) != expected_plan_sha256:
        raise ValueError("Execution plan differs from supplied parent freeze hash")
    plan = json.loads(plan_path.read_text())
    design = validate_execution_plan(plan, root=root)
    output = fresh_private_run(root, root / plan["output_directory"])
    receipt = {"stage": "authorized_post_freeze_execution", "plan_sha256": expected_plan_sha256,
               "script_sha256": sha256(Path(__file__)), "runner_script_sha256": plan["runner_script_sha256"],
               "prepared_summary_sha256": plan["prepared_summary_sha256"],
               "started_at_utc": datetime.now(timezone.utc).isoformat(), "method_contract": execution_contract(),
               "runtime": runtime_versions(), "source_pins": SOURCE_PINS, "ready_for_inference": True}
    exclusive_text(output / "execution-start.json", json.dumps(receipt, indent=2) + "\n")
    try:
        cache = output / "released-log2cpm.f64"
        cache_audit = normalize_count_cache(root / RAW / "MacroMap_raw_expression.txt.gz", cache,
                                           design["feature_ids"], design["sample_ids"], chunk_rows=256)
        frozen_raw_sha = SOURCE_PINS[RAW + "/MacroMap_raw_expression.txt.gz"]["sha256"]
        if cache_audit["source_sha256"] != frozen_raw_sha:
            raise ValueError("Normalized raw source differs from frozen pin; no matching or fit permitted")
        if sha256(plan_path) != expected_plan_sha256:
            raise ValueError("Execution plan changed during normalization")
        validate_execution_plan(plan, root=root)
        normalized = np.memmap(cache, dtype="<f8", mode="r", shape=(len(design["feature_ids"]), len(design["samples"])))
        summaries = execute_fixed_context(normalized, design, output)
        del normalized
        if sha256(plan_path) != expected_plan_sha256:
            raise ValueError("Execution plan changed; no complete publication")
        validate_execution_plan(plan, root=root)
        if sha256(cache) != cache_audit["cache_sha256"]:
            raise ValueError("Normalized cache changed; no complete publication")
        complete = {**receipt, "stage": "execution_complete_all_prespecified_statuses_retained",
                    "completed_at_utc": datetime.now(timezone.utc).isoformat(), "cache": cache_audit,
                    "endpoints": summaries,
                    "artifact_sha256": {str(path.relative_to(output)): sha256(path) for path in sorted(output.rglob("*")) if path.is_file()}}
        exclusive_text(output / "execution-complete.json", json.dumps(original._json_safe(complete), indent=2, allow_nan=False) + "\n")
        return complete
    except BaseException as exc:
        exclusive_text(output / "execution-failed.json", json.dumps({**receipt, "stage": "execution_failed_no_complete_endpoint",
                       "error_type": type(exc).__name__, "error": str(exc)}, indent=2) + "\n")
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / "work/macromap-design/preparation")
    args = parser.parse_args()
    summary = prepare(args.root, args.output)
    print(json.dumps({key: summary[key] for key in ("stage", "ready_for_inference", "released_genes", "samples",
                                                  "mapped_lines", "protocol_design", "fold_gate_status")}, indent=2))


if __name__ == "__main__":
    main()
