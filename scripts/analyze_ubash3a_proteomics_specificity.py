#!/usr/bin/env python3
"""Frozen whole-sequence peptide specificity and deposited-search metadata audit.

This script does not search mass spectra. It requires checksum-pinned complete
reference files and does not substitute missing references with empty inputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import html
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ubash3a-proteomics-specificity"
PLAN = ROOT / "config/ubash3a-proteomics-specificity-plan.json"
LOCK = ROOT / "config/ubash3a-proteomics-reference-lock.json"
OUT = ROOT / "data/derived/ubash3a-proteomics-specificity"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta(path):
    """Yield original full header and upper-case sequence; reject broken records."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    header, parts = None, []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(parts).upper()
                    if not sequence or not re.fullmatch(r"[A-Z*]+", sequence):
                        raise ValueError(f"Invalid sequence: {header}")
                    yield header, sequence
                header, parts = line[1:], []
                if not header:
                    raise ValueError("Empty FASTA header")
            else:
                if header is None:
                    raise ValueError("Sequence precedes FASTA header")
                parts.append(line)
        if header is not None:
            sequence = "".join(parts).upper()
            if not sequence or not re.fullmatch(r"[A-Z*]+", sequence):
                raise ValueError(f"Invalid sequence: {header}")
            yield header, sequence


def boundary(sequence, position):
    """Trypsin cleavage boundary at zero-based gap, including protein termini."""
    return position in (0, len(sequence)) or (
        sequence[position - 1] in "KR" and sequence[position] != "P"
    )


def match_positions(sequence, peptide, il_equivalent=False):
    if il_equivalent:
        sequence, peptide = sequence.replace("I", "L"), peptide.replace("I", "L")
    start = 0
    while (start := sequence.find(peptide, start)) >= 0:
        yield start
        start += 1  # Keep overlapping occurrences.


def match_record(sequence, peptide, start):
    end = start + len(peptide)
    left, right = boundary(sequence, start), boundary(sequence, end)
    missed = sum(boundary(sequence, gap) for gap in range(start + 1, end))
    return {
        "start_1based": start + 1,
        "end_1based_inclusive": end,
        "matched_sequence": sequence[start:end],
        "preceding_residue": sequence[start - 1] if start else "-",
        "following_residue": sequence[end] if end < len(sequence) else "-",
        "n_tryptic": left,
        "c_tryptic": right,
        "internal_tryptic_cleavages": missed,
        "fits_frozen_candidate_digest": left and right and missed <= 1 and 7 <= len(peptide) <= 35,
    }


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError("Explicit columns required for empty table")
        fields = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extract_settings(prefix):
    # PD stores processing-tree settings as escaped text, not XML attributes.
    found = re.findall(r"<search_summary>(.*?)</search_summary>", prefix, re.S)
    if len(found) != 1:
        raise ValueError(f"Expected one complete settings element, found {len(found)}")
    settings = html.unescape(found[0])
    names = re.findall(r"^\s*([^\r\n]*\.raw)\s*$", settings, re.M | re.I)
    raw_names = [name.strip().replace("/", "\\").rsplit("\\", 1)[-1] for name in names]
    if not raw_names or len(set(raw_names)) != len(raw_names):
        raise ValueError("Missing or repeated raw-file references")
    experiment = re.search(r"^Search name:\s*(.+)$", settings, re.M)
    if not experiment:
        raise ValueError("Missing experiment name")
    return experiment.group(1).strip(), raw_names, settings


def audit_metadata(output):
    inventory_path = ROOT / "data/derived/ubash3a-proteomics-readiness/PXD000376-file-inventory.csv"
    inventory = list(csv.DictReader(inventory_path.open()))
    raw_index = {row["file_name"]: row for row in inventory if row["file_name"].lower().endswith(".raw")}
    case_index = {}
    for name in raw_index:
        case_index.setdefault(name.casefold(), []).append(name)
    mapping, groups = [], []
    for row in inventory:
        if not row["file_name"].endswith(".pep.xml"):
            continue
        path = RAW / (row["file_name"] + ".header-prefix")
        experiment, names, settings = extract_settings(path.read_text())
        databases = re.findall(r"^Protein Database:\s*(.+)$", settings, re.M)
        enzymes = re.findall(r"^Enzyme Name:\s*(.+)$", settings, re.M)
        precursor = re.findall(r"^Precursor Mass Tolerance:\s*(.+)$", settings, re.M)
        fragment = re.findall(r"^Fragment Mass Tolerance:\s*(.+)$", settings, re.M)
        (output / "search-settings").mkdir(parents=True, exist_ok=True)
        (output / "search-settings" / (row["file_name"] + ".txt")).write_text(settings + "\n")
        matched, size = 0, 0
        for original in names:
            aliases = case_index.get(original.casefold(), [])
            if len(aliases) > 1:
                raise ValueError(f"Ambiguous archive filename case: {original}")
            name = original if original in raw_index else (aliases[0] if aliases else "")
            status = "exact" if original in raw_index else ("case_only_alias" if name else "missing_from_archive")
            if name:
                matched += 1
                size += int(raw_index[name]["size_bytes"])
            mapping.append({
                "experiment": experiment, "search_file": row["file_name"],
                "source_raw_name": original, "archive_raw_name": name,
                "mapping_status": status, "archive_bytes": raw_index[name]["size_bytes"] if name else "",
                "settings_source_sha256": sha256(path),
            })
        groups.append({
            "experiment": experiment, "search_file": row["file_name"],
            "referenced_raw_files": len(names), "matched_raw_files": matched,
            "missing_raw_files": len(names) - matched, "matched_archive_bytes": size,
            "database_names_as_published": " | ".join(databases),
            "enzymes_as_published": " | ".join(enzymes),
            "precursor_tolerances_as_published": " | ".join(precursor),
            "fragment_tolerances_as_published": " | ".join(fragment),
            "full_settings_sha256": hashlib.sha256(settings.encode()).hexdigest(),
            "sample_scope": "source experiment label; not an independently verified donor identity",
        })
    if len(groups) != 12:
        raise ValueError(f"Expected all 12 deposited search headers, found {len(groups)}")
    write_csv(output / "experiment-raw-mapping.csv", mapping)
    write_csv(output / "experiment-search-metadata.csv", groups)
    selected = json.loads((ROOT / "config/ubash3a-proteomics-spectrum-selection.json").read_text())
    selected_rows = [row for row in mapping if row["experiment"] == selected["experiment_name"]]
    assert {row["archive_raw_name"] for row in selected_rows} == {row["file_name"] for row in selected["files"]}
    assert sum(int(row["archive_bytes"]) for row in selected_rows) == selected["download_bytes"]
    mapped_names = {row["archive_raw_name"] for row in mapping if row["archive_raw_name"]}
    return {"search_headers": len(groups), "file_references": len(mapping),
            "archive_raw_files": len(raw_index), "unique_mapped_archive_raw_files": len(mapped_names),
            "archive_raw_files_without_header_mapping": sorted(set(raw_index) - mapped_names),
            "unmatched_references": [row for row in mapping if not row["archive_raw_name"]],
            "case_only_aliases": sum(row["mapping_status"] == "case_only_alias" for row in mapping)}


def analyze(output):
    output.mkdir(parents=True, exist_ok=True)
    plan, lock = json.loads(PLAN.read_text()), json.loads(LOCK.read_text())
    assert lock["locked_before_candidate_matching"] is True
    candidates = plan["candidates"]
    assert len(candidates) == 4 and len(set(candidates)) == 4
    for path, expected in plan["input_pins"].items():
        assert sha256(ROOT / path) == expected, path
    matches, summaries, parent = [], [], Counter()
    for source in lock["sources"]:
        path = ROOT / source["local_path"]
        assert sha256(path) == source["sha256"], path
        records, residues, nonstandard = 0, 0, Counter()
        seen = set()
        for header, sequence in fasta(path):
            identifier = header.split()[0]
            if identifier in seen:
                raise ValueError(f"Repeated source identifier: {source['set_id']} {identifier}")
            seen.add(identifier)
            records += 1
            residues += len(sequence)
            nonstandard.update(c for c in sequence if c not in STANDARD)
            identity = hashlib.sha256(sequence.encode()).hexdigest()
            for peptide in candidates:
                for comparison in plan["comparisons"]:
                    for start in match_positions(sequence, peptide, comparison == "I_L_equivalent"):
                        record = {"candidate": peptide, "comparison": comparison,
                                  "set_id": source["set_id"], "source_group": source["source_group"],
                                  "record_id": identifier, "header": header,
                                  "sequence_sha256": identity, "protein_length": len(sequence),
                                  **match_record(sequence, peptide, start)}
                        matches.append(record)
                        if source["source_group"] == "conditional_T29" and comparison == "exact":
                            assert record["fits_frozen_candidate_digest"]
                            parent[peptide] += 1
        assert records == source["records"], (source["set_id"], records)
        summaries.append({"set_id": source["set_id"], "source_group": source["source_group"],
                          "records": records, "residues": residues,
                          "nonstandard_residue_counts_json": json.dumps(dict(sorted(nonstandard.items()))),
                          "local_path": source["local_path"], "sha256": source["sha256"]})
    assert parent == Counter({p: 2 for p in candidates}), parent
    matches.sort(key=lambda row: (row["candidate"], row["comparison"], row["set_id"], row["record_id"], row["start_1based"]))
    write_csv(output / "all-sequence-matches.csv", matches)
    write_csv(output / "reference-summary.csv", summaries)
    candidate_rows = []
    for peptide in candidates:
        row = {"candidate": peptide, "length": len(peptide)}
        for comparison in plan["comparisons"]:
            subset = [m for m in matches if m["candidate"] == peptide and m["comparison"] == comparison]
            for group in ("human_reference", "contaminant", "previous_models", "conditional_T29"):
                group_rows = [m for m in subset if m["source_group"] == group]
                row[f"{comparison}_{group}_records"] = len({(m["set_id"], m["record_id"]) for m in group_rows})
                row[f"{comparison}_{group}_tryptic_records"] = len({(m["set_id"], m["record_id"]) for m in group_rows if m["fits_frozen_candidate_digest"]})
        row["interpretation"] = (
            "has retrieved background sequence competitors" if any(
                m["candidate"] == peptide and m["source_group"] in {"human_reference", "contaminant"}
                for m in matches)
            else "no literal or I/L-equivalent match in retrieved backgrounds; not universal uniqueness or protein detection"
        )
        row["nested_in_other_candidate"] = " | ".join(other for other in candidates if peptide != other and peptide in other)
        candidate_rows.append(row)
    write_csv(output / "candidate-specificity.csv", candidate_rows)
    result = {"analysis": "sequence specificity only; no new spectrum identifications",
              "plan_sha256": sha256(PLAN), "reference_lock_sha256": sha256(LOCK),
              "reference_sets": len(summaries), "records": sum(x["records"] for x in summaries),
              "match_rows": len(matches), "candidate_count": len(candidates),
              "all_parent_checks_pass": True,
              "scope_limitations": ["Unknown/ambiguous source residues are not expanded into possible sequences",
                                    "Unannotated proteins, all population variants and other species are not exhaustive",
                                    "Nested/overlapping peptide strings are not independent protein evidence",
                                    "Conditional T29 translation does not establish RNA or protein existence"]}
    (output / "specificity-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = audit_metadata(args.output)
    (args.output / "metadata-summary.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"metadata": metadata, "specificity": None if args.metadata_only else analyze(args.output)}, indent=2))
