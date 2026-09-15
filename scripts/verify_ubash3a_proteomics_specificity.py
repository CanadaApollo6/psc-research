#!/usr/bin/env python3
"""Independent SeqIO/regular-expression checks of sequence and XML-header results."""
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-proteomics-specificity"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plan = json.loads((ROOT / "config/ubash3a-proteomics-specificity-plan.json").read_text())
    lock = json.loads((ROOT / "config/ubash3a-proteomics-reference-lock.json").read_text())
    observed = list(csv.DictReader((OUT / "all-sequence-matches.csv").open()))
    expected = []
    source_counts = {}
    # No imports from the primary matching, FASTA or boundary implementation.
    for source in lock["sources"]:
        path = ROOT / source["local_path"]
        assert sha(path) == source["sha256"]
        opener = gzip.open if path.suffix == ".gz" else open
        records = 0
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                records += 1
                sequence = str(record.seq).upper()
                gaps = {0, len(sequence)} | {m.start() + 1 for m in re.finditer(r"[KR](?!P)", sequence)}
                for peptide in plan["candidates"]:
                    for comparison in plan["comparisons"]:
                        pattern = peptide if comparison == "exact" else "".join("[IL]" if aa in "IL" else aa for aa in peptide)
                        for match in re.finditer("(?=(" + pattern + "))", sequence):
                            start, end = match.start(), match.start() + len(peptide)
                            missed = len([gap for gap in gaps if start < gap < end])
                            expected.append({"candidate": peptide, "comparison": comparison,
                                             "set_id": source["set_id"], "source_group": source["source_group"],
                                             "record_id": record.id, "header": record.description,
                                             "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                             "protein_length": str(len(sequence)),
                                             "start_1based": str(start + 1), "end_1based_inclusive": str(end),
                                             "matched_sequence": sequence[start:end],
                                             "preceding_residue": sequence[start - 1] if start else "-",
                                             "following_residue": sequence[end] if end < len(sequence) else "-",
                                             "n_tryptic": str(start in gaps), "c_tryptic": str(end in gaps),
                                             "internal_tryptic_cleavages": str(missed),
                                             "fits_frozen_candidate_digest": str(start in gaps and end in gaps and missed <= 1 and 7 <= len(peptide) <= 35)})
        assert records == source["records"]
        source_counts[source["set_id"]] = records
    key = lambda row: tuple(row[field] for field in sorted(row))
    assert sorted(observed, key=key) == sorted(expected, key=key)
    primary_mapping = list(csv.DictReader((OUT / "experiment-raw-mapping.csv").open()))
    xml_mapping = []
    for path in sorted((ROOT / "data/raw/ubash3a-proteomics-specificity").glob("*.pep.xml.header-prefix")):
        prefix = path.read_text()
        # The bounded prefix contains the complete settings element; its outer
        # analysis wrapper can extend beyond 64 KiB into other result metadata.
        start = prefix.index("<search_summary>")
        end = prefix.index("</search_summary>") + len("</search_summary>")
        tree = ET.fromstring(prefix[start:end])
        settings = tree.text
        name = next(line.split(":", 1)[1].strip() for line in settings.splitlines() if line.startswith("Search name:"))
        for line in settings.splitlines():
            if line.strip().lower().endswith(".raw"):
                xml_mapping.append((name, line.strip().replace("/", "\\").split("\\")[-1]))
    assert Counter(xml_mapping) == Counter((row["experiment"], row["source_raw_name"]) for row in primary_mapping)
    result = {"all_checks_pass": True, "methods": ["Biopython SeqIO versus primary streaming FASTA parser",
              "Overlapping regular expressions versus substring search", "Independent cleavage-gap set versus per-position boundary function",
              "XML tree parsing versus primary escaped-text extraction"],
              "reference_counts": source_counts, "sequence_match_rows": len(expected),
              "sequence_match_fields_compared": len(expected) * len(expected[0]),
              "xml_raw_file_references": len(xml_mapping),
              "inputs": {str(path.relative_to(ROOT)): sha(path) for path in (OUT / "all-sequence-matches.csv", OUT / "experiment-raw-mapping.csv", ROOT / "config/ubash3a-proteomics-reference-lock.json")}}
    (ROOT / "reports/ubash3a-proteomics-specificity-verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
