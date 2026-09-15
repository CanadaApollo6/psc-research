#!/usr/bin/env python3
"""Independent b/y fragment accounting for every reported fixed-target assignment.

This does not replace Comet competition or estimate novel-peptide FDR. It uses
the original centroid peaks, a fixed 0.02-Da matching window, and base b/y ions
without neutral losses. Its counts are not expected to equal Comet scores.
"""
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from Bio.SeqUtils import molecular_weight

from analyze_ubash3a_proteomics_specificity import sha256, write_csv
from qualify_ubash3a_spectra import decode_array, terms, NS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-proteomics-spectra"
PROTON = 1.007276466621
WATER = 18.0105646837
ISOTOPE = 1.0033548378


def modification_masses(sequence, encoded):
    masses = defaultdict(float)
    static_c = set()
    if encoded not in ("", "-"):
        for item in encoded.split(","):
            fields = item.split("_")
            assert len(fields) in (3, 4), item
            position, kind, mass = int(fields[0]) - 1, fields[1], float(fields[2])
            assert kind in {"S", "V"}
            if len(fields) == 4:
                assert fields[3] in {"N", "C"}
                position = -1 if fields[3] == "N" else len(sequence)
            if kind == "S" and 0 <= position < len(sequence) and sequence[position] == "C":
                assert abs(mass - 57.021464) < 0.00001
                static_c.add(position)
            masses[position] += mass
    for position, residue in enumerate(sequence):
        if residue == "C" and position not in static_c:
            masses[position] += 57.021464
    return dict(masses)


def peptide_mass(sequence):
    return molecular_weight(sequence, seq_type="protein", monoisotopic=True)


def fragments(sequence, modifications, precursor_charge):
    max_charge = max(1, min(3, precursor_charge - 1))
    ions = []
    for gap in range(1, len(sequence)):
        b = peptide_mass(sequence[:gap]) - WATER + sum(m for p, m in modifications.items() if p < gap)
        y = peptide_mass(sequence[gap:]) + sum(m for p, m in modifications.items() if p >= gap)
        for charge in range(1, max_charge + 1):
            ions.extend([("b", gap, charge, (b + charge * PROTON) / charge),
                         ("y", len(sequence) - gap, charge, (y + charge * PROTON) / charge)])
    return ions


def longest_run(values):
    longest = current = 0
    previous = None
    for value in sorted(set(values)):
        current = current + 1 if previous is not None and value == previous + 1 else 1
        longest = max(longest, current)
        previous = value
    return longest


def fragment_metrics(row, mz, intensities):
    sequence = row["plain_peptide"]
    modifications = modification_masses(sequence, row["modifications"])
    neutral = peptide_mass(sequence) + sum(modifications.values())
    assert abs(neutral - float(row["calc_neutral_mass"])) < 0.001, (sequence, neutral, row["calc_neutral_mass"])
    match_rows, observed_indices = [], set()
    for series, ordinal, charge, theoretical in fragments(sequence, modifications, int(row["charge"])):
        insertion = int(np.searchsorted(mz, theoretical))
        nearby = [index for index in (insertion - 1, insertion) if 0 <= index < len(mz)]
        closest = min(nearby, key=lambda index: abs(float(mz[index]) - theoretical)) if nearby else None
        matched = closest is not None and abs(float(mz[closest]) - theoretical) <= 0.02
        if matched:
            observed_indices.add(closest)
        match_rows.append({"series": series, "ordinal": ordinal, "charge": charge,
                           "theoretical_mz": theoretical, "matched_within_0_02_Da": matched,
                           "observed_mz": float(mz[closest]) if matched else "",
                           "mass_error_Da": float(mz[closest]) - theoretical if matched else "",
                           "observed_intensity": float(intensities[closest]) if matched else ""})
    matched = [row for row in match_rows if row["matched_within_0_02_Da"]]
    b = {row["ordinal"] for row in matched if row["series"] == "b"}
    y = {row["ordinal"] for row in matched if row["series"] == "y"}
    cleavage_gaps = b | {len(sequence) - ordinal for ordinal in y}
    denominator = float(intensities.sum())
    isotopes = [(float(row["exp_neutral_mass"]) - neutral - shift * ISOTOPE) / neutral * 1e6 for shift in (0, 1, 2)]
    best_shift = min(range(3), key=lambda shift: abs(isotopes[shift]))
    return {"independent_neutral_mass": neutral, "mass_agrees_with_Comet_within_0_001_Da": True,
            "nearest_allowed_isotope_shift": best_shift, "precursor_error_ppm_at_nearest_shift": isotopes[best_shift],
            "independent_b_ordinals_matched": len(b), "independent_y_ordinals_matched": len(y),
            "independent_cleavage_gaps_supported": len(cleavage_gaps), "possible_cleavage_gaps": len(sequence) - 1,
            "longest_consecutive_b_run": longest_run(b), "longest_consecutive_y_run": longest_run(y),
            "distinct_centroid_peaks_matched": len(observed_indices),
            "centroid_intensity_fraction_matched": sum(float(intensities[index]) for index in observed_indices) / denominator if denominator else 0,
            "fragment_match_window_Da": 0.02, "neutral_loss_ions_included_in_independent_count": False}, match_rows


def requested_spectra(path, wanted):
    found = {}
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] != "spectrum":
            continue
        scan = element.attrib["id"].rsplit("scan=", 1)[-1]
        if scan in wanted:
            arrays = {}
            for child in element.findall("m:binaryDataArrayList/m:binaryDataArray", NS):
                cv = terms(child)
                if "MS:1000514" in cv:
                    arrays["mz"] = decode_array(child, int(element.attrib["defaultArrayLength"]))
                if "MS:1000515" in cv:
                    arrays["intensity"] = decode_array(child, int(element.attrib["defaultArrayLength"]))
            assert scan not in found and set(arrays) == {"mz", "intensity"}
            found[scan] = arrays
        element.clear()
    assert set(found) == set(wanted)
    return found


def main():
    nominations_path = OUT / "candidate-all-reported-assignments.csv"
    nominations = list(csv.DictReader(nominations_path.open()))
    grouped = defaultdict(list)
    for row in nominations:
        grouped[row["file"]].append(row)
    reviews, matches = [], []
    candidate_set = {row["candidate"].replace("I", "L") for row in nominations}
    for file, candidates in grouped.items():
        mzml = ROOT / "work/ubash3a-proteomics-spectra/converted" / (file + ".mzML")
        spectra = requested_spectra(mzml, {row["scan"] for row in candidates})
        with mzml.with_suffix(".txt").open() as handle:
            next(handle)
            all_rows = list(csv.DictReader(handle, delimiter="\t"))
        for candidate in candidates:
            data = spectra[candidate["scan"]]
            relevant = [row for row in all_rows if row["scan"] == candidate["comet_scan"] and row["charge"] == candidate["charge"]]
            competitors = [row for row in relevant if row["plain_peptide"].replace("I", "L") not in candidate_set and not any(p.startswith("DECOY_") for p in row["protein"].split(","))]
            competitor = min(competitors, key=lambda row: int(row["num"])) if competitors else None
            for role, assignment in [("candidate", candidate)] + ([("best_exported_noncandidate_target", competitor)] if competitor else []):
                metrics, ions = fragment_metrics(assignment, data["mz"], data["intensity"])
                key = {"candidate": candidate["candidate"], "file": file, "scan": candidate["scan"],
                       "comet_scan": candidate["comet_scan"],
                       "candidate_rank": candidate["num"], "candidate_query_hit_index": candidate["query_hit_index"], "role": role,
                       "assignment_sequence": assignment["plain_peptide"], "assignment_modified_peptide": assignment["modified_peptide"],
                       "assignment_charge": assignment["charge"], "assignment_rank": assignment["num"]}
                reviews.append({**key, "comet_expectation": assignment["e-value"], "comet_xcorr": assignment["xcorr"],
                                "comet_matched_ions": assignment["ions_matched"], "comet_total_ions": assignment["ions_total"],
                                "protein_ids": assignment["protein"], "primary_candidate_win": candidate["is_primary_scan_winner"],
                                "candidate_global_spectrum_q": candidate["candidate_global_spectrum_q_value"],
                                "exported_noncandidate_target_available": competitor is not None, **metrics,
                                "interpretation": "Fragment accounting only; does not establish novel-peptide error control or protein existence"})
                matches.extend({**key, **ion} for ion in ions)
    fields = list(reviews[0]) if reviews else ["candidate", "file", "scan", "role", "assignment_sequence", "interpretation"]
    write_csv(OUT / "candidate-fragment-review.csv", reviews, fields)
    full_path = ROOT / "work/ubash3a-proteomics-spectra/analysis/candidate-fragment-ion-matches.csv"
    write_csv(full_path, matches, list(matches[0]) if matches else ["candidate", "file", "scan", "series", "ordinal", "charge", "theoretical_mz", "matched_within_0_02_Da"])
    report = {"all_reported_candidate_assignments_reviewed": True, "candidate_assignment_count": len(nominations),
              "distinct_candidate_spectra": len({(row["file"], row["scan"]) for row in nominations}),
              "candidate_and_competitor_review_rows": len(reviews), "fragment_ion_rows": len(matches),
              "method": "Independent Biopython monoisotopic peptide masses, exact exported modifications, b/y ions through min(3,z-1), nearest centroid within fixed 0.02 Da; no neutral losses or learned scoring",
              "limitation": "Independent base-ion counts are not a replacement for full-database competition; only first five Comet assignments were exported. Missing noncandidate competitor is retained as unavailable.",
              "biological_detection_claim": False,
              "files": {str(path.relative_to(ROOT)): sha256(path) for path in [nominations_path, OUT / "candidate-fragment-review.csv", full_path]}}
    (ROOT / "reports/ubash3a-proteomics-fragment-review.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
