#!/usr/bin/env python3
"""Generated-spectrum qualification of the pinned Comet binary.

No measured spectrum or T29 candidate is used. These are software controls,
not analytical sensitivity, specificity, biological evidence or FDR estimates.
"""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess

from Bio.SeqUtils import molecular_weight

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-spectra/synthetic-controls"
PROTON = 1.007276466621
WATER = 18.0105646837


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def neutral_mass(sequence):
    return molecular_weight(sequence, seq_type="protein", monoisotopic=True)


def generated_spectrum(sequence, modifications, charge=2):
    """Generate b/y ions independently using Biopython whole-peptide masses."""
    neutral = neutral_mass(sequence) + sum(modifications.values())
    ions = []
    for gap in range(1, len(sequence)):
        bmass = neutral_mass(sequence[:gap]) - WATER + sum(v for p, v in modifications.items() if p < gap)
        ymass = neutral_mass(sequence[gap:]) + sum(v for p, v in modifications.items() if p >= gap)
        ions.extend([(bmass + PROTON, 10000.0), (ymass + PROTON, 10000.0)])
    # Uninformative low-intensity ions ensure the minimum peak-count rule is met.
    ions.extend((70.1234 + 43.137 * index, 1.0) for index in range(5))
    return (neutral + charge * PROTON) / charge, sorted(ions), neutral


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    tool = json.loads((ROOT / "config/ubash3a-proteomics-tool-install.json").read_text())["sources"][1]
    binary = ROOT / tool["installed_binary"]
    assert sha(binary) == tool["binary_sha256"]
    cases = [
        {"name": "unmodified", "sequence": "PEPTIDER", "modifications": {}, "decoy": False},
        {"name": "oxidation_M", "sequence": "MPEPTIDER", "modifications": {0: 15.994915}, "decoy": False},
        {"name": "phospho_S", "sequence": "GASSTVQK", "modifications": {2: 79.966331}, "decoy": False},
        {"name": "pyro_N_terminal_Q", "sequence": "QPEPTIDER", "modifications": {0: -17.026549}, "decoy": False},
        {"name": "acetyl_protein_N", "sequence": "MPEPTIDER", "modifications": {-1: 42.010565}, "decoy": False},
        {"name": "fixed_carbamidomethyl_C", "sequence": "CPEPTIDER", "modifications": {0: 57.021464}, "decoy": False},
        {"name": "internal_decoy", "sequence": "PEPTIDER"[:-1][::-1] + "R", "modifications": {}, "decoy": True},
    ]
    database = WORK / "controls.fasta"
    unique_sequences = sorted({case["sequence"] for case in cases if not case["decoy"]})
    database.write_text("".join(f">CONTROL_{index}\n{sequence}\n" for index, sequence in enumerate(unique_sequences)))
    mgf = WORK / "controls.mgf"
    with mgf.open("w") as handle:
        for scan, case in enumerate(cases, start=1):
            precursor, ions, neutral = generated_spectrum(case["sequence"], case["modifications"])
            case.update(scan=scan, neutral_mass=neutral)
            handle.write(f"BEGIN IONS\nTITLE=synthetic_{case['name']}\nSCANS={scan}\nPEPMASS={precursor:.8f}\nCHARGE=2+\nRTINSECONDS={scan * 10}\n")
            handle.writelines(f"{mz:.8f} {intensity:.1f}\n" for mz, intensity in ions)
            handle.write("END IONS\n\n")
        handle.write("BEGIN IONS\nTITLE=outside_frozen_mass_range\nSCANS=8\nPEPMASS=9999\nCHARGE=2+\n")
        handle.writelines(f"{100 + 50 * i}.1234 1000\n" for i in range(12))
        handle.write("END IONS\n")
    base = (ROOT / "config/ubash3a-proteomics-comet.params").read_text()
    parameters = WORK / "controls.params"
    parameters.write_text(re.sub(r"^database_name\s*=.*$", "database_name = " + str(database), base, flags=re.M))
    expected = WORK / "expected.json"
    expected.write_text(json.dumps(cases, indent=2) + "\n")
    began = datetime.datetime.now(datetime.timezone.utc).isoformat()
    args = [str(binary), "-P" + str(parameters), str(mgf)]
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    (WORK / "comet.stdout.txt").write_text(result.stdout)
    (WORK / "comet.stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(result.stderr + result.stdout)
    output = mgf.with_suffix(".txt")
    with output.open() as handle:
        version_line = next(handle).strip()
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rank1 = {int(row["scan"]): row for row in rows if row["num"] == "1"}
    checks = []
    for case in cases:
        row = rank1[case["scan"]]
        sequence = row["plain_peptide"]
        neutral_column = "calc_neutral_mass"
        observed_decoy = row["protein"].startswith("DECOY_")
        checks.append({"name": case["name"], "scan": case["scan"], "expected_peptide": case["sequence"],
                       "observed_peptide": sequence, "modified_peptide": row["modified_peptide"],
                       "peptide_matches_with_IL_equivalence": sequence.replace("I", "L") == case["sequence"].replace("I", "L"),
                       "expected_decoy": case["decoy"], "observed_decoy": observed_decoy,
                       "decoy_class_matches": observed_decoy == case["decoy"],
                       "calculated_neutral_mass": float(row[neutral_column]), "expected_neutral_mass": case["neutral_mass"],
                       "mass_matches": abs(float(row[neutral_column]) - case["neutral_mass"]) < 0.001,
                       "xcorr": float(row["xcorr"]), "matched_ions": int(row["ions_matched"]),
                       "total_ions": int(row["ions_total"]), "expectation_value": float(row["e-value"])})
    assert 8 not in rank1, "Out-of-range precursor unexpectedly searched"
    for check in checks:
        assert check["peptide_matches_with_IL_equivalence"] and check["decoy_class_matches"] and check["mass_matches"], check
        assert check["matched_ions"] >= 10, check
    report = {"started_at_utc": began, "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "binary_sha256": sha(binary), "args": args, "exit_code": result.returncode,
              "output_version": version_line, "checks": checks, "outside_mass_range_excluded": True,
              "all_checks_pass": True, "public_or_candidate_spectra_used": False,
              "qualification_only_not_biological_or_FDR_validation": True,
              "files": {str(path.relative_to(ROOT)): sha(path) for path in (database, mgf, parameters, expected, output, WORK / "comet.stdout.txt", WORK / "comet.stderr.txt")}}
    destination = ROOT / "reports/ubash3a-proteomics-comet-qualification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
