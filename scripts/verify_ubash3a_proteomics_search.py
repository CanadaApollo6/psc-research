#!/usr/bin/env python3
"""Independent pepXML/pandas checks of the completed spectrum analysis."""
import csv
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from Bio import SeqIO
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-proteomics-spectra"
NS = {"p": "http://regis-web.systemsbiology.net/pepXML"}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_database():
    lock = json.loads((ROOT / "config/ubash3a-proteomics-reference-lock.json").read_text())
    db = json.loads((ROOT / "config/ubash3a-proteomics-search-database.json").read_text())
    db_path = ROOT / db["database_path"]
    assert sha(db_path) == db["database_sha256"]
    actual = {}
    for record in SeqIO.parse(db_path, "fasta"):
        sequence = str(record.seq)
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        assert record.id == "P_" + sequence_hash[:24]
        assert record.id not in actual
        actual[record.id] = sequence_hash
    expected_pairs = set()
    expected_db = set()
    for source in lock["sources"]:
        path = ROOT / source["local_path"]
        assert sha(path) == source["sha256"]
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequence_hash = hashlib.sha256(str(record.seq).upper().encode()).hexdigest()
                identifier = "P_" + sequence_hash[:24]
                assert actual[identifier] == sequence_hash
                expected_pairs.add((source["set_id"], record.id, identifier, sequence_hash, record.description, source["source_group"]))
                expected_db.add(identifier)
    assert set(actual) == expected_db
    observed_pairs = set()
    with gzip.open(ROOT / db["membership_path"], "rt") as handle:
        for line in handle:
            row = json.loads(line)
            key = (row["set_id"], row["source_record_id"], row["database_id"], row["sequence_sha256"], row["source_header"], row["source_group"])
            assert key not in observed_pairs
            observed_pairs.add(key)
    assert expected_pairs == observed_pairs
    assert len(actual) == db["unique_exact_protein_sequences"] and len(expected_pairs) == db["source_records"]
    return {"all_source_records_preserved": len(expected_pairs), "unique_sequences_verified": len(actual)}


def parse_pepxml(path, file):
    rows, raw_unprefixed_decoy_alternatives = [], 0
    for _, query in ET.iterparse(path, events=("end",)):
        if query.tag.rsplit("}", 1)[-1] != "spectrum_query":
            continue
        comet_scan, charge = int(query.attrib["start_scan"]), int(query.attrib["assumed_charge"])
        native_match = re.search(r"(?:^|\s)scan=([0-9]+)(?=\s|$)", query.attrib["spectrumNativeID"])
        assert native_match is not None
        scan = int(native_match.group(1))
        assert query.attrib["start_scan"] == query.attrib["end_scan"]
        for hit_index, hit in enumerate(query.findall("p:search_result/p:search_hit", NS), 1):
            scores = {score.attrib["name"]: float(score.attrib["value"]) for score in hit.findall("p:search_score", NS)}
            protein = hit.attrib["protein"]
            references = [protein]
            for alt in hit.findall("p:alternative_protein", NS):
                ref = alt.attrib["protein"]
                # Preserve XML exactly: alternative names lose decoy prefixes,
                # including decoys attached to a target parent. PIN supplies
                # their independently exported labels later.
                if protein.startswith("DECOY_") and not ref.startswith("DECOY_"):
                    raw_unprefixed_decoy_alternatives += 1
                references.append(ref)
            rows.append({"file": file, "scan": scan, "comet_scan": comet_scan,
                         "charge": charge, "rank": int(hit.attrib["hit_rank"]), "query_hit_index": hit_index,
                         "peptide": hit.attrib["peptide"], "expect": scores["expect"], "xcorr": scores["xcorr"],
                         "exp_mass": float(query.attrib["precursor_neutral_mass"]), "calc_mass": float(hit.attrib["calc_neutral_pep_mass"]),
                         "matched_ions": int(hit.attrib["num_matched_ions"]), "total_ions": int(hit.attrib["tot_num_ions"]),
                         "protein_refs": tuple(sorted(references)), "hit_decoy": any(ref.startswith("DECOY_") for ref in references)})
        query.clear()
    return rows, raw_unprefixed_decoy_alternatives


def parse_pin(path):
    rows = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        protein_index = header.index("Proteins")
        assert protein_index == len(header) - 1
        for values in reader:
            if not values:
                continue
            assert len(values) > protein_index
            fields = dict(zip(header[:protein_index], values[:protein_index]))
            match = re.search(r"_([0-9]+)_([0-9]+)_([0-9]+)$", fields["SpecId"])
            assert match is not None
            key = tuple(map(int, match.groups()))
            assert key[0] == int(fields["ScanNr"]) and key not in rows
            references = tuple(sorted(values[protein_index:]))
            rows[key] = {"protein_refs": references,
                         "hit_decoy": any(ref.startswith("DECOY_") for ref in references),
                         "PIN_reported_label": int(fields["Label"]),
                         "modified_peptide": fields["Peptide"],
                         "MH_exp_mass": float(fields["ExpMass"]),
                         "MH_calc_mass": float(fields["CalcMass"]),
                         "xcorr": float(fields["Xcorr"])}
    return rows


def independent_q(frame):
    counts = frame.groupby("expect", sort=True)["is_decoy"].agg(["count", "sum"])
    cumulative_d = counts["sum"].cumsum()
    cumulative_t = (counts["count"] - counts["sum"]).cumsum()
    rates = ((cumulative_d + 1) / cumulative_t.replace(0, np.nan)).fillna(1).clip(upper=1)
    qvalues = rates.iloc[::-1].cummin().iloc[::-1]
    return frame.assign(independent_q=frame["expect"].map(qvalues))


def check_search_outputs():
    execution_path = ROOT / "reports/ubash3a-proteomics-search-execution.json"
    execution = json.loads(execution_path.read_text())
    assert execution["complete_experiment"] and len(execution["files"]) == 14
    rows, alternative_count, scalar_comparisons = [], 0, 0
    omitted_prefix_count, mixed_class_rows = 0, 0
    for file in execution["files"]:
        paths = file["output_sha256"]
        xml_path = next(ROOT / path for path in paths if path.endswith(".pep.xml"))
        txt_path = next(ROOT / path for path in paths if path.endswith(".txt"))
        pin_path = next(ROOT / path for path in paths if path.endswith(".pin"))
        assert sha(xml_path) == paths[str(xml_path.relative_to(ROOT))]
        assert sha(txt_path) == paths[str(txt_path.relative_to(ROOT))]
        assert sha(pin_path) == paths[str(pin_path.relative_to(ROOT))]
        xml_rows, n_alt = parse_pepxml(xml_path, txt_path.stem)
        alternative_count += n_alt
        pin_rows = parse_pin(pin_path)
        with txt_path.open() as handle:
            next(handle)
            txt, hit_counts = {}, {}
            for row in csv.DictReader(handle, delimiter="\t"):
                key = int(row["scan"]), int(row["charge"])
                hit_counts[key] = hit_counts.get(key, 0) + 1
                txt[(*key, hit_counts[key])] = row
        assert len(xml_rows) == len(txt) == len(pin_rows)
        for row in xml_rows:
            key = row["comet_scan"], row["charge"], row["query_hit_index"]
            other, pin = txt[key], pin_rows[key]
            assert Counter(ref.removeprefix("DECOY_") for ref in row["protein_refs"]) == Counter(ref.removeprefix("DECOY_") for ref in pin["protein_refs"])
            omitted = sum(ref.startswith("DECOY_") for ref in pin["protein_refs"]) - sum(ref.startswith("DECOY_") for ref in row["protein_refs"])
            assert omitted >= 0
            omitted_prefix_count += omitted
            mixed_class_rows += pin["hit_decoy"] and any(not ref.startswith("DECOY_") for ref in pin["protein_refs"])
            row["protein_refs"] = pin["protein_refs"]
            row["hit_decoy"] = pin["hit_decoy"]
            assert pin["modified_peptide"] == other["modified_peptide"]
            assert abs(pin["MH_exp_mass"] - row["exp_mass"] - 1.007276466621) <= 0.000002
            assert abs(pin["MH_calc_mass"] - row["calc_mass"] - 1.007276466621) <= 0.000002
            assert abs(pin["xcorr"] - row["xcorr"]) <= 0.000051
            assert row["rank"] == int(other["num"])
            assert row["peptide"] == other["plain_peptide"]
            assert row["protein_refs"] == tuple(sorted(other["protein"].split(",")))
            for key, field in [("expect", "e-value"), ("xcorr", "xcorr"), ("exp_mass", "exp_neutral_mass"), ("calc_mass", "calc_neutral_mass"), ("matched_ions", "ions_matched"), ("total_ions", "ions_total")]:
                assert row[key] == float(other[field]), (txt_path.name, key, row[key], other[field])
            scalar_comparisons += 8
        rows.extend(xml_rows)
    frame = pd.DataFrame(rows)
    query_keys = ["file", "scan", "charge"]
    rank_one = frame.loc[frame["rank"] == 1].copy()
    first = rank_one.sort_values(query_keys + ["hit_decoy", "peptide", "query_hit_index"],
                                 ascending=[True, True, True, False, True, True]).drop_duplicates(query_keys).copy()
    assert not first.duplicated(query_keys).any()
    joined = frame.merge(first[query_keys + ["xcorr"]].rename(columns={"xcorr": "winner_xcorr"}), on=query_keys, validate="many_to_one")
    decoy_tie = joined.loc[joined["xcorr"] == joined["winner_xcorr"]].groupby(query_keys)["hit_decoy"].any().rename("is_decoy")
    first = first.merge(decoy_tie, on=query_keys, validate="one_to_one")
    # Sort implements the same declared tie rule using independent tabular ops.
    winners = first.sort_values(["file", "scan", "expect", "is_decoy", "xcorr", "peptide", "charge"],
                               ascending=[True, True, True, False, False, True, True]).drop_duplicates(["file", "scan"]).copy()
    winners["peptide_IL"] = winners["peptide"].str.replace("I", "L", regex=False)
    winners = independent_q(winners)
    peptides = winners.sort_values(["peptide_IL", "expect", "is_decoy", "file", "scan"],
                                   ascending=[True, True, False, True, True]).drop_duplicates("peptide_IL").copy()
    peptides = independent_q(peptides)
    original_path = ROOT / "work/ubash3a-proteomics-spectra/analysis/all-spectrum-winners.csv"
    primary = pd.read_csv(original_path)
    assert len(primary) == len(winners)
    paired = winners.merge(primary, on=["file", "scan"], suffixes=("_xml", "_primary"), validate="one_to_one")
    assert len(paired) == len(primary)
    assert (paired["charge_xml"] == paired["charge_primary"]).all()
    assert (paired["comet_scan_xml"] == paired["comet_scan_primary"]).all()
    assert (paired["query_hit_index_xml"] == paired["query_hit_index_primary"]).all()
    assert (paired["peptide"] == paired["plain_peptide"]).all()
    assert (paired["is_decoy_xml"] == paired["is_decoy_primary"]).all()
    assert np.allclose(paired["independent_q"], paired["q_value"], rtol=0, atol=1e-12)
    peptide_lookup = peptides.set_index("peptide_IL")["independent_q"]
    assert np.allclose(primary["peptide_IL"].map(peptide_lookup), primary["peptide_q_value"], rtol=0, atol=1e-12)
    candidates = json.loads((ROOT / "config/ubash3a-proteomics-specificity-plan.json").read_text())["candidates"]
    expected = set()
    for candidate in candidates:
        subset = frame.loc[frame["peptide"].str.replace("I", "L", regex=False) == candidate.replace("I", "L")]
        expected.update((candidate, row.file, row.scan, row.charge, row.rank, row.query_hit_index) for row in subset.itertuples())
    with (OUT / "candidate-all-reported-assignments.csv").open() as handle:
        actual = {(row["candidate"], row["file"], int(row["scan"]), int(row["charge"]), int(row["num"]), int(row["query_hit_index"])) for row in csv.DictReader(handle)}
    assert actual == expected
    return {"all_pepxml_queries_parsed_through_EOF": True, "TXT_XML_assignment_rows": len(frame),
            "TXT_XML_scalar_comparisons": scalar_comparisons,
            "XML_unprefixed_alternatives_under_decoy_parents": alternative_count,
            "XML_omitted_decoy_prefixes_recovered_from_independent_PIN_output": omitted_prefix_count,
            "mixed_target_decoy_membership_assignment_rows": mixed_class_rows,
            "PIN_protein_memberships_modified_peptides_and_mass_score_fields_checked": len(frame),
            "independent_spectrum_winners": len(winners), "independent_unique_IL_peptides": len(peptides),
            "queries_with_multiple_rank_one_assignments": int((rank_one.groupby(query_keys).size() > 1).sum()),
            "scan_identity": "Original instrument scan from explicit spectrumNativeID; exported Comet scan and hit position retained separately",
            "all_primary_spectrum_and_peptide_q_values_agree_within": 1e-12,
            "all_candidate_assignment_keys_agree": True, "candidate_assignments_verified": len(expected),
            "independent_spectrum_q01_targets": int(((~winners["is_decoy"]) & (winners["independent_q"] <= .01)).sum()),
            "source_execution_sha256": sha(execution_path)}


if __name__ == "__main__":
    database = check_database()
    search = check_search_outputs()
    report = {"all_checks_pass": True, "database": database, "search": search,
              "independence_boundary": "XML supplies independently parsed peptides/scores/native identities and PIN supplies complete target-decoy membership labels missing from XML alternatives; these are compared with TXT and used in a separate pandas q-value implementation. Same spectra, same search engine and same frozen database. This is computational consistency, not independent biological replication or novel-class FDR calibration."}
    (ROOT / "reports/ubash3a-proteomics-search-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
