#!/usr/bin/env python3
"""Summarize the complete frozen search with explicit target-decoy competition."""
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from analyze_ubash3a_proteomics_specificity import sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-proteomics-spectra"


def proteins(value):
    return value.split(",") if value else []


def has_decoy(value):
    return any(protein.startswith("DECOY_") for protein in proteins(value))


def query_scan_map(path):
    """Use explicit native IDs; Comet's exported scan can be an mzML ordinal."""
    result, reverse = {}, {}
    for _, query in ET.iterparse(path, events=("end",)):
        if query.tag.rsplit("}", 1)[-1] != "spectrum_query":
            continue
        native = re.findall(r"(?:^|\s)scan=(\d+)(?:$|\s)", query.attrib["spectrumNativeID"])
        assert len(native) == 1
        original, exported = native[0], query.attrib["start_scan"]
        assert exported == query.attrib["end_scan"]
        key = exported, query.attrib["assumed_charge"]
        assert key not in result
        assert reverse.setdefault(original, exported) == exported
        result[key] = original
        query.clear()
    return result


def best_per_scan(rows):
    """Only rank-one query winners contribute; retain conservative reported ties."""
    queries = defaultdict(list)
    for row in rows:
        queries[(row["file"], int(row["scan"]), int(row["charge"]))].append(row)
    scans = defaultdict(list)
    for key, query in queries.items():
        first = [row for row in query if int(row["num"]) == 1]
        assert first, key
        # The exported rank is dense XCorr rank: multiple assignments can all
        # have rank 1. Keep one deterministic representative, decoy first.
        winner = dict(min(first, key=lambda row: (
            not has_decoy(row["protein"]), row["plain_peptide"],
            int(row.get("query_hit_index", 1)))))
        # Comet ranks within each spectrum by XCorr. Reported equal top scores
        # containing a decoy remain conservative decoy wins for error estimation.
        tied = [row for row in query if float(row["xcorr"]) == float(winner["xcorr"])]
        winner["is_decoy"] = any(has_decoy(row["protein"]) for row in tied)
        winner["reported_top_score_tie_count"] = len(tied)
        winner["rank_one_assignment_count"] = len(first)
        winner["target_decoy_top_score_tie"] = winner["is_decoy"] and any(not has_decoy(row["protein"]) for row in tied)
        scans[key[:2]].append(winner)
    winners = []
    for key, query_winners in sorted(scans.items()):
        best_score = min(float(row["e-value"]) for row in query_winners)
        tied = [row for row in query_winners if float(row["e-value"]) == best_score]
        chosen = dict(sorted(tied, key=lambda row: (not row["is_decoy"], -float(row["xcorr"]), row["plain_peptide"], int(row["charge"])))[0])
        chosen["is_decoy"] = any(row["is_decoy"] for row in tied)
        chosen["charge_queries_for_scan"] = len(query_winners)
        chosen["score"] = best_score
        chosen["peptide_IL"] = chosen["plain_peptide"].replace("I", "L")
        winners.append(chosen)
    return winners


def q_values(rows):
    """Tie-grouped conservative (D+1)/T q-values; return copies in score order."""
    ordered = sorted((dict(row) for row in rows), key=lambda row: (row["score"], row["file"], int(row["scan"])))
    groups = []
    target, decoy, index = 0, 0, 0
    while index < len(ordered):
        end = index
        while end < len(ordered) and ordered[end]["score"] == ordered[index]["score"]:
            decoy += bool(ordered[end]["is_decoy"])
            target += not ordered[end]["is_decoy"]
            end += 1
        rate = min(1.0, (decoy + 1) / target) if target else 1.0
        groups.append((index, end, rate, target, decoy))
        index = end
    running = 1.0
    for start, end, rate, target, decoy in reversed(groups):
        running = min(running, rate)
        for row in ordered[start:end]:
            row.update(q_value=running, threshold_fdr=rate, targets_at_threshold=target, decoys_at_threshold=decoy)
    return ordered


def unique_peptide_winners(winners):
    grouped = defaultdict(list)
    for row in winners:
        grouped[row["peptide_IL"]].append(row)
    selected = []
    for peptide, rows in grouped.items():
        best = min(row["score"] for row in rows)
        tied = [row for row in rows if row["score"] == best]
        chosen = dict(sorted(tied, key=lambda row: (not row["is_decoy"], row["file"], int(row["scan"])))[0])
        chosen["is_decoy"] = any(row["is_decoy"] for row in tied)
        chosen["spectrum_wins_for_peptide"] = len(rows)
        chosen["files_for_peptide"] = len({row["file"] for row in rows})
        selected.append(chosen)
    return selected


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    execution_path = ROOT / "reports/ubash3a-proteomics-search-execution.json"
    execution = json.loads(execution_path.read_text())
    assert execution["complete_experiment"] and execution["file_count"] == 14
    plan = json.loads((ROOT / "config/ubash3a-proteomics-specificity-plan.json").read_text())
    candidates = plan["candidates"]
    gene_members = list(csv.DictReader((ROOT / "data/derived/ubash3a-proteomics-specificity/ubash3a-search-database-members.csv").open()))
    gene_ids = {row["database_id"] for row in gene_members}
    human_gene_ids = {row["database_id"] for row in gene_members if row["source_group"] == "human_reference"}
    t29_ids = {row["database_id"] for row in gene_members if row["source_group"] == "conditional_T29"}
    all_rows, source_counts, loaded_counts = [], {}, {}
    for file in execution["files"]:
        paths = [path for path in file["output_sha256"] if path.endswith(".txt")]
        assert len(paths) == 1
        path = ROOT / paths[0]
        assert sha256(path) == file["output_sha256"][paths[0]]
        logs = [name for name in file["log_sha256"] if name.endswith(".stdout.txt")]
        assert len(logs) == 1
        log_path = ROOT / logs[0]
        assert sha256(log_path) == file["log_sha256"][logs[0]]
        batches = [int(value) for value in re.findall(r"Load spectra:\s*(\d+)", log_path.read_text())]
        assert batches, log_path
        loaded_counts[path.stem] = sum(batches)
        with path.open() as handle:
            assert next(handle).startswith("CometVersion 2026.02 rev. 2 (6edec91)")
            rows = list(csv.DictReader(handle, delimiter="\t"))
        xml_names = [name for name in file["output_sha256"] if name.endswith(".pep.xml")]
        assert len(xml_names) == 1
        xml_path = ROOT / xml_names[0]
        assert sha256(xml_path) == file["output_sha256"][xml_names[0]]
        native_scans = query_scan_map(xml_path)
        hit_counts = Counter()
        for row in rows:
            row["file"] = path.stem
            key = row["scan"], row["charge"]
            hit_counts[key] += 1
            row["query_hit_index"] = hit_counts[key]
            row["comet_scan"] = row["scan"]
            row["scan"] = native_scans[key]
            assert all(math.isfinite(float(row[field])) for field in ("e-value", "xcorr", "exp_neutral_mass", "calc_neutral_mass"))
            assert float(row["e-value"]) >= 0
        all_rows.extend(rows)
        source_counts[path.stem] = len(rows)
    rank_one_counts = Counter((row["file"], row["scan"], row["charge"]) for row in all_rows if row["num"] == "1")
    for row in all_rows:
        row["rank_one_assignments_in_query"] = rank_one_counts[(row["file"], row["scan"], row["charge"])]
    winners = q_values(best_per_scan(all_rows))
    peptide_winners = q_values(unique_peptide_winners(winners))
    peptide_q = {row["peptide_IL"]: row["q_value"] for row in peptide_winners}
    by_scan = {(row["file"], row["scan"]): row for row in winners}
    for row in winners:
        row["peptide_q_value"] = peptide_q[row["peptide_IL"]]
        row["passes_global_spectrum_q01"] = not row["is_decoy"] and row["q_value"] <= 0.01
        row["passes_global_peptide_q01"] = not row["is_decoy"] and row["peptide_q_value"] <= 0.01
        row["nonstandard_residue_sequence"] = any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in row["plain_peptide"])
    # Full background PSM/peptide tables remain locally reproducible large outputs.
    local = ROOT / "work/ubash3a-proteomics-spectra/analysis"
    local.mkdir(parents=True, exist_ok=True)
    write_csv(local / "all-spectrum-winners.csv", winners)
    write_csv(local / "all-peptide-winners.csv", peptide_winners)
    gene_rows = []
    for row in winners:
        refs = proteins(row["protein"])
        target_refs = [ref.removeprefix("DECOY_") for ref in refs]
        if set(target_refs) & gene_ids:
            gene_rows.append({**row, "all_reported_protein_refs_in_UBASH3A_set": set(target_refs) <= gene_ids,
                             "has_human_reference_UBASH3A_protein_member": bool(set(target_refs) & human_gene_ids),
                             "has_conditional_T29_protein_member": bool(set(target_refs) & t29_ids),
                             "interpretation": "Gene-associated peptide assignment; source isoform and translation not established by shared sequence"})
    gene_fields = list(winners[0]) + ["all_reported_protein_refs_in_UBASH3A_set", "has_human_reference_UBASH3A_protein_member", "has_conditional_T29_protein_member", "interpretation"]
    write_csv(OUT / "ubash3a-spectrum-winners.csv", gene_rows, gene_fields)
    nominations = []
    for row in all_rows:
        matching = [p for p in candidates if p.replace("I", "L") == row["plain_peptide"].replace("I", "L")]
        for candidate in matching:
            winner = by_scan[(row["file"], row["scan"])]
            exact_win = row["num"] == "1" and row["charge"] == winner["charge"] and row["query_hit_index"] == winner["query_hit_index"]
            nominations.append({"candidate": candidate, **row,
                                "is_primary_scan_winner": exact_win,
                                "scan_winner_peptide": winner["plain_peptide"], "scan_winner_modified_peptide": winner["modified_peptide"],
                                "scan_winner_proteins": winner["protein"], "scan_winner_expectation": winner["e-value"],
                                "scan_winner_decoy": winner["is_decoy"], "scan_winner_q_value": winner["q_value"],
                                "candidate_global_spectrum_q_value": winner["q_value"] if exact_win else "",
                                "candidate_global_peptide_q_value": peptide_q.get(candidate.replace("I", "L"), ""),
                                "novel_class_q_value": "not_estimable", "protein_detection_claim": False})
    nomination_fields = ["candidate"] + list(all_rows[0]) + ["is_primary_scan_winner", "scan_winner_peptide", "scan_winner_modified_peptide", "scan_winner_proteins", "scan_winner_expectation", "scan_winner_decoy", "scan_winner_q_value", "candidate_global_spectrum_q_value", "candidate_global_peptide_q_value", "novel_class_q_value", "protein_detection_claim"]
    write_csv(OUT / "candidate-all-reported-assignments.csv", nominations, nomination_fields)
    candidate_rows = []
    for candidate in candidates:
        rows = [r for r in nominations if r["candidate"] == candidate]
        winning = [r for r in rows if r["is_primary_scan_winner"] and not r["scan_winner_decoy"]]
        passing = [r for r in winning if r["candidate_global_spectrum_q_value"] <= 0.01]
        candidate_rows.append({"candidate": candidate, "reported_assignments_any_rank": len(rows),
                               "distinct_spectra_any_rank": len({(r["file"], r["scan"]) for r in rows}),
                               "reported_rank_one_assignments": sum(r["num"] == "1" for r in rows),
                               "rank_one_assignments_in_tied_queries": sum(r["num"] == "1" and r["rank_one_assignments_in_query"] > 1 for r in rows),
                               "primary_target_spectrum_wins": len(winning), "global_spectrum_q01_wins": len(passing),
                               "minimum_reported_expectation": min((float(r["e-value"]) for r in rows), default=""),
                               "minimum_winner_spectrum_q": min((r["candidate_global_spectrum_q_value"] for r in winning), default=""),
                               "novel_class_q_value": "not_estimable", "protein_detection_claim": False,
                               "interpretation": "Reported assignment requires specificity/fragment/competitor review; global q is insufficient for novel detection" if rows else "Not among the five reported competitors per query; not evidence of biological absence"})
    write_csv(OUT / "candidate-search-summary.csv", candidate_rows)
    conversion = json.loads((ROOT / "reports/ubash3a-proteomics-conversion.json").read_text())
    file_rows = []
    for source in conversion["files"]:
        name = Path(source["mzml_path"]).stem
        rows = [r for r in winners if r["file"] == name]
        passed = [r for r in rows if r["passes_global_spectrum_q01"]]
        file_rows.append({"file": name, "converted_MS2": source["spectra"], "Comet_loaded_spectrum_queries": loaded_counts[name], "reported_assignments": source_counts[name],
                          "spectra_with_reported_winner": len(rows), "MS2_without_reported_winner": source["spectra"] - len(rows),
                          "target_spectrum_wins": sum(not r["is_decoy"] for r in rows),
                          "decoy_spectrum_wins": sum(r["is_decoy"] for r in rows),
                          "passing_global_spectrum_q01": len(passed),
                          "unique_IL_peptides_among_passing_spectra": len({r["peptide_IL"] for r in passed}),
                          "UBASH3A_associated_passing_spectra": sum(r["passes_global_spectrum_q01"] for r in gene_rows if r["file"] == name)})
        assert file_rows[-1]["MS2_without_reported_winner"] >= 0
    write_csv(OUT / "file-search-summary.csv", file_rows)
    result = {"complete_experiment": True, "files": 14, "converted_MS2": conversion["total_MS2_spectra"],
              "Comet_loaded_spectrum_queries": sum(loaded_counts.values()),
              "all_reported_assignments": len(all_rows), "spectrum_winners": len(winners),
              "target_spectrum_wins": sum(not r["is_decoy"] for r in winners),
              "decoy_spectrum_wins": sum(r["is_decoy"] for r in winners),
              "spectrum_q01_target_wins": sum(r["passes_global_spectrum_q01"] for r in winners),
              "unique_IL_peptide_q01_targets": sum(not r["is_decoy"] and r["q_value"] <= .01 for r in peptide_winners),
              "UBASH3A_associated_spectrum_q01_wins": sum(r["passes_global_spectrum_q01"] for r in gene_rows),
              "UBASH3A_associated_unique_peptides_in_passing_spectra": len({r["peptide_IL"] for r in gene_rows if r["passes_global_spectrum_q01"]}),
              "ordinary_reference_UBASH3A_spectrum_q01_wins": sum(r["passes_global_spectrum_q01"] and r["has_human_reference_UBASH3A_protein_member"] for r in gene_rows),
              "ordinary_reference_UBASH3A_unique_peptides_passing_both_q01_levels": len({r["peptide_IL"] for r in gene_rows if r["passes_global_spectrum_q01"] and r["passes_global_peptide_q01"] and r["has_human_reference_UBASH3A_protein_member"]}),
              "candidate_reported_assignments": len(nominations), "candidate_summary": candidate_rows,
              "ranking": "best rank-one query winner per file/scan, ranked by expectation value; conservative target-decoy ties",
              "error_estimator": "tie-grouped (D+1)/T and monotone minimum; spectrum and unique unmodified I/L-equivalent peptide summaries",
              "novel_class_error_rate": "Not reliably estimable from four overlapping/nested targets",
              "technical_scope": "Only first five assignments per query are exported; their dense XCorr ranks can repeat. A non-reported target is not an exhaustive targeted spectrum test",
              "counting_scope": "Original scan IDs come from the explicit pepXML spectrumNativeID, with the exported Comet scan preserved separately. Converted MS2 scans, engine-loaded spectrum queries, and scans with reported winners are separate. Mass/peak filters or absent competitors can leave converted scans without reported winners; multiple charge queries and tied rank-one assignments are collapsed for FDR.",
              "measurement_scope": "One complete public experiment; repeated/pooled runs are not independent donors; no protein-absence inference",
              "search_execution_sha256": sha256(execution_path),
              "large_output_hashes": {str(path.relative_to(ROOT)): sha256(path) for path in (local / "all-spectrum-winners.csv", local / "all-peptide-winners.csv")}}
    (OUT / "search-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
