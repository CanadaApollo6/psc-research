"""Classify the frozen TRAILS raw-read audit, retaining partial and failed support."""

import argparse
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path

from Bio.Align import PairwiseAligner
import pysam

from analyze_ubash3a_upf1_reads import covered, exclusions, structure

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-raw-read-audit"
OUT = ROOT / "data/derived/ubash3a-raw-read-audit"
PLAN = ROOT / "config/ubash3a-raw-read-plan.json"
NORMAL, PLUS29 = (42434954, 42437488), (42434983, 42437488)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path, rows, empty_fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else empty_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def genomic_transcript_strand(read):
    """Minimap2's ts tag is query-relative, unlike the desired genomic strand."""
    if not read.has_tag("ts"):
        return "unknown"
    tag = read.get_tag("ts")
    if tag not in ("+", "-"):
        return "unknown"
    return ("-" if tag == "+" else "+") if read.is_reverse else tag


def event_of(chain):
    return "plus29" if PLUS29 in chain else "normal" if NORMAL in chain else "other"


def alignment_identity(read):
    denominator = sum(n for op, n in read.cigartuples or [] if op in (0, 1, 2, 7, 8))
    return 1 - read.get_tag("NM") / denominator if denominator and read.has_tag("NM") else None


def end_run(sequence, base, right=True):
    sequence = sequence.upper()
    return len(sequence) - len(sequence.rstrip(base) if right else sequence.lstrip(base))


def classify(read, query_has_supplementary, anchor_min=20):
    chain, anchors, blocks = structure(read.reference_start, read.cigartuples or [])
    event = event_of(chain)
    reasons = exclusions(read.flag, read.mapping_quality)
    if query_has_supplementary or read.has_tag("SA"):
        reasons.append("query_has_supplementary_or_SA")
    strand = genomic_transcript_strand(read)
    if strand == "-":
        reasons.append("inferred_negative_genomic_transcript_strand")
    if event != "other":
        index = chain.index(PLUS29 if event == "plus29" else NORMAL)
        if min(anchors[index]) < anchor_min:
            reasons.append(f"target_anchor_below_{anchor_min}")
    else:
        index = None
    return {"chain": chain, "anchors": anchors, "blocks": blocks, "event": event,
            "index": index, "exclusions": reasons, "strand": strand,
            "qualified": event != "other" and not reasons}


def terminal_bases(blocks, exon):
    a, b = exon
    return sum(max(0, min(y, b) - max(x, a) + 1) for x, y in blocks)


def reference_matches(info, refs, anchor_min):
    full, downstream, branches = [], [], []
    if not info["qualified"]:
        return full, downstream, branches
    chain, anchors, index = info["chain"], info["anchors"], info["index"]
    tail = chain[index+1:]
    for ref in refs:
        if ref["event"] != info["event"]:
            continue
        target = PLUS29 if info["event"] == "plus29" else NORMAL
        suffix = ref["chain"][ref["chain"].index(target)+1:]
        terminal_ok = terminal_bases(info["blocks"], ref["exons"][-1]) >= 20
        if tail == suffix and all(min(a) >= anchor_min for a in anchors[index+1:]) and terminal_ok:
            downstream.append(ref["target_label"])
        if chain == ref["chain"] and all(min(a) >= anchor_min for a in anchors) and terminal_ok and covered(info["blocks"], 42403946, 42403948):
            full.append(ref["target_label"])
    discriminators = {"boundary_B": (42442596, 42447057),
                      "canonical_downstream_C_U": (42442596, 42443312),
                      "alternative_ending_T": (42437580, 42437855)}
    for label, junction in discriminators.items():
        if junction in tail and min(anchors[chain.index(junction)]) >= anchor_min:
            # Require each intervening junction to have qualifying support too.
            stop = chain.index(junction)
            if all(min(a) >= anchor_min for a in anchors[index+1:stop+1]):
                branches.append(label)
    return full, downstream, branches


def local_template_scores(read, genome):
    """Independent affine-gap comparison at shared flanks; no calibrated probability."""
    pairs = [(q, r+1) for q, r in read.get_aligned_pairs(matches_only=True)
             if q is not None and r is not None]
    left = [(q, r) for q, r in pairs if 42434904 <= r <= 42434954]
    right = [(q, r) for q, r in pairs if 42437488 <= r <= 42437537]
    if len(left) < 10 or len(right) < 10 or not read.query_sequence:
        return {"dp_status": "insufficient_shared_flanks", "dp_normal_score": None,
                "dp_plus29_score": None, "dp_plus29_minus_normal": None,
                "dp_query_length": None, "dp_flank_start": None, "dp_flank_end": None}
    q0, r0 = min(left, key=lambda p: p[1])
    q1, r1 = max(right, key=lambda p: p[1])
    query = read.query_sequence[q0:q1+1].upper()
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                             open_gap_score=-5, extend_gap_score=-1)
    templates = {"normal": genome.fetch("21", r0-1, NORMAL[0]).upper() + genome.fetch("21", NORMAL[1]-1, r1).upper(),
                 "plus29": genome.fetch("21", r0-1, PLUS29[0]).upper() + genome.fetch("21", PLUS29[1]-1, r1).upper()}
    scores = {key: aligner.score(template, query) for key, template in templates.items()}
    return {"dp_status": "descriptive_global_affine_gap_scores_reference_A_templates",
            "dp_normal_score": scores["normal"], "dp_plus29_score": scores["plus29"],
            "dp_plus29_minus_normal": scores["plus29"]-scores["normal"],
            "dp_query_length": len(query), "dp_flank_start": r0, "dp_flank_end": r1}


def read_refs():
    path = ROOT / "data/derived/ubash3a-cd4-experiment/assay-targets.csv"
    refs = list(csv.DictReader(path.open()))
    for ref in refs:
        ref["chain"] = [tuple(x) for x in json.loads(ref["junction_chain_exon_end_next_start"])]
        ref["exons"] = [tuple(x) for x in json.loads(ref["exons_1based_inclusive"])]
    return refs


def run(output=OUT):
    plan = json.loads(PLAN.read_text())
    descriptor_path = ROOT / "config/ubash3a-raw-read-coding-descriptors.json"
    descriptors = json.loads(descriptor_path.read_text())
    assert descriptors["raw_plan_sha256"] == sha(PLAN)
    assert descriptors["T29_consequences_sha256"] == sha(ROOT / "data/derived/ubash3a-alternative-ending/consequences.csv")
    for path, expected in plan["input_pins"].items():
        assert sha(ROOT / path) == expected
    refs = read_refs()
    observations, summaries, comparisons, architecture_counts = [], [], [], []
    provenance = {}
    with pysam.FastaFile(str(WORK / "GRCh38.primary.fa")) as genome:
        for sample in plan["sample_ids"]:
            folder = WORK / sample
            metadata = json.loads((folder / "alignment.json").read_text())
            assert metadata["plan_sha256"] == sha(PLAN)
            for path, expected in metadata["output_hashes"].items():
                assert sha(folder / path) == expected
            provenance[sample] = metadata
            sensitivity = {}
            with pysam.AlignmentFile(str(folder / "target-sensitivity.bam"), "rb") as bam:
                for query, records in itertools.groupby(bam, key=lambda r: r.query_name):
                    group = list(records)
                    has_supplement = any(r.is_supplementary or r.has_tag("SA") for r in group)
                    sensitivity[query] = []
                    for read in group:
                        if read.is_unmapped or read.reference_name != "21":
                            continue
                        info = classify(read, has_supplement)
                        if not read.is_secondary and not read.is_supplementary:
                            sensitivity[query].append(info)
            local = []
            with pysam.AlignmentFile(str(folder / "target-groups.bam"), "rb") as bam:
                for query, records in itertools.groupby(bam, key=lambda r: r.query_name):
                    group = list(records)
                    has_supplement = any(r.is_supplementary or r.has_tag("SA") for r in group)
                    primary = next(r for r in group if not r.is_secondary and not r.is_supplementary)
                    raw_sequence = primary.get_forward_sequence() or ""
                    sequence_hash = hashlib.sha256(raw_sequence.upper().encode()).hexdigest()
                    competing = [{"reference": r.reference_name, "start_1based": r.reference_start+1,
                                  "end_1based": r.reference_end, "flag": r.flag, "mapq": r.mapping_quality,
                                  "score": r.get_tag("AS") if r.has_tag("AS") else None}
                                 for r in group if not r.is_unmapped and (r.is_secondary or r.is_supplementary)]
                    for alignment_index, read in enumerate(group):
                        if read.is_unmapped or read.reference_name != "21" or read.reference_start >= 42447684 or read.reference_end < 42403447:
                            continue
                        info = classify(read, has_supplement)
                        info10 = classify(read, has_supplement, 10)
                        full, downstream, branches = reference_matches(info, refs, 20)
                        full10, downstream10, branches10 = reference_matches(info10, refs, 10)
                        terminal_window = genome.fetch("21", read.reference_end, read.reference_end+20).upper()
                        scores = local_template_scores(read, genome) if info["event"] != "other" else local_template_scores_empty()
                        changed = sensitivity.get(query, [])
                        consistent = any(s["event"] == info["event"] and s["qualified"] for s in changed) if info["event"] != "other" else False
                        row = {"sample": sample, "query_name": query, "alignment_index_in_query": alignment_index,
                               "flag": read.flag, "mapq": read.mapping_quality, "raw_ts_tag": read.get_tag("ts") if read.has_tag("ts") else "",
                               "alignment_strand": "-" if read.is_reverse else "+", "inferred_genomic_transcript_strand": info["strand"],
                               "start_1based": read.reference_start+1, "end_1based": read.reference_end,
                               "cigar": read.cigarstring, "alignment_identity": alignment_identity(read),
                               "query_length": len(raw_sequence), "soft_clipped_bases": sum(n for op,n in read.cigartuples or [] if op == 4),
                               "event": info["event"], "qualified_anchor20": info["qualified"], "qualified_anchor10": info10["qualified"],
                               "exclusions_anchor20": ";".join(info["exclusions"]),
                               "junctions_json": json.dumps(info["chain"], separators=(",", ":")),
                               "anchors_json": json.dumps(info["anchors"], separators=(",", ":")),
                               "full_splice_chain_and_start_matches": ";".join(full), "complete_downstream_chain_matches": ";".join(downstream),
                               "discriminating_branch_matches": ";".join(branches),
                               "full_chain_matches_anchor10": ";".join(full10), "downstream_matches_anchor10": ";".join(downstream10),
                               "branch_matches_anchor10": ";".join(branches10),
                               "shared_start_coordinates_aligned": covered(info["blocks"], *descriptors["T29_assumed_start_1based"]),
                               "T29_stop_coordinates_aligned": covered(info["blocks"], *descriptors["T29_conditional_stop_1based"]),
                               "previous_plus29_stop_coordinates_aligned": covered(info["blocks"], *descriptors["previous_plus29_conditional_stop_1based"]),
                               "T29_catalog_end_coordinate_aligned": covered(info["blocks"], descriptors["T29_catalog_end_1based"], descriptors["T29_catalog_end_1based"]),
                               "sensitivity_same_event_qualified": consistent,
                               "sensitivity_primary_events": ";".join(s["event"] for s in changed),
                               "sequence_sha256": sequence_hash, "competing_alignments_json": json.dumps(competing, separators=(",", ":")),
                               "raw_left_A_run": end_run(raw_sequence, "A", False), "raw_left_T_run": end_run(raw_sequence, "T", False),
                               "raw_right_A_run": end_run(raw_sequence, "A"), "raw_right_T_run": end_run(raw_sequence, "T"),
                               "genomic_next20_A_count": terminal_window.count("A"),
                               "genomic_next20_A_run_max": max([len(s) for s in terminal_window.replace("C", " ").replace("G", " ").replace("T", " ").replace("N", " ").split()] or [0]),
                               **scores}
                        local.append(row)
            observations.extend(local)
            summary = {"sample": sample, **metadata["stats"], "locus_alignment_rows": len(local),
                       "normal_exact_all_alignments": sum(r["event"] == "normal" for r in local),
                       "plus29_exact_all_alignments": sum(r["event"] == "plus29" for r in local)}
            for event in ["normal", "plus29"]:
                for anchor in [20, 10]:
                    subset = [r for r in local if r["event"] == event and r[f"qualified_anchor{anchor}"]]
                    assert len({r["query_name"] for r in subset}) == len(subset)
                    summary[f"{event}_qualified_anchor{anchor}"] = len(subset)
                    summary[f"{event}_distinct_exact_sequences_anchor{anchor}"] = len({r["sequence_sha256"] for r in subset})
                summary[f"{event}_anchor20_also_sensitivity_qualified"] = sum(r["event"] == event and r["qualified_anchor20"] and r["sensitivity_same_event_qualified"] for r in local)
            summaries.append(summary)
            for ref in refs:
                for anchor in [20, 10]:
                    full_key = "full_splice_chain_and_start_matches" if anchor == 20 else "full_chain_matches_anchor10"
                    downstream_key = "complete_downstream_chain_matches" if anchor == 20 else "downstream_matches_anchor10"
                    matches = [r for r in local if ref["target_label"] in r[full_key].split(";")]
                    architecture_counts.append({"sample": sample, "target_label": ref["target_label"], "source_model": ref["source_model"],
                                                "anchor_min": anchor, "full_splice_chain_and_start_reads": len(matches),
                                                "full_chain_distinct_exact_sequences": len({r["sequence_sha256"] for r in matches}),
                                                "complete_downstream_chain_reads": sum(ref["target_label"] in r[downstream_key].split(";") for r in local)})
            for row in local:
                if row["event"] != "other" and not row["flag"] & (256 | 2048):
                    comparisons.append({key: row[key] for key in ["sample", "query_name", "event", "qualified_anchor20", "qualified_anchor10", "sensitivity_same_event_qualified", "sensitivity_primary_events", "dp_status", "dp_normal_score", "dp_plus29_score", "dp_plus29_minus_normal"]})
    output.mkdir(parents=True, exist_ok=True)
    files = {"raw-locus-alignments.csv": (observations, ["sample", "query_name", "event"]),
             "raw-read-summary.csv": (summaries, ["sample"]),
             "raw-architecture-counts.csv": (architecture_counts, ["sample", "target_label", "anchor_min"]),
             "raw-alignment-sensitivity.csv": (comparisons, ["sample", "query_name", "event"])}
    for name, (rows, empty_fields) in files.items():
        write_csv(output / name, rows, empty_fields)
    audit = {"plan_sha256": sha(PLAN), "script_sha256": sha(Path(__file__)),
             "additive_coding_descriptors_sha256": sha(descriptor_path),
             "shared_CIGAR_parser_sha256": sha(ROOT / "scripts/analyze_ubash3a_upf1_reads.py"),
             "samples": provenance, "files_sha256": {name: sha(output / name) for name in files},
             "interpretation": "Amplified cDNA reads from one donor; exact sequences are not original-molecule identities; no inferred genotype, RNA survival, or protein detection."}
    (output / "raw-read-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"summaries": summaries, "architecture_counts": architecture_counts}, indent=2))


def local_template_scores_empty():
    return {"dp_status": "not_a_target_junction", "dp_normal_score": None,
            "dp_plus29_score": None, "dp_plus29_minus_normal": None,
            "dp_query_length": None, "dp_flank_start": None, "dp_flank_end": None}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    run(parser.parse_args().output_dir)
