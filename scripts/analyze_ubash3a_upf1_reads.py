"""Classify the prospectively specified UBASH3A long-read locus extractions."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-hypothesis-test"
ALIGNED = {0, 7, 8}
REFERENCE = {0, 2, 3, 7, 8}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data, schema):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]) if data else schema, lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def structure(start0, cigar):
    """Return exact splice coordinates and actually aligned one-based intervals."""
    position = start0
    chain, anchors, blocks = [], [], []
    for i, (op, length) in enumerate(cigar):
        if not 0 <= op <= 8 or length <= 0:
            raise ValueError("Unsupported or invalid CIGAR operation")
        if op == 3:
            chain.append((position, position + length + 1))
            sides = []
            for step in (-1, 1):
                j, anchor = i + step, 0
                while 0 <= j < len(cigar) and cigar[j][0] in ALIGNED:
                    anchor += cigar[j][1]
                    j += step
                sides.append(anchor)
            anchors.append(tuple(sides))
        if op in ALIGNED:
            blocks.append((position + 1, position + length))
        if op in REFERENCE:
            position += length
    return chain, anchors, blocks


def covered(blocks, start, end):
    return all(any(a <= p <= b for a, b in blocks) for p in range(start, end + 1))


def exclusions(flag, mapq):
    labels = [(4, "unmapped"), (256, "secondary"), (2048, "supplementary"),
              (512, "QC_failed"), (1024, "duplicate_flag")]
    result = [label for mask, label in labels if flag & mask]
    if mapq < 20:
        result.append("MAPQ_below_20")
    return result


def run(output=OUT):
    plan_path = ROOT / "config/ubash3a-upf1-read-test-plan.json"
    plan = json.loads(plan_path.read_text())
    manifest = json.loads((ROOT / "config/ubash3a-upf1-read-sources.json").read_text())
    assert manifest["plan_sha256"] == sha(plan_path)
    ref_path = OUT / "reference-junction-chains.csv"
    annotation_audit = json.loads((OUT / "annotation-audit.json").read_text())
    assert sha(ref_path) == annotation_audit["files_sha256"][ref_path.name]
    refs = list(csv.DictReader(ref_path.open()))
    for ref in refs:
        ref["chain"] = [tuple(x) for x in json.loads(ref["full_junction_chain_exon_end_next_start"])]
        ref["suffix"] = [tuple(x) for x in json.loads(ref["junction_chain_after_target"])]
        ref["coding"] = int(ref["coding_start_genomic_1based"])
    canonical_suffix = next(r["suffix"] for r in refs if r["transcript"] == "ENST00000319294.11")
    boundary_suffix = next(r["suffix"] for r in refs if r["transcript"] == "ENST00000398367.2")
    for ref in refs:
        ref["architecture"] = "boundary" if ref["suffix"] == boundary_suffix else "canonical_downstream" if ref["suffix"] == canonical_suffix else "other_reference"
    samples = {r["sample_id"]: r for r in manifest["samples"]}
    assert set(samples) == {s["sample_id"] for s in plan["samples"]}
    observations, summaries, per_reference, match_rows, group_counts = [], [], [], [], []
    for sample in plan["samples"]:
        source = samples[sample["sample_id"]]
        for name, pin in source["local_files"].items():
            assert sha(ROOT / name) == pin["sha256"]
        local = []
        with pysam.AlignmentFile(str(ROOT / source["local_bam"]), "rb") as bam:
            for index, read in enumerate(bam.fetch(until_eof=True), 1):
                chain, anchors, blocks = structure(read.reference_start, read.cigartuples or [])
                excluded = exclusions(read.flag, read.mapping_quality)
                event = "plus29" if (42434983, 42437488) in chain else "normal" if (42434954, 42437488) in chain else "other"
                target = (42434983 if event == "plus29" else 42434954, 42437488)
                target_index = chain.index(target) if event != "other" else None
                event_qualified = not excluded and target_index is not None and min(anchors[target_index]) >= 20
                full, downstream, architectures = [], [], set()
                for ref in refs:
                    if not event_qualified or event != ref["event"]:
                        continue
                    suffix_match = chain[target_index + 1:] == ref["suffix"] and all(min(x) >= 20 for x in anchors[target_index + 1:])
                    full_match = chain == ref["chain"] and all(min(x) >= 20 for x in anchors) and covered(blocks, ref["coding"], ref["coding"])
                    if event == "plus29":
                        full_match &= covered(blocks, 42442545, 42442547)
                    if full_match:
                        full.append(ref["transcript"])
                        architectures.add(ref["architecture"])
                    if suffix_match:
                        downstream.append(ref["transcript"])
                    if full_match or suffix_match:
                        match_rows.append({"sample_id": sample["sample_id"], "read_index": index, "query_name": read.query_name,
                                           "event": event, "reference_transcript": ref["transcript"], "full_chain": full_match,
                                           "downstream_chain": suffix_match, "reference_coding_start_aligned": covered(blocks, ref["coding"], ref["coding"]),
                                           "predicted_stop_aligned": covered(blocks, 42442545, 42442547)})
                assert len(architectures) <= 1
                row = {"sample_id": sample["sample_id"], "platform": sample["platform"], "condition": sample["condition"],
                       "read_index": index, "query_name": read.query_name, "flag": read.flag, "mapq": read.mapping_quality,
                       "alignment_strand": "-" if read.is_reverse else "+", "start_1based": read.reference_start + 1,
                       "end_1based": read.reference_end, "cigar": read.cigarstring, "exclusion_reasons": ";".join(excluded),
                       "event_from_coordinates": event, "event_passes_quality_and_anchor": event_qualified,
                       "junctions_json": json.dumps(chain, separators=(",", ":")), "junction_anchors_json": json.dumps(anchors, separators=(",", ":")),
                       "full_chain_reference_matches": ";".join(full), "downstream_chain_reference_matches": ";".join(downstream),
                       "full_chain_architecture": next(iter(architectures), "")}
                local.append(row)
                observations.append(row)
        assert len(local) == source["all_locus_alignments"]
        summaries.append({"sample_id": sample["sample_id"], "platform": sample["platform"], "condition": sample["condition"],
                          "replicate_label": sample["replicate_label"], "all_locus_alignments": len(local),
                          "quality_eligible_alignments": sum(not r["exclusion_reasons"] for r in local),
                          "exact_normal_coordinate_alignments": sum(r["event_from_coordinates"] == "normal" for r in local),
                          "exact_plus29_coordinate_alignments": sum(r["event_from_coordinates"] == "plus29" for r in local),
                          "qualified_normal_event_reads": sum(r["event_passes_quality_and_anchor"] and r["event_from_coordinates"] == "normal" for r in local),
                          "qualified_plus29_event_reads": sum(r["event_passes_quality_and_anchor"] and r["event_from_coordinates"] == "plus29" for r in local),
                          "qualified_full_reference_reads": sum(bool(r["full_chain_reference_matches"]) for r in local)})
        for ref in refs:
            per_reference.append({"sample_id": sample["sample_id"], "platform": sample["platform"], "condition": sample["condition"],
                                  "reference_transcript": ref["transcript"], "event": ref["event"], "architecture": ref["architecture"],
                                  "full_chain_read_count": sum(ref["transcript"] in r["full_chain_reference_matches"].split(";") and r["event_from_coordinates"] == ref["event"] for r in local),
                                  "downstream_chain_read_count": sum(ref["transcript"] in r["downstream_chain_reference_matches"].split(";") and r["event_from_coordinates"] == ref["event"] for r in local)})
        for architecture in ["boundary", "canonical_downstream", "other_reference"]:
            for event in ["normal", "plus29"]:
                group_counts.append({"sample_id": sample["sample_id"], "platform": sample["platform"], "condition": sample["condition"],
                                     "architecture": architecture, "event": event,
                                     "full_chain_read_count": sum(r["full_chain_architecture"] == architecture and r["event_from_coordinates"] == event for r in local)})
    contrasts = []
    for platform in ["PacBio", "ONT_dRNA"]:
        changes = {}
        for architecture in ["boundary", "canonical_downstream"]:
            counts = {f"{condition}_{event}": sum(r["full_chain_read_count"] for r in group_counts if r["platform"] == platform and r["architecture"] == architecture and r["condition"] == condition and r["event"] == event)
                      for condition in ["control", "UPF1_depleted_12h"] for event in ["normal", "plus29"]}
            estimable = all(v > 0 for v in counts.values())
            effect = math.log2((counts["UPF1_depleted_12h_plus29"] / counts["UPF1_depleted_12h_normal"]) / (counts["control_plus29"] / counts["control_normal"])) if estimable else None
            changes[architecture] = effect
            contrasts.append({"platform": platform, "architecture": architecture, **counts, "log2_change_retained_to_normal_odds": effect,
                              "status": "descriptive_estimate" if estimable else "not_estimable", "reason": "" if estimable else "One or more required full-chain counts are zero; no pseudocount used."})
        difference = changes["boundary"] - changes["canonical_downstream"] if all(v is not None for v in changes.values()) else None
        for row in contrasts:
            if row["platform"] == platform:
                row["boundary_minus_canonical_difference"] = difference
    output.mkdir(parents=True, exist_ok=True)
    write(output / "upf1-locus-alignments.csv", observations, ["sample_id", "query_name", "event_from_coordinates"])
    write(output / "upf1-read-summary.csv", summaries, ["sample_id"])
    write(output / "upf1-read-reference-matches.csv", match_rows, ["sample_id", "read_index", "query_name", "event", "reference_transcript", "full_chain", "downstream_chain"])
    write(output / "upf1-read-reference-counts.csv", per_reference, ["sample_id"])
    write(output / "upf1-read-architecture-counts.csv", group_counts, ["sample_id"])
    write(output / "upf1-read-contrasts.csv", contrasts, ["platform"])
    names = ["upf1-locus-alignments.csv", "upf1-read-summary.csv", "upf1-read-reference-matches.csv", "upf1-read-reference-counts.csv", "upf1-read-architecture-counts.csv", "upf1-read-contrasts.csv"]
    audit = {"plan_sha256": sha(plan_path), "read_source_manifest_sha256": sha(ROOT / "config/ubash3a-upf1-read-sources.json"),
             "script_sha256": sha(Path(__file__)), "samples": len(samples), "all_locus_alignments": len(observations),
             "reference_count_rows_including_zeros": len(per_reference), "architecture_count_rows_including_zeros": len(group_counts),
             "files_sha256": {name: sha(output / name) for name in names},
             "interpretation": "Read-level source audit and descriptive conditional contrasts; no expression, decay or protein inference from a zero or unavailable contrast."}
    (output / "upf1-read-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"samples": len(samples), "all_locus_alignments": len(observations), "contrast_statuses": [r["status"] for r in contrasts]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    run(parser.parse_args().output_dir)
