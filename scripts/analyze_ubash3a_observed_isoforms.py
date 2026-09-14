"""Audit deposited transcript structures for the fixed UBASH3A +29 hypothesis.

Catalog entries are counted as models, never as molecules or independent donors.
No NMD treatment response is estimated by this annotation-only script.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-hypothesis-test-plan.json"
INPUTS = ROOT / "config/ubash3a-test-annotation-inputs.json"
OUT = ROOT / "data/derived/ubash3a-hypothesis-test"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attributes(text):
    """Read GTF or GFF3 attributes without guessing coordinate conventions."""
    result = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field and (" " not in field or field.index("=") < field.index(" ")):
            key, value = field.split("=", 1)
            value = unquote(value)
        else:
            match = re.fullmatch(r'(\S+)\s+(?:"([^"]*)"|([^\s";]+))', field)
            if match is None:
                raise ValueError("Unsupported annotation attribute: " + field)
            key, quoted, unquoted = match.groups()
            value = quoted if quoted is not None else unquoted
        if key in result and result[key] != value:
            raise ValueError("Conflicting duplicate annotation attribute: " + key)
        result[key] = value
    return result


def transcript_ids(kind, attrs):
    if "transcript_id" in attrs:
        return [attrs["transcript_id"]]
    if kind in {"transcript", "mRNA"} and "ID" in attrs:
        return [attrs["ID"]]
    if kind == "exon" and "Parent" in attrs:
        return attrs["Parent"].split(",")
    return []


def chromosome(value):
    return value[3:] if value.startswith("chr") else value


def text_lines(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        yield from stream


def records(path):
    for number, line in enumerate(text_lines(path), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"Expected nine annotation fields at {path.name}:{number}")
        if fields[2] not in {"exon", "transcript", "mRNA"}:
            continue
        start, end = int(fields[3]), int(fields[4])
        if start < 1 or end < start:
            raise ValueError("Invalid one-based annotation interval")
        yield number, fields, start, end


def collect_locus(path, gene_start, gene_end, symbol="UBASH3A", gene_id="ENSG00000160185"):
    selected, total_exon_rows = set(), 0
    for _, fields, start, end in records(path):
        total_exon_rows += fields[2] == "exon"
        coordinate_hit = chromosome(fields[0]) == "21" and start <= gene_end and end >= gene_start
        label_candidate = symbol in fields[8] or gene_id in fields[8]
        if not coordinate_hit and not label_candidate:
            continue
        attrs = attributes(fields[8])
        label_hit = any(v == symbol or v.split(".")[0] == gene_id for k, v in attrs.items() if k in {"gene_name", "gene_id", "gene", "gene_symbol"})
        if coordinate_hit or label_hit:
            selected.update(transcript_ids(fields[2], attrs))
    if total_exon_rows == 0:
        raise ValueError("Source contains no exon records")
    collected = {tid: {"exons": [], "chromosomes": set(), "strands": set(), "source_labels": set(), "attributes": [], "source_line_numbers": [], "coordinate_selected": False, "gene_label_selected": False} for tid in selected}
    for number, fields, start, end in records(path):
        # Most records are on other chromosomes; candidate IDs with explicit
        # UBASH3A labels on other chromosomes are still recovered below.
        if chromosome(fields[0]) != "21" and symbol not in fields[8] and gene_id not in fields[8] and not any(tid in fields[8] for tid in selected):
            continue
        attrs = attributes(fields[8])
        for tid in transcript_ids(fields[2], attrs):
            if tid not in collected:
                continue
            item = collected[tid]
            item["chromosomes"].add(chromosome(fields[0]))
            item["strands"].add(fields[6])
            item["source_labels"].add(fields[1])
            item["source_line_numbers"].append(number)
            if attrs not in item["attributes"]:
                item["attributes"].append(attrs)
            item["coordinate_selected"] |= chromosome(fields[0]) == "21" and start <= gene_end and end >= gene_start
            item["gene_label_selected"] |= any(v == symbol or v.split(".")[0] == gene_id for k, v in attrs.items() if k in {"gene_name", "gene_id", "gene", "gene_symbol"})
            if fields[2] == "exon":
                item["exons"].append((start, end))
    return collected, total_exon_rows


def junction_chain(exons):
    ordered = sorted(set(exons))
    if not ordered:
        raise ValueError("Transcript has no exons")
    if any(a[1] >= b[0] for a, b in zip(ordered, ordered[1:])):
        raise ValueError("Transcript has overlapping exon intervals")
    return tuple((a[1], b[0]) for a, b in zip(ordered, ordered[1:]))


def interval_contains(exons, position):
    return any(start <= position <= end for start, end in exons)


def event_and_suffix(chain, plan):
    target = plan["downstream_exon_start_1based"]
    normal = (plan["normal_upstream_exon_end_1based"], target)
    retained = (plan["retained_upstream_exon_end_1based"], target)
    hits = [(i, "normal" if pair == normal else "plus29") for i, pair in enumerate(chain) if pair in {normal, retained}]
    if len(hits) > 1:
        raise ValueError("Transcript repeats the target junction")
    if not hits:
        return "other", None, ()
    index, event = hits[0]
    return event, index, chain[index + 1:]


def reference_models(plan):
    gene = json.loads((ROOT / "data/raw/ubash3a-consequences/ensembl-gene.json").read_text())
    previous = json.loads((ROOT / "config/ubash3a-consequence-plan.json").read_text())
    wanted = {t["id"] for t in previous["expected_transcripts"]}
    refs = []
    for t in sorted(gene["Transcript"], key=lambda t: t["id"]):
        if t["id"] not in wanted:
            continue
        for event in ("normal", "plus29"):
            exons = sorted((e["start"], e["end"] + (29 if event == "plus29" and e["end"] == plan["normal_upstream_exon_end_1based"] else 0)) for e in t["Exon"])
            chain = junction_chain(exons)
            _, index, suffix = event_and_suffix(chain, plan)
            refs.append({"transcript": t["id"] + "." + str(t["version"]), "event": event, "biotype": t["biotype"], "coding_start_genomic": t["Translation"]["start"], "exons": exons, "chain": chain, "suffix": suffix, "target_index": index})
    return gene, refs


def write_csv(path, rows, fallback_fields=None):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or fallback_fields
    if not fields:
        raise ValueError("Missing output schema")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output_dir=OUT):
    plan = json.loads(PLAN.read_text())
    inputs = json.loads(INPUTS.read_text())
    if sha(ROOT / plan["protocol"]) != plan["protocol_sha256"] or sha(PLAN) != inputs["plan_sha256"]:
        raise ValueError("Frozen plan or protocol changed")
    for path, expected in plan["reference_pins"].items():
        if sha(ROOT / path) != expected:
            raise ValueError("Prior pinned file changed: " + path)
    gene, refs = reference_models(plan)
    observed, exon_rows, matches, summaries, reference_rows = [], [], [], [], []
    for ref in refs:
        reference_rows.append({"transcript": ref["transcript"], "event": ref["event"], "biotype": ref["biotype"], "coding_start_genomic_1based": ref["coding_start_genomic"], "first_exon_start_1based": ref["exons"][0][0], "final_exon_end_1based": ref["exons"][-1][1], "exon_count": len(ref["exons"]), "all_exons_1based_inclusive": json.dumps(ref["exons"], separators=(",", ":")), "full_junction_chain_exon_end_next_start": json.dumps(ref["chain"], separators=(",", ":")), "junction_chain_after_target": json.dumps(ref["suffix"], separators=(",", ":"))})
    for source in inputs["annotations"]:
        path = ROOT / source["path"]
        if sha(path) != source["sha256"]:
            raise ValueError("Annotation source changed: " + source["source_id"])
        models, total_exons = collect_locus(path, gene["start"], gene["end"])
        local = []
        for tid, item in sorted(models.items()):
            exons = sorted(set(item["exons"]))
            eligible = item["chromosomes"] == {"21"} and item["strands"] == {"+"} and bool(exons)
            chain = junction_chain(exons) if eligible else ()
            event, index, suffix = event_and_suffix(chain, plan) if eligible else ("ineligible_context", None, ())
            same_chain, same_downstream = [], []
            for ref in refs:
                if event != ref["event"] or not eligible:
                    continue
                full = chain == ref["chain"]
                downstream = suffix == ref["suffix"]
                if full:
                    same_chain.append(ref["transcript"])
                if downstream:
                    same_downstream.append(ref["transcript"])
                if full or downstream:
                    matches.append({"source_id": source["source_id"], "source_transcript_id": tid, "event": event, "reference_transcript": ref["transcript"], "full_internal_chain_matches": full, "downstream_chain_matches": downstream, "reference_coding_start_is_exonic": interval_contains(exons, ref["coding_start_genomic"]), "retained_stop_is_fully_exonic": all(interval_contains(exons, p) for p in range(plan["first_retained_stop_genomic_interval_1based"][0], plan["first_retained_stop_genomic_interval_1based"][1] + 1)), "first_exon_start_difference_nt": exons[0][0] - ref["exons"][0][0], "last_exon_end_difference_nt": exons[-1][1] - ref["exons"][-1][1], "exact_all_exon_coordinates": exons == ref["exons"]})
            row = {"source_id": source["source_id"], "family": source["family"], "context": source["context"], "source_transcript_id": tid, "chromosomes": ";".join(sorted(item["chromosomes"])), "strands": ";".join(sorted(item["strands"])), "coordinate_selected": item["coordinate_selected"], "gene_label_selected": item["gene_label_selected"], "eligible_chr21_positive": eligible, "event": event, "exon_count": len(exons), "duplicate_exon_rows": len(item["exons"]) - len(exons), "first_exon_start_1based": exons[0][0] if exons else None, "last_exon_end_1based": exons[-1][1] if exons else None, "full_chain_reference_matches": ";".join(same_chain), "downstream_chain_reference_matches": ";".join(same_downstream), "target_acceptor_donor_ends_1based": ";".join(str(a) for a, b in chain if b == plan["downstream_exon_start_1based"]), "junction_chain_exon_end_next_start": json.dumps(chain, separators=(",", ":")), "source_feature_labels": ";".join(sorted(item["source_labels"])), "source_attributes_json": json.dumps(item["attributes"], separators=(",", ":")), "source_line_numbers": ";".join(map(str, item["source_line_numbers"])), "support_unit": source["support_unit"], "interpretation_boundary": source["interpretation_boundary"]}
            observed.append(row)
            local.append(row)
            exon_rows.extend({"source_id": source["source_id"], "source_transcript_id": tid, "exon_number_in_genomic_order": i + 1, "start_1based": start, "end_1based": end} for i, (start, end) in enumerate(exons))
        summary = {"source_id": source["source_id"], "family": source["family"], "context": source["context"], "total_source_exon_rows": total_exons, "locus_or_label_models": len(local), "eligible_positive_locus_models": sum(r["eligible_chr21_positive"] for r in local), "normal_junction_models": sum(r["event"] == "normal" for r in local), "plus29_junction_models": sum(r["event"] == "plus29" for r in local), "plus29_full_chain_matches": sum(r["event"] == "plus29" and bool(r["full_chain_reference_matches"]) for r in local), "plus29_boundary_downstream_matches": sum(r["event"] == "plus29" and "ENST00000398367.2" in r["downstream_chain_reference_matches"] for r in local), "other_or_ineligible_models": sum(r["event"] not in {"normal", "plus29"} for r in local), "support_unit": source["support_unit"], "interpretation_boundary": source["interpretation_boundary"]}
        summaries.append(summary)
        print(json.dumps({k: summary[k] for k in ("source_id", "locus_or_label_models", "normal_junction_models", "plus29_junction_models")}), flush=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "reference-junction-chains.csv", reference_rows)
    write_csv(output_dir / "observed-transcript-models.csv", observed)
    write_csv(output_dir / "observed-exons.csv", exon_rows)
    write_csv(output_dir / "reference-chain-matches.csv", matches, ["source_id", "source_transcript_id", "event", "reference_transcript"])
    write_csv(output_dir / "annotation-summary.csv", summaries)
    audit = {"analysis_id": plan["plan_id"], "plan_sha256": sha(PLAN), "input_specification_sha256": sha(INPUTS), "script_sha256": sha(Path(__file__)), "sources": len(summaries), "observed_model_rows": len(observed), "exon_rows": len(exon_rows), "reference_match_rows": len(matches), "prior_reference_pins_unchanged": True, "files_sha256": {p.name: sha(p) for p in sorted(output_dir.iterdir()) if p.is_file() and p.name in {"reference-junction-chains.csv", "observed-transcript-models.csv", "observed-exons.csv", "reference-chain-matches.csv", "annotation-summary.csv"}}, "scope": "Deposited annotation structures only. Counts are transcript-model counts, not molecules, donors, protein detections or NMD-response effects."}
    (output_dir / "annotation-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    run(args.output_dir)
