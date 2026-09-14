"""Independently check fixed UBASH3A outputs using genomic exon extension and Bio.Seq.

This script deliberately does not import the primary reconstruction. It edits a
genomic reference before splicing, extends the selected exon interval, maps CDS
coordinates by interval lengths, and translates with Biopython's standard code.
"""

import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ubash3a-consequences"
OUT = ROOT / "data/derived/ubash3a-consequences"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def read_csv(name):
    with (OUT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def run():
    comparisons = 0

    def equal(actual, expected, label):
        nonlocal comparisons
        comparisons += 1
        if actual != expected:
            raise ValueError(f"Independent check failed ({label}): {actual!r} != {expected!r}")

    def fields(row, values):
        for key, value in values.items():
            equal(row[key], "" if value is None else str(value), row.get("scenario_id", row.get("transcript_id", "")) + ":" + key)

    plan = read_json(ROOT / "config/ubash3a-consequence-plan.json")
    manifest = read_json(ROOT / "config/ubash3a-consequence-sources.json")
    audit = read_json(OUT / "audit.json")
    equal(digest(ROOT / "config/ubash3a-consequence-plan.json"), manifest["plan_sha256"], "source plan freeze")
    equal(digest(ROOT / plan["protocol"]), plan["protocol_sha256"], "protocol freeze")
    equal(digest(ROOT / "config/ubash3a-consequence-sources.json"), audit["source_manifest_sha256"], "manifest pin")
    equal(digest(ROOT / "scripts/analyze_ubash3a_consequences.py"), audit["script_sha256"], "primary script pin")
    for path, pin in plan["reference_pins"].items():
        equal(digest(ROOT / path), pin, "prior artifact:" + path)
    for source in manifest["sources"]:
        path = ROOT / source["local_cache_path"]
        equal(digest(path), source["sha256"], "source hash:" + source["source_id"])
        equal(path.stat().st_size, source["bytes"], "source size:" + source["source_id"])
    for filename, pin in audit["files_sha256"].items():
        equal(digest(OUT / filename), pin, "output pin:" + filename)

    gene = read_json(RAW / "ensembl-gene.json")
    equal(gene, read_json(ROOT / "data/raw/UBASH3A-transcripts-GRCh38.json"), "fresh versus preserved expanded gene annotation")
    genome = read_json(RAW / "ensembl-gene-sequence.json")["seq"].upper()
    equal(len(genome), gene["end"] - gene["start"] + 1, "genomic length")
    equal(gene["strand"], 1, "gene strand")
    variant_offset = plan["variant"]["position_1based"] - gene["start"]
    equal(genome[variant_offset], "A", "genomic reference allele")
    inventory = {r["transcript_id"]: r for r in read_csv("transcript-inventory.csv")}
    equal(set(inventory), {t["id"] for t in gene["Transcript"]}, "complete transcript inventory")
    selected = {}
    for transcript in gene["Transcript"]:
        tid = transcript["id"]
        exons = sorted(transcript["Exon"], key=lambda e: e["start"])
        candidates = [i for i, (a, b) in enumerate(zip(exons, exons[1:])) if (a["end"], b["start"]) == (plan["canonical_upstream_exon_end_1based"], plan["downstream_exon_start_1based"])]
        equal(len(candidates) <= 1, True, "unique target junction:" + tid)
        fields(inventory[tid], {"version": transcript["version"], "display_name": transcript["display_name"], "biotype": transcript["biotype"], "is_canonical": transcript["is_canonical"], "selected": bool(candidates), "reason": "exact_consecutive_exon_boundaries" if candidates else "lacks_exact_consecutive_exon_boundaries", "exon_count": len(exons), "reference_cdna_length": sum(e["end"] - e["start"] + 1 for e in exons), "annotated_protein_length": transcript.get("Translation", {}).get("length"), "target_upstream_exon_number": candidates[0] + 1 if candidates else None})
        if candidates:
            selected[tid] = (transcript, exons, candidates[0])
    equal(set(selected), {t["id"] for t in plan["expected_transcripts"]}, "selected transcript set")
    for expected in plan["expected_transcripts"]:
        transcript = selected[expected["id"]][0]
        equal((transcript["version"], transcript["biotype"]), (expected["version"], expected["biotype"]), "frozen version and biotype:" + expected["id"])

    fasta = {name: SeqIO.to_dict(SeqIO.parse(OUT / name, "fasta")) for name in ("transcripts.fasta", "orfs-with-stop.fasta", "proteins.fasta", "retained-segments.fasta")}
    rows = {r["scenario_id"]: r for r in read_csv("consequences.csv")}
    exon_rows = read_csv("exon-coordinates.csv")
    alignment = read_csv("protein-reference-alignment.csv")
    features = read_csv("protein-features.csv")
    expected_keys = {t["id"] + "." + str(t["version"]) + "__" + s["key"] for t, _, _ in selected.values() for s in plan["scenarios"]}
    equal(set(rows), expected_keys, "scenario table keys")
    for name in ("transcripts.fasta", "orfs-with-stop.fasta", "proteins.fasta"):
        equal(set(fasta[name]), expected_keys, "FASTA keys:" + name)
    equal(set(fasta["retained-segments.fasta"]), {"plus29_A", "plus29_C"}, "retained FASTA keys")
    constructed, boundary_details = {}, []
    reference_checks = []
    for tid, (transcript, source_exons, target_index) in sorted(selected.items()):
        reference_peptide = read_json(RAW / (tid + "-protein.json"))["seq"]
        original_lengths = [e["end"] - e["start"] + 1 for e in source_exons]
        insertion_boundary = sum(original_lengths[:target_index + 1])
        for scenario in plan["scenarios"]:
            key = tid + "." + str(transcript["version"]) + "__" + scenario["key"]
            retained = scenario["splice"] == "plus29"
            # Mutate the genomic allele first, then splice the extended intervals.
            allele_genome = genome[:variant_offset] + scenario["allele"] + genome[variant_offset + 1:]
            intervals = [(e["start"], e["end"] + (29 if retained and i == target_index else 0)) for i, e in enumerate(source_exons)]
            lengths = [end - start + 1 for start, end in intervals]
            boundaries = list(itertools.accumulate(lengths))
            cdna = "".join(allele_genome[start - gene["start"]:end - gene["start"] + 1] for start, end in intervals)
            genomic_map = [p for start, end in intervals for p in range(start, end + 1)]
            cds_start = sum(end - start + 1 for start, end in intervals if end < transcript["Translation"]["start"])
            start_exon = next((start, end) for start, end in intervals if start <= transcript["Translation"]["start"] <= end)
            cds_start += transcript["Translation"]["start"] - start_exon[0]
            equal(cdna[cds_start:cds_start + 3], "ATG", "annotated start:" + key)
            coding = cdna[cds_start:]
            translated = str(Seq(coding[:len(coding) // 3 * 3]).translate(table=1))
            equal("*" in translated, True, "first stop available:" + key)
            peptide = translated.split("*", 1)[0]
            stop_start = cds_start + len(peptide) * 3  # zero based
            stop_end = stop_start + 3  # exclusive, also one-based inclusive end
            orf = cdna[cds_start:stop_end]
            equal(str(fasta["transcripts.fasta"][key].seq), cdna, "independent genomic splicing:" + key)
            equal(str(fasta["proteins.fasta"][key].seq), peptide, "Biopython translation:" + key)
            equal(str(fasta["orfs-with-stop.fasta"][key].seq), orf, "complete ORF:" + key)
            if not retained:
                equal(cdna, read_json(RAW / (tid + "-cdna.json"))["seq"], "API cDNA:" + key)
                equal(orf, read_json(RAW / (tid + "-cds.json"))["seq"], "API CDS including stop:" + key)
                equal(peptide, reference_peptide, "API protein:" + key)
                equal(genomic_map[stop_end - 1], transcript["Translation"]["end"], "annotated CDS end includes stop:" + key)
                equal(len(peptide), transcript["Translation"]["length"], "lookup protein length:" + key)
                if scenario["key"] == "normal_A":
                    reference_checks.append({"transcript": tid + "." + str(transcript["version"]), "cdna_cds_protein_and_lookup_match": True})
            prefix = len(list(itertools.takewhile(lambda pair: pair[0] == pair[1], zip(reference_peptide, peptide))))
            stop_exon_index = next(i for i, boundary in enumerate(boundaries) if stop_start < boundary)
            last_distance = boundaries[-2] - stop_end
            next_junction = next((b for b in boundaries[:-1] if b >= stop_end), None)
            kept_segment = allele_genome[plan["retained_interval_1based_inclusive"][0] - gene["start"]:plan["retained_interval_1based_inclusive"][1] - gene["start"] + 1]
            if retained:
                equal(str(fasta["retained-segments.fasta"][scenario["key"]].seq), kept_segment, "retained genomic segment:" + key)
            fields(rows[key], {"transcript_id": tid, "transcript_version": transcript["version"], "biotype": transcript["biotype"], "is_canonical": transcript["is_canonical"], "scenario": scenario["key"], "allele": scenario["allele"], "splice": scenario["splice"], "consequence": "coding_29nt_retention_frameshift" if retained else "normal_splicing_variant_removed", "cdna_length_nt": len(cdna), "added_nt": len(cdna) - sum(original_lengths), "retained_segment": kept_segment if retained else "", "retention_after_reference_cdna_base_1based": insertion_boundary if retained else None, "retention_after_reference_coding_base_1based": insertion_boundary - cds_start if retained else None, "variant_cdna_position_1based": genomic_map.index(plan["variant"]["position_1based"]) + 1 if retained else None, "cds_start_cdna_1based": cds_start + 1, "reference_protein_length_aa": len(reference_peptide), "protein_length_aa": len(peptide), "protein_changed": peptide != reference_peptide, "first_changed_amino_acid_1based": prefix + 1 if retained else None, "unchanged_prefix_aa": prefix, "novel_peptide_tail": peptide[prefix:] if retained else "", "reference_suffix_replaced_or_removed_aa": len(reference_peptide) - prefix, "stop_codon": cdna[stop_start:stop_end], "stop_start_cdna_1based": stop_start + 1, "stop_end_cdna_1based": stop_end, "stop_start_coding_1based": stop_start - cds_start + 1, "stop_genomic_positions_1based": ";".join(map(str, genomic_map[stop_start:stop_end])), "no_stop_in_transcript": False, "cdna_sha256": hashlib.sha256(cdna.encode()).hexdigest(), "protein_sha256": hashlib.sha256(peptide.encode()).hexdigest(), "nmd_status": "positional_features_only", "final_junction_cdna_boundary": boundaries[-2], "stop_end_to_final_junction_nt": last_distance, "stop_end_to_next_junction_nt": next_junction - stop_end if next_junction is not None else None, "stop_in_final_exon": stop_exon_index == len(intervals) - 1, "stop_spans_junction": stop_end > boundaries[stop_exon_index], "stop_exon_number": stop_exon_index + 1, "stop_exon_length_nt": lengths[stop_exon_index], "stop_exon_gt_400_nt": lengths[stop_exon_index] > 400, "stop_start_within_first_150_coding_nt": stop_start - cds_start < 150, "nmd_55nt_flag": last_distance > 55, "nmd_50nt_flag": last_distance > 50})
            per_exon = [r for r in exon_rows if r["scenario_id"] == key]
            equal(len(per_exon), len(intervals), "exon row count:" + key)
            for i, row in enumerate(per_exon):
                fields(row, {"exon_number": i + 1, "exon_id": source_exons[i]["id"], "genomic_start_1based": intervals[i][0], "genomic_end_1based": intervals[i][1], "cdna_start_0based": boundaries[i] - lengths[i], "cdna_end_0based_exclusive": boundaries[i]})
            constructed[key] = (cdna, peptide, genomic_map, cds_start)
            boundary_details.append({"scenario_id": key, "stop_start_cdna_1based": stop_start + 1, "stop_end_cdna_1based": stop_end, "final_junction_cdna_boundary": boundaries[-2], "stop_end_to_final_junction_nt": last_distance, "stop_first_base_to_final_junction_nt": boundaries[-2] - (stop_start + 1), "primary_gt55_from_stop_end": last_distance > 55, "sensitivity_gt50_from_stop_end": last_distance > 50, "sensitivity_gt50_from_stop_first_base": boundaries[-2] - (stop_start + 1) > 50})

    uniprot = read_json(RAW / "uniprot-entry.json")
    canonical = [key for key in constructed if key.endswith("__normal_A") and rows[key]["is_canonical"] == "1"]
    equal(len(canonical), 1, "unique canonical transcript")
    _, anchor_peptide, anchor_positions, anchor_start = constructed[canonical[0]]
    equal(anchor_peptide, uniprot["sequence"]["value"], "canonical exact UniProt anchor")
    # Index each genomic nucleotide by reference residue AND codon offset.
    per_base = {p: (aa_index + 1, offset) for aa_index in range(len(anchor_peptide)) for offset, p in enumerate(anchor_positions[anchor_start + aa_index * 3:anchor_start + aa_index * 3 + 3])}
    checked_mapping = {}
    for key, (_, peptide, positions, start) in constructed.items():
        mapping_rows = [r for r in alignment if r["scenario_id"] == key]
        equal(len(mapping_rows), len(peptide), "alignment row count:" + key)
        mapped = []
        for i, row in enumerate(mapping_rows):
            codon = positions[start + 3 * i:start + 3 * i + 3]
            origins = [per_base.get(p) for p in codon]
            reference_index = origins[0][0] if all(origin is not None for origin in origins) and [origin[1] for origin in origins] == [0, 1, 2] and len({origin[0] for origin in origins}) == 1 else None
            if reference_index is not None:
                equal(peptide[i], anchor_peptide[reference_index - 1], "projected amino acid:" + key)
                mapped.append((i + 1, reference_index))
            fields(row, {"peptide_position_1based": i + 1, "amino_acid": peptide[i], "genomic_codon_positions_1based": ";".join(map(str, codon)), "uniprot_position_1based": reference_index, "status": "exact_genomic_codon_and_amino_acid" if reference_index is not None else "no_reference_codon_match"})
        checked_mapping[key] = mapped
        for feature in [f for f in uniprot["features"] if f["type"] in ("Domain", "Region", "Active site", "Binding site")]:
            lo, hi = feature["location"]["start"]["value"], feature["location"]["end"]["value"]
            hits = [i for i, reference_index in mapped if lo <= reference_index <= hi]
            intact = len(hits) == hi - lo + 1 and all(b == a + 1 for a, b in zip(hits, hits[1:]))
            matches = [r for r in features if r["scenario_id"] == key and r["feature_name"] == feature.get("description", "") and r["uniprot_start_1based"] == str(lo)]
            equal(len(matches), 1, "unique feature row:" + key)
            fields(matches[0], {"transcript_id": rows[key]["transcript_id"], "scenario": rows[key]["scenario"], "uniprot_accession": uniprot["primaryAccession"], "feature_type": feature["type"], "uniprot_end_1based": hi, "reference_length_aa": hi - lo + 1, "reference_residues_present": len(hits), "peptide_first_mapped_residue": min(hits) if hits else None, "peptide_last_mapped_residue": max(hits) if hits else None, "sequence_status": "intact_reference_sequence" if intact else ("partial_reference_sequence" if hits else "reference_sequence_absent")})
            equal(json.loads(matches[0]["evidence"]), feature.get("evidences", []), "feature source evidence:" + key)
    for transcript, _, _ in selected.values():
        versioned = transcript["id"] + "." + str(transcript["version"])
        equal(constructed[versioned + "__normal_A"][:2], constructed[versioned + "__normal_C"][:2], "normal alleles identical:" + versioned)
        a, c = constructed[versioned + "__plus29_A"], constructed[versioned + "__plus29_C"]
        equal([p for x, y, p in zip(a[0], c[0], a[2]) if x != y], [plan["variant"]["position_1based"]], "only retained variant differs:" + versioned)
        peptide_differences = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a[1], c[1])) if x != y]
        equal(peptide_differences, [(int(rows[versioned + "__plus29_A"]["first_changed_amino_acid_1based"]), "M", "L")], "retained peptide allele difference:" + versioned)

    report = {"verification_id": "ubash3a-consequences-independent-v1", "verified_at_utc": datetime.now(timezone.utc).isoformat(), "status": "all_checks_passed", "scalar_comparisons": comparisons, "methods": {"primary": "Insert 29 nt into reconstructed normal cDNA; literal standard codon table", "independent": "Replace genomic allele, extend only the chosen genomic exon interval, concatenate exons; Biopython translation", "feature_check": "Per-genomic-nucleotide reference residue and codon-offset mapping", "shared_inputs": "Both implementations use the same pinned source annotation and sequences; this verifies computation, not the biological truth or completeness of the annotations.", "biopython_version": Bio.__version__}, "reference_checks": reference_checks, "counts": {"source_files": len(manifest["sources"]), "transcripts_inventoried": len(inventory), "selected_transcripts": len(selected), "scenarios": len(rows), "aligned_peptide_residues": len(alignment), "feature_rows": len(features), "exon_rows": len(exon_rows)}, "stop_coordinate_sensitivity": boundary_details, "boundary_interpretation": "ENST00000398367.2 retained stop is 49 nt from its final nucleotide, or 51 nt from its first nucleotide, to the final exon boundary. It is below the frozen >55 nt flag under both conventions; the >50 nt sensitivity classification changes with codon convention. This is a boundary case, not established NMD escape.", "plan_sha256": digest(ROOT / "config/ubash3a-consequence-plan.json"), "primary_audit_sha256": digest(OUT / "audit.json"), "checker_sha256": digest(Path(__file__))}
    destination = ROOT / "reports/ubash3a-consequences-independent-check.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "scalar_comparisons": comparisons, "counts": report["counts"], "report": str(destination)}, indent=2))


if __name__ == "__main__":
    run()
