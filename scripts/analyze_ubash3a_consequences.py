"""Reconstruct the fixed UBASH3A +29-nt splice scenarios from pinned sequences.

The primary implementation uses a literal standard genetic code. The independent
checker uses Biopython and reconstructs the altered genomic exon intervals anew.
No expression, clinical outcome, or splice-frequency data are analyzed here.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import platform

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ubash3a-consequences"
OUT = ROOT / "data/derived/ubash3a-consequences"
PLAN_PATH = ROOT / "config/ubash3a-consequence-plan.json"
SOURCES_PATH = ROOT / "config/ubash3a-consequence-sources.json"
STOPS = {"TAA", "TAG", "TGA"}
STANDARD_CODE = dict(zip(
    ("".join(c) for c in itertools.product("TCAG", repeat=3)),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    strict=True,
))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def ordered_exons(transcript):
    if transcript["strand"] != 1:
        raise ValueError("This fixed genomic reconstruction requires positive strand")
    exons = sorted(transcript["Exon"], key=lambda e: e["start"])
    if not exons or any(e["strand"] != 1 or e["end"] < e["start"] for e in exons):
        raise ValueError("Invalid exon coordinates or strand")
    if any(a["end"] >= b["start"] for a, b in zip(exons, exons[1:])):
        raise ValueError("Overlapping exon intervals")
    return exons


def reconstruct_reference(genome, region_start, exons):
    """Return cDNA, a nucleotide-level genomic map and exon cDNA intervals."""
    sequence, coordinates, blocks = [], [], []
    offset = 0
    for i, exon in enumerate(exons, 1):
        start, end = exon["start"] - region_start, exon["end"] - region_start + 1
        if start < 0 or end > len(genome) or end <= start:
            raise ValueError("Exon outside captured genomic sequence")
        part = genome[start:end]
        sequence.append(part)
        coordinates.extend(range(exon["start"], exon["end"] + 1))
        blocks.append({"exon_number": i, "exon_id": exon["id"], "genomic_start_1based": exon["start"], "genomic_end_1based": exon["end"], "cdna_start_0based": offset, "cdna_end_0based_exclusive": offset + len(part)})
        offset += len(part)
    cdna = "".join(sequence)
    if set(cdna) - set("ACGT"):
        raise ValueError("Ambiguous cDNA base; exact translation is unavailable")
    return cdna, coordinates, blocks


def translate_from_start(cdna, start):
    if start < 0 or start + 3 > len(cdna):
        raise ValueError("Invalid translation start")
    if cdna[start:start + 3] != "ATG":
        raise ValueError("The annotated start is not ATG")
    peptide = []
    for offset in range(start, len(cdna) - 2, 3):
        codon = cdna[offset:offset + 3]
        if codon not in STANDARD_CODE:
            raise ValueError("Ambiguous translated codon")
        aa = STANDARD_CODE[codon]
        if aa == "*":
            return {"peptide": "".join(peptide), "stop_start_0based": offset, "stop_end_0based_exclusive": offset + 3, "stop_codon": codon, "orf_with_stop": cdna[start:offset + 3], "no_stop_in_transcript": False}
        peptide.append(aa)
    return {"peptide": "".join(peptide), "stop_start_0based": None, "stop_end_0based_exclusive": None, "stop_codon": None, "orf_with_stop": cdna[start:start + len(peptide) * 3], "no_stop_in_transcript": True}


def reference_coding_bounds(cdna, coordinates, translation, source_cds, source_protein):
    try:
        start = coordinates.index(translation["start"])
        end = coordinates.index(translation["end"]) + 1
    except ValueError as exc:
        raise ValueError("Annotated CDS endpoint is not exonic") from exc
    if cdna[start:end] == source_cds:
        convention = "lookup_translation_bounds_equal_API_CDS"
    elif cdna[start:end + 3] == source_cds and source_cds[-3:] in STOPS:
        convention = "API_CDS_adds_stop_after_lookup_translation_bounds"
    else:
        raise ValueError("Genomic translation bounds do not reproduce the API CDS")
    observed = translate_from_start(cdna, start)
    if observed["peptide"] != source_protein or observed["no_stop_in_transcript"]:
        raise ValueError("Reference translation does not reproduce source protein and terminal stop")
    if observed["orf_with_stop"] != source_cds:
        raise ValueError("API CDS differs from the complete start-to-first-stop ORF")
    if len(source_protein) != translation["length"]:
        raise ValueError("Reference protein length disagrees with lookup annotation")
    return start, observed, convention


def insert_retained_segment(cdna, coordinates, blocks, insert_at, segment, genomic_start):
    if insert_at not in [b["cdna_end_0based_exclusive"] for b in blocks[:-1]]:
        raise ValueError("Retention position must be an internal exon junction")
    if not segment or set(segment) - set("ACGT"):
        raise ValueError("Invalid retained sequence")
    n = len(segment)
    upstream = next(i for i, b in enumerate(blocks[:-1]) if b["cdna_end_0based_exclusive"] == insert_at)
    if genomic_start != blocks[upstream]["genomic_end_1based"] + 1 or genomic_start + n >= blocks[upstream + 1]["genomic_start_1based"]:
        raise ValueError("Retained interval must start at the donor and leave an intron")
    altered_blocks = []
    for block in blocks:
        b = dict(block)
        if b["cdna_end_0based_exclusive"] == insert_at:
            b["cdna_end_0based_exclusive"] += n
            b["genomic_end_1based"] += n
        elif b["cdna_start_0based"] >= insert_at:
            b["cdna_start_0based"] += n
            b["cdna_end_0based_exclusive"] += n
        altered_blocks.append(b)
    return cdna[:insert_at] + segment + cdna[insert_at:], coordinates[:insert_at] + list(range(genomic_start, genomic_start + n)) + coordinates[insert_at:], altered_blocks


def common_prefix_length(a, b):
    for i, (left, right) in enumerate(zip(a, b)):
        if left != right:
            return i
    return min(len(a), len(b))


def nmd_features(translation, start, blocks, settings):
    stop_start, stop_end = translation["stop_start_0based"], translation["stop_end_0based_exclusive"]
    junctions = [b["cdna_end_0based_exclusive"] for b in blocks[:-1]]
    if stop_start is None:
        return {"nmd_status": "unavailable_no_stop", "final_junction_cdna_boundary": junctions[-1] if junctions else None}
    exon = next(b for b in blocks if b["cdna_start_0based"] <= stop_start < b["cdna_end_0based_exclusive"])
    downstream = [j for j in junctions if j >= stop_end]
    last_distance = junctions[-1] - stop_end if junctions else None
    next_distance = min(downstream) - stop_end if downstream else None
    long_length = exon["cdna_end_0based_exclusive"] - exon["cdna_start_0based"]
    return {
        "nmd_status": "positional_features_only",
        "final_junction_cdna_boundary": junctions[-1] if junctions else None,
        "stop_end_to_final_junction_nt": last_distance,
        "stop_end_to_next_junction_nt": next_distance,
        "stop_in_final_exon": exon["exon_number"] == len(blocks),
        "stop_spans_junction": stop_end > exon["cdna_end_0based_exclusive"],
        "stop_exon_number": exon["exon_number"],
        "stop_exon_length_nt": long_length,
        "stop_exon_gt_400_nt": long_length > settings["long_stop_exon_gt_nt"],
        "stop_start_within_first_150_coding_nt": stop_start - start < settings["start_proximity_coding_nt"],
        "nmd_55nt_flag": last_distance is not None and last_distance > settings["primary_stop_end_to_final_junction_gt_nt"],
        "nmd_50nt_flag": last_distance is not None and last_distance > settings["sensitivity_stop_end_to_final_junction_gt_nt"],
    }


def codon_alignment(peptide, coordinates, start, anchor_codons, uniprot_sequence):
    rows = []
    for i, aa in enumerate(peptide):
        key = tuple(coordinates[start + 3 * i:start + 3 * i + 3])
        if len(key) != 3:
            raise ValueError("Incomplete amino-acid coordinate map")
        reference_position = anchor_codons.get(key)
        if reference_position is not None and uniprot_sequence[reference_position - 1] != aa:
            raise ValueError("Matched reference genomic codon has a different amino acid")
        rows.append({"peptide_position_1based": i + 1, "amino_acid": aa, "genomic_codon_positions_1based": ";".join(map(str, key)), "uniprot_position_1based": reference_position, "status": "exact_genomic_codon_and_amino_acid" if reference_position is not None else "no_reference_codon_match"})
    mapped = [r["uniprot_position_1based"] for r in rows if r["uniprot_position_1based"] is not None]
    if mapped != sorted(set(mapped)):
        raise ValueError("Reference projection is ambiguous or out of order")
    return rows


def summarize_feature(feature, alignment):
    start, end = feature["location"]["start"], feature["location"]["end"]
    if start.get("modifier") != "EXACT" or end.get("modifier") != "EXACT":
        raise ValueError("Feature bounds are not exact")
    lo, hi = start["value"], end["value"]
    hits = [r for r in alignment if r["uniprot_position_1based"] is not None and lo <= r["uniprot_position_1based"] <= hi]
    positions = [r["peptide_position_1based"] for r in hits]
    complete = len(hits) == hi - lo + 1 and positions == list(range(positions[0], positions[0] + len(positions))) if positions else False
    status = "intact_reference_sequence" if complete else ("partial_reference_sequence" if hits else "reference_sequence_absent")
    return {"feature_type": feature["type"], "feature_name": feature.get("description", ""), "uniprot_start_1based": lo, "uniprot_end_1based": hi, "reference_length_aa": hi - lo + 1, "reference_residues_present": len(hits), "peptide_first_mapped_residue": min(positions) if positions else None, "peptide_last_mapped_residue": max(positions) if positions else None, "sequence_status": status, "evidence": json.dumps(feature.get("evidences", []), separators=(",", ":")), "interpretation_boundary": "Sequence preservation or interruption; not an assay of folding or function"}


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, records):
    with path.open("w") as stream:
        for name, description, sequence in records:
            stream.write(f">{name} {description}\n")
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start + 80] + "\n")


def run(output_dir=OUT):
    plan, manifest = read_json(PLAN_PATH), read_json(SOURCES_PATH)
    if sha256(ROOT / plan["protocol"]) != plan["protocol_sha256"] or sha256(PLAN_PATH) != manifest["plan_sha256"]:
        raise ValueError("Analysis plan or protocol changed")
    for path, expected in plan["reference_pins"].items():
        if sha256(ROOT / path) != expected:
            raise ValueError(f"Earlier source or analysis changed: {path}")
    for source in manifest["sources"]:
        path = ROOT / source["local_cache_path"]
        if sha256(path) != source["sha256"] or path.stat().st_size != source["bytes"]:
            raise ValueError(f"Source integrity failed: {source['source_id']}")
    source_ids = {s["source_id"] for s in manifest["sources"]}
    if not {"ensembl-gene", "ensembl-gene-sequence", "ensembl-variant", "uniprot-entry"} <= source_ids:
        raise ValueError("A required input lacks a source pin")
    gene = read_json(RAW / "ensembl-gene.json")
    genome = read_json(RAW / "ensembl-gene-sequence.json")["seq"].upper()
    if len(genome) != gene["end"] - gene["start"] + 1 or set(genome) - set("ACGT"):
        raise ValueError("Invalid full-gene reference sequence")
    variant = plan["variant"]
    if genome[variant["position_1based"] - gene["start"]] != variant["reference"]:
        raise ValueError("Variant reference allele mismatch")
    variation_record = read_json(RAW / "ensembl-variant.json")
    matches = [m for m in variation_record["mappings"] if m["assembly_name"] == "GRCh38" and m["seq_region_name"] == "21" and m["start"] == variant["position_1based"] and m["end"] == variant["position_1based"]]
    if not any(m["allele_string"].split("/")[0] == "A" and "C" in m["allele_string"].split("/") for m in matches):
        raise ValueError("Current variant record does not confirm the A/C coordinate")
    selected, inventory, references, reference_checks = [], [], {}, []
    for transcript in sorted(gene["Transcript"], key=lambda t: t["id"]):
        exons = ordered_exons(transcript)
        indices = [i for i in range(len(exons) - 1) if exons[i]["end"] == plan["canonical_upstream_exon_end_1based"] and exons[i + 1]["start"] == plan["downstream_exon_start_1based"]]
        if len(indices) > 1:
            raise ValueError("Repeated target junction within a transcript")
        tid = transcript["id"]
        inventory.append({"transcript_id": tid, "version": transcript["version"], "display_name": transcript["display_name"], "biotype": transcript["biotype"], "is_canonical": transcript["is_canonical"], "selected": bool(indices), "reason": "exact_consecutive_exon_boundaries" if indices else "lacks_exact_consecutive_exon_boundaries", "exon_count": len(exons), "reference_cdna_length": transcript["length"], "annotated_protein_length": transcript.get("Translation", {}).get("length"), "target_upstream_exon_number": indices[0] + 1 if indices else None})
        if not indices:
            continue
        expected = next((x for x in plan["expected_transcripts"] if x["id"] == tid), None)
        if expected is None or transcript["version"] != expected["version"] or transcript["biotype"] != expected["biotype"]:
            raise ValueError(f"Unexpected qualifying transcript or version: {tid}")
        for suffix in ("cdna", "cds", "protein"):
            if tid + "-" + suffix not in source_ids:
                raise ValueError(f"Unpinned sequence input: {tid}-{suffix}")
        cdna, positions, blocks = reconstruct_reference(genome, gene["start"], exons)
        source_cdna = read_json(RAW / (tid + "-cdna.json"))["seq"]
        source_cds = read_json(RAW / (tid + "-cds.json"))["seq"]
        source_protein = read_json(RAW / (tid + "-protein.json"))["seq"]
        if cdna != source_cdna or len(cdna) != transcript["length"]:
            raise ValueError(f"Genomic exon concatenation differs from source cDNA: {tid}")
        start, translated, convention = reference_coding_bounds(cdna, positions, transcript["Translation"], source_cds, source_protein)
        reference_checks.append({"transcript_id": tid, "version": transcript["version"], "cdna_matches_genomic_exons": True, "cds_matches_genomic_annotation": True, "protein_matches_translation": True, "cds_includes_terminal_stop": True, "cds_length_nt": len(source_cds), "protein_length_aa": len(source_protein), "cds_coordinate_convention": convention})
        references[tid] = {"transcript": transcript, "cdna": cdna, "positions": positions, "blocks": blocks, "start": start, "translation": translated, "insert_at": blocks[indices[0]]["cdna_end_0based_exclusive"]}
        selected.append(tid)
    if set(selected) != {t["id"] for t in plan["expected_transcripts"]}:
        raise ValueError("Selected transcript set changed")
    uniprot = read_json(RAW / "uniprot-entry.json")
    uniprot_sequence = uniprot["sequence"]["value"]
    exact_anchors = [tid for tid, ref in references.items() if ref["translation"]["peptide"] == uniprot_sequence]
    canonical = [tid for tid in exact_anchors if references[tid]["transcript"]["is_canonical"]]
    if len(canonical) != 1:
        raise ValueError("No unique canonical Ensembl reference exactly matching UniProt")
    anchor = references[canonical[0]]
    anchor_codons = {tuple(anchor["positions"][anchor["start"] + 3 * i:anchor["start"] + 3 * i + 3]): i + 1 for i in range(len(uniprot_sequence))}
    if len(anchor_codons) != len(uniprot_sequence):
        raise ValueError("Non-unique reference genomic codons")
    features = [f for f in uniprot.get("features", []) if f["type"] in ("Domain", "Region", "Active site", "Binding site")]
    segment_start, segment_end = plan["retained_interval_1based_inclusive"]
    segment_a = genome[segment_start - gene["start"]:segment_end - gene["start"] + 1]
    local_variant = variant["position_1based"] - segment_start
    if len(segment_a) != 29 or segment_a[local_variant] != "A":
        raise ValueError("Retained segment coordinates or reference allele do not match")
    segment_c = segment_a[:local_variant] + "C" + segment_a[local_variant + 1:]
    results, exon_rows, feature_rows, alignment_rows = [], [], [], []
    cdna_fasta, orf_fasta, peptide_fasta = [], [], []
    scenario_cache = {}
    for tid in selected:
        ref = references[tid]
        transcript = ref["transcript"]
        versioned = tid + "." + str(transcript["version"])
        for scenario in plan["scenarios"]:
            key = versioned + "__" + scenario["key"]
            retained = scenario["splice"] == "plus29"
            segment = segment_a if scenario["allele"] == "A" else segment_c
            if retained:
                cdna, positions, blocks = insert_retained_segment(ref["cdna"], ref["positions"], ref["blocks"], ref["insert_at"], segment, segment_start)
                start = ref["start"] + (len(segment) if ref["insert_at"] <= ref["start"] else 0)
            else:
                cdna, positions, blocks, start = ref["cdna"], ref["positions"], ref["blocks"], ref["start"]
            translated = translate_from_start(cdna, start)
            peptide, baseline = translated["peptide"], ref["translation"]["peptide"]
            prefix = common_prefix_length(peptide, baseline)
            changed = peptide != baseline
            if not retained:
                consequence = "normal_splicing_variant_removed"
            elif ref["insert_at"] <= ref["start"]:
                consequence = "five_prime_UTR_retention"
            elif ref["insert_at"] >= ref["translation"]["stop_end_0based_exclusive"]:
                consequence = "three_prime_UTR_retention"
            else:
                consequence = "coding_29nt_retention_frameshift"
            stop_start, stop_end = translated["stop_start_0based"], translated["stop_end_0based_exclusive"]
            row = {"scenario_id": key, "transcript_id": tid, "transcript_version": transcript["version"], "biotype": transcript["biotype"], "is_canonical": transcript["is_canonical"], "scenario": scenario["key"], "allele": scenario["allele"], "splice": scenario["splice"], "consequence": consequence, "cdna_length_nt": len(cdna), "added_nt": len(cdna) - len(ref["cdna"]), "retained_segment": segment if retained else "", "retention_after_reference_cdna_base_1based": ref["insert_at"] if retained else None, "retention_after_reference_coding_base_1based": ref["insert_at"] - ref["start"] if retained else None, "variant_cdna_position_1based": ref["insert_at"] + local_variant + 1 if retained else None, "cds_start_cdna_1based": start + 1, "reference_protein_length_aa": len(baseline), "protein_length_aa": len(peptide), "protein_changed": changed, "first_changed_amino_acid_1based": prefix + 1 if changed else None, "unchanged_prefix_aa": prefix, "novel_peptide_tail": peptide[prefix:] if changed else "", "reference_suffix_replaced_or_removed_aa": len(baseline) - prefix if changed else 0, "stop_codon": translated["stop_codon"], "stop_start_cdna_1based": stop_start + 1 if stop_start is not None else None, "stop_end_cdna_1based": stop_end, "stop_start_coding_1based": stop_start - start + 1 if stop_start is not None else None, "stop_genomic_positions_1based": ";".join(map(str, positions[stop_start:stop_end])) if stop_start is not None else "", "no_stop_in_transcript": translated["no_stop_in_transcript"], "cdna_sha256": hashlib.sha256(cdna.encode()).hexdigest(), "protein_sha256": hashlib.sha256(peptide.encode()).hexdigest()}
            row.update(nmd_features(translated, start, blocks, plan["nmd_features"]))
            results.append(row)
            scenario_cache[(tid, scenario["key"])] = (cdna, peptide)
            description = f"splice={scenario['splice']} allele={scenario['allele']} conditional_sequence_only"
            cdna_fasta.append((key, "mature_RNA_in_DNA_alphabet " + description, cdna))
            orf_fasta.append((key, "annotated_start_to_first_stop_inclusive " + description, translated["orf_with_stop"]))
            peptide_fasta.append((key, "terminal_stop_omitted " + description, peptide))
            exon_rows.extend({"scenario_id": key, **block} for block in blocks)
            alignment = codon_alignment(peptide, positions, start, anchor_codons, uniprot_sequence)
            alignment_rows.extend({"scenario_id": key, **r} for r in alignment)
            feature_rows.extend({"scenario_id": key, "transcript_id": tid, "scenario": scenario["key"], "uniprot_accession": uniprot["primaryAccession"], **summarize_feature(feature, alignment)} for feature in features)
    for tid in selected:
        if scenario_cache[(tid, "normal_A")] != scenario_cache[(tid, "normal_C")]:
            raise ValueError("Normal splicing unexpectedly retains an intronic allele difference")
        a, c = scenario_cache[(tid, "plus29_A")][0], scenario_cache[(tid, "plus29_C")][0]
        mismatches = [i for i, (x, y) in enumerate(zip(a, c)) if x != y]
        if len(a) != len(c) or mismatches != [references[tid]["insert_at"] + local_variant]:
            raise ValueError("Retained allele scenarios differ outside the target variant")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "transcript-inventory.csv", inventory)
    write_csv(output_dir / "consequences.csv", results)
    write_csv(output_dir / "exon-coordinates.csv", exon_rows)
    write_csv(output_dir / "protein-features.csv", feature_rows)
    write_csv(output_dir / "protein-reference-alignment.csv", alignment_rows)
    write_fasta(output_dir / "transcripts.fasta", cdna_fasta)
    write_fasta(output_dir / "orfs-with-stop.fasta", orf_fasta)
    write_fasta(output_dir / "proteins.fasta", peptide_fasta)
    write_fasta(output_dir / "retained-segments.fasta", [("plus29_A", "GRCh38_21:42434955-42434983_reference", segment_a), ("plus29_C", "GRCh38_21:42434955-42434983_rs1893592_A_to_C", segment_c)])
    audit = {"analysis_id": plan["plan_id"], "plan_sha256": sha256(PLAN_PATH), "source_manifest_sha256": sha256(SOURCES_PATH), "script_sha256": sha256(Path(__file__)), "reference_checks": reference_checks, "source_count": len(manifest["sources"]), "gene_transcripts_inventoried": len(inventory), "selected_transcripts": selected, "scenario_count": len(results), "retained_segments": {"A": segment_a, "C": segment_c}, "uniprot": {"accession": uniprot["primaryAccession"], "audit": uniprot["entryAudit"], "canonical_length_aa": len(uniprot_sequence), "exact_Ensembl_matches": exact_anchors, "anchor_transcript": canonical[0], "projection_method": "Unique exact genomic codon triplets plus identical amino acids; unmatched codons remain unmapped."}, "output_counts": {"protein_features": len(feature_rows), "protein_alignment_residues": len(alignment_rows), "exon_scenarios": len(exon_rows)}, "runtime_python": platform.python_version(), "files_sha256": {p.name: sha256(p) for p in sorted(output_dir.iterdir()) if p.is_file() and p.name != "audit.json"}, "scope": "Deterministic sequence consequences; NMD features are conditional flags, not measured degradation, protein abundance, allele effects or treatment outcomes."}
    (output_dir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    audit = run(args.output_dir)
    print(json.dumps({"selected_transcripts": audit["selected_transcripts"], "scenarios": audit["scenario_count"], "reference_checks_passed": len(audit["reference_checks"]), "output_directory": str(args.output_dir)}, indent=2))


if __name__ == "__main__":
    main()
