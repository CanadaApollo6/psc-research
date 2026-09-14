"""Reconstruct the fixed T29 catalog RNA and conditional translation products."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

import analyze_ubash3a_consequences as prior

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-alternative-ending-plan.json"
OUT = ROOT / "data/derived/ubash3a-alternative-ending"


def digest_peptides(sequence, missed_options=(0, 1), lo=7, hi=35):
    cuts = [0] + [i+1 for i, aa in enumerate(sequence) if aa in "KR" and (i+1 == len(sequence) or sequence[i+1] != "P")]
    if cuts[-1] != len(sequence):
        cuts.append(len(sequence))
    for i, start in enumerate(cuts[:-1]):
        for missed in missed_options:
            j = i+missed+1
            if j >= len(cuts):
                continue
            end = cuts[j]
            if lo <= end-start <= hi:
                yield start, end, missed, sequence[start:end]


def run(output=OUT):
    plan = json.loads(PLAN.read_text())
    for path, pin in plan["input_pins"].items():
        assert prior.sha256(ROOT / path) == pin, path
    targets = list(csv.DictReader((ROOT / "data/derived/ubash3a-cd4-experiment/assay-targets.csv").open()))
    target = next(t for t in targets if t["target_label"] == plan["candidate_label"])
    assert target["source_model"] == plan["candidate_source_model"]
    intervals = json.loads(target["exons_1based_inclusive"])
    exons = [{"start": a, "end": b, "id": f"T29_exon_{i}"} for i,(a,b) in enumerate(intervals,1)]
    gene = json.loads((ROOT / "data/raw/ubash3a-consequences/ensembl-gene-sequence.json").read_text())["seq"].upper()
    start = plan["gene_reference_start_1based"]
    cdna, coordinates, blocks = prior.reconstruct_reference(gene, start, exons)
    coding_start = coordinates.index(plan["conditional_coding_start_1based"])
    variant_index = coordinates.index(plan["allele_position_1based"])
    assert cdna[coding_start:coding_start+3] == "ATG"
    assert cdna[variant_index] == "A"
    old_proteins = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / "data/derived/ubash3a-consequences/proteins.fasta", "fasta")}
    assert len(old_proteins) == 20
    old_map = list(csv.DictReader((ROOT / "data/derived/ubash3a-consequences/protein-reference-alignment.csv").open()))
    anchor_codons = {tuple(map(int,r["genomic_codon_positions_1based"].split(";"))): int(r["uniprot_position_1based"])
                     for r in old_map if r["scenario_id"] == "ENST00000319294.11__normal_A" and r["uniprot_position_1based"]}
    uniprot = json.loads((ROOT / "data/raw/ubash3a-consequences/uniprot-entry.json").read_text())
    uniprot_sequence = uniprot["sequence"]["value"]
    settings = json.loads((ROOT / "config/ubash3a-consequence-plan.json").read_text())["nmd_features"]
    rnas, orfs, proteins = [], [], []
    consequences, comparisons, codons, features, peptides, checks = [], [], [], [], [], []
    independent_coordinates = [position for a,b in intervals for position in range(a,b+1)]
    assert independent_coordinates == coordinates
    independent_cdna = "".join(gene[position-start] for position in independent_coordinates)
    assert independent_cdna == cdna
    for allele in plan["hypothetical_alleles"]:
        scenario = f"T29__{allele}"
        rna = cdna[:variant_index] + allele + cdna[variant_index+1:]
        translation = prior.translate_from_start(rna, coding_start)
        peptide = translation["peptide"]
        independent_rna = "".join(allele if p == plan["allele_position_1based"] else gene[p-start] for p in independent_coordinates)
        whole_codons = independent_rna[coding_start:coding_start+(len(independent_rna)-coding_start)//3*3]
        independent_translation = str(Seq(whole_codons).translate())
        independent_peptide = independent_translation.split("*",1)[0]
        assert independent_rna == rna and independent_peptide == peptide
        if "*" in independent_translation:
            first_stop = coding_start + 3*independent_translation.index("*")
            assert first_stop == translation["stop_start_0based"]
        else:
            assert translation["no_stop_in_transcript"]
        mapping = prior.codon_alignment(peptide, coordinates, coding_start, anchor_codons, uniprot_sequence)
        codons.extend({"scenario": scenario, **r} for r in mapping)
        for feature in uniprot.get("features", []):
            if feature["type"] in ("Domain", "Region", "Active site", "Binding site"):
                features.append({"scenario": scenario, **prior.summarize_feature(feature,mapping)})
        stop = translation["stop_start_0based"]
        stop_positions = coordinates[stop:stop+3] if stop is not None else []
        consequence = {"scenario": scenario, "source_model": target["source_model"], "hypothetical_allele": allele,
                       "translation_start_status": "assumed_shared_reference_ATG_not_measured_for_T29",
                       "cdna_length_nt": len(rna), "coding_start_cdna_1based": coding_start+1,
                       "protein_length_aa": len(peptide), "stop_codon": translation["stop_codon"],
                       "stop_cdna_start_1based": stop+1 if stop is not None else None,
                       "stop_genomic_positions_1based": ";".join(map(str,stop_positions)),
                       "no_stop_in_transcript": translation["no_stop_in_transcript"],
                       "rna_sha256": hashlib.sha256(rna.encode()).hexdigest(),
                       "protein_sha256": hashlib.sha256(peptide.encode()).hexdigest(),
                       **prior.nmd_features(translation,coding_start,blocks,settings)}
        consequences.append(consequence)
        for name, previous in old_proteins.items():
            prefix = prior.common_prefix_length(previous,peptide)
            comparisons.append({"scenario": scenario, "previous_scenario": name, "previous_length_aa": len(previous),
                                "candidate_length_aa": len(peptide), "unchanged_prefix_aa": prefix,
                                "identical_full_protein": peptide==previous,
                                "candidate_tail_after_first_difference": peptide[prefix:]})
        for begin,end,missed,sequence in digest_peptides(peptide,tuple(plan["missed_cleavages"]),*plan["tryptic_peptide_lengths"]):
            mass_sequence = sequence.replace("I","L")
            matches = [name for name,seq in old_proteins.items() if mass_sequence in seq.replace("I","L")]
            peptides.append({"scenario": scenario, "peptide": sequence, "start_aa_1based": begin+1,
                             "end_aa_1based": end, "missed_cleavages": missed,
                             "IL_equivalent_matches_among_previous_products": ";".join(matches),
                             "absent_from_previous_twenty_products": not matches,
                             "human_proteome_specificity": "not_yet_tested", "mass_spectrum_detection": "not_tested"})
        checks.append({"scenario": scenario, "independent_RNA_equal": True, "independent_Biopython_translation_equal": True,
                       "independent_first_stop_equal": True, "bases_checked": len(rna), "amino_acids_checked": len(peptide)})
        rnas.append((scenario,"reference_exon_reconstruction_hypothetical_allele_no_measured_ends",rna))
        orfs.append((scenario,"assumed_shared_start_first_stop_included",translation["orf_with_stop"]))
        proteins.append((scenario,"conditional_translation_not_observed_protein",peptide))
    assert [i for i,(a,b) in enumerate(zip(rnas[0][2],rnas[1][2])) if a!=b] == [variant_index]
    output.mkdir(parents=True, exist_ok=True)
    for name,rows in [("consequences.csv",consequences),("previous-product-comparisons.csv",comparisons),
                      ("protein-codon-map.csv",codons),("protein-features.csv",features),("tryptic-candidates.csv",peptides),
                      ("exon-coordinates.csv",blocks)]:
        prior.write_csv(output/name,rows)
    for name,records in [("transcripts.fasta",rnas),("orfs-with-stop.fasta",orfs),("proteins.fasta",proteins)]:
        prior.write_fasta(output/name,records)
    primary_names = ["consequences.csv","previous-product-comparisons.csv","protein-codon-map.csv","protein-features.csv",
                     "tryptic-candidates.csv","exon-coordinates.csv","transcripts.fasta","orfs-with-stop.fasta","proteins.fasta"]
    record = {"plan_sha256": prior.sha256(PLAN), "script_sha256": prior.sha256(Path(__file__)),
              "prior_translation_library_sha256": prior.sha256(ROOT/"scripts/analyze_ubash3a_consequences.py"),
              "verification": checks, "allele_RNAs_differ_only_at_fixed_position": True,
              "files_sha256": {name:prior.sha256(output/name) for name in primary_names},
              "distinct_candidate_peptides_absent_from_prior_products_IL_equivalent": len({r['peptide'].replace('I','L') for r in peptides if r['absent_from_previous_twenty_products']}),
              "biological_boundary": "Conditional sequence results. No observed translation initiation, protein, decay, donor genotype, or human-proteome-specific peptide."}
    (output/"audit.json").write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({"consequences":consequences,"verification":checks,
                      "peptides_absent_from_prior_products":record["distinct_candidate_peptides_absent_from_prior_products_IL_equivalent"]},indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=OUT)
    run(parser.parse_args().output_dir)
