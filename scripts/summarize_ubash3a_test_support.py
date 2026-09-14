"""Preserve source support values and prepare RNA/protein discriminators."""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-hypothesis-test"
CD4 = ["NaiveCD4", "Th1", "Th2", "Th17", "Tfh", "Fra1", "Fra2.aTreg", "Fra3", "LAG3Treg", "MemoryCD4", "Thx"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, rows):
    assert rows
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output=OUT):
    output.mkdir(parents=True, exist_ok=True)
    parent = json.loads((ROOT / "config/ubash3a-hypothesis-test-plan.json").read_text())
    for name, pin in parent["reference_pins"].items():
        assert sha(ROOT / name) == pin
    sources = json.loads((ROOT / "config/ubash3a-hypothesis-test-sources.json").read_text())["sources"]
    by_id = {s["source_id"]: s for s in sources}
    for key in ["trails-isoform_counts_matrix-tsv", "trails-isoform-info"]:
        assert sha(ROOT / by_id[key]["local_cache_path"]) == by_id[key]["sha256"]
    models = {r["source_transcript_id"]: r for r in csv.DictReader((OUT / "observed-transcript-models.csv").open()) if r["source_id"] == "trails-TRAILS-gtf-gz"}
    count_rows, count_ids, source_matrix_rows, plus_summaries = [], set(), 0, []
    with (ROOT / by_id["trails-isoform_counts_matrix-tsv"]["local_cache_path"]).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames[1:]
        assert columns[:11] == CD4 and len(columns) == 29
        for row in reader:
            source_matrix_rows += 1
            tid, _, gene = row["ids"].partition("_ENSG")
            if tid not in models:
                continue
            assert tid not in count_ids and gene
            count_ids.add(tid)
            nums = {cell: Decimal(row[cell]) for cell in columns}
            assert all(v.is_finite() and v >= 0 for v in nums.values())
            for cell in columns:
                count_rows.append({"source_transcript_id": tid, "source_matrix_id": row["ids"], "event": models[tid]["event"],
                                   "eligible_chr21_positive": models[tid]["eligible_chr21_positive"], "cell_subset": cell,
                                   "primary_CD4_context": cell in CD4, "source_count_value": row[cell]})
            if models[tid]["event"] == "plus29":
                plus_summaries.append({"source_transcript_id": tid, "sum_source_count_units_all_29_subsets": str(sum(nums.values())),
                                       "sum_source_count_units_CD4_subsets": str(sum(nums[c] for c in CD4)),
                                       "CD4_subsets_with_positive_source_value": sum(nums[c] > 0 for c in CD4),
                                       "CD4_subsets": len(CD4), "unique_donors": 1})
    info_rows, info_ids = [], set()
    with gzip.open(ROOT / by_id["trails-isoform-info"]["local_cache_path"], "rt") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            tid = row["isoform"]
            if tid in models:
                assert tid not in info_ids
                info_ids.add(tid)
                assert int(row["start"]) + 1 == int(models[tid]["first_exon_start_1based"])
                assert int(row["end"]) == int(models[tid]["last_exon_end_1based"])
                assert row["strand"] == models[tid]["strands"]
                info_rows.append(row)
    write(output / "trails-locus-source-counts.csv", count_rows)
    write(output / "trails-locus-source-annotations.csv", info_rows)
    protein_groups = defaultdict(list)
    for record in SeqIO.parse(ROOT / "data/derived/ubash3a-consequences/proteins.fasta", "fasta"):
        if "__plus29_" in record.id:
            protein_groups[str(record.seq)].append(record.id)
    protein_rows = [{"protein_length_aa": len(seq), "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                     "number_of_reference_RNAs_encoding_this_sequence": len(ids), "reference_scenarios": ";".join(sorted(ids)),
                     "interpretation": "Conditional translation products from the prior reconstruction; a peptide cannot identify which of these RNAs produced it."}
                    for seq, ids in sorted(protein_groups.items(), key=lambda x: (len(x[0]), x[1]))]
    write(output / "protein-sequence-equivalence.csv", protein_rows)
    old_sources = json.loads((ROOT / "config/ubash3a-consequence-sources.json").read_text())["sources"]
    genome_source = next(s for s in old_sources if s["file"] == "ensembl-gene-sequence.json")
    assert sha(ROOT / genome_source["local_cache_path"]) == genome_source["sha256"]
    genome_record = json.loads((ROOT / genome_source["local_cache_path"]).read_text())
    genome = genome_record["seq"].upper()
    gene_start = int(genome_record["id"].split(":")[3])
    assert genome[42434957 - gene_start] == "A"
    refs = list(csv.DictReader((OUT / "reference-junction-chains.csv").open()))
    candidates = []
    for label, tid in [("canonical_downstream", "ENST00000319294.11"), ("boundary", "ENST00000398367.2"), ("other_NMD_reference", "ENST00000635325.1")]:
        ref = next(r for r in refs if r["transcript"] == tid and r["event"] == "plus29")
        candidates.append({"label": label, "source_model": tid, "exons": json.loads(ref["all_exons_1based_inclusive"]), "evidence": "prior conditional reference"})
    for model in models.values():
        if model["event"] == "plus29":
            exons = [(int(r["start_1based"]), int(r["end_1based"])) for r in csv.DictReader((OUT / "observed-exons.csv").open()) if r["source_id"] == "trails-TRAILS-gtf-gz" and r["source_transcript_id"] == model["source_transcript_id"]]
            candidates.append({"label": "TRAILS_catalog_candidate", "source_model": model["source_transcript_id"], "exons": exons, "evidence": "deposited catalog model; individual molecules not checked"})
    for item in candidates:
        item["chain"] = [(a[1], b[0]) for a, b in zip(item["exons"], item["exons"][1:])]
        item["target"] = item["chain"].index((42434983, 42437488))
        item["suffix"] = item["chain"][item["target"]:]
    windows, fasta = [], []
    for item in candidates:
        # Find the first downstream junction that distinguishes this architecture.
        offset = next(i for i in range(1, len(item["suffix"])) if all(other is item or other["suffix"][:i + 1] != item["suffix"][:i + 1] for other in candidates))
        donor_index = item["target"] + offset
        for left in [20, 39]:
            blocks = [(item["exons"][item["target"]][1] - left + 1, item["exons"][item["target"]][1])]
            blocks += [tuple(x) for x in item["exons"][item["target"] + 1:donor_index + 1]]
            blocks += [(item["exons"][donor_index + 1][0], item["exons"][donor_index + 1][0] + 19)]
            positions = [p for a, b in blocks for p in range(a, b + 1)]
            allele_covered = 42434957 in positions
            for allele in ["A", "C"] if allele_covered else ["not_covered"]:
                bases = [genome[p - gene_start] for p in positions]
                if allele_covered:
                    bases[positions.index(42434957)] = allele
                seq = "".join(bases)
                name = f"{item['label']}__left{left}__allele_{allele}"
                windows.append({"window_id": name, "source_model": item["source_model"], "evidence": item["evidence"],
                                "left_anchor_nt": left, "right_anchor_nt": 20, "window_nt": len(seq),
                                "rs1893592_covered": allele_covered, "constructed_allele": allele,
                                "discriminating_junction_exon_end_next_start": json.dumps(item["chain"][donor_index]),
                                "genomic_blocks_1based_inclusive": json.dumps(blocks, separators=(",", ":")),
                                "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                                "interpretation": "Reference-genome cDNA planning window; not an observed genotype or validated primer design. Does not determine the full upstream RNA chain."})
                fasta.append(">" + name + " constructed_cDNA_window\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)))
    write(output / "diagnostic-RNA-windows.csv", windows)
    (output / "diagnostic-RNA-windows.fasta").write_text("\n".join(fasta) + "\n")
    names = ["trails-locus-source-counts.csv", "trails-locus-source-annotations.csv", "protein-sequence-equivalence.csv", "diagnostic-RNA-windows.csv", "diagnostic-RNA-windows.fasta"]
    audit = {"source_matrix_rows": source_matrix_rows, "selected_locus_models": len(models),
             "count_rows": len(count_rows), "count_ids_matched": len(count_ids), "count_ids_missing": sorted(set(models) - count_ids),
             "info_ids_matched": len(info_ids), "info_ids_missing": sorted(set(models) - info_ids),
             "source_info_start_convention": "All matched source starts equal GTF first exon start minus one; source ends agree exactly.",
             "source_count_unit_boundary": "Original values from the deposited isoform_counts_matrix.tsv are preserved. Source FLAIR command requests a count table and also --tpm; the deposited table export/normalization step is not fully documented. Report source count units, not unique molecules, independent donors or validated full-length reads per subset.",
             "plus29_source_support": plus_summaries, "protein_equivalence_groups": len(protein_rows), "diagnostic_windows": len(windows),
             "genome_source_sha256": genome_source["sha256"], "script_sha256": sha(Path(__file__)),
             "files_sha256": {name: sha(output / name) for name in names}}
    (output / "support-and-discriminators-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: audit[k] for k in ["count_ids_matched", "info_ids_matched", "protein_equivalence_groups", "diagnostic_windows", "plus29_source_support"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    run(parser.parse_args().output_dir)
