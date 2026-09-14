"""Check source support, sequence discriminators and the sparse BAM result."""

import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import struct

import pysam

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/ubash3a-hypothesis-test"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(path):
    data, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            assert name not in data
            data[name] = ""
        elif line.strip():
            data[name] += line.strip()
    return data


def index_stats(path):
    data = path.read_bytes()
    assert data[:4] == b"BAI\x01"
    nref = struct.unpack_from("<i", data, 4)[0]
    offset, counts = 8, []
    for _ in range(nref):
        nbins = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        mapped = None
        for _ in range(nbins):
            number, nchunks = struct.unpack_from("<Ii", data, offset)
            offset += 8
            if number == 37450:
                assert nchunks == 2
                mapped = struct.unpack_from("<Q", data, offset + 16)[0]
            offset += 16 * nchunks
        nintervals = struct.unpack_from("<i", data, offset)[0]
        offset += 4 + 8 * nintervals
        counts.append(mapped)
    assert len(data) - offset in {0, 8}
    return counts


def main():
    source_path = ROOT / "data/raw/ubash3a-hypothesis-test/trails-isoform_counts_matrix.tsv"
    requested = list(csv.DictReader((OUT / "trails-locus-source-counts.csv").open()))
    wanted = {r["source_matrix_id"] for r in requested}
    source_values = {}
    with source_path.open() as source:
        header = next(source).rstrip("\n").split("\t")
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in wanted:
                assert fields[0] not in source_values
                source_values[fields[0]] = dict(zip(header[1:], fields[1:], strict=True))
    for row in requested:
        assert row["source_count_value"] == source_values[row["source_matrix_id"]][row["cell_subset"]]
    positive = [r for r in requested if r["event"] == "plus29"]
    support = {"all_subsets_source_units": str(sum(Fraction(r["source_count_value"]) for r in positive)),
               "CD4_source_units": str(sum(Fraction(r["source_count_value"]) for r in positive if r["primary_CD4_context"] == "True")),
               "positive_CD4_subsets": sum(Fraction(r["source_count_value"]) > 0 for r in positive if r["primary_CD4_context"] == "True")}
    primary_support = json.loads((OUT / "support-and-discriminators-audit.json").read_text())["plus29_source_support"][0]
    assert support["all_subsets_source_units"] == primary_support["sum_source_count_units_all_29_subsets"]
    assert support["CD4_source_units"] == primary_support["sum_source_count_units_CD4_subsets"]
    assert support["positive_CD4_subsets"] == primary_support["CD4_subsets_with_positive_source_value"]
    proteins = {k: v for k, v in fasta(ROOT / "data/derived/ubash3a-consequences/proteins.fasta").items() if "__plus29_" in k}
    protein_groups = list(csv.DictReader((OUT / "protein-sequence-equivalence.csv").open()))
    assigned = set()
    for group in protein_groups:
        names = group["reference_scenarios"].split(";")
        seq = proteins[names[0]]
        assert all(proteins[name] == seq for name in names)
        assert {name for name, candidate in proteins.items() if candidate == seq} == set(names)
        assert len(seq) == int(group["protein_length_aa"])
        assert hashlib.sha256(seq.encode()).hexdigest() == group["sequence_sha256"]
        assert not assigned.intersection(names)
        assigned.update(names)
    assert assigned == set(proteins)
    windows = fasta(OUT / "diagnostic-RNA-windows.fasta")
    rnarefs = fasta(ROOT / "data/derived/ubash3a-consequences/transcripts.fasta")
    genome = json.loads((ROOT / "data/raw/ubash3a-consequences/ensembl-gene-sequence.json").read_text())
    gene_origin = int(genome["id"].split(":")[3])
    window_rows = list(csv.DictReader((OUT / "diagnostic-RNA-windows.csv").open()))
    for row in window_rows:
        seq = windows[row["window_id"]]
        assert len(seq) == int(row["window_nt"])
        assert hashlib.sha256(seq.encode()).hexdigest() == row["sequence_sha256"]
        alleles = ["A", "C"] if row["constructed_allele"] == "not_covered" else [row["constructed_allele"]]
        if row["source_model"].startswith("ENST"):
            for allele in alleles:
                assert seq in rnarefs[row["source_model"] + "__plus29_" + allele]
        else:
            exons = [(int(r["start_1based"]), int(r["end_1based"])) for r in csv.DictReader((OUT / "observed-exons.csv").open()) if r["source_id"] == "trails-TRAILS-gtf-gz" and r["source_transcript_id"] == row["source_model"]]
            coords = [p for a, b in exons for p in range(a, b + 1)]
            for allele in alleles:
                full = "".join(allele if p == 42434957 else genome["seq"][p - gene_origin].upper() for p in coords)
                assert seq in full
        intervals = json.loads(row["genomic_blocks_1based_inclusive"])
        assert any(a <= 42434957 <= b for a, b in intervals) == (row["rs1893592_covered"] == "True")
        assert (intervals[0][1], intervals[1][0]) == (42434983, 42437488)
        assert [intervals[-2][1], intervals[-1][0]] == json.loads(row["discriminating_junction_exon_end_next_start"])
    reads = json.loads((ROOT / "config/ubash3a-upf1-read-sources.json").read_text())
    reported = list(csv.DictReader((OUT / "upf1-locus-alignments.csv").open()))
    sample_summaries = {r["sample_id"]: r for r in csv.DictReader((OUT / "upf1-read-summary.csv").open())}
    bam_audits, alignment_checks = [], 0
    for sample in reads["samples"]:
        stats = index_stats(ROOT / sample["source_index_path"])
        header_path = next(ROOT / n for n in sample["local_files"] if n.endswith(".header.sam"))
        sq = [dict(field.split(":", 1) for field in line.split("\t")[1:]) for line in header_path.read_text().splitlines() if line.startswith("@SQ\t")]
        assert len(stats) == len(sq)
        chrom_index = next(i for i, row in enumerate(sq) if row["SN"].removeprefix("chr") == "21")
        assert int(sq[chrom_index]["LN"]) == 46709983
        assert stats[chrom_index] is not None and stats[chrom_index] > 0
        sam_rows = [line.split("\t") for line in pysam.samtools.view(str(ROOT / sample["local_bam"])).splitlines()]
        assert len(sam_rows) == sample["all_locus_alignments"] == int(sample_summaries[sample["sample_id"]]["all_locus_alignments"])
        for i, fields in enumerate(sam_rows, 1):
            actual = next(r for r in reported if r["sample_id"] == sample["sample_id"] and int(r["read_index"]) == i)
            pos = int(fields[3])
            junctions = []
            for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
                length = int(length)
                if operation == "N":
                    junctions.append([pos - 1, pos + length])
                if operation in "MDN=X":
                    pos += length
            assert actual["query_name"] == fields[0] and int(actual["flag"]) == int(fields[1])
            assert actual["cigar"] == fields[5] and int(actual["mapq"]) == int(fields[4])
            assert int(actual["end_1based"]) == pos - 1
            assert json.loads(actual["junctions_json"]) == junctions
            # All actual retained source alignments are unspliced in this audit.
            assert not junctions and actual["event_from_coordinates"] == "other"
            assert actual["full_chain_reference_matches"] == actual["downstream_chain_reference_matches"] == ""
            alignment_checks += 1
        bam_audits.append({"sample_id": sample["sample_id"], "original_index_chr21_mapped_alignment_count": stats[chrom_index],
                           "UBASH3A_subset_alignments": len(sam_rows), "GRCh38_chr21_length_checked": True})
    assert alignment_checks == len(reported)
    result = {"status": "pass", "script_sha256": sha(Path(__file__)), "source_count_values_exactly_checked": len(requested),
              "plus29_support_independent_fraction_sum": support, "protein_scenarios_checked": len(proteins),
              "protein_groups_checked": len(protein_groups), "RNA_windows_checked": len(windows),
              "RNA_verification": "Reference windows occur in both appropriate prior reconstructed RNAs; new catalog windows occur in independently assembled genomic exon sequences. Alleles remain hypothetical.",
              "BAM_verification": "Independent binary BAI metadata parse and samtools SAM-text CIGAR parse; source indexes contain chromosome-21 mappings in every library. Original full BAMs were not downloaded or hashed.",
              "BAM_samples": bam_audits, "actual_alignments_checked": alignment_checks}
    (ROOT / "reports/ubash3a-test-support-independent-check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "pass", "count_values": len(requested), "RNA_windows": len(windows), "BAM_samples": len(bam_audits)}))


if __name__ == "__main__":
    main()
