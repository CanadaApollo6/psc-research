#!/usr/bin/env python3
"""Qualify complete sources before matching; build a reproducible search database."""
import argparse
from collections import Counter
import datetime
import gzip
import hashlib
import json
from pathlib import Path

from analyze_ubash3a_proteomics_specificity import boundary, fasta, match_positions, sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/ubash3a-proteomics-specificity"
RAW = ROOT / "data/raw/ubash3a-proteomics-specificity"
OUT = ROOT / "data/derived/ubash3a-proteomics-specificity"
LOCK = ROOT / "config/ubash3a-proteomics-reference-lock.json"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def lock_references():
    assert not LOCK.exists(), "Do not replace an existing prospective reference lock"
    download = json.loads((WORK / "uniprot-pagination-completion.json").read_text())
    assert download["complete"] and download["records"] == 232837 and download["page_count"] == 466
    assert download["release"] == "2026_03"
    for page in download["pages"]:
        assert sha256(ROOT / page["local_path"]) == page["sha256"]
        headers = {k.lower(): v for k, v in page["response"]["headers"].items()}
        assert headers["x-uniprot-release"] == "2026_03" and int(headers["x-total-results"]) == 232837
    uniprot = ROOT / download["local_path"]
    assert sha256(uniprot) == download["sha256"]
    canonical, isoforms, reviewed, unreviewed = 0, 0, 0, 0
    for header, sequence in fasta(uniprot):
        fields = header.split()[0].split("|")
        assert fields[0] in {"sp", "tr"} and len(fields) == 3
        is_isoform = "-" in fields[1]
        isoforms += is_isoform
        canonical += not is_isoform
        reviewed += fields[0] == "sp"
        unreviewed += fields[0] == "tr"
    assert (canonical, isoforms) == (210706, 22131)
    md5_path = RAW / "gencode.v50.MD5SUMS"
    checksum_line = [line for line in md5_path.read_text().splitlines() if line.split()[-1] == "gencode.v50.pc_translations.fa.gz"]
    assert len(checksum_line) == 1
    gencode = RAW / "gencode.v50.pc_translations.fa.gz"
    md5 = hashlib.md5(gencode.read_bytes()).hexdigest()
    assert md5 == checksum_line[0].split()[0] == "509a89a4c3bb8c594baac50d4d57e589"
    specs = [
        ("uniprot_2026_03_all_human_isoforms", "human_reference", uniprot, 232837,
         {"release": "2026_03", "canonical_records": canonical, "isoform_records": isoforms,
          "reviewed_records_including_isoforms": reviewed, "unreviewed_records": unreviewed,
          "release_date": download["release_date"], "source_url": "https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606&format=fasta&includeIsoform=true&size=500",
          "completeness": "All 466 returned pages, stable release and total headers, unique record IDs and terminal no-next-link"}),
        ("gencode_v50_comprehensive_translations", "human_reference", gencode, 382428,
         {"release": "50", "genome_build": "GRCh38.p14", "publisher_md5": md5,
          "source_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz",
          "scope": "Comprehensive PC, NMD, nonstop, IG, TR, polymorphic pseudogene and protein-coding LoF translations; chromosomes, scaffolds, patches and alternatives"}),
        ("gpm_crap_20190304", "contaminant", RAW / "crap.fasta", 116,
         {"source_url": "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta", "publisher_ftp_mdtm": "20190304195801", "publisher_bytes": 42092}),
        ("previous_twenty_models", "previous_models", ROOT / "data/derived/ubash3a-consequences/proteins.fasta", 20,
         {"scope": "Five reference models x normal/plus29 x hypothetical A/C; conditional proteins, not measurements"}),
        ("T29_conditional_alleles", "conditional_T29", ROOT / "data/derived/ubash3a-alternative-ending/proteins.fasta", 2,
         {"scope": "Two conditional translations under the shared-start assumption; protein existence unconfirmed"}),
    ]
    sources = []
    for set_id, group, path, expected, metadata in specs:
        ids = set()
        records, residues, unusual = 0, 0, Counter()
        for header, sequence in fasta(path):
            identifier = header.split()[0]
            assert identifier not in ids, (set_id, identifier)
            ids.add(identifier)
            records += 1
            residues += len(sequence)
            unusual.update(c for c in sequence if c not in "ACDEFGHIKLMNPQRSTVWY")
        assert records == expected
        sources.append({"set_id": set_id, "source_group": group, "local_path": str(path.relative_to(ROOT)),
                        "records": records, "residues": residues, "nonstandard_residues": dict(unusual),
                        "size_bytes": path.stat().st_size, "sha256": sha256(path), **metadata})
    record = {"locked_at_utc": now(), "locked_before_candidate_matching": True,
              "sources": sources, "total_records": sum(row["records"] for row in sources),
              "scope_complete_for_retrieved_sources": True, "uniprot_full_page_receipt": download,
              "qualifications": "No exhaustive population-variant or unannotated-protein coverage is implied"}
    LOCK.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"locked": True, "records": record["total_records"], "sources": [{"set": row["set_id"], "records": row["records"]} for row in sources]}), flush=True)


def build_database():
    lock = json.loads(LOCK.read_text())
    assert json.loads((ROOT / "reports/ubash3a-proteomics-specificity-verification.json").read_text())["all_checks_pass"]
    candidates = json.loads((ROOT / "config/ubash3a-proteomics-specificity-plan.json").read_text())["candidates"]
    destination = ROOT / "work/ubash3a-proteomics-spectra"
    destination.mkdir(parents=True, exist_ok=True)
    database = destination / "search-background.fasta"
    members = destination / "search-background-membership.jsonl.gz"
    sequences, gene_rows, collisions = {}, [], []
    record_count = 0
    with database.open("w") as output, members.open("wb") as compressed:
        with gzip.GzipFile(filename="", fileobj=compressed, mode="wb", mtime=0) as zipped:
            for source in lock["sources"]:
                path = ROOT / source["local_path"]
                assert sha256(path) == source["sha256"]
                for header, sequence in fasta(path):
                    record_count += 1
                    identity = hashlib.sha256(sequence.encode()).hexdigest()
                    identifier = "P_" + identity[:24]
                    if identifier in sequences:
                        assert sequences[identifier] == sequence
                    else:
                        sequences[identifier] = sequence
                        output.write(">" + identifier + "\n")
                        output.writelines(sequence[pos:pos + 80] + "\n" for pos in range(0, len(sequence), 80))
                    record = {"database_id": identifier, "set_id": source["set_id"], "source_group": source["source_group"],
                              "source_record_id": header.split()[0], "source_header": header,
                              "sequence_sha256": identity, "protein_length": len(sequence)}
                    zipped.write((json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode())
                    if source["source_group"] in {"previous_models", "conditional_T29"} or "GN=UBASH3A " in header + " " or "|UBASH3A|" in header:
                        gene_rows.append(record)
    for identifier, sequence in sequences.items():
        for peptide in candidates:
            partner = peptide[:-1][::-1] + peptide[-1]
            for start in match_positions(sequence, partner, True):
                end = start + len(partner)
                # Comet also treats explicit stop characters as boundaries.
                n_term = boundary(sequence, start) or (start > 0 and sequence[start - 1] == "*")
                c_term = boundary(sequence, end) or (end < len(sequence) and sequence[end] == "*")
                missed = sum(boundary(sequence, gap) for gap in range(start + 1, end))
                collisions.append({"candidate": peptide, "decoy_generating_sequence": partner,
                                   "matched_partner_with_IL_equivalence": sequence[start:end], "database_id": identifier,
                                   "start_1based": start + 1, "end_1based_inclusive": end,
                                   "full_tryptic": n_term and c_term, "missed_cleavages": missed,
                                   "eligible_under_search_digestion": n_term and c_term and missed <= 2})
    write_csv(OUT / "ubash3a-search-database-members.csv", gene_rows)
    collision_columns = ["candidate", "decoy_generating_sequence", "matched_partner_with_IL_equivalence", "database_id", "start_1based", "end_1based_inclusive", "full_tryptic", "missed_cleavages", "eligible_under_search_digestion"]
    write_csv(OUT / "candidate-decoy-collisions.csv", collisions, collision_columns)
    parameter = ROOT / "config/ubash3a-proteomics-comet.params"
    plan = ROOT / "config/ubash3a-proteomics-search-plan.json"
    assert sha256(parameter) == json.loads(plan.read_text())["parameter_sha256"]
    result = {"built_at_utc": now(), "before_any_public_spectrum_search": True,
              "source_records": record_count, "unique_exact_protein_sequences": len(sequences),
              "removed_exact_duplicates_with_memberships_preserved": record_count - len(sequences),
              "database_path": str(database.relative_to(ROOT)), "database_sha256": sha256(database),
              "database_bytes": database.stat().st_size,
              "membership_path": str(members.relative_to(ROOT)), "membership_sha256": sha256(members),
              "reference_lock_sha256": sha256(LOCK), "search_plan_sha256": sha256(plan), "parameter_sha256": sha256(parameter),
              "decoys": "Comet internal digestion-aware peptide reversal with final residue retained; concatenated competition; classic FASTA mode",
              "candidate_decoy_partner_matches": len(collisions),
              "candidate_eligible_decoy_collisions": sum(row["eligible_under_search_digestion"] for row in collisions),
              "source_sequence_transformations": "Exact full-sequence deduplication only; upper-case FASTA normalization, no ambiguous-residue deletion or stop removal",
              "nonstandard_sequence_scope": "Comet recognizes explicit stops as boundaries; unknown-residue-containing peptides are not evidence for a concrete amino-acid sequence"}
    (ROOT / "config/ubash3a-proteomics-search-database.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["lock", "database"])
    args = parser.parse_args()
    lock_references() if args.action == "lock" else build_database()
