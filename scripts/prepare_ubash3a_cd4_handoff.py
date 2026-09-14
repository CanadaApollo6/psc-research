"""Package an unexecuted laboratory proposal using previously verified references.

This does not analyze experimental data or generate validated primers.
Run from any directory with the repository's Python; standard library only.
"""

import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-cd4-experiment-design.json"
OUT = ROOT / "data/derived/ubash3a-cd4-experiment"
ZIP = ROOT / "reports/ubash3a-cd4-lab-handoff.zip"
RECORD = ROOT / "reports/ubash3a-cd4-experiment-preparation.json"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(name, fields, rows):
    path = OUT / name
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def main():
    plan = json.loads(PLAN.read_text())
    for item in plan["input_files"]:
        assert digest((ROOT / item["path"]).read_bytes()) == item["sha256"], item["path"]
    assert not plan["bench_experiment_started"]
    assert not plan["biological_measurements_collected"]
    assert not plan["confirmatory_protocol_frozen"]
    assert all(value is None for value in plan["confirmatory_freeze_required"].values())
    OUT.mkdir(parents=True, exist_ok=True)
    references = read_csv(ROOT / "data/derived/ubash3a-hypothesis-test/reference-junction-chains.csv")
    aliases = {
        "ENST00000398367.2": "B", "ENST00000319294.11": "C",
        "ENST00000291535.12": "U", "ENST00000635325.1": "R", "ENST00001131803.1": "C2",
    }
    roles = {"B": "primary_boundary", "C": "primary_canonical", "U": "secondary_upstream_matched",
             "R": "exploratory_preexisting_NMD_annotation", "C2": "exploratory_other_reference"}
    targets = []
    for row in references:
        alias = aliases[row["transcript"]]
        exons = json.loads(row["all_exons_1based_inclusive"])
        chain = [[left[1], right[0]] for left, right in zip(exons, exons[1:])]
        assert chain == json.loads(row["full_junction_chain_exon_end_next_start"])
        event = row["event"]
        expected = plan["plus29_junction"] if event == "plus29" else plan["normal_junction"]
        assert expected in chain
        targets.append({
            "target_label": alias + ("29" if event == "plus29" else "0"),
            "source_model": row["transcript"], "role": roles[alias],
            "evidence_level": "conditional_plus29_reference" if event == "plus29" else "reference_annotation",
            "assembly": "GRCh38", "chromosome": "21", "strand": "+", "event": event,
            "exons_1based_inclusive": json.dumps(exons, separators=(",", ":")),
            "junction_chain_exon_end_next_start": json.dumps(chain, separators=(",", ":")),
            "reference_coding_start_1based": row["coding_start_genomic_1based"],
            "reference_biotype": row["biotype"],
            "lab_validation_status": "not_performed", "new_experiment_donor_count": "",
            "interpretation": "Full original-molecule linkage required; local +29 alone is insufficient; not a validated primer or observed donor genotype.",
        })
    candidates = [row for row in read_csv(ROOT / "data/derived/ubash3a-hypothesis-test/observed-transcript-models.csv")
                  if row["source_transcript_id"] == plan["exploratory_trails_model"]]
    assert len(candidates) == 1
    candidate = candidates[0]
    exon_rows = [row for row in read_csv(ROOT / "data/derived/ubash3a-hypothesis-test/observed-exons.csv")
                 if row["source_id"] == candidate["source_id"] and row["source_transcript_id"] == candidate["source_transcript_id"]]
    exons = sorted([[int(row["start_1based"]), int(row["end_1based"])] for row in exon_rows])
    chain = [[left[1], right[0]] for left, right in zip(exons, exons[1:])]
    assert len(exons) == 11 and exons[-1] == [42437855, 42440816]
    assert chain == json.loads(candidate["junction_chain_exon_end_next_start"])
    targets.append({
        "target_label": "T29", "source_model": candidate["source_transcript_id"], "role": "exploratory_alternative_ending",
        "evidence_level": "published_catalog_model_prior_audit_no_new_molecules", "assembly": "GRCh38", "chromosome": "21",
        "strand": "+", "event": "plus29", "exons_1based_inclusive": json.dumps(exons, separators=(",", ":")),
        "junction_chain_exon_end_next_start": json.dumps(chain, separators=(",", ":")),
        "reference_coding_start_1based": "", "reference_biotype": "", "lab_validation_status": "not_performed",
        "new_experiment_donor_count": "", "interpretation": "Exploratory catalog candidate; full coding identity and ends require experimental confirmation. Cannot replace B29.",
    })
    assert len(targets) == len({row["target_label"] for row in targets}) == 11
    generated = [write_csv("assay-targets.csv", list(targets[0]), targets)]
    schemas = {
        "sample-manifest-template.csv": "sample_id,donor_code,study_phase,cell_state,culture_id,condition,perturbation_id,exposure_duration_h,chase_time_h,viable_cell_count,purity_fraction,RNA_preparation_id,RT_preparation_id,library_id,batch_id,raw_data_file,raw_data_sha256,QC_status,exclusion_reason",
        "molecule-evidence-template.csv": "sample_id,RT_preparation_id,library_id,original_molecule_family_id,raw_read_ids,assembly,strand,observed_exons_1based_inclusive,target_event,full_chain_assignment,upstream_coding_linkage,three_prime_end_evidence,five_prime_end_evidence,rs1893592_RNA_base,age_label_evidence,assignment_confidence,QC_status,exclusion_reason",
        "decay-measurements-template.csv": "sample_id,donor_code,cell_state,condition,target_label,chase_time_h,measurement_method,measurement_unit,measured_abundance,measurement_standard_error,age_label_method,labeled_abundance_or_probability,recovery_factor,viable_cell_count,division_measurement,limit_of_detection,limit_of_quantification,molecule_assignment_version,QC_status,missing_or_failure_reason",
    }
    for name, header in schemas.items():
        generated.append(write_csv(name, header.split(","), []))
    source_rows = [
        ("TRAILS", "Long-read sequencing for 29 immune cell subsets reveals disease-linked isoforms", "https://www.nature.com/articles/s41467-024-48615-4", "Prior immune-cell catalog context; not a CD4 NMD intervention."),
        ("Penter_2024", "Integrative genotyping of cancer and immune phenotypes by long-read sequencing", "https://www.nature.com/articles/s41467-023-44137-7", "Targeted immune-transcript enrichment precedent; does not validate this UBASH3A assay."),
        ("Boehm_2025", "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome", "https://pubmed.ncbi.nlm.nih.gov/40934927/", "Engineered human cell-line perturbation; auxin conditions do not transfer to unmodified CD4 cells."),
        ("SMG1_2025", "Identification of nonsense-mediated decay inhibitors that alter the tumor immune landscape", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11832170/", "SMG1 inhibition precedent in tumor models; primary-CD4 qualification required."),
        ("Herzog_2017", "Thiol-linked alkylation of RNA to assess expression dynamics", "https://www.nature.com/articles/nmeth.4435", "Metabolic labeling principle; default short or 3-prime reads do not resolve the proposed full RNAs."),
        ("CD8_metabolic_labeling", "Meta-unstable mRNAs in activated CD8+ T cells are defined by interlinked AU-rich elements and m6A mRNA methylation", "https://www.nature.com/articles/s41467-025-67762-w", "Human CD8-cell turnover-method precedent; not CD4 or this UBASH3A endpoint."),
        ("nano_ID_2020", "Native molecule sequencing by nano-ID reveals synthesis and stability of RNA isoforms", "https://pmc.ncbi.nlm.nih.gov/articles/PMC7545145/", "Isoform-resolved metabolic sequencing precedent in K562 cells; not a validated targeted CD4 workflow."),
    ]
    generated.append(write_csv("method-sources.csv", ["source_id", "title", "url", "consulted_date", "application_boundary"],
                               [dict(zip(["source_id", "title", "url", "consulted_date", "application_boundary"], [key, title, url, "2026-09-14", boundary]))
                                for key, title, url, boundary in source_rows]))
    members = {
        "START-HERE.md": (ROOT / "docs/ubash3a-cd4-experiment-brief.md").read_bytes(),
        "ubash3a-cd4-analysis-plan.md": (ROOT / "docs/ubash3a-cd4-analysis-plan.md").read_bytes(),
        "proposed-design.json": PLAN.read_bytes(),
        "conditional-reference-transcripts.fasta": (ROOT / "data/derived/ubash3a-consequences/transcripts.fasta").read_bytes(),
        "reference-junction-chains.csv": (ROOT / "data/derived/ubash3a-hypothesis-test/reference-junction-chains.csv").read_bytes(),
        "diagnostic-RNA-windows.csv": (ROOT / "data/derived/ubash3a-hypothesis-test/diagnostic-RNA-windows.csv").read_bytes(),
        "diagnostic-RNA-windows.fasta": (ROOT / "data/derived/ubash3a-hypothesis-test/diagnostic-RNA-windows.fasta").read_bytes(),
        "protein-sequence-equivalence.csv": (ROOT / "data/derived/ubash3a-hypothesis-test/protein-sequence-equivalence.csv").read_bytes(),
    }
    members.update({path.name: path.read_bytes() for path in generated})
    members["PACKAGE-README.txt"] = (
        "UBASH3A CD4 laboratory feasibility packet, 2026-09-14.\n"
        "Read START-HERE.md, then ubash3a-cd4-analysis-plan.md.\n"
        "DESIGN ONLY. No bench experiment, donor measurements, validated primers, orders or contacts.\n"
        "All three measurement templates are intentionally header-only.\n"
        "Reference FASTA sequences are conditional planning inputs, not experimentally observed full RNAs.\n"
        "proposed-design.json has unset laboratory fields; it is not a frozen confirmatory protocol.\n"
        "Input paths in that file identify the source research repository, not files to load from this ZIP.\n"
        "Published references in method-sources.csv describe methodological precedents, not validation here.\n"
        "Included reference chains have 10 rows; assay-targets.csv adds the exploratory TRAILS model.\n"
        "PACKAGE-MANIFEST.json lists SHA-256 hashes for every other member.\n"
    ).encode()
    member_hashes = {name: digest(data) for name, data in sorted(members.items())}
    members["PACKAGE-MANIFEST.json"] = (json.dumps({"status": "design_only_no_experiment", "files": member_hashes}, indent=2) + "\n").encode()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(members.items()):
            item = zipfile.ZipInfo(name, date_time=(2026, 9, 14, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, data)
    with zipfile.ZipFile(ZIP) as archive:
        assert archive.testzip() is None
        assert set(archive.namelist()) == set(members)
        for name, data in members.items():
            assert archive.read(name) == data
        for name in schemas:
            assert list(csv.DictReader(io.StringIO(archive.read(name).decode()))) == []
    record = {
        "status": "preparation_verified_experiment_not_started", "date": "2026-09-14",
        "design_sha256": digest(PLAN.read_bytes()), "builder_sha256": digest(Path(__file__).read_bytes()),
        "prior_input_files_verified": len(plan["input_files"]), "target_structures": len(targets),
        "empty_measurement_templates": len(schemas), "method_sources": len(source_rows),
        "package_member_count": len(members), "package_sha256": digest(ZIP.read_bytes()),
        "package_bytes": ZIP.stat().st_size, "package_members_verified": True,
        "new_biological_measurements": 0, "laboratory_access_confirmed": False,
        "protocol_frozen": False, "validated_primers": False,
        "generated_file_sha256": {str(path.relative_to(ROOT)): digest(path.read_bytes()) for path in generated},
        "package_contents_sha256": member_hashes,
    }
    RECORD.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({key: record[key] for key in ["status", "prior_input_files_verified", "target_structures", "empty_measurement_templates", "package_member_count", "package_bytes", "new_biological_measurements"]}, indent=2))


if __name__ == "__main__":
    main()
