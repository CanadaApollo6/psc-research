"""Independently verify metadata joins, source totals, and raw-read classification."""

import collections
import csv
import hashlib
import itertools
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio.Seq import Seq
import pysam

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/ubash3a-raw-read-audit"
OUT = ROOT / "data/derived/ubash3a-raw-read-audit"
WORK = ROOT / "work/ubash3a-raw-read-audit"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def independently_walk_cigar(read):
    """Use textual CIGAR and one-based positions, independently of primary parser."""
    tokens = [(int(n),op) for n,op in re.findall(r"(\d+)([MIDNSHP=X])", read.cigarstring or "")]
    position = read.reference_start+1
    junctions, anchors, blocks = [], [], []
    for i,(length,op) in enumerate(tokens):
        if op == "N":
            junctions.append((position-1,position+length))
            sides = []
            for step in [-1,1]:
                cursor, count = i+step,0
                while 0<=cursor<len(tokens) and tokens[cursor][1] in "M=X":
                    count += tokens[cursor][0]
                    cursor += step
                sides.append(count)
            anchors.append(sides)
        if op in "M=X":
            blocks.append((position,position+length-1))
        if op in "MDN=X":
            position += length
    return junctions,anchors,blocks


def affine_score(reference, query):
    """Independent three-state dynamic programming, using the fixed score scheme."""
    negative = float("-inf")
    size = len(query)+1
    previous_match = [negative]*size
    previous_insert = [negative]*size
    previous_delete = [negative]*size
    previous_match[0] = 0
    for j in range(1,size):
        previous_delete[j] = -5-(j-1)
    for i,base in enumerate(reference,1):
        current_match = [negative]*size
        current_insert = [negative]*size
        current_delete = [negative]*size
        current_insert[0] = -5-(i-1)
        for j,other in enumerate(query,1):
            current_match[j] = max(previous_match[j-1],previous_insert[j-1],previous_delete[j-1])+(2 if base==other else -3)
            current_insert[j] = max(previous_match[j]-5,previous_insert[j]-1,previous_delete[j]-5)
            current_delete[j] = max(current_match[j-1]-5,current_insert[j-1]-5,current_delete[j-1]-1)
        previous_match,previous_insert,previous_delete = current_match,current_insert,current_delete
    return max(previous_match[-1],previous_insert[-1],previous_delete[-1])


def metadata_checks():
    sources = json.loads((ROOT / "config/ubash3a-raw-read-metadata-sources.json").read_text())["sources"]
    for source in sources:
        assert sha(ROOT/source["path"]) == source["sha256"]
    experiments = {r.attrib["accession"]:r for r in ET.parse(RAW/"DRA016285.experiment.xml").getroot().findall("EXPERIMENT")}
    samples = {r.attrib["accession"]:r for r in ET.parse(RAW/"DRA016285.sample.xml").getroot().findall("SAMPLE")}
    runs = {r.find("EXPERIMENT_REF").attrib["accession"]:r.attrib["accession"]
            for r in ET.parse(RAW/"DRA016285.run.xml").getroot().findall("RUN")}
    table = list(csv.DictReader((OUT/"sample-run-inventory.csv").open()))
    assert len(table)==len(experiments)==len(samples)==len(runs)==29
    assert {r["experiment"] for r in table}==set(experiments)
    for row in table:
        exp = experiments[row["experiment"]]
        sid = exp.find("DESIGN/SAMPLE_DESCRIPTOR").attrib["accession"]
        sample_name = next(a.findtext("VALUE") for a in samples[sid].findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE") if a.findtext("TAG")=="sample_name")
        assert row["sample"]==sid and row["run"]==runs[row["experiment"]]
        assert row["source_sample_name"]==sample_name
        assert row["source_library_name"]==exp.findtext("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME")
        assert row["platform"]==exp.findtext(".//INSTRUMENT_MODEL")
        assert row["layout"]==list(exp.find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT"))[0].tag
    selected = json.loads((ROOT/"config/ubash3a-raw-read-sample-selection.json").read_text())["initial_samples"]
    for row in selected:
        assert row in table
    assert {r["source_library_name"] for r in selected}=={"NCD4","MemCD4"}
    return {"all_experiment_sample_run_joins":len(table),"selected_libraries":len(selected),"metadata_source_hashes":len(sources)}


def run():
    plan_path = ROOT/"config/ubash3a-raw-read-plan.json"
    plan = json.loads(plan_path.read_text())
    audit = json.loads((OUT/"raw-read-audit.json").read_text())
    assert audit["plan_sha256"]==sha(plan_path)
    descriptor_path=ROOT/"config/ubash3a-raw-read-coding-descriptors.json"
    assert audit["additive_coding_descriptors_sha256"]==sha(descriptor_path)
    descriptors=json.loads(descriptor_path.read_text())
    for name,pin in plan["input_pins"].items():
        assert sha(ROOT/name)==pin
    for name,pin in audit["files_sha256"].items():
        assert sha(OUT/name)==pin
    metadata = metadata_checks()
    totals = json.loads((ROOT/"config/ubash3a-raw-read-single-run-total-sources.json").read_text())["sources"]
    expected_totals = {}
    for source in totals:
        assert sha(ROOT/source["local_path"])==source["sha256"]
        records = json.loads((ROOT/source["local_path"]).read_text())
        assert len(records)==1
        expected_totals[source["run"]]=records[0]
    selected = json.loads((ROOT/"config/ubash3a-raw-read-sample-selection.json").read_text())["initial_samples"]
    accession = {r["source_library_name"]:r["run"] for r in selected}
    table = list(csv.DictReader((OUT/"raw-locus-alignments.csv").open()))
    refs = list(csv.DictReader((ROOT/"data/derived/ubash3a-cd4-experiment/assay-targets.csv").open()))
    lookup = {(r["sample"],r["query_name"],int(r["alignment_index_in_query"])):r for r in table}
    assert len(lookup)==len(table)
    verified,dp_verified,denominators = 0,0,[]
    with pysam.FastaFile(str(WORK/"GRCh38.primary.fa")) as genome:
        for sample in plan["sample_ids"]:
            source = audit["samples"][sample]
            stats = source["stats"]
            total = expected_totals[accession[sample]]
            assert stats["input_reads"]==int(total["read_count"])
            assert stats["input_bases"]==int(total["base_count"])
            assert stats["output_query_groups"]==stats["input_reads"]
            denominators.append({"sample":sample,"reads":stats["input_reads"],"bases":stats["input_bases"],"source_totals_match":True})
            normal_dp_count = 0
            with pysam.AlignmentFile(str(WORK/sample/"target-groups.bam"),"rb") as bam:
                for query,records in itertools.groupby(bam,key=lambda r:r.query_name):
                    group = list(records)
                    primary = next(r for r in group if not r.flag & (256|2048))
                    sequence = str(Seq(primary.query_sequence).reverse_complement()) if primary.flag & 16 else primary.query_sequence
                    sequence_hash = hashlib.sha256(sequence.upper().encode()).hexdigest()
                    has_supplementary = any(r.flag & 2048 or r.has_tag("SA") for r in group)
                    for index,read in enumerate(group):
                        if read.is_unmapped or read.reference_name!="21" or read.reference_start>=42447684 or read.reference_end<42403447:
                            continue
                        row = lookup[(sample,query,index)]
                        junctions,anchors,blocks = independently_walk_cigar(read)
                        assert [list(j) for j in junctions]==json.loads(row["junctions_json"])
                        assert anchors==json.loads(row["anchors_json"])
                        event = "other"
                        target = None
                        for name,j in plan["junctions_1based_exonic_endpoints"].items():
                            if tuple(j) in junctions:
                                event,target=name,tuple(j)
                        assert event==row["event"]
                        strand = read.get_tag("ts") if read.has_tag("ts") else "unknown"
                        if read.flag & 16 and strand!="unknown":
                            strand = "+" if strand=="-" else "-"
                        assert strand==row["inferred_genomic_transcript_strand"]
                        assert sequence_hash==row["sequence_sha256"]
                        coverage_fields={"shared_start_coordinates_aligned":descriptors["T29_assumed_start_1based"],
                                         "T29_stop_coordinates_aligned":descriptors["T29_conditional_stop_1based"],
                                         "previous_plus29_stop_coordinates_aligned":descriptors["previous_plus29_conditional_stop_1based"],
                                         "T29_catalog_end_coordinate_aligned":[descriptors["T29_catalog_end_1based"]]*2}
                        for field,(lo,hi) in coverage_fields.items():
                            expected_covered=all(any(a<=base<=b for a,b in blocks) for base in range(lo,hi+1))
                            assert expected_covered==(row[field]=="True")
                        for anchor in [20,10]:
                            eligible = (target is not None and not read.flag & (4|256|512|1024|2048)
                                        and read.mapping_quality>=20 and not has_supplementary and strand!="-"
                                        and min(anchors[junctions.index(target)])>=anchor)
                            assert eligible==(row[f"qualified_anchor{anchor}"]=="True")
                            full_matches,downstream_matches=[],[]
                            if eligible:
                                for ref in refs:
                                    if ref["event"]!=event:
                                        continue
                                    ref_chain=[tuple(j) for j in json.loads(ref["junction_chain_exon_end_next_start"])]
                                    last_start,last_end=json.loads(ref["exons_1based_inclusive"])[-1]
                                    last_bases=sum(max(0,min(b,last_end)-max(a,last_start)+1) for a,b in blocks)
                                    tail_index=junctions.index(target)+1
                                    expected_tail=ref_chain[ref_chain.index(target)+1:]
                                    if junctions[tail_index:]==expected_tail and all(min(a)>=anchor for a in anchors[tail_index:]) and last_bases>=20:
                                        downstream_matches.append(ref["target_label"])
                                    start_covered=all(any(a<=base<=b for a,b in blocks) for base in [42403946,42403947,42403948])
                                    if junctions==ref_chain and all(min(a)>=anchor for a in anchors) and start_covered and last_bases>=20:
                                        full_matches.append(ref["target_label"])
                            full_key="full_splice_chain_and_start_matches" if anchor==20 else "full_chain_matches_anchor10"
                            downstream_key="complete_downstream_chain_matches" if anchor==20 else "downstream_matches_anchor10"
                            assert ";".join(full_matches)==row[full_key]
                            assert ";".join(downstream_matches)==row[downstream_key]
                        verify_dp = (event=="plus29" or event=="normal" and normal_dp_count<10)
                        if verify_dp and row["dp_flank_start"] and not read.flag & (256|2048):
                            lo,hi=int(row["dp_flank_start"]),int(row["dp_flank_end"])
                            pairs={r+1:q for q,r in read.get_aligned_pairs(matches_only=True) if q is not None and r is not None}
                            observed=read.query_sequence[pairs[lo]:pairs[hi]+1].upper()
                            for kind,donor in [("normal",42434954),("plus29",42434983)]:
                                template=genome.fetch("21",lo-1,donor).upper()+genome.fetch("21",42437487,hi).upper()
                                assert affine_score(template,observed)==float(row[f"dp_{kind}_score"])
                            dp_verified+=1
                            normal_dp_count+=int(event=="normal")
                        verified+=1
            with pysam.AlignmentFile(str(WORK/sample/"target-groups.sorted.bam"),"rb") as bam:
                count=sum(1 for _ in bam.fetch("21",42403446,42447684))
                assert count==sum(r["sample"]==sample for r in table)
    assert verified==len(table)
    summaries=list(csv.DictReader((OUT/"raw-read-summary.csv").open()))
    for summary in summaries:
        local=[r for r in table if r["sample"]==summary["sample"]]
        assert len(local)==int(summary["locus_alignment_rows"])
        for event in ["normal","plus29"]:
            assert sum(r["event"]==event for r in local)==int(summary[f"{event}_exact_all_alignments"])
            for anchor in [20,10]:
                qualifying=[r for r in local if r["event"]==event and r[f"qualified_anchor{anchor}"]=="True"]
                assert len(qualifying)==int(summary[f"{event}_qualified_anchor{anchor}"])
                assert len({r["sequence_sha256"] for r in qualifying})==int(summary[f"{event}_distinct_exact_sequences_anchor{anchor}"])
    result={"plan_sha256":sha(plan_path),"audit_sha256":sha(OUT/"raw-read-audit.json"),
            "verification_script_sha256":sha(Path(__file__)),"metadata_checks":metadata,
            "independent_source_totals":denominators,"independent_alignment_rows":verified,
            "independent_affine_DP_reads":dp_verified,
            "DP_sampling_rule":"All primary exact +29 candidates with shared flanks, plus first ten normal candidates per library in source order",
            "all_checks_passed":True}
    (ROOT/"reports/ubash3a-raw-read-verification.json").write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    run()
