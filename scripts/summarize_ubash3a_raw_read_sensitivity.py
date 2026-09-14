"""Summarize all nominated-read alignments under both frozen mapping settings."""

import collections
import csv
import hashlib
import json
from pathlib import Path

import pysam
from analyze_ubash3a_upf1_reads import structure

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/"work/ubash3a-raw-read-audit"
OUT=ROOT/"data/derived/ubash3a-raw-read-audit"


def main():
    rows=[]
    events=[]
    pins={}
    for sample in ["NCD4","MemCD4"]:
        metadata=json.loads((WORK/sample/"alignment.json").read_text())
        for setting,name in [("primary","target-groups.bam"),("splice_flank_no","target-sensitivity.bam")]:
            path=WORK/sample/name
            h=hashlib.sha256(path.read_bytes()).hexdigest()
            assert h==metadata["output_hashes"][name]
            pins[str(path.relative_to(ROOT))]=h
            counts=collections.Counter()
            all_events=collections.Counter()
            with pysam.AlignmentFile(str(path),"rb") as bam:
                for read in bam:
                    if read.is_unmapped or read.reference_name!="21":
                        continue
                    chain,_,_=structure(read.reference_start,read.cigartuples or [])
                    kind="secondary_or_supplementary" if read.flag & (256|2048) else "primary"
                    for donor,acceptor in chain:
                        if (donor,acceptor) in [(42434954,42437488),(42434983,42437488),(42437580,42437855)]:
                            label="normal" if donor==42434954 else "plus29" if donor==42434983 else "T_specific_downstream"
                            all_events[(label,kind)]+=1
                        if kind=="primary" and (acceptor==42437488 or donor==42437580):
                            counts[(donor,acceptor)]+=1
            for (donor,acceptor),number in sorted(counts.items()):
                rows.append({"sample":sample,"setting":setting,"donor_exon_end_1based":donor,
                             "acceptor_exon_start_1based":acceptor,"primary_alignment_count":number,
                             "quality_filter":"none_beyond_primary_mapped_alignment"})
            for label in ["normal","plus29","T_specific_downstream"]:
                for kind in ["primary","secondary_or_supplementary"]:
                    events.append({"sample":sample,"setting":setting,"event":label,"alignment_class":kind,
                                   "alignment_count":all_events[(label,kind)]})
    for name,data in [("raw-junction-neighborhood.csv",rows),("raw-all-setting-event-counts.csv",events)]:
        with (OUT/name).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=data[0].keys(),lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    audit={"input_alignment_sha256":pins,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "output_sha256":{name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in ["raw-junction-neighborhood.csv","raw-all-setting-event-counts.csv"]},
           "scope":"All output alignments for nominated target queries; the sensitivity run cannot find queries absent from primary-setting nomination. Counts are alignment-coordinate descriptors, not donor or molecule counts."}
    (OUT/"raw-sensitivity-summary-audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(events,indent=2))


if __name__=="__main__":
    main()
