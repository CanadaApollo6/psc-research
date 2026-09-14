"""Rebuild the complete PXD000376 file inventory from pinned metadata pages."""

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data/raw/ubash3a-proteomics-readiness"
OUT=ROOT/"data/derived/ubash3a-proteomics-readiness"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output=OUT):
    manifest=json.loads((ROOT/"config/ubash3a-proteomics-readiness-sources.json").read_text())
    for source in manifest["sources"]:
        if source.get("local_path"):
            assert sha(ROOT/source["local_path"])==source["sha256"]
    for field in ["source_API_documentation","study_design_source"]:
        source=manifest.get(field,{})
        if source.get("local_path"):
            assert sha(ROOT/source["local_path"])==source["sha256"]
    records=[r for page in range(3) for r in json.loads((RAW/f"PXD000376-files-page{page}.json").read_text())]
    expected=json.loads((RAW/"PXD000376-file-count.json").read_text())
    assert len(records)==expected==273
    assert len({r["accession"] for r in records})==len(records)
    assert all("PXD000376" in r["projectAccessions"] for r in records)
    rows=[]
    for r in records:
        assert isinstance(r["fileSizeBytes"],int) and r["fileSizeBytes"]>=0
        rows.append({"project":"PXD000376","file_accession":r["accession"],"file_name":r["fileName"],
                     "category":r["fileCategory"]["value"],"size_bytes":r["fileSizeBytes"],
                     "publisher_checksum_as_returned":r.get("checksum",""),
                     "public_locations_json":json.dumps(r["publicFileLocations"],separators=(",",":")),
                     "sample_assignment_status":"not_yet_confirmed_from_metadata_or_methods",
                     "retrieval_status":"metadata_only_not_downloaded"})
    output.mkdir(parents=True,exist_ok=True)
    with (output/"PXD000376-file-inventory.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys(),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts=collections.Counter(r["category"] for r in rows)
    bytes_by_category={c:sum(r["size_bytes"] for r in rows if r["category"]==c) for c in counts}
    extensions=collections.Counter(Path(r["file_name"]).suffix.lower() for r in rows)
    audit={"project":"PXD000376","complete_file_count_endpoint":expected,"unique_file_records":len(rows),
           "metadata_pages":[100,100,73],"files_by_category":dict(counts),"bytes_by_category":bytes_by_category,
           "total_listed_bytes":sum(r["size_bytes"] for r in rows),"files_by_extension":dict(extensions),
           "files_with_nonempty_publisher_checksum":sum(bool(r["publisher_checksum_as_returned"]) for r in rows),
           "inventory_sha256":sha(output/"PXD000376-file-inventory.csv"),"spectrum_files_downloaded":0,
           "mass_spectra_searched":0,"human_proteome_specificity_tested":False,
           "source_metadata_warning":"Project protocol says CD45RO depletion of resting memory CD4 cells; retain verbatim and reconcile with paper/sample mapping before assigning files or comparisons."}
    (output/"PXD000376-readiness-audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps({"files":len(rows),"bytes":audit["total_listed_bytes"],"spectra_searched":0}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=OUT)
    run(parser.parse_args().output_dir)
