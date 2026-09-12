"""Verify the pinned public ETS2 release and inventory its benchmark inputs.

Does not execute source code or calculate expression/program/drug scores.
"""

import csv
import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFIXES = (
    "RNA-seq/RNAseq_CRISPR/", "RNA-seq/RNAseq_MEKi/",
    "RNA-seq/RNAseq_overexpression/", "RNA-seq/RNAseq_biopsies/",
    "Spatial_transcriptomics/",
)


def digest(path, algorithm="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def main():
    manifest = json.loads((ROOT / "config/ets2-benchmark-sources.json").read_text())
    directory = ROOT / manifest["source_directory"]
    for source in manifest["sources"]:
        path = directory / source["file"]
        if path.stat().st_size != source["bytes"] or digest(path) != source["sha256"]:
            raise ValueError(f"Pinned source mismatch: {source['file']}")
    archive = directory / "ets2-public-release.zip"
    if digest(archive, "md5") != manifest["published_zip_md5"]:
        raise ValueError("Publisher archive MD5 mismatch")
    rows = []
    with zipfile.ZipFile(archive) as stream:
        base = stream.namelist()[0]
        names = set(stream.namelist())
        for item in stream.infolist():
            path = item.filename[len(base):]
            if item.is_dir() or not path.startswith(PREFIXES) or path.endswith(".Rhistory"):
                continue
            content = stream.read(item)
            rows.append({
                "member": path, "bytes": item.file_size,
                "sha256": hashlib.sha256(content).hexdigest(),
                "role": "author analysis code" if path.endswith((".r", ".R")) else
                        "public released data; endpoint semantics need review before analysis",
                "raw_perturbation_count_matrix_present_in_same_folder": str(
                    base + path.rsplit("/", 1)[0] + "/raw_counts.txt" in names
                ).lower(),
            })
    target = ROOT / "data/derived/ets2-benchmark-inventory.csv"
    with target.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Verified three sources and inventoried {len(rows)} public release members.")


if __name__ == "__main__":
    main()
