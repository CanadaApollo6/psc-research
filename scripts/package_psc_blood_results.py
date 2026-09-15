#!/usr/bin/env python3
"""Losslessly package large aggregate RNA CSVs; frozen original bytes stay local."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    folder = ROOT / "data/derived/psc-blood-rnaseq"
    manifest = json.loads((folder / "rnaseq-output-manifest.json").read_text())
    packaged = []
    for name in ["rnaseq-feature-universe.csv", "rnaseq-primary.csv", "rnaseq-leave-one-plate-out.csv"]:
        source = folder / name
        original_hash = digest(source)
        if original_hash != manifest["outputs"][name]["sha256"]:
            raise ValueError("Cannot package a changed frozen output")
        dest = folder / (name + ".gz")
        with source.open("rb") as inp, dest.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
                shutil.copyfileobj(inp, out)
        h = hashlib.sha256()
        with gzip.open(dest, "rb") as inp:
            for block in iter(lambda: inp.read(1024 * 1024), b""):
                h.update(block)
        if h.hexdigest() != original_hash:
            raise ValueError("Compression round trip changed the aggregate results")
        packaged.append({"file":dest.name,"bytes":dest.stat().st_size,"sha256":digest(dest),"uncompressed_file":name,"uncompressed_bytes":source.stat().st_size,"uncompressed_sha256":original_hash})
    receipt = {"transformation":"gzip only; exact uncompressed SHA-256 verified; no scientific values or frozen files changed", "files":packaged}
    (folder / "compressed-artifacts.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__ == "__main__":
    main()
