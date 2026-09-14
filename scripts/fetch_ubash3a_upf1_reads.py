"""Capture only the fixed UBASH3A interval from eight public indexed BAMs."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pysam
import requests

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-upf1-read-test-plan.json"
RAW = ROOT / "data/raw/ubash3a-hypothesis-test"
MANIFEST = ROOT / "config/ubash3a-upf1-read-sources.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def metadata(url):
    response = requests.head(url, allow_redirects=True, timeout=30)
    response.raise_for_status()
    return {"url": response.url, "status": response.status_code,
            "bytes": int(response.headers["Content-Length"]),
            "etag": response.headers.get("ETag"),
            "last_modified": response.headers.get("Last-Modified"),
            "accept_ranges": response.headers.get("Accept-Ranges")}


def fetch(sample, index_source, interval, old=None, offline=False, cache=RAW, expected_metadata=None):
    stem = sample["sample_id"]
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / (stem + ".UBASH3A.bam")
    header = cache / (stem + ".header.sam")
    prefix = cache / (stem + ".prefix.bgzf")
    if old:
        missing = []
        for name, pin in old["local_files"].items():
            path = ROOT / name
            if not path.exists():
                missing.append(path)
                continue
            if path.stat().st_size != pin["bytes"] or sha(path) != pin["sha256"]:
                raise ValueError("Changed locus cache: " + name)
        if missing:
            if offline:
                raise FileNotFoundError(missing[0])
            staging = ROOT / "work/ubash3a-hypothesis-test/rehydrate" / stem
            fresh = fetch(sample, index_source, interval, cache=staging,
                          expected_metadata=old["source_metadata_before"])
            by_name = {Path(name).name: (ROOT / name, pin) for name, pin in fresh["local_files"].items()}
            for name, pin in old["local_files"].items():
                if by_name[Path(name).name][1] != pin:
                    raise ValueError("Re-extracted subset differs from its historical pin")
            (staging / "restoration-receipt.json").write_text(json.dumps(fresh, indent=2) + "\n")
            for path in missing:
                temporary = path.with_suffix(path.suffix + ".rehydrate.part")
                if temporary.exists():
                    raise ValueError("An incomplete restoration needs inspection")
                shutil.copyfile(by_name[path.name][0], temporary)
                temporary.replace(path)
        return old
    if offline:
        raise FileNotFoundError(target)
    if any(p.exists() for p in (target, header, prefix, target.with_suffix(".bam.part"))):
        raise ValueError("Unpinned partial extraction requires inspection: " + stem)
    index = ROOT / index_source["local_cache_path"]
    if sha(index) != index_source["sha256"]:
        raise ValueError("Pinned index changed")
    start_time = now()
    before = metadata(sample["bam_url"])
    if expected_metadata and any(before[key] != expected_metadata[key] for key in ("bytes", "etag", "last_modified")):
        raise ValueError("Remote BAM identity differs from the historical extraction source")
    with requests.get(sample["bam_url"], headers={"Range": "bytes=0-65535"},
                      stream=True, timeout=(15, 30)) as response:
        if response.status_code != 206 or response.headers.get("Content-Range") != f"bytes 0-65535/{before['bytes']}":
            raise ValueError("Server did not honor the bounded byte-range probe")
        blocks, size = [], 0
        for block in response.iter_content(16384):
            size += len(block)
            if size > 65536:
                raise ValueError("Prefix probe exceeded its byte range")
            blocks.append(block)
    body = b"".join(blocks)
    if len(body) != 65536 or body[:2] != b"\x1f\x8b":
        raise ValueError("Invalid compressed BAM prefix")
    prefix.write_bytes(body)
    temporary = target.with_suffix(".bam.part")
    count = 0
    with pysam.AlignmentFile(sample["bam_url"], "rb", index_filename=str(index),
                             require_index=True) as source:
        # Chromosome 21 distinguishes the claimed GRCh38 assembly from GRCh37.
        chrom = "chr21" if "chr21" in source.references else "21"
        if source.get_reference_length(chrom) != 46709983:
            raise ValueError("BAM chromosome-21 length does not match GRCh38")
        header.write_text(source.text)
        with pysam.AlignmentFile(str(temporary), "wb", template=source) as out:
            for read in source.fetch(chrom, interval[0] - 1, interval[1]):
                out.write(read)
                count += 1
                if count > 100000 or temporary.stat().st_size > 100_000_000:
                    raise ValueError("Locus extraction exceeded its safety bound")
    after = metadata(sample["bam_url"])
    if any(before[key] != after[key] for key in ("bytes", "etag", "last_modified")):
        raise ValueError("Source identity changed during extraction")
    temporary.replace(target)
    pysam.index(str(target))
    with pysam.AlignmentFile(str(target), "rb") as saved:
        if sum(1 for _ in saved.fetch(until_eof=True)) != count:
            raise ValueError("Extracted BAM readback count mismatch")
    paths = [target, Path(str(target) + ".bai"), header, prefix]
    return {**sample, "retrieval_started_at_utc": start_time, "retrieved_at_utc": now(),
            "source_metadata_before": before, "source_metadata_after": after,
            "source_index_sha256": sha(index), "source_index_path": str(index.relative_to(ROOT)),
            "region_1based_inclusive": [chrom, *interval], "all_locus_alignments": count,
            "local_bam": str(target.relative_to(ROOT)),
            "local_files": {str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": sha(p)} for p in paths},
            "representation": "pysam/htslib indexed interval extraction; all overlapping alignments retained with the source header. Hashes identify the subset, prefix, header and indexes, not the entire remote BAM."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    plan = json.loads(PLAN.read_text())
    if sha(ROOT / "config/ubash3a-hypothesis-test-plan.json") != plan["parent_plan_sha256"]:
        raise ValueError("Parent plan changed")
    if sha(RAW / "E-MTAB-14725.sdrf.txt") != plan["sdrf_sha256"]:
        raise ValueError("Sample metadata changed")
    parent = json.loads((ROOT / "config/ubash3a-hypothesis-test-plan.json").read_text())
    gene_path = ROOT / "data/raw/ubash3a-consequences/ensembl-gene.json"
    if sha(gene_path) != parent["reference_pins"][str(gene_path.relative_to(ROOT))]:
        raise ValueError("Locus coordinates changed")
    gene = json.loads(gene_path.read_text())
    sources = {s["source_id"]: s for s in json.loads((ROOT / "config/ubash3a-hypothesis-test-sources.json").read_text())["sources"]}
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {
        "plan_sha256": sha(PLAN), "created_at_utc": now(), "fetch_script_sha256_at_start": sha(Path(__file__)),
        "pysam_version": pysam.__version__, "htslib_version": pysam.__samtools_version__,
        "samples": [], "failed_attempts": []}
    if manifest["plan_sha256"] != sha(PLAN):
        raise ValueError("Read-test plan changed")
    existing = {s["sample_id"]: s for s in manifest["samples"]}
    failures = 0
    for sample in plan["samples"]:
        try:
            row = fetch(sample, sources[sample["sample_id"] + "-bai"],
                        [gene["start"], gene["end"]], existing.get(sample["sample_id"]), args.offline)
            existing[sample["sample_id"]] = row
            print(json.dumps({"sample": sample["sample_id"], "status": "verified", "locus_alignments": row["all_locus_alignments"]}), flush=True)
        except Exception as exc:
            failures += 1
            error = {"sample_id": sample["sample_id"], "at_utc": now(), "error": str(exc)}
            manifest["failed_attempts"].append(error)
            print(json.dumps(error), flush=True)
        manifest["samples"] = sorted(existing.values(), key=lambda s: s["sample_id"])
        if not args.offline:
            MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
