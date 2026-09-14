"""Run the frozen two-library, whole-primary-assembly TRAILS alignment audit."""

import argparse
import collections
import gzip
import hashlib
import itertools
import json
import shutil
import subprocess
import threading
from datetime import datetime, timezone
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "config/ubash3a-raw-read-plan.json"
WORK = ROOT / "work/ubash3a-raw-read-audit"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def check_plan():
    plan = json.loads(PLAN.read_text())
    assert plan["frozen_before_raw_target_inspection"]
    for name, expected in plan["input_pins"].items():
        assert sha(ROOT / name) == expected, name
    version = subprocess.check_output([str(ROOT / plan["aligner"]), "--version"], text=True).strip()
    assert version == plan["aligner_version"]
    return plan


def prepare_reference(plan):
    source = next(r for r in plan["downloads"] if r["id"] == "GRCh38_primary")
    archive = ROOT / source["local_path"]
    receipt = json.loads(archive.with_name(archive.name + ".receipt.json").read_text())
    assert sha(archive) == receipt["sha256"]
    WORK.mkdir(parents=True, exist_ok=True)
    fasta = WORK / "GRCh38.primary.fa"
    index = WORK / "GRCh38.primary.splice.mmi"
    record = WORK / "reference.json"
    if record.exists():
        old = json.loads(record.read_text())
        assert old["source_sha256"] == receipt["sha256"]
        assert sha(fasta) == old["fasta_sha256"]
        assert sha(index) == old["minimap_index_sha256"]
        return fasta, index
    temporary = fasta.with_suffix(".fa.part")
    with gzip.open(archive, "rb") as compressed, temporary.open("wb") as output:
        shutil.copyfileobj(compressed, output, length=8 * 1024 * 1024)
    temporary.replace(fasta)
    pysam.faidx(str(fasta))
    with pysam.FastaFile(str(fasta)) as reference:
        seq = json.loads((ROOT / "data/raw/ubash3a-consequences/ensembl-gene-sequence.json").read_text())["seq"].upper()
        assert reference.fetch("21", 42403446, 42447684).upper() == seq
        sequence_lengths = dict(zip(reference.references, reference.lengths, strict=True))
        assert sum(sequence_lengths.values()) < 4_000_000_000, "Index must be one part"
        assert {str(i) for i in range(1, 23)} | {"X", "Y", "MT"} <= set(sequence_lengths)
    command = [str(ROOT / plan["aligner"]), *plan["index_options"], "-d", str(index), str(fasta)]
    with (WORK / "index.log").open("w") as log:
        subprocess.run(command, stderr=log, check=True)
    dump(record, {"source_url": source["url"], "source_sha256": receipt["sha256"],
                  "fasta_sha256": sha(fasta), "minimap_index_sha256": sha(index),
                  "fasta_index_sha256": sha(Path(str(fasta) + ".fai")),
                  "sequence_lengths": sequence_lengths, "total_bases": sum(sequence_lengths.values()),
                  "one_part_index_required": True, "gene_reference_exact_match": True,
                  "index_command": command, "decompression_integrity": "gzip_CRC_and_length_passed",
                  "completed_at_utc": datetime.now(timezone.utc).isoformat()})
    print(json.dumps({"stage": "reference_ready", "sequences": len(sequence_lengths),
                      "bases": sum(sequence_lengths.values())}), flush=True)
    return fasta, index


def fastq_records(stream):
    """Validate the archive's four-line FASTQ representation without changing bases."""
    while True:
        header = stream.readline()
        if not header:
            break
        sequence, separator, quality = (stream.readline() for _ in range(3))
        if not header.startswith(b"@") or not separator.startswith(b"+") or not quality:
            raise ValueError("Invalid or wrapped FASTQ record; do not silently truncate")
        seq, qual = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
        if not seq or len(seq) != len(qual):
            raise ValueError("FASTQ sequence/quality length mismatch")
        if set(seq.upper()) - set(b"ACGTRYWSKMBDHVN") or min(qual) < 33 or max(qual) > 126:
            raise ValueError("Invalid FASTQ alphabet or quality")
        name = header[1:].split()[0]
        if not name:
            raise ValueError("Empty FASTQ read name")
        yield header + sequence + separator + quality, len(seq)


def feed_fastq(source, target, stats, errors, log_path):
    try:
        with log_path.open("w") as log:
            process = subprocess.Popen(["bzip2", "-dc", str(source)], stdout=subprocess.PIPE, stderr=log)
            try:
                for record, length in fastq_records(process.stdout):
                    target.write(record)
                    stats["input_reads"] += 1
                    stats["input_bases"] += length
                    stats["shortest_read"] = min(stats.get("shortest_read", length), length)
                    stats["longest_read"] = max(stats.get("longest_read", length), length)
                    if stats["input_reads"] % 250000 == 0:
                        print(json.dumps({"stage": "mapping_input", "sample": stats["sample"],
                                          "reads": stats["input_reads"], "bases": stats["input_bases"]}), flush=True)
                if process.wait() != 0:
                    raise RuntimeError("bzip2 CRC/decompression failed")
                stats["complete_fastq_and_bzip2_integrity"] = True
            finally:
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                    process.wait()
    except BaseException as exc:
        errors.append(repr(exc))
    finally:
        target.close()


def locus_overlap(read, plan):
    left, right = plan["locus_interval_1based_inclusive"]
    return (not read.is_unmapped and read.reference_name == plan["chromosome"]
            and read.reference_start < right and read.reference_end >= left)


def map_sample(plan, sample, index):
    source = next(r for r in plan["downloads"] if r["id"] == sample)
    fastq = ROOT / source["local_path"]
    receipt = json.loads(fastq.with_name(fastq.name + ".receipt.json").read_text())
    assert sha(fastq) == receipt["sha256"]
    folder = WORK / sample
    folder.mkdir(parents=True, exist_ok=True)
    completed = folder / "alignment.json"
    if completed.exists():
        old = json.loads(completed.read_text())
        assert old["plan_sha256"] == sha(PLAN)
        assert old["input_sha256"] == receipt["sha256"]
        for name, expected in old["output_hashes"].items():
            assert sha(folder / name) == expected
        return old
    stats = collections.Counter(sample=sample)
    errors = []
    target_bam, target_fastq = folder / "target-groups.bam", folder / "target-reads.fastq.gz"
    command = [str(ROOT / plan["aligner"]), *plan["primary_alignment_options"], str(index), "-"]
    with (folder / "minimap-primary.log").open("w") as log:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log)
        feeder = threading.Thread(target=feed_fastq,
                                  args=(fastq, process.stdin, stats, errors, folder / "decompression.log"))
        feeder.start()
        try:
            with pysam.AlignmentFile(process.stdout, "r") as sam, pysam.AlignmentFile(str(target_bam), "wb", header=sam.header) as out, gzip.open(target_fastq, "wt") as selected:
                for query, records in itertools.groupby(sam, key=lambda r: r.query_name):
                    group = list(records)
                    stats["output_query_groups"] += 1
                    primaries = [r for r in group if not r.is_secondary and not r.is_supplementary]
                    if len(primaries) != 1:
                        raise ValueError(f"Expected one primary/unmapped record per unique read: {query}")
                    primary = primaries[0]
                    stats["mapped_primary_reads"] += int(not primary.is_unmapped)
                    stats["unmapped_reads"] += int(primary.is_unmapped)
                    stats["all_output_alignments"] += len(group)
                    if not any(locus_overlap(r, plan) for r in group):
                        continue
                    stats["target_query_groups"] += 1
                    stats["target_all_locus_and_competing_alignments"] += len(group)
                    for read in group:
                        out.write(read)
                    # The primary SAM record retains the entire query, including soft clips.
                    if not primary.query_sequence or any(op == 5 for op, _ in primary.cigartuples or []):
                        raise ValueError("Primary record does not retain the full source read")
                    sequence, quality = primary.get_forward_sequence(), primary.get_forward_qualities()
                    if quality is None or len(quality) != len(sequence):
                        raise ValueError("Target source-read qualities unavailable")
                    selected.write(f"@{query}\n{sequence}\n+\n{''.join(chr(q+33) for q in quality)}\n")
        except BaseException:
            process.terminate()
            raise
        finally:
            process.stdout.close()
            feeder.join()
        returncode = process.wait()
    if errors or returncode:
        raise RuntimeError(f"Incomplete alignment: {errors}; minimap2 exit={returncode}")
    assert stats["input_reads"] == stats["output_query_groups"]
    assert stats["complete_fastq_and_bzip2_integrity"]
    sorted_bam = folder / "target-groups.sorted.bam"
    pysam.sort("-m", "256M", "-@", "2", "-o", str(sorted_bam), str(target_bam))
    pysam.index(str(sorted_bam))
    sensitivity = folder / "target-sensitivity.bam"
    second = [str(ROOT / plan["aligner"]), *plan["primary_alignment_options"],
              *plan["sensitivity_extra_options"], str(index), str(target_fastq)]
    with (folder / "minimap-sensitivity.log").open("w") as log:
        process = subprocess.Popen(second, stdout=subprocess.PIPE, stderr=log)
        with pysam.AlignmentFile(process.stdout, "r") as sam, pysam.AlignmentFile(str(sensitivity), "wb", header=sam.header) as out:
            for read in sam:
                out.write(read)
        process.stdout.close()
        if process.wait():
            raise RuntimeError("Target-read sensitivity alignment failed")
    result = {"sample": sample, "stats": dict(stats), "plan_sha256": sha(PLAN),
              "script_sha256": sha(Path(__file__)), "input_sha256": receipt["sha256"],
              "command": command, "sensitivity_command": second,
              "sensitivity_scope": "only reads nominated by a primary-setting locus or competing alignment",
              "completed_at_utc": datetime.now(timezone.utc).isoformat(),
              "output_hashes": {p.name: sha(p) for p in [target_bam, target_fastq, sorted_bam, Path(str(sorted_bam)+".bai"), sensitivity]}}
    dump(completed, result)
    print(json.dumps({"stage": "sample_complete", **dict(stats)}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--sample", choices=["NCD4", "MemCD4"])
    args = parser.parse_args()
    plan = check_plan()
    _, index = prepare_reference(plan)
    if not args.reference_only:
        for sample in [args.sample] if args.sample else plan["sample_ids"]:
            map_sample(plan, sample, index)


if __name__ == "__main__":
    main()
