"""Shared coordinate, provenance and deterministic-output helpers for coloc."""

import bisect
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + '.part')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def gzip_bytes(plain):
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode='wb', filename='', mtime=0) as stream:
        stream.write(plain)
    return output.getvalue()


class ChainMap:
    """Map single bases through UCSC chain blocks; never guess gap mappings."""

    def __init__(self, text):
        self.blocks = {}
        header = None
        for line in text.splitlines():
            if not line:
                header = None
                continue
            fields = line.split()
            if fields[0] == 'chain':
                if len(fields) != 13 or fields[4] != '+':
                    raise ValueError('Unsupported chain header')
                header = fields
                source_pos, dest_pos = int(fields[5]), int(fields[10])
                continue
            if header is None or len(fields) not in (1, 3):
                raise ValueError('Malformed chain block')
            size = int(fields[0])
            if size <= 0:
                raise ValueError('Empty chain block')
            block = (source_pos, source_pos + size, header[7], dest_pos, int(header[8]), header[9], header[12])
            self.blocks.setdefault(header[2], []).append(block)
            source_pos += size
            dest_pos += size
            if len(fields) == 3:
                source_pos += int(fields[1])
                dest_pos += int(fields[2])
            elif source_pos != int(header[6]) or dest_pos != int(header[11]):
                raise ValueError('Chain block lengths disagree with header')
        self.starts = {}
        self.max_ends = {}
        for chrom, blocks in self.blocks.items():
            blocks.sort()
            self.starts[chrom] = [b[0] for b in blocks]
            maximum, ends = 0, []
            for block in blocks:
                maximum = max(maximum, block[1]); ends.append(maximum)
            self.max_ends[chrom] = ends

    @classmethod
    def from_file(cls, path):
        return cls(gzip.decompress(Path(path).read_bytes()).decode())

    def map_base(self, chrom, position0):
        blocks = self.blocks.get(chrom, [])
        index = bisect.bisect_right(self.starts.get(chrom, []), position0) - 1
        hits = set()
        while index >= 0 and self.max_ends[chrom][index] > position0:
            start, end, dest, dest_start, dest_size, strand, _ = blocks[index]
            if start <= position0 < end:
                coordinate = dest_start + position0 - start
                if strand == '-':
                    coordinate = dest_size - coordinate - 1
                hits.add((dest, coordinate, strand))
            index -= 1
        return sorted(hits)

    def enclosing_destination(self, chrom, start0, end0):
        endpoints = []
        for start, end, dest, dest_start, dest_size, strand, _ in self.blocks.get(chrom, []):
            lo, hi = max(start, start0), min(end, end0)
            if lo >= hi or dest != chrom:
                continue
            for base in (lo, hi - 1):
                value = dest_start + base - start
                endpoints.append(value if strand == '+' else dest_size - value - 1)
        if not endpoints:
            raise ValueError('No same-chromosome chain coverage for region')
        return min(endpoints), max(endpoints) + 1


def reciprocal_base(forward, reverse, chrom, position0):
    hits = forward.map_base(chrom, position0)
    if len(hits) != 1:
        return None
    dest, coordinate, strand = hits[0]
    if dest != chrom or reverse.map_base(dest, coordinate) != [(chrom, position0, strand)]:
        return None
    return hits[0]


def complement(base):
    return base.translate(str.maketrans('ACGT', 'TGCA'))


def verify_plan():
    path = ROOT / 'config/coloc-plan.json'
    plan = json.loads(path.read_text())
    for filename, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / filename).read_bytes()) != digest:
            raise ValueError('Frozen coloc input changed: ' + filename)
    return plan
