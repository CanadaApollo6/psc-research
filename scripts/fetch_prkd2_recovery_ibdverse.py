"""Read four prespecified chromosome members from the public IBDverse ZIP.

Uses exact byte ranges, full raw-DEFLATE EOF and central-directory CRC checks.
Preserves the complete PRKD2 source vector; does not compute associations.
"""

import argparse
import csv
import hashlib
import json
import re
import struct
import zlib
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha
from fetch_coverage_repair_parquet import RangeReader


class AuthorRangeReader(RangeReader):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.directory = ROOT / 'data/raw/prkd2-recovery-zip-ranges'
        self.directory.mkdir(exist_ok=True)


def variant_position(identifier):
    fields = re.split('[_:]', identifier)
    if len(fields) >= 2 and fields[0].removeprefix('chr') == '19' and fields[1].isdigit():
        return int(fields[1])
    return None


def scan_member(blocks, gene, positions, expected_crc, expected_size):
    decoder, crc, count, size = zlib.decompressobj(-15), 0, 0, 0
    pending, header, names = b'', None, []
    gene_lines, point_lines = [], []
    digest = hashlib.sha256()

    def consume(raw):
        nonlocal header, names, count
        if header is None:
            header = raw.rstrip(b'\r')
            names = header.decode().split('\t')
            if 'phenotype_id' not in names or 'variant_id' not in names:
                raise ValueError('Unexpected nominal TSV identifiers')
            return
        count += 1
        fields = raw.rstrip(b'\r').decode().split('\t')
        if len(fields) != len(names):
            raise ValueError('Malformed nominal TSV row')
        if fields[names.index('phenotype_id')].split('.')[0] in (gene, 'PRKD2'):
            gene_lines.append(raw.rstrip(b'\r'))
        if variant_position(fields[names.index('variant_id')]) in positions:
            point_lines.append(raw.rstrip(b'\r'))
        if len(gene_lines) > 100000 or len(point_lines) > 100000 or count > 20000000:
            raise ValueError('Bounded member row cap exceeded')

    for block in blocks:
        output = decoder.decompress(block)
        crc = zlib.crc32(output, crc)
        size += len(output)
        digest.update(output)
        pending += output
        pieces = pending.split(b'\n')
        pending = pieces.pop()
        for raw in pieces:
            consume(raw)
    tail = decoder.flush()
    if tail:
        raise ValueError('Unexpected buffered deflate output')
    if pending:
        consume(pending)
    if not decoder.eof or decoder.unused_data or crc != expected_crc or size != expected_size:
        raise ValueError('ZIP member EOF, size or CRC mismatch')
    return {'source_rows': count, 'source_bytes': size, 'source_sha256': digest.hexdigest(),
            'crc32': crc, 'deflate_eof_verified': True, 'header': names}, header, gene_lines, point_lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    plan_path = ROOT / 'config/prkd2-recovery-author-plan.json'
    plan = json.loads(plan_path.read_text())
    for filename, digest in plan['frozen_files'].items():
        if sha((ROOT / filename).read_bytes()) != digest:
            raise ValueError('Changed author source selection input')
    metadata = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-zip-metadata.json').read_text())
    source = metadata['queries'][0]
    path = ROOT / 'config/prkd2-recovery-ibdverse-queries.json'
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'reader_sha256': sha(Path(__file__).read_bytes())}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'queries': []}
    if any(record[k] != v for k, v in pins.items()):
        raise ValueError('Changed author query reader or plan')
    targets = list(csv.DictReader(open(ROOT / 'data/derived/prkd2-recovery-targets.csv')))
    positions = {int(row[key]) for row in targets for key in ('position_grch37', 'position')}
    for selected in plan['ibdverse']:
        query = next((q for q in record['queries'] if q['id'] == selected['id']), None)
        if query is None:
            query = {**selected, **{k: source[k] for k in ('url', 'remote_size', 'etag', 'last_modified')},
                     'ranges': [], 'started_at_utc': datetime.now(timezone.utc).isoformat()}
            record['queries'].append(query)
        reader = AuthorRangeReader(query, record, path, args.offline)
        member = selected['member']
        offset, length = int(member['local_header_offset']), int(member['compressed_bytes'])
        if length > plan['transfer_cap_per_context_bytes']:
            raise ValueError('Prespecified member exceeds per-context cap')
        header = reader.fetch(offset, 30)
        fields = struct.unpack('<4s5H3L2H', header)
        if fields[0] != b'PK\x03\x04' or fields[3] != 8 or fields[2] & 1:
            raise ValueError('Unexpected or encrypted ZIP member')
        name_len, extra_len = fields[-2:]
        name = reader.fetch(offset + 30, name_len + extra_len)[:name_len].decode()
        if name != member['filename']:
            raise ValueError('Local ZIP member identity differs from central directory')
        start = offset + 30 + name_len + extra_len
        parts = []
        for relative in range(0, length, 8388608):
            size = min(8388608, length - relative)
            reader.fetch(start + relative, size)
            parts.append((start + relative, size))
            if not args.offline:
                print(json.dumps({'id': query['id'], 'verified_compressed_bytes': relative + size, 'member_bytes': length}), flush=True)
        result, header, genes, points = scan_member(
            (reader.fetch(a, b) for a, b in parts), plan['gene_id'], positions,
            int(member['crc32']), int(member['uncompressed_bytes']))
        result['outputs'] = []
        for label, lines in [('PRKD2', genes), ('target-positions', points)]:
            plain = header + b'\n' + b'\n'.join(lines) + (b'\n' if lines else b'')
            data = gzip_bytes(plain)
            output = ROOT / f"data/derived/prkd2-recovery-queries/{query['id']}-{label}.tsv.gz"
            if output.exists() and output.read_bytes() != data:
                raise ValueError('Changed author source extraction')
            output.write_bytes(data)
            result['outputs'].append({'file': str(output.relative_to(ROOT)), 'rows': len(lines), 'sha256': sha(data), 'uncompressed_sha256': sha(plain)})
        result['transfer_bytes'] = sum(r['bytes'] for r in query['ranges'])
        if 'result' in query and query['result'] != result:
            raise ValueError('Author source replay differs')
        query['result'] = result
        save_json(path, record)
        print(json.dumps({'id': query['id'], 'status': 'complete', **result}), flush=True)


if __name__ == '__main__':
    main()
