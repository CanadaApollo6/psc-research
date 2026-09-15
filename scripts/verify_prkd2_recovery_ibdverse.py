"""Independent standard-ZIP extraction using only verified cached byte ranges."""

import bisect
import csv
import gzip
import hashlib
import io
import json
import zipfile
from functools import lru_cache
from pathlib import Path

from coloc_common import ROOT, save_json, sha


class CachedArchive(io.RawIOBase):
    def __init__(self, records, size):
        self.records = sorted(records, key=lambda r: int(r['range'].split('-')[0]))
        self.starts = [int(r['range'].split('-')[0]) for r in self.records]
        self.size, self.position = size, 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        if self.position < 0:
            raise ValueError('Negative cached archive seek')
        return self.position

    @lru_cache(maxsize=4)
    def block(self, i):
        r = self.records[i]
        data = (ROOT / r['file']).read_bytes()
        if len(data) != r['bytes'] or sha(data) != r['sha256']:
            raise ValueError('Changed cached ZIP range')
        return data

    def read(self, length=-1):
        if length < 0:
            length = self.size - self.position
        end = min(self.size, self.position + length)
        out = []
        while self.position < end:
            i = bisect.bisect_right(self.starts, self.position) - 1
            if i < 0:
                raise ValueError('Uncached ZIP range')
            start = self.starts[i]
            data = self.block(i)
            stop = min(end, start + len(data))
            if stop <= self.position:
                raise ValueError('Gap in cached ZIP range')
            out.append(data[self.position - start:stop - start])
            self.position = stop
        return b''.join(out)

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)


def main():
    q = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-queries.json').read_text())
    metadata = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-zip-metadata.json').read_text())
    source = metadata['queries'][0]
    ranges = source['ranges'] + [r for s in q['queries'] for r in s['ranges']]
    reader = CachedArchive(ranges, source['remote_size'])
    with open(ROOT / 'data/derived/prkd2-recovery-targets.csv') as stream:
        targets = list(csv.DictReader(stream))
    positions = {x[k].encode() for x in targets for k in ('position', 'position_grch37')}
    outcomes = []
    with zipfile.ZipFile(reader) as archive:
        for spec in q['queries']:
            gene, points, total = [], [], 0
            digest = hashlib.sha256()
            info = archive.getinfo(spec['member']['filename'])
            with archive.open(info) as stream:
                header = stream.readline()
                digest.update(header)
                names = header.rstrip(b'\r\n').split(b'\t')
                gi, vi = names.index(b'phenotype_id'), names.index(b'variant_id')
                for line in stream:
                    digest.update(line)
                    total += 1
                    fields = line.rstrip(b'\r\n').split(b'\t')
                    if len(fields) != len(names):
                        raise ValueError('Malformed independent source row')
                    if fields[gi].split(b'.')[0] == b'ENSG00000105287':
                        gene.append(line.rstrip(b'\r\n'))
                    variant = fields[vi].split(b':')
                    if len(variant) == 4 and variant[0] == b'chr19' and variant[1] in positions:
                        points.append(line.rstrip(b'\r\n'))
            expected = spec['result']
            if digest.hexdigest() != expected['source_sha256'] or total != expected['source_rows']:
                raise ValueError('Independent entire-member scan differs')
            for rows, result in zip([gene, points], expected['outputs']):
                plain = header.rstrip(b'\r\n') + b'\n' + b'\n'.join(rows) + (b'\n' if rows else b'')
                if plain != gzip.decompress((ROOT / result['file']).read_bytes()) or len(rows) != result['rows']:
                    raise ValueError('Independent retained ZIP rows differ')
            outcome = {'id': spec['id'], 'source_rows': total, 'PRKD2_rows': len(gene),
                       'target_position_rows': len(points), 'full_source_sha256': digest.hexdigest(),
                       'standard_zip_CRC_and_EOF_verified': True, 'exact_extraction_match': True}
            outcomes.append(outcome)
            print(json.dumps(outcome), flush=True)
    save_json(ROOT / 'reports/prkd2-recovery-ibdverse-verification.json', {
        'all_checks_pass': True, 'verifier_sha256': sha(Path(__file__).read_bytes()),
        'method': 'Standard-library ZIP reader and independent byte-field parser; primary reader uses explicit raw-DEFLATE blocks.',
        'contexts': outcomes})


if __name__ == '__main__':
    main()
