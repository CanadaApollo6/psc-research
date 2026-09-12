"""Bounded HTTP-range reader for the three frozen r8 OneK1K Parquet files.

Every transferred byte range is cached and hashed. Position statistics and
identifier columns select row groups. Metadata mode does not read effects.
"""

import argparse
import io
import json
import re
import time
import urllib.request
from datetime import datetime, timezone

import pyarrow
import pyarrow.parquet as pq

from coloc_common import ROOT, gzip_bytes, save_json, sha


class RangeReader(io.RawIOBase):
    def __init__(self, spec, manifest, manifest_path, offline):
        self.spec, self.manifest, self.manifest_path = spec, manifest, manifest_path
        self.offline, self.position = offline, 0
        self.directory = ROOT / 'data/raw/repair-parquet-ranges'
        self.directory.mkdir(exist_ok=True)
        if 'remote_size' not in spec:
            self.fetch(0, 8)

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.spec['remote_size'] + offset
        if self.position < 0: raise ValueError('Negative range')
        return self.position

    def fetch(self, start, length):
        if length <= 0: return b''
        identifier = f'{start}-{start + length - 1}'
        previous = next((r for r in self.spec['ranges'] if r['range'] == identifier), None)
        if previous:
            data = (ROOT / previous['file']).read_bytes()
            if len(data) != length or sha(data) != previous['sha256']: raise ValueError('Changed cached byte range')
            return data
        if self.offline: raise ValueError('Uncached Parquet byte range: ' + identifier)
        if sum(r['bytes'] for r in self.spec['ranges']) + length > 300000000:
            raise ValueError('Frozen 300 MB per-dataset transfer cap exceeded')
        request = urllib.request.Request(self.spec['url'], headers={'Range': 'bytes=' + identifier, 'User-Agent': 'PSC-public-research/0.6'})
        with urllib.request.urlopen(request, timeout=60) as response:
            match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
            if response.status != 206 or not match or tuple(map(int, match.groups()[:2])) != (start, start + length - 1):
                raise ValueError('Server did not honor exact byte range')
            size = int(match[3])
            identity = {'remote_size': size, 'etag': response.headers.get('ETag'), 'last_modified': response.headers.get('Last-Modified')}
            if any(k in self.spec and self.spec[k] != v for k, v in identity.items()):
                raise ValueError('Remote Parquet identity changed')
            data = response.read(length + 1)
            if len(data) != length: raise ValueError('Truncated or oversized byte-range response')
        self.spec.update(identity)
        path = self.directory / (self.spec['id'] + '-' + identifier + '.data')
        path.write_bytes(data)
        self.spec['ranges'].append({'range': identifier, 'file': str(path.relative_to(ROOT)), 'bytes': len(data),
                                    'sha256': sha(data), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()})
        save_json(self.manifest_path, self.manifest)
        time.sleep(2)
        return data

    def read(self, length=-1):
        if length < 0: length = self.spec['remote_size'] - self.position
        length = min(length, self.spec['remote_size'] - self.position)
        data = self.fetch(self.position, length)
        self.position += len(data)
        return data

    def readinto(self, buffer):
        data = self.read(len(buffer)); buffer[:len(data)] = data
        return len(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata-only', action='store_true')
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    plan_path = ROOT / 'config/coverage-repair-alternative-plan.json'
    plan = json.loads(plan_path.read_text())
    access_path = ROOT / 'config/coverage-repair-alternative-parquet-query-plan.json'
    access = json.loads(access_path.read_text())
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest: raise ValueError('Frozen alternative input changed')
    path = ROOT / 'config/coverage-repair-parquet-queries.json'
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'access_plan_sha256': sha(access_path.read_bytes()),
            'reader_sha256': sha((ROOT / 'scripts/fetch_coverage_repair_parquet.py').read_bytes()), 'pyarrow_version': pyarrow.__version__}
    initial_path = ROOT / 'config/coverage-repair-parquet-metadata-attempt.json'
    if sha(initial_path.read_bytes()) != access['metadata_attempt_sha256']: raise ValueError('Changed original metadata attempt')
    initial = json.loads(initial_path.read_text())
    if [s['url'] for s in initial['queries']] != [s['url'] for s in access['sources']]: raise ValueError('Changed source URL')
    manifest = json.loads(path.read_text()) if path.exists() else {**pins, 'queries': initial['queries']}
    if any(manifest.get(k) != v for k, v in pins.items()): raise ValueError('Frozen Parquet reader or plan changed')
    for spec in manifest['queries']:
        reader = RangeReader(spec, manifest, path, args.offline)
        parquet = pq.ParquetFile(reader)
        schema = parquet.schema_arrow.names
        gene_column = 'gene_id' if 'gene_id' in schema else 'molecular_trait_id'
        if gene_column not in schema: raise ValueError('No gene or trait identifier')
        position_index = schema.index('position')
        candidate = []
        for i in range(parquet.metadata.num_row_groups):
            statistics = parquet.metadata.row_group(i).column(position_index).statistics
            if not statistics or not statistics.has_min_max: raise ValueError('Missing row-group position bounds; do not scan whole archive')
            lo, hi = statistics.min, statistics.max
            if lo <= plan['end_grch38_exclusive'] and hi > plan['start_grch38_0based']:
                candidate.append(i)
        selected = []
        for i in candidate:
            identifiers = parquet.read_row_group(i, columns=['chromosome', gene_column], use_threads=False).to_pandas()
            keep = identifiers.chromosome.eq(plan['chromosome'].removeprefix('chr')) & identifiers[gene_column].str.split('.').str[0].eq(plan['gene_id'])
            if keep.any(): selected.append(i)
        spec['metadata'] = {'schema': str(parquet.schema_arrow), 'total_rows': parquet.metadata.num_rows,
                            'total_row_groups': parquet.metadata.num_row_groups, 'candidate_gene_row_groups': selected,
                            'identifier_only_row_groups': candidate,
                            'identifier_only_rows': sum(parquet.metadata.row_group(i).num_rows for i in candidate),
                            'candidate_rows': sum(parquet.metadata.row_group(i).num_rows for i in selected)}
        save_json(path, manifest)
        if args.metadata_only:
            print(json.dumps({'id': spec['id'], **spec['metadata']}), flush=True)
            continue
        if spec['metadata']['candidate_rows'] > 1000000: raise ValueError('Frozen nominal row cap exceeded')
        frame = parquet.read_row_groups(selected, use_threads=False).to_pandas()
        frame = frame[frame[gene_column].str.split('.').str[0].eq(plan['gene_id']) & frame.chromosome.eq(plan['chromosome'].removeprefix('chr'))]
        # Preserve full source gene rows, then explicitly record regional filtering.
        gene_rows = len(frame)
        frame = frame[frame.position.gt(plan['start_grch38_0based']) & frame.position.le(plan['end_grch38_exclusive'])]
        file = ROOT / f"data/derived/coverage-repair-queries/{spec['id']}-alternative.tsv.gz"
        plain = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
        data = gzip_bytes(plain)
        if file.exists() and file.read_bytes() != data: raise ValueError('Changed extracted alternative rows')
        file.write_bytes(data)
        result = {'file': str(file.relative_to(ROOT)), 'sha256': sha(data), 'uncompressed_sha256': sha(plain),
                  'gene_rows_before_region_filter': gene_rows, 'rows': len(frame), 'transfer_bytes': sum(r['bytes'] for r in spec['ranges'])}
        if 'result' in spec and spec['result'] != result: raise ValueError('Parquet replay differs')
        spec['result'] = result
        save_json(path, manifest)
        print(json.dumps({'id': spec['id'], **result}), flush=True)


if __name__ == '__main__': main()
