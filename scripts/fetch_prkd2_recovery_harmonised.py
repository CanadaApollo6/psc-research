"""Audit the complete public harmonised PSC copy for missing allele identities."""

import concurrent.futures
import csv
import gzip
import hashlib
import json
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha


def main():
    sources = json.loads((ROOT / 'config/prkd2-recovery-extra-sources.json').read_text())['sources']
    head = next(s for s in sources if s['id'] == 'psc-harmonised-head')
    plan_path = ROOT / 'config/prkd2-recovery-harmonised-plan.json'
    plan = json.loads(plan_path.read_text())
    size = int(head['headers']['Content-Length'])
    if size != plan['source_size'] or size > 500000000:
        raise ValueError('PSC harmonised archive size exceeds plan')
    path = ROOT / 'config/prkd2-recovery-harmonised-query.json'
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'reader_sha256': sha(Path(__file__).read_bytes()), 'source_head': head}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'ranges': [], 'started_at_utc': datetime.now(timezone.utc).isoformat()}
    if any(record[k] != v for k, v in pins.items()):
        raise ValueError('Changed harmonised archive plan')
    if record.get('complete'):
        if sha((ROOT / record['output']['file']).read_bytes()) != record['output']['sha256']:
            raise ValueError('Changed harmonised regional rows')
        print(json.dumps({'status': 'verified-cache'}))
        return
    directory = ROOT / 'data/raw/prkd2-recovery-harmonised-ranges'
    directory.mkdir(exist_ok=True)
    archive = ROOT / 'data/raw/prkd2-recovery-GCST004030-harmonised.tsv.gz'
    if not archive.exists():
        lock = threading.Lock()
        last = [0.0]
        starts = list(range(0, size, 8388608))

        def download(start):
            previous = next((r for r in record['ranges'] if r['start'] == start), None)
            if previous:
                if sha((ROOT / previous['file']).read_bytes()) != previous['sha256']:
                    raise ValueError('Changed cached range')
                return previous
            end = min(size - 1, start + 8388607)
            with lock:
                time.sleep(max(0, 2 - (time.monotonic() - last[0])))
                last[0] = time.monotonic()
            req = urllib.request.Request(head['url'], headers={'Range': f'bytes={start}-{end}', 'User-Agent': 'PSC-public-research/PRKD2-recovery-v1'})
            with urllib.request.urlopen(req, timeout=60) as response:
                if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{size}' or response.headers.get('ETag') != head['headers']['ETag']:
                    raise ValueError('Changed remote range identity')
                data = response.read(end - start + 2)
            if len(data) != end - start + 1:
                raise ValueError('Incomplete range')
            file = directory / f'{start}-{end}.data'
            file.write_bytes(data)
            return {'start': start, 'end': end, 'file': str(file.relative_to(ROOT)), 'bytes': len(data),
                    'sha256': sha(data), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            for result in pool.map(download, starts):
                if not any(r['start'] == result['start'] for r in record['ranges']):
                    record['ranges'].append(result)
                    save_json(path, record)
                if len(record['ranges']) % 8 == 0:
                    print(json.dumps({'verified_bytes': sum(r['bytes'] for r in record['ranges']), 'total': size}), flush=True)
        record['ranges'].sort(key=lambda r: r['start'])
        if [r['start'] for r in record['ranges']] != starts:
            raise ValueError('Incomplete byte coverage')
        with archive.with_suffix('.assembling').open('wb') as stream:
            for result in record['ranges']:
                data = (ROOT / result['file']).read_bytes()
                if sha(data) != result['sha256']:
                    raise ValueError('Changed range during assembly')
                stream.write(data)
        archive.with_suffix('.assembling').replace(archive)
    digest, md5 = hashlib.sha256(), hashlib.md5()
    with archive.open('rb') as stream:
        while block := stream.read(8388608):
            digest.update(block)
            md5.update(block)
    if archive.stat().st_size != size or md5.hexdigest() != plan['publisher_md5']:
        raise ValueError('Publisher checksum mismatch')
    lines, count = [], 0
    with gzip.open(archive, 'rt', newline='') as stream:
        header = stream.readline().rstrip('\r\n')
        names = header.split('\t')
        for line in stream:
            count += 1
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != len(names):
                raise ValueError('Malformed harmonised row')
            row = dict(zip(names, fields))
            in_union = any(row[chrom] == '19' and row[pos].isdigit() and plan['start_union_0based'] < int(row[pos]) <= plan['end_union_exclusive']
                           for chrom, pos in [('chromosome', 'base_pair_location'), ('hm_chrom', 'hm_pos')])
            is_control = row['variant_id'] in plan['control_rsids'] or row['hm_rsid'] in plan['control_rsids']
            if in_union or is_control:
                lines.append(line.rstrip('\r\n'))
            if count > 20000000 or len(lines) > 50000:
                raise ValueError('Frozen row cap exceeded')
    plain = (header + '\n' + '\n'.join(lines) + '\n').encode()
    output = ROOT / 'data/derived/prkd2-recovery-queries/PSC-harmonised-union-region.tsv.gz'
    data = gzip_bytes(plain)
    output.write_bytes(data)
    record.update({'complete': True, 'completed_at_utc': datetime.now(timezone.utc).isoformat(), 'source_rows': count,
                   'gzip_eof_verified': True, 'archive_sha256': digest.hexdigest(), 'publisher_md5_verified': md5.hexdigest(),
                   'output': {'file': str(output.relative_to(ROOT)), 'sha256': sha(data), 'uncompressed_sha256': sha(plain), 'rows': len(lines)}})
    save_json(path, record)
    print(json.dumps({'status': 'complete', 'total_source_rows': count, 'retained_union_rows': len(lines)}), flush=True)


if __name__ == '__main__':
    main()
