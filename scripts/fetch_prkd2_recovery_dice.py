"""Audit the entire original public DICE classical-monocyte association file."""

import gzip
import hashlib
import io
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha


def main():
    plan_path = ROOT / 'config/prkd2-recovery-plan.json'
    plan = json.loads(plan_path.read_text())
    sources = json.loads((ROOT / 'config/prkd2-recovery-sources.json').read_text())['sources']
    head = next(s for s in sources if s['id'] == 'dice-original-monocytes-head')
    size = int(head['headers']['Content-Length'])
    if size > 4000000000:
        raise ValueError('Original DICE archive exceeds frozen 4 GB cap')
    record_path = ROOT / 'config/prkd2-recovery-dice-stream.json'
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'source_head': head,
            'reader_sha256': sha(Path(__file__).read_bytes())}
    if record_path.exists():
        record = json.loads(record_path.read_text())
        if any(record[k] != v for k, v in pins.items()):
            raise ValueError('Changed stream inputs')
        for result in record['outputs']:
            if sha((ROOT / result['file']).read_bytes()) != result['sha256']:
                raise ValueError('Changed retained association rows')
        print(json.dumps({'status': 'verified-cache', 'rows': record['total_rows']}))
        return
    archive = ROOT / 'data/raw/prkd2-recovery-original-DICE-MONOCYTES.vcf.gz'
    started = datetime.now(timezone.utc).isoformat()
    if not archive.exists():
        part = archive.with_suffix('.part')
        if part.exists():
            raise ValueError('Partial download preserved; explicitly audit before resume')
        request = urllib.request.Request(head['url'], headers={'User-Agent': 'PSC-public-research/PRKD2-recovery-v1'})
        with urllib.request.urlopen(request, timeout=60) as response, part.open('wb') as out:
            if response.status != 200 or int(response.headers['Content-Length']) != size or response.headers.get('ETag') != head['headers']['ETag']:
                raise ValueError('Original archive identity changed')
            count = 0
            while block := response.read(8 * 1024 * 1024):
                count += len(block)
                if count > size:
                    raise ValueError('Oversized archive')
                out.write(block)
                if count % (128 * 1024 * 1024) == 0:
                    print(json.dumps({'downloaded': count, 'total': size}), flush=True)
        if count != size:
            raise ValueError('Truncated archive')
        part.replace(archive)
    if archive.stat().st_size != size:
        raise ValueError('Cached archive size mismatch')
    digest = hashlib.sha256()
    with archive.open('rb') as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    headers, gene_rows, target_rows = [], [], []
    gene_marker = ('Gene=' + plan['gene_id']).encode()
    target_marker = ('\t' + plan['primary_missing_rsid'] + '\t').encode()
    total_rows, uncompressed_bytes = 0, 0
    uncompressed_hash = hashlib.sha256()
    with gzip.open(archive, 'rb') as stream:
        for line in stream:
            uncompressed_bytes += len(line)
            if uncompressed_bytes > 50000000000:
                raise ValueError('Frozen 50 GB uncompressed cap exceeded')
            uncompressed_hash.update(line)
            if line.startswith(b'#'):
                headers.append(line)
                continue
            total_rows += 1
            if gene_marker in line:
                info = dict(x.split('=', 1) for x in line.decode().rstrip().split('\t')[7].split(';'))
                if info['Gene'].split('.')[0] == plan['gene_id']:
                    gene_rows.append(line)
            if target_marker in line:
                target_rows.append(line)
            if total_rows % 10000000 == 0:
                print(json.dumps({'scanned_rows': total_rows, 'gene_rows': len(gene_rows), 'target_rows': len(target_rows)}), flush=True)
    outputs = []
    directory = ROOT / 'data/derived/prkd2-recovery-queries'
    directory.mkdir(exist_ok=True)
    for name, rows in [('original-DICE-PRKD2', gene_rows), ('original-DICE-rs112445263-all-genes', target_rows)]:
        plain = b''.join(headers + rows)
        data = gzip_bytes(plain)
        path = directory / (name + '.vcf.gz')
        path.write_bytes(data)
        outputs.append({'file': str(path.relative_to(ROOT)), 'rows': len(rows), 'sha256': sha(data), 'uncompressed_sha256': sha(plain)})
    save_json(record_path, {**pins, 'started_at_utc': started, 'completed_at_utc': datetime.now(timezone.utc).isoformat(),
                           'archive': str(archive.relative_to(ROOT)), 'compressed_bytes': size, 'compressed_sha256': digest.hexdigest(),
                           'uncompressed_bytes': uncompressed_bytes, 'uncompressed_sha256': uncompressed_hash.hexdigest(),
                           'total_rows': total_rows, 'gzip_eof_verified': True, 'outputs': outputs,
                           'build': 'GRCh37.p19 per original DICE download page',
                           'limitations': 'Source has rounded Beta/Statistic and no explicit SE, per-variant N, or full component vectors. Presence does not establish suitability for colocalisation.'})
    print(json.dumps({'status': 'complete', 'rows': total_rows, 'gene_rows': len(gene_rows), 'target_rows': len(target_rows)}), flush=True)


if __name__ == '__main__':
    main()
