"""Read the original public BLUEPRINT monocyte gene-QTL archive through EOF."""

import gzip
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, gzip_bytes, save_json, sha


def main():
    plan_path = ROOT / 'config/prkd2-recovery-plan.json'
    plan = json.loads(plan_path.read_text())
    sources = json.loads((ROOT / 'config/prkd2-recovery-sources.json').read_text())['sources']
    head = next(s for s in sources if s['id'] == 'blueprint-original-monocytes-head')
    size = int(head['headers']['Content-Length'])
    if size > 2000000000:
        raise ValueError('BLUEPRINT archive exceeds 2 GB acquisition cap')
    record_path = ROOT / 'config/prkd2-recovery-blueprint-stream.json'
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
    archive = ROOT / 'data/raw/prkd2-recovery-original-BLUEPRINT-monocytes.txt.gz'
    started = datetime.now(timezone.utc).isoformat()
    if not archive.exists():
        part = archive.with_suffix('.part')
        if part.exists():
            raise ValueError('Partial archive preserved; explicitly audit before resume')
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
    gene_rows, target_rows = [], []
    gene_marker = (' ' + plan['gene_id'] + '.').encode()
    target_marker = (' ' + plan['primary_missing_rsid'] + ' ').encode()
    total_rows, uncompressed_bytes = 0, 0
    uncompressed_hash = hashlib.sha256()
    with gzip.open(archive, 'rb') as stream:
        for line in stream:
            uncompressed_bytes += len(line)
            if uncompressed_bytes > 20000000000:
                raise ValueError('20 GB uncompressed cap exceeded')
            uncompressed_hash.update(line)
            total_rows += 1
            if gene_marker in line:
                fields = line.decode().split()
                if len(fields) != 9 or fields[2].split('.')[0] != plan['gene_id']:
                    raise ValueError('Unexpected original BLUEPRINT gene row')
                gene_rows.append(line)
            if target_marker in line:
                target_rows.append(line)
            if total_rows % 10000000 == 0:
                print(json.dumps({'scanned_rows': total_rows, 'gene_rows': len(gene_rows), 'target_rows': len(target_rows)}), flush=True)
    outputs = []
    directory = ROOT / 'data/derived/prkd2-recovery-queries'
    directory.mkdir(exist_ok=True)
    for name, rows in [('original-BLUEPRINT-PRKD2', gene_rows), ('original-BLUEPRINT-rs112445263-all-genes', target_rows)]:
        plain = b''.join(rows)
        data = gzip_bytes(plain)
        path = directory / (name + '.txt.gz')
        path.write_bytes(data)
        outputs.append({'file': str(path.relative_to(ROOT)), 'rows': len(rows), 'sha256': sha(data), 'uncompressed_sha256': sha(plain)})
    save_json(record_path, {**pins, 'started_at_utc': started, 'completed_at_utc': datetime.now(timezone.utc).isoformat(),
                           'archive': str(archive.relative_to(ROOT)), 'compressed_bytes': size, 'compressed_sha256': digest.hexdigest(),
                           'uncompressed_bytes': uncompressed_bytes, 'uncompressed_sha256': uncompressed_hash.hexdigest(),
                           'total_rows': total_rows, 'gzip_eof_verified': True, 'outputs': outputs,
                           'columns': ['chr:pos_ref_alt', 'rsid', 'phenotypeID', 'p.value', 'beta', 'Bonferroni.p.value', 'FDR', 'alt_allele_frequency', 'std.error_of_beta'],
                           'column_source_id': 'blueprint-qtl-readme',
                           'limitations': 'Original cohort and analysis; not an independent replication or interchangeable Catalogue component input. Genome build, exact edits, per-gene N and effect convention require reconciliation.'})
    print(json.dumps({'status': 'complete', 'rows': total_rows, 'gene_rows': len(gene_rows), 'target_rows': len(target_rows)}), flush=True)


if __name__ == '__main__':
    main()
