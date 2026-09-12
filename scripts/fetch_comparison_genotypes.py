"""Fetch bounded indexed public VCF subsets needed for the fixed matching rule."""

import argparse
import csv
import gzip
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    if args.offline and args.refresh:
        parser.error('Cannot combine offline and refresh')
    plan_path = ROOT / 'data/derived/psc-comparator-query-plan.json'
    target_path = ROOT / 'data/derived/psc-comparator-genotype-targets.csv'
    plans = json.loads(plan_path.read_text())
    targets = list(csv.DictReader(io.StringIO(target_path.read_text())))
    source_manifest = json.loads((ROOT / 'config/psc-comparator-sources.json').read_text())
    run_path = ROOT / 'config/psc-comparator-genotype-queries.json'
    run = json.loads(run_path.read_text()) if run_path.exists() else {
        'plan_sha256': hashlib.sha256(plan_path.read_bytes()).hexdigest(),
        'targets_sha256': hashlib.sha256(target_path.read_bytes()).hexdigest(),
        'pysam_version': pysam.__version__, 'queries': []}
    if run['plan_sha256'] != hashlib.sha256(plan_path.read_bytes()).hexdigest() or run['targets_sha256'] != hashlib.sha256(target_path.read_bytes()).hexdigest():
        raise ValueError('Genotype query inputs changed')
    if pysam.__version__ != run['pysam_version']:
        raise ValueError('Use the recorded pysam version')
    for plan in plans:
        label = plan['archive_region']
        previous = next((q for q in run['queries'] if q['archive_region'] == label), None)
        path = ROOT / f'data/raw/match-genotypes-{label}.vcf.gz'
        if previous and path.exists():
            if hashlib.sha256(path.read_bytes()).hexdigest() != previous['sha256']:
                raise ValueError('Genotype subset checksum changed')
            if not args.refresh:
                print(json.dumps({'region': label, 'status': 'verified-cache'}), flush=True)
                continue
        if args.offline:
            raise FileNotFoundError(path)
        index = next(s for s in source_manifest['sources'] if s['id'] == plan['index_source_id'])
        index_path = ROOT / 'data/raw' / index['file']
        if hashlib.sha256(index_path.read_bytes()).hexdigest() != index['sha256']:
            raise ValueError('Genotype index changed')
        keys = {(r['chromosome'], int(r['position_grch37_1based']), r['reference_bases'], r['alternate_bases']) for r in targets if r['archive_region'] == label}
        selected, seen, scanned, uncompressed = [], set(), 0, 0
        with pysam.TabixFile(plan['url'], index=str(index_path)) as table:
            headers = list(table.header)
            for line in table.fetch(plan['chromosome'], plan['start_grch37_0based'], plan['end_grch37_0based_exclusive']):
                scanned += 1
                uncompressed += len(line)
                if scanned > 75000 or uncompressed > 768 * 1024**2:
                    raise ValueError('Regional VCF query exceeded its bound')
                fields = line.split('\t', 5)
                key = fields[0], int(fields[1]), fields[3], fields[4]
                if key in keys:
                    if key in seen:
                        raise ValueError('Duplicate exact coordinate/allele VCF record')
                    seen.add(key)
                    selected.append(line)
        data = ('\n'.join(headers + selected) + '\n').encode()
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, filename='', mode='wb', mtime=0) as stream:
            stream.write(data)
        compressed = buffer.getvalue()
        record = {**plan, 'retrieved_at_utc': datetime.now(timezone.utc).isoformat(), 'file': str(path.relative_to(ROOT)),
                  'sha256': hashlib.sha256(compressed).hexdigest(), 'uncompressed_sha256': hashlib.sha256(data).hexdigest(),
                  'index_sha256': index['sha256'], 'regional_rows_scanned': scanned, 'regional_bytes_scanned': uncompressed,
                  'target_coordinate_allele_keys': len(keys), 'exact_target_rows_retained': len(selected), 'bytes': len(compressed)}
        if previous and record['uncompressed_sha256'] != previous['uncompressed_sha256']:
            raise ValueError('Remote genotype subset changed; existing result retained')
        if not previous:
            if path.exists():
                raise ValueError('Unpinned genotype file already exists')
            path.write_bytes(compressed)
            run['queries'].append(record)
            run_path.write_text(json.dumps(run, indent=2) + '\n')
        print(json.dumps({'region': label, 'scanned_rows': scanned, 'retained_rows': len(selected), 'target_keys': len(keys)}), flush=True)


if __name__ == '__main__':
    main()
