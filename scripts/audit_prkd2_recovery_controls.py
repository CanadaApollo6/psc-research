"""Verify the true rs313839 control and the fixed public-reference site checks."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pysam

from coloc_common import ROOT, gzip_bytes, save_json, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['controls', 'reference'])
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    plan_path = ROOT / 'config/prkd2-recovery-followup-plan.json'
    plan = json.loads(plan_path.read_text())
    parent = json.loads((ROOT / 'config/prkd2-recovery-plan.json').read_text())
    sources = json.loads((ROOT / 'config/prkd2-recovery-sources.json').read_text())['sources']
    path = ROOT / ('config/prkd2-recovery-' + args.stage + '-audit.json')
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'reader_sha256': sha(Path(__file__).read_bytes())}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'queries': []}
    if any(record[k] != v for k, v in pins.items()):
        raise ValueError('Changed control audit plan')
    if args.stage == 'controls':
        target = next(t for t in plan['targets'] if t['rsid'] == 'rs313839')
        registry = {s['id']: s for s in sources}
        queries = []
        for ds in parent['cohorts']:
            if ds['release'] != 'r7' or ds['prior_analysis']:
                continue
            dsid = ds['dataset_id']
            index = registry[dsid + '-index']
            queries.append({'id': dsid, 'url': index['url'].removesuffix('.tbi'), 'index': index,
                            'targets': [target], 'chromosome': '19'})
    else:
        extra = json.loads((ROOT / 'config/prkd2-recovery-extra-sources.json').read_text())['sources']
        index = next(s for s in extra if s['id'] == 'reference-30x-20201028-chr19-index')
        queries = [{'id': '1000G-30x-20201028', 'url': index['url'].removesuffix('.tbi'), 'index': index,
                    'targets': plan['targets'], 'chromosome': 'chr19'}]
    directory = ROOT / 'data/derived/prkd2-recovery-queries'
    directory.mkdir(exist_ok=True)
    for query in queries:
        previous = next((q for q in record['queries'] if q['id'] == query['id']), None)
        if previous:
            if sha((ROOT / previous['file']).read_bytes()) != previous['sha256']:
                raise ValueError('Changed site audit result')
            continue
        if args.offline:
            raise ValueError('Missing site audit cache')
        index = query['index']
        if sha((ROOT / index['file']).read_bytes()) != index['sha256']:
            raise ValueError('Changed public index')
        counts, lines, stream_hash = {}, [], hashlib.sha256()
        with pysam.TabixFile(query['url'], index=str(ROOT / index['file'])) as table:
            if args.stage == 'reference':
                header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
            else:
                prefix = registry[query['id'] + '-prefix']
                import zlib
                header = zlib.decompressobj(31).decompress((ROOT / prefix['file']).read_bytes()).decode().splitlines()[0]
            for target in query['targets']:
                count = 0
                for line in table.fetch(query['chromosome'], target['position'] - 1, target['position']):
                    count += 1
                    if count > 100000:
                        raise ValueError('Site query row cap exceeded')
                    stream_hash.update((line + '\n').encode())
                    lines.append('\t'.join(line.split('\t')[:8]) if args.stage == 'reference' else line)
                counts[str(target['position'])] = count
                time.sleep(2)
        plain = (header + '\n' + '\n'.join(lines) + ('\n' if lines else '')).encode()
        data = gzip_bytes(plain)
        output = directory / (query['id'] + '-' + args.stage + '.tsv.gz')
        output.write_bytes(data)
        result = {'id': query['id'], 'url': query['url'], 'index_sha256': index['sha256'],
                  'file': str(output.relative_to(ROOT)), 'sha256': sha(data), 'uncompressed_sha256': sha(plain),
                  'all_returned_rows_sha256': stream_hash.hexdigest(), 'all_gene_or_site_counts': counts,
                  'rows': len(lines), 'completed_at_utc': datetime.now(timezone.utc).isoformat()}
        record['queries'].append(result)
        save_json(path, record)
        print(json.dumps({'id': query['id'], 'rows': len(lines), 'counts': counts}), flush=True)


if __name__ == '__main__':
    main()
