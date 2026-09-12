"""Query all 90 frozen CD4 dataset/variant positions with bounded indexed reads."""

import argparse
import gzip
import hashlib
import io
import json
import time
import zlib
from datetime import datetime, timezone
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.part')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def verify_inputs():
    plan = json.loads((ROOT / 'config/cd4-rna-plan.json').read_text())
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Input preceding the CD4 plan changed: ' + name)
    if sha((ROOT / 'data/derived/cd4-rna-planned-pairs.csv').read_bytes()) != plan['planned_pairs_sha256']:
        raise ValueError('Planned pair universe changed')
    sources = json.loads((ROOT / 'config/cd4-rna-sources.json').read_text())['sources']
    for source in sources:
        if not source.get('sha256') or sha((ROOT / 'data/raw' / source['file']).read_bytes()) != source['sha256']:
            raise ValueError('CD4 source is not pinned or changed: ' + source['id'])
    return plan, {s['id']: s for s in sources}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--offline', action='store_true')
    group.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    plan, sources = verify_inputs()
    identities = ['config/cd4-rna-plan.json', 'config/cd4-rna-sources.json', 'data/derived/cd4-rna-planned-pairs.csv']
    lock_payload = {'frozen_file_sha256': {p: sha((ROOT / p).read_bytes()) for p in identities}, 'pysam_version': pysam.__version__}
    lock_path = ROOT / 'config/cd4-rna-input-lock.json'
    if lock_path.exists():
        previous = json.loads(lock_path.read_text()); previous.pop('frozen_at_utc')
        if previous != lock_payload:
            raise ValueError('Existing CD4 input lock differs')
    elif args.offline:
        raise ValueError('No prior CD4 input lock exists')
    else:
        write_json(lock_path, {**lock_payload, 'frozen_at_utc': datetime.now(timezone.utc).isoformat()})
    run_path = ROOT / 'config/cd4-rna-queries.json'
    run = json.loads(run_path.read_text()) if run_path.exists() else {
        'input_lock_sha256': sha(lock_path.read_bytes()), 'reader_sha256': sha(Path(__file__).read_bytes()),
        'started_at_utc': datetime.now(timezone.utc).isoformat(), 'queries': [], 'status': 'running'}
    if run['input_lock_sha256'] != sha(lock_path.read_bytes()) or run['reader_sha256'] != sha(Path(__file__).read_bytes()):
        raise ValueError('Query inputs or reader changed; preserve the existing run')
    output = ROOT / 'data/derived/cd4-rna-query-rows'
    output.mkdir(exist_ok=True)
    bounds = plan['transfer_bounds']
    for dataset in plan['datasets']:
        dataset_id = dataset['dataset_id']
        index = sources[dataset_id + '-index']
        header_source = sources[dataset_id + '-header']
        header = zlib.decompressobj(31).decompress((ROOT / 'data/raw' / header_source['file']).read_bytes()).decode().splitlines()[0]
        required = {'variant', 'gene_id', 'molecular_trait_id', 'chromosome', 'position', 'ref', 'alt', 'beta', 'se', 'pvalue', 'an', 'ac'}
        if not required <= set(header.split('\t')):
            raise ValueError('Source header lacks required gene-level fields')
        for variant in plan['variants']:
            identifier = dataset_id + '-' + variant['rsid']
            previous = next((q for q in run['queries'] if q['query_id'] == identifier), None)
            if previous and previous['status'] == 'complete':
                data = (ROOT / previous['file']).read_bytes()
                if sha(data) != previous['sha256'] or sha(gzip.decompress(data)) != previous['uncompressed_sha256']:
                    raise ValueError('Cached nominal QTL rows changed')
                if not args.refresh:
                    continue
            if args.offline:
                raise ValueError('Missing completed query: ' + identifier)
            chromosome = variant['chromosome'].removeprefix('chr')
            position = variant['position_grch38_1based']
            record = {'query_id': identifier, 'dataset_id': dataset_id, 'rsid': variant['rsid'], 'url': dataset['url'],
                      'chromosome': chromosome, 'start0': position - 1, 'end0': position,
                      'queried_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'running',
                      'index_sha256': index['sha256'], 'header_sha256': header_source['sha256']}
            try:
                queried, nbytes = [], 0
                with pysam.TabixFile(dataset['url'], index=str(ROOT / 'data/raw' / index['file'])) as table:
                    contig_present = chromosome in table.contigs
                    iterator = table.fetch(chromosome, position-1, position) if contig_present else []
                    for line in iterator:
                        nbytes += len(line)
                        if len(queried) >= bounds['maximum_rows_per_position'] or nbytes > bounds['maximum_uncompressed_bytes_per_position']:
                            raise ValueError('Position query exceeded its transfer bound')
                        fields = dict(zip(header.split('\t'), line.split('\t')))
                        if len(line.split('\t')) != len(header.split('\t')) or fields['chromosome'] != chromosome or int(fields['position']) != position:
                            raise ValueError('Indexed query returned an inconsistent row')
                        queried.append(line)
                plain = (header + '\n' + '\n'.join(queried) + ('\n' if queried else '')).encode()
                compressed = io.BytesIO()
                with gzip.GzipFile(fileobj=compressed, filename='', mode='wb', mtime=0) as stream:
                    stream.write(plain)
                data = compressed.getvalue()
                path = output / (identifier + '.tsv.gz')
                record.update(status='complete', rows=len(queried), contig_present=contig_present,
                              file=str(path.relative_to(ROOT)), sha256=sha(data), uncompressed_sha256=sha(plain), bytes=len(data))
                if previous and previous['status'] == 'complete':
                    if previous['uncompressed_sha256'] != record['uncompressed_sha256']:
                        raise ValueError('Remote nominal QTL response changed; existing result retained')
                    print(json.dumps({'query_id': identifier, 'status': 'remote-matches-cache', 'rows': len(queried)}), flush=True)
                    time.sleep(bounds['minimum_seconds_between_queries'])
                    continue
                path.write_bytes(data)
            except Exception as exc:
                if previous and previous['status'] == 'complete':
                    raise
                record.update(status='failed', error_type=type(exc).__name__, error_detail=str(exc))
            run['queries'] = [q for q in run['queries'] if q['query_id'] != identifier] + [record]
            write_json(run_path, run)
            print(json.dumps({'query_id': identifier, 'status': record['status'], 'rows': record.get('rows')}), flush=True)
            time.sleep(bounds['minimum_seconds_between_queries'])
    complete = sum(q['status'] == 'complete' for q in run['queries'])
    if not args.offline and not args.refresh:
        if run['status'] != 'complete':
            run.update(status='complete' if complete == bounds['queries'] else 'incomplete', finished_at_utc=datetime.now(timezone.utc).isoformat())
            write_json(run_path, run)
    print(json.dumps({'status': 'verified-cache' if args.offline else run['status'], 'complete_queries': complete,
                      'planned_queries': bounds['queries'], 'rows': sum(q.get('rows', 0) for q in run['queries'])}))


if __name__ == '__main__':
    main()
