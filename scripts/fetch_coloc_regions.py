"""Acquire complete frozen regions as bounded public aggregate/reference subsets."""

import argparse
import csv
import gzip
import io
import json
import time
import urllib.request
import zlib
from datetime import datetime, timezone

import pysam

from coloc_common import ROOT, gzip_bytes, save_json, sha, verify_plan


def now():
    return datetime.now(timezone.utc).isoformat()


def read_gwas_range(data, offset):
    """Only complete interior lines; preserve their original absolute offsets."""
    lines = data.splitlines(keepends=True)
    records = []
    position = offset
    for index, line in enumerate(lines):
        start = position
        position += len(line)
        if index == 0 or index == len(lines)-1 or line.startswith(b'#'):
            continue
        fields = line.decode('ascii').split()
        if len(fields) != 13 or not fields[0].isdigit():
            raise ValueError('Unexpected complete GWAS row')
        records.append(((int(fields[0]), int(fields[2])), start, line.decode().rstrip('\r\n')))
    if any(a[0] > b[0] for a, b in zip(records, records[1:])):
        raise ValueError('GWAS source is not coordinate sorted')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['gwas', 'qtl', 'reference'])
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    plan = verify_plan()
    sources_path = ROOT / 'config/coloc-sources.json'
    sources = {s['id']: s for s in json.loads(sources_path.read_text())['sources']}
    for source in sources.values():
        if not source.get('sha256') or sha((ROOT/'data/raw'/source['file']).read_bytes()) != source['sha256']:
            raise ValueError('Supporting source is not pinned: ' + source['id'])
    path = ROOT / ('config/coloc-' + args.stage + '-queries.json')
    frozen = {'plan_sha256': sha((ROOT/'config/coloc-plan.json').read_bytes()),
              'source_manifest_sha256': sha(sources_path.read_bytes()),
              'reader_sha256': sha(__import__('pathlib').Path(__file__).read_bytes()), 'pysam_version': pysam.__version__}
    run = json.loads(path.read_text()) if path.exists() else {**frozen, 'started_at_utc': now(), 'queries': [], 'ranges': []}
    if any(run.get(k) != v for k, v in frozen.items()):
        raise ValueError('Reader or frozen query inputs changed')
    bounds = plan['transfer_bounds']
    output = ROOT/'data/derived/coloc-query-rows'
    output.mkdir(exist_ok=True)

    def complete(identifier):
        prior = next((q for q in run['queries'] if q['id'] == identifier), None)
        if prior:
            data = (ROOT/prior['file']).read_bytes()
            if sha(data) != prior['sha256'] or sha(gzip.decompress(data)) != prior['uncompressed_sha256']:
                raise ValueError('Regional cache checksum changed')
            return True
        if args.offline:
            raise ValueError('Missing completed region: ' + identifier)
        return False

    def retain(identifier, plain, filename, fields):
        data = gzip_bytes(plain)
        target = ROOT/filename
        if target.exists() and target.read_bytes() != data:
            raise ValueError('Unpinned regional file already differs')
        target.write_bytes(data)
        run['queries'].append({'id': identifier, 'file': filename, 'sha256': sha(data),
                               'uncompressed_sha256': sha(plain), 'bytes': len(data), 'finished_at_utc': now(), **fields})
        save_json(path, run)
        print(json.dumps({'id': identifier, 'rows': fields.get('rows'), 'status': 'complete'}), flush=True)

    def gwas_fetch(start, end):
        total = plan['gwas']['source_bytes']
        start, end = max(0, int(start)), min(total-1, int(end))
        filename = f'coloc-gwas-bytes-{start}-{end}.txt'
        known = next((r for r in run['ranges'] if r['file'] == filename), None)
        cache = ROOT/'data/raw'/filename
        if known:
            data = cache.read_bytes()
            if sha(data) != known['sha256']:
                raise ValueError('GWAS byte-range cache changed')
        else:
            if len(run['ranges']) >= bounds['maximum_gwas_range_requests'] or end-start+1 > bounds['maximum_gwas_range_bytes']:
                raise ValueError('GWAS region request cap exceeded')
            request = urllib.request.Request(plan['gwas']['url'], headers={'User-Agent': 'PSC-public-research/coloc-v01', 'Range': f'bytes={start}-{end}'})
            with urllib.request.urlopen(request, timeout=45) as response:
                if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{total}':
                    raise ValueError('GWAS range or source size changed')
                data = response.read(end-start+2)
            if len(data) != end-start+1:
                raise ValueError('Incomplete or oversized GWAS range')
            cache.write_bytes(data)
            run['ranges'].append({'file': filename, 'start': start, 'end': end, 'sha256': sha(data), 'bytes': len(data), 'retrieved_at_utc': now()})
            save_json(path, run)
        return read_gwas_range(data, start)

    def gwas_boundary(target):
        low, high = 0, plan['gwas']['source_bytes']-65536
        for _ in range(28):
            mid = (low+high)//2
            rows = gwas_fetch(mid, mid+65535)
            if rows[0][0] <= target <= rows[-1][0]:
                return next(r[1] for r in rows if r[0] >= target)
            if high-low <= 65536:
                rows = gwas_fetch(low-65536, high+131072)
                eligible = [r[1] for r in rows if r[0] >= target]
                if not eligible or target < rows[0][0]:
                    raise ValueError('Cannot bracket the GWAS boundary')
                return eligible[0]
            if target < rows[0][0]: high = mid
            else: low = mid
        raise ValueError('GWAS boundary search did not converge')

    if args.stage == 'gwas':
        header = next(line for line in (ROOT/'data/raw/GCST004030-original-header.txt').read_text().splitlines() if line.startswith('#chr '))
        for region in plan['regions']:
            identifier = 'GWAS-' + region['gene_name']
            if complete(identifier): continue
            chrom = int(region['chromosome'][3:])
            lo = (chrom, region['start_grch37_0based']+1)
            hi = (chrom, region['end_grch37_exclusive'])
            start, end = gwas_boundary(lo), gwas_boundary(hi)
            rows = gwas_fetch(start-65536, end+65536)
            if not (rows[0][0] < lo and rows[-1][0] > hi):
                raise ValueError('Final GWAS region lacks both flanking completeness witnesses')
            selected = [r[2] for r in rows if lo <= r[0] <= hi]
            plain = (header+'\n'+'\n'.join(selected)+'\n').encode()
            retain(identifier, plain, str((output/(identifier+'.txt.gz')).relative_to(ROOT)),
                   {'rows': len(selected), 'region': region, 'source_url': plan['gwas']['url'],
                    'bracketing_first_coordinate': rows[0][0], 'bracketing_last_coordinate': rows[-1][0]})
    elif args.stage == 'qtl':
        for dataset in plan['datasets']:
            ds = dataset['dataset_id']
            header = zlib.decompressobj(31).decompress((ROOT/'data/raw'/sources[ds+'-header']['file']).read_bytes()).decode().splitlines()[0]
            names = header.split('\t'); gene_column = names.index('gene_id')
            chrom_column, pos_column = names.index('chromosome'), names.index('position')
            for region in plan['regions']:
                identifier = ds+'-'+region['gene_name']
                if complete(identifier): continue
                chrom = region['chromosome'][3:]
                selected, scanned, size = [], 0, 0
                from hashlib import sha256
                stream_hash = sha256()
                with pysam.TabixFile(dataset['url'], index=str(ROOT/'data/raw'/sources[ds+'-index']['file'])) as table:
                    for line in table.fetch(chrom, region['start_grch38_0based'], region['end_grch38_exclusive']):
                        scanned += 1; size += len(line)+1; stream_hash.update((line+'\n').encode())
                        if scanned > bounds['maximum_nominal_rows_per_region'] or size > bounds['maximum_nominal_uncompressed_bytes_per_region']:
                            raise ValueError('Nominal regional query cap exceeded')
                        fields = line.split('\t')
                        if len(fields) != len(names) or fields[chrom_column] != chrom or not region['start_grch38_0based'] < int(fields[pos_column]) <= region['end_grch38_exclusive']:
                            raise ValueError('Unexpected nominal regional row')
                        if fields[gene_column].split('.')[0] == region['gene_id']:
                            selected.append(line)
                plain = (header+'\n'+'\n'.join(selected)+('\n' if selected else '')).encode()
                retain(identifier, plain, str((output/(identifier+'.tsv.gz')).relative_to(ROOT)),
                       {'rows': len(selected), 'region': region, 'dataset_id': ds, 'source_url': dataset['url'],
                        'all_gene_rows_scanned': scanned, 'all_gene_bytes_scanned': size, 'all_gene_stream_sha256': stream_hash.hexdigest()})
                time.sleep(bounds['minimum_seconds_between_indexed_queries'])
    else:
        with (ROOT/'data/raw'/sources['kgp-phase3-panel']['file']).open() as stream:
            eur = {r['sample'] for r in csv.DictReader(stream, delimiter='\t') if r['super_pop'] == 'EUR'}
        if len(eur) != 503:
            raise ValueError('Unexpected EUR donor registry')
        for region in plan['regions']:
            identifier = 'EUR-' + region['gene_name']
            if complete(identifier): continue
            chrom = region['chromosome'][3:]
            index = sources['match-kgp-chr'+chrom+'-index']
            url = index['url'].removesuffix('.tbi')
            selected, scanned, size = [], 0, 0
            with pysam.TabixFile(url, index=str(ROOT/'data/raw'/index['file'])) as table:
                header = next(line for line in table.header if line.startswith('#CHROM')).split('\t')
                columns = [i for i, name in enumerate(header) if name in eur]
                if len(columns) != 503: raise ValueError('Reference sample IDs are incomplete')
                names = [header[i] for i in columns]
                for line in table.fetch(chrom, region['start_grch37_0based'], region['end_grch37_exclusive']):
                    scanned += 1; size += len(line)+1
                    if scanned > bounds['maximum_reference_rows_per_region'] or size > bounds['maximum_reference_uncompressed_bytes_per_region']:
                        raise ValueError('Reference regional query cap exceeded')
                    fields = line.split('\t')
                    if len(fields) != len(header): raise ValueError('Malformed public VCF row')
                    if len(fields[3]) != 1 or len(fields[4]) != 1 or fields[3] not in 'ACGT' or fields[4] not in 'ACGT': continue
                    if fields[0] != chrom or not region['start_grch37_0based'] < int(fields[1]) <= region['end_grch37_exclusive']:
                        raise ValueError('Unexpected reference regional coordinate')
                    gt_index = fields[8].split(':').index('GT')
                    dosages = []
                    for column in columns:
                        gt = fields[column].split(':')[gt_index].replace('|', '/').split('/')
                        if gt == ['.', '.']: dosages.append('-1')
                        elif len(gt) == 2 and all(a in ('0','1') for a in gt): dosages.append(str(sum(map(int,gt))))
                        else: raise ValueError('Unsupported reference genotype')
                    selected.append('\t'.join(fields[:5]+dosages))
            plain = ('chromosome\tposition_grch37\trsid\tref\talt\t'+'\t'.join(names)+'\n'+'\n'.join(selected)+'\n').encode()
            retain(identifier, plain, 'data/raw/coloc-'+identifier+'-dosages.tsv.gz',
                   {'rows': len(selected), 'region': region, 'source_url': url, 'donors': 503,
                    'all_reference_rows_scanned': scanned, 'all_reference_bytes_scanned': size, 'source_build': 'GRCh37'})
            time.sleep(bounds['minimum_seconds_between_indexed_queries'])
    expected = 20 if args.stage == 'qtl' else 2
    if len(run['queries']) != expected:
        raise ValueError('Incomplete planned query set')
    if 'finished_at_utc' not in run:
        run['finished_at_utc'] = now(); save_json(path, run)
    print(json.dumps({'stage': args.stage, 'complete_regions': len(run['queries']), 'status': 'verified-cache' if args.offline else 'complete'}), flush=True)


if __name__ == '__main__':
    main()
