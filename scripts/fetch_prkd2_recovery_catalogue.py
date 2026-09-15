"""Screen every frozen monocyte dataset without selecting by target effects."""

import argparse
import csv
import gzip
import io
import json
import time
import zlib
from datetime import datetime, timezone
from pathlib import Path

import pyarrow.parquet as pq
import pysam

from coloc_common import ROOT, gzip_bytes, save_json, sha
from fetch_coverage_repair_parquet import RangeReader
from prkd2_recovery_sources import fetch


class RecoveryRangeReader(RangeReader):
    def __init__(self, spec, manifest, manifest_path, offline):
        self.spec, self.manifest, self.manifest_path = spec, manifest, manifest_path
        self.offline, self.position = offline, 0
        self.directory = ROOT / 'data/raw/prkd2-recovery-parquet-ranges'
        self.directory.mkdir(exist_ok=True)
        if 'remote_size' not in spec:
            self.fetch(0, 8)


def exact_target(row, gene):
    return (row.get('gene_id', '').split('.')[0] == gene and
            str(row.get('chromosome', '')).removeprefix('chr') == '19' and
            int(row['position']) == 46700099 and row.get('ref') == 'C' and row.get('alt') == 'A')


def retain(identifier, header, lines):
    directory = ROOT / 'data/derived/prkd2-recovery-queries'
    directory.mkdir(exist_ok=True)
    plain = (header + '\n' + '\n'.join(lines) + ('\n' if lines else '')).encode()
    path = directory / (identifier + '.tsv.gz')
    data = gzip_bytes(plain)
    if path.exists() and path.read_bytes() != data:
        raise ValueError('Changed extracted rows')
    path.write_bytes(data)
    return {'file': str(path.relative_to(ROOT)), 'rows': len(lines), 'sha256': sha(data), 'uncompressed_sha256': sha(plain)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('release', choices=['r7', 'r8_beta'])
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    plan_path = ROOT / 'config/prkd2-recovery-plan.json'
    plan = json.loads(plan_path.read_text())
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest:
            raise ValueError('Changed frozen input: ' + name)
    path = ROOT / ('config/prkd2-recovery-' + args.release + '-queries.json')
    pins = {'plan_sha256': sha(plan_path.read_bytes()), 'reader_sha256': sha(Path(__file__).read_bytes()),
            'range_dependency_sha256': sha((ROOT / 'scripts/fetch_coverage_repair_parquet.py').read_bytes())}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'queries': []}
    if any(record.get(k) != v for k, v in pins.items()):
        raise ValueError('Changed query plan or reader')
    positions = sorted({int(r['matched_edit'].split(':')[1]) for r in plan['prior_missing_evidence']} | {46702450})
    for ds in plan['cohorts']:
        if ds['release'] != args.release or ds['prior_analysis']:
            continue
        identifier = ds['dataset_id']
        spec = next((s for s in record['queries'] if s['id'] == identifier), None)
        if spec and spec.get('complete'):
            for item in spec.get('outputs', []):
                if sha((ROOT / item['file']).read_bytes()) != item['sha256']:
                    raise ValueError('Changed query output')
            continue
        if spec is None:
            spec = {'id': identifier, 'dataset': ds, 'started_at_utc': datetime.now(timezone.utc).isoformat(), 'outputs': []}
            record['queries'].append(spec)
            save_json(path, record)
        try:
            if args.release == 'r7':
                url = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/' + ds['study_id'] + '/' + identifier + '/' + identifier + '.all.tsv.gz'
                spec['url'] = url
                index = fetch(identifier + '-index', url + '.tbi', limit=10000000)
                prefix = fetch(identifier + '-prefix', url, byte_range='bytes=0-65535', limit=65536)
                if index.get('error') or prefix.get('error'):
                    raise ValueError('Archive index/header access failed; not biological absence')
                header = zlib.decompressobj(31).decompress((ROOT / prefix['file']).read_bytes()).decode().splitlines()[0]
                names = header.split('\t')
                if not {'gene_id', 'position', 'ref', 'alt'}.issubset(names):
                    raise ValueError('Unresolved gene/allele schema')
                lines = []
                point_counts = {}
                with pysam.TabixFile(url, index=str(ROOT / index['file'])) as table:
                    for position in positions:
                        point_lines = list(table.fetch('19', position - 1, position))
                        if len(point_lines) > plan['query_rules']['r7_point_row_cap']:
                            raise ValueError('Frozen point row cap exceeded')
                        point_counts[str(position)] = len(point_lines)
                        lines.extend(point_lines)
                        time.sleep(2)
                    spec['outputs'].append(retain(identifier + '-points-all-genes', header, lines))
                    rows = list(csv.DictReader(io.StringIO(header + '\n' + '\n'.join(lines)), delimiter='\t'))
                    spec['exact_primary_target_gene_rows'] = sum(exact_target(r, plan['gene_id']) for r in rows)
                    spec['all_gene_point_counts'] = point_counts
                    if spec['exact_primary_target_gene_rows']:
                        selected, count = [], 0
                        for line in table.fetch('19', plan['start_grch38_0based'], plan['end_grch38_exclusive']):
                            count += 1
                            if count > plan['query_rules']['r7_region_row_cap']:
                                raise ValueError('Frozen regional row cap exceeded')
                            fields = line.split('\t')
                            if len(fields) != len(names):
                                raise ValueError('Malformed regional row')
                            if fields[names.index('gene_id')].split('.')[0] == plan['gene_id']:
                                selected.append(line)
                        spec['outputs'].append(retain(identifier + '-PRKD2-region', header, selected))
                        spec['regional_all_gene_rows'] = count
                spec['status'] = 'candidate-for-coverage-audit' if spec['exact_primary_target_gene_rows'] else 'missing-required-PSC-variant'
            else:
                spec.setdefault('url', 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/sumstats/' + ds['study_id'] + '/' + identifier + '/' + identifier + '.all.parquet')
                spec.setdefault('ranges', [])
                reader = RecoveryRangeReader(spec, record, path, args.offline)
                parquet = pq.ParquetFile(reader)
                names = parquet.schema_arrow.names
                gene_column = 'gene_id' if 'gene_id' in names else 'molecular_trait_id'
                if gene_column not in names or 'position' not in names:
                    raise ValueError('Unresolved Parquet schema')
                candidate, selected = [], []
                for i in range(parquet.metadata.num_row_groups):
                    statistics = parquet.metadata.row_group(i).column(names.index('position')).statistics
                    if not statistics or not statistics.has_min_max:
                        raise ValueError('No position bounds; whole-archive scan not authorized by plan')
                    if statistics.min <= plan['end_grch38_exclusive'] and statistics.max > plan['start_grch38_0based']:
                        candidate.append(i)
                for i in candidate:
                    ids = parquet.read_row_group(i, columns=['chromosome', gene_column], use_threads=False).to_pandas()
                    keep = ids.chromosome.astype(str).str.removeprefix('chr').eq('19') & ids[gene_column].str.split('.').str[0].eq(plan['gene_id'])
                    if keep.any():
                        selected.append(i)
                candidate_rows = sum(parquet.metadata.row_group(i).num_rows for i in selected)
                spec['metadata'] = {'schema': str(parquet.schema_arrow), 'total_rows': parquet.metadata.num_rows,
                                    'total_row_groups': parquet.metadata.num_row_groups, 'identifier_row_groups': candidate,
                                    'selected_row_groups': selected, 'selected_candidate_rows': candidate_rows}
                save_json(path, record)
                if candidate_rows > plan['query_rules']['r8_candidate_row_cap']:
                    raise ValueError('Frozen nominal row cap exceeded')
                if selected:
                    frame = parquet.read_row_groups(selected, use_threads=False).to_pandas()
                    frame = frame[frame[gene_column].str.split('.').str[0].eq(plan['gene_id']) & frame.chromosome.astype(str).str.removeprefix('chr').eq('19')]
                    spec['full_gene_rows'] = len(frame)
                    frame = frame[frame.position.gt(plan['start_grch38_0based']) & frame.position.le(plan['end_grch38_exclusive'])]
                else:
                    frame = parquet.read_row_groups([], use_threads=False).to_pandas()
                    spec['full_gene_rows'] = 0
                plain = frame.to_csv(sep='\t', index=False, lineterminator='\n')
                head, *body = plain.splitlines()
                spec['outputs'].append(retain(identifier + '-PRKD2-region', head, body))
                rows = frame.to_dict('records')
                spec['exact_primary_target_gene_rows'] = sum(exact_target(r, plan['gene_id']) for r in rows)
                spec['status'] = 'candidate-for-coverage-audit' if spec['exact_primary_target_gene_rows'] else 'missing-required-PSC-variant'
                spec['transfer_bytes'] = sum(r['bytes'] for r in spec['ranges'])
            spec['complete'] = True
            spec['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
        except Exception as exc:
            spec['error'] = str(exc)
            spec['status'] = 'access-or-format-failure'
            spec['complete'] = False
        save_json(path, record)
        print(json.dumps({'dataset': identifier, 'status': spec['status'], 'target_rows': spec.get('exact_primary_target_gene_rows'), 'error': spec.get('error')}), flush=True)
        time.sleep(2)


if __name__ == '__main__':
    main()
