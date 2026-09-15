"""Pin full reference alleles and two-build sequence for the PRKD2 extension."""

import csv
import gzip
import hashlib
import io
import json
import urllib.request
from collections import Counter
from datetime import datetime, timezone

import pysam

from prkd2_common import ROOT, gzip_bytes, save_json, sha, verify_plan


def digest_file(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    original = verify_plan()
    plan_path = ROOT / 'config/prkd2-expanded-plan.json'
    plan = json.loads(plan_path.read_text())
    assert plan['original_plan_sha256'] == sha((ROOT / 'config/prkd2-plan.json').read_bytes())
    output = ROOT / 'config/prkd2-expanded-sources.json'
    if output.exists():
        record = json.loads(output.read_text())
        for item in record['files']:
            assert digest_file(ROOT / item['file']) == item['sha256'], item['file']
        print(json.dumps({'status': 'verified-cache', 'files': len(record['files'])}))
        return
    record = {'plan_sha256': sha(plan_path.read_bytes()), 'started_at_utc': datetime.now(timezone.utc).isoformat(), 'files': []}
    base = 'https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/'
    for name, cap in [('Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz', 50000000), ('CHECKSUMS', 2000000)]:
        path = ROOT / ('data/raw/prkd2-' + name)
        if not path.exists():
            with urllib.request.urlopen(base + name, timeout=60) as response:
                data = response.read(cap + 1)
                if response.status != 200 or len(data) > cap:
                    raise ValueError('Unexpected reference download')
                path.write_bytes(data)
        record['files'].append({'file': str(path.relative_to(ROOT)), 'url': base + name,
                                'sha256': digest_file(path), 'bytes': path.stat().st_size})
    sequence37 = ''.join(line.strip() for line in gzip.decompress((ROOT / 'data/raw/prkd2-Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz').read_bytes()).decode().splitlines() if not line.startswith('>'))
    if len(sequence37) != 59128983 or set(sequence37) - set('ACGTN'):
        raise ValueError('Unexpected GRCh37 chromosome 19')
    region = original['regions'][0]
    padding = plan['bounds']['reference_padding_bases']
    starts = {'GRCh37': region['start_grch37_0based'] - padding, 'GRCh38': region['start_grch38_0based'] - padding}
    ends = {'GRCh37': region['end_grch37_exclusive'] + padding, 'GRCh38': region['end_grch38_exclusive'] + padding}
    reference38 = ROOT / 'work/ubash3a-raw-read-audit/GRCh38.primary.fa'
    with pysam.FastaFile(str(reference38)) as reference:
        if reference.get_reference_length('19') != 58617616:
            raise ValueError('Unexpected GRCh38 chromosome 19')
        sequence38 = reference.fetch('19', starts['GRCh38'], ends['GRCh38'])
    # Whole FASTA hash is recorded; the upstream compressed archive remains pinned separately.
    record['local_GRCh38_reference'] = {'file': str(reference38.relative_to(ROOT)), 'sha256': digest_file(reference38),
        'compressed_source_receipt': 'data/raw/ubash3a-raw-read-audit/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.receipt.json',
        'compressed_source_receipt_sha256': sha((ROOT / 'data/raw/ubash3a-raw-read-audit/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.receipt.json').read_bytes())}
    for build, sequence in [('GRCh37', sequence37[starts['GRCh37']:ends['GRCh37']]), ('GRCh38', sequence38)]:
        if len(sequence) != ends[build] - starts[build] or set(sequence) - set('ACGTN'):
            raise ValueError('Reference extraction failed')
        path = ROOT / ('data/raw/prkd2-' + build + '-region.txt')
        path.write_text(sequence + '\n')
        record['files'].append({'file': str(path.relative_to(ROOT)), 'sha256': digest_file(path), 'bytes': path.stat().st_size,
            'build': build, 'chromosome': 'chr19', 'start_1based': starts[build] + 1, 'end_1based': ends[build],
            'sequence_sha256': sha(sequence.encode()), 'length': len(sequence)})
    print(json.dumps({'status': 'two-build-sequences-pinned'}), flush=True)
    source = next(s for s in json.loads((ROOT / 'config/prkd2-sources.json').read_text())['sources'] if s['id'] == 'match-kgp-chr19-index')
    url = source['url'].removesuffix('.tbi')
    panel = ROOT / 'data/raw/junction-kgp-phase3-panel.data'
    eur = {r['sample'] for r in csv.DictReader(panel.open(), delimiter='\t') if r['super_pop'] == 'EUR'}
    assert len(eur) == 503
    selected, scanned, total_bytes = [], 0, 0
    statuses = Counter()
    digest = hashlib.sha256()
    with pysam.TabixFile(url, index=str(ROOT / 'data/raw' / source['file'])) as table:
        header = next(line for line in table.header if line.startswith('#CHROM')).split('\t')
        columns = [i for i, name in enumerate(header) if name in eur]
        assert len(columns) == 503
        for line in table.fetch('19', starts['GRCh37'], ends['GRCh37']):
            scanned += 1
            total_bytes += len(line) + 1
            digest.update((line + '\n').encode())
            if scanned > plan['bounds']['maximum_reference_vcf_rows'] or total_bytes > plan['bounds']['maximum_reference_vcf_uncompressed_bytes']:
                raise ValueError('Reference query exceeds frozen bounds')
            fields = line.split('\t')
            if len(fields) != len(header) or fields[0] != '19' or not starts['GRCh37'] < int(fields[1]) <= ends['GRCh37']:
                raise ValueError('Malformed or out-of-region public variant')
            if ',' in fields[4]:
                statuses['multiallelic'] += 1
                continue
            if not fields[3] or not fields[4] or set(fields[3] + fields[4]) - set('ACGT'):
                statuses['symbolic_or_non_DNA'] += 1
                continue
            if fields[3] == fields[4]:
                raise ValueError('Identical REF and ALT in reference')
            gt_index = fields[8].split(':').index('GT')
            values = []
            for column in columns:
                gt = fields[column].split(':')[gt_index].replace('|', '/').split('/')
                if gt == ['.', '.']:
                    values.append('-1')
                elif len(gt) == 2 and all(a in ('0', '1') for a in gt):
                    values.append(str(sum(map(int, gt))))
                else:
                    raise ValueError('Unexpected biallelic reference genotype')
            statuses['retained_biallelic'] += 1
            selected.append('\t'.join(fields[:5] + values))
    plain = ('chromosome\tposition_grch37\trsid\tref\talt\t' + '\t'.join(header[i] for i in columns) + '\n' + '\n'.join(selected) + '\n').encode()
    path = ROOT / 'data/raw/prkd2-expanded-EUR-dosages.tsv.gz'
    path.write_bytes(gzip_bytes(plain))
    record['files'].append({'file': str(path.relative_to(ROOT)), 'url': url, 'sha256': digest_file(path), 'bytes': path.stat().st_size,
        'uncompressed_sha256': sha(plain), 'all_reference_rows_scanned': scanned, 'all_reference_bytes_scanned': total_bytes,
        'all_reference_stream_sha256': digest.hexdigest(), 'status_counts': dict(statuses),
        'source_index_sha256': source['sha256'], 'start_0based': starts['GRCh37'], 'end_exclusive': ends['GRCh37'],
        'donors': 503, 'source_build': 'GRCh37', 'sample_bearing_file_excluded_from_Git': True})
    record.update(finished_at_utc=datetime.now(timezone.utc).isoformat(), fetcher_sha256=sha(__import__('pathlib').Path(__file__).read_bytes()))
    save_json(output, record)
    print(json.dumps({'status': 'complete', 'files': len(record['files']), 'reference_rows': len(selected), 'statuses': dict(statuses)}), flush=True)


if __name__ == '__main__':
    main()
