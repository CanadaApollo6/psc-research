"""Retrieve complete fixed peak files and evaluate every fixed interval."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

import requests

from prkd2_mechanism_common import ROOT, RAW, PLAN, now, sha, save_new, write_csv


def overlaps(start0, end0, query_start0, query_end0):
    return start0 < query_end0 and query_start0 < end0


def download(file, offline=False):
    folder = RAW / 'peaks'
    folder.mkdir(parents=True, exist_ok=True)
    accession = file['file_accession']
    path = folder / (accession + '.bed.gz')
    receipt_path = folder / (accession + '-receipt.json')
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if sha(path) != receipt['sha256'] or path.stat().st_size != file['bytes']:
            raise ValueError('Cached peak file changed')
        return path, receipt
    if path.exists():
        raise ValueError('Unpinned peak file already exists')
    if offline:
        raise FileNotFoundError(path)
    started = now()
    part = path.with_suffix(path.suffix + '.part')
    md5, h256, total = hashlib.md5(), hashlib.sha256(), 0
    # Binary data download, not an ENCODE metadata API request. All metadata
    # selection and audit calls went through the required ENCODE skill.
    with requests.get(file['download_url'], headers={'User-Agent': 'PSC-public-research/1.0'}, timeout=45, stream=True) as response:
        response.raise_for_status()
        with part.open('wb') as stream:
            for chunk in response.iter_content(131072):
                total += len(chunk)
                if total > file['bytes']:
                    raise ValueError('Peak transfer exceeds exact selected size')
                stream.write(chunk)
                md5.update(chunk)
                h256.update(chunk)
        if total != file['bytes'] or md5.hexdigest() != file['publisher_md5']:
            raise ValueError('Peak file size or publisher MD5 mismatch')
        safe_headers = {k: response.headers[k] for k in ['Content-Type','Content-Length','ETag','Last-Modified'] if k in response.headers}
    part.replace(path)
    receipt = {'file_accession': accession, 'experiment': file['experiment'], 'url': file['download_url'],
               'started_at_utc': started, 'finished_at_utc': now(), 'bytes': total,
               'publisher_md5': file['publisher_md5'], 'computed_md5': md5.hexdigest(),
               'sha256': h256.hexdigest(), 'response_headers': safe_headers, 'file': str(path.relative_to(ROOT))}
    save_new(receipt_path, receipt)
    print(json.dumps({'peak_file': accession, 'bytes': total, 'publisher_md5_matches': True}), flush=True)
    return path, receipt


def parse_file(path, file, start0, end0):
    regional = []
    count, chr19_count = 0, 0
    with gzip.open(path, 'rt') as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) < 3:
                raise ValueError('Malformed BED row')
            chrom, begin, finish = fields[0], int(fields[1]), int(fields[2])
            if begin < 0 or finish <= begin:
                raise ValueError('Invalid BED coordinates')
            count += 1
            if chrom != 'chr19':
                continue
            chr19_count += 1
            if not overlaps(begin, finish, start0, end0):
                continue
            regional.append({'experiment': file['experiment'], 'file_accession': file['file_accession'], 'group': file['group'],
                             'source_line': line_number, 'chromosome': chrom, 'start0': begin, 'end0': finish,
                             'source_fields_json': json.dumps(fields, separators=(',', ':'))})
    return regional, {'total_source_rows': count, 'chr19_source_rows': chr19_count, 'regional_rows': len(regional), 'gzip_eof_verified': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--output-directory', type=Path, default=ROOT/'data/derived/prkd2-mechanism-encode')
    args = parser.parse_args()
    plan = json.loads(PLAN.read_text())
    inputs_path = ROOT/'config/prkd2-mechanism-inputs.json'
    inputs = json.loads(inputs_path.read_text())
    qualification_path = ROOT/'config/prkd2-mechanism-encode-audit-resolution.json'
    qualification = json.loads(qualification_path.read_text())
    if qualification['plan_sha256'] != sha(PLAN) or not all(r['qualified'] for r in qualification['selected']):
        raise ValueError('Full selected-file audit qualification required')
    for name, expected in qualification['source_sha256'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Qualified metadata changed')
    start0, end0 = inputs['variants'][0]['interval_start0'], inputs['variants'][0]['interval_end0']
    points = [{'label': v['rsid'], 'kind': v['role'], 'point1': v['position_grch38_1based']} for v in inputs['variants']]
    points.append({'label': 'PRKD2_fixed_QTL_TSS', 'kind': 'fixed TSS context', 'point1': plan['region']['tss_grch38_1based']})
    all_regional, summaries, receipts, hits = [], [], [], []
    for file in qualification['selected']:
        path, receipt = download(file, args.offline)
        regional, accounting = parse_file(path, file, start0, end0)
        all_regional.extend(regional)
        receipts.append({**receipt, **accounting})
        for point in points:
            pos0 = point['point1'] - 1
            exact = [r for r in regional if overlaps(r['start0'], r['end0'], pos0, pos0+1)]
            nearby = [r for r in regional if overlaps(r['start0'], r['end0'], pos0-250, pos0+251)]
            distances = [0 if r['start0'] <= pos0 < r['end0'] else r['start0'] - pos0 if pos0 < r['start0'] else pos0-(r['end0']-1) for r in regional]
            summaries.append({'label': point['label'], 'kind': point['kind'], 'position1': point['point1'],
                              'experiment': file['experiment'], 'file_accession': file['file_accession'], 'group': file['group'],
                              'exact_base_peak_count': len(exact), 'within_250bp_peak_count': len(nearby),
                              'nearest_regional_peak_distance_bp': min(distances) if distances else None,
                              'experiment_audit_levels': ';'.join(file['experiment_audit_levels']),
                              'file_audit_levels': ';'.join(file['file_audit_levels']),
                              'source_donor_ids': ';'.join(file['donor_ids']),
                              'interpretation': 'Observed regulatory context, not an allele effect or a proven gene connection'})
            for peak in nearby:
                hits.append({'label': point['label'], 'exact_base_overlap': peak in exact, **peak})
    out = args.output_directory
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out/'peak-overlap-summary.csv', summaries)
    if all_regional: write_csv(out/'regional-peaks.csv', all_regional)
    if hits: write_csv(out/'point-nearby-peaks.csv', hits)
    report = {'plan_sha256': sha(PLAN), 'input_sha256': sha(inputs_path), 'qualification_sha256': sha(qualification_path),
              'files': receipts, 'total_compressed_bytes': sum(r['bytes'] for r in receipts),
              'total_source_rows': sum(r['total_source_rows'] for r in receipts),
              'coordinate_convention': 'BED zero-based half-open; edited base [position1-1,position1), flank [position1-251,position1+250)',
              'points': points, 'summary_rows': len(summaries), 'regional_peaks': len(all_regional),
              'source_scope': 'Fixed eight experiments; three donor identifiers with extensive assay reuse. No allele-specific assay or direct gene-link measurement.',
              'analysis_script_sha256': sha(Path(__file__))}
    (out/'peak-analysis.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'files': len(receipts), 'compressed_bytes': report['total_compressed_bytes'],
                      'source_rows': report['total_source_rows'], 'regional_peaks': len(all_regional), 'summary_rows': len(summaries)}))


if __name__ == '__main__':
    main()
