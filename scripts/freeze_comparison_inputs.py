"""Lock a fully matched, source-verified comparison before any prediction call."""

import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(name):
    return list(csv.DictReader(io.StringIO((ROOT / name).read_text())))


def main():
    protocol = json.loads((ROOT / 'config/psc-controlled-comparison.json').read_text())
    original = json.loads((ROOT / 'config/psc-controlled-comparison-lock.json').read_text())
    if digest(ROOT / 'config/psc-controlled-comparison.json') != original['protocol_sha256']:
        raise ValueError('Original protocol changed')
    for name in ('config/finemap-sources.json', 'config/psc-comparator-sources.json'):
        for source in json.loads((ROOT / name).read_text())['sources']:
            if digest(ROOT / 'data/raw' / source['file']) != source['sha256']:
                raise ValueError('Source checksum changed: ' + source['id'])
    statuses = csv_rows('data/derived/psc-comparator-anchor-status.csv')
    covariates = {r['rsid']: r for r in csv_rows('data/derived/psc-comparator-covariates.csv')}
    windows = {w['archive_region']: w for w in csv_rows('data/derived/psc-comparison-windows.csv')}
    genes = csv_rows('data/derived/psc-comparison-gene-universe.csv')
    metadata_record = json.loads((ROOT / 'config/psc-comparator-model-metadata.json').read_text())
    metadata_path = ROOT / 'data/raw' / metadata_record['file']
    if digest(metadata_path) != metadata_record['sha256']:
        raise ValueError('Model metadata changed')
    tracks = [r for r in csv.DictReader(io.StringIO(metadata_path.read_text())) if r['ontology_curie'] == protocol['model']['primary_ontology']]
    grouped = {}
    for row in tracks:
        grouped.setdefault(row['name'], []).append(row['strand'])
    if not grouped or any(sorted(strands) not in (['+', '-'], ['.']) for strands in grouped.values()):
        raise ValueError('Unexpected primary track/strand inventory')
    inputs, regional_genes = [], {}
    for status in statuses:
        if status['status'] != 'matched':
            continue
        label = status['archive_region']
        universe = [g for g in genes if g['archive_region'] == label and g['present_in_model_transcript_tss_universe'] == 'True']
        if not universe or len({g['gene_id'] for g in universe}) != len(universe):
            raise ValueError('Missing or duplicate fixed gene universe')
        regional_genes[label] = universe
        for order, rsid in enumerate((status['anchor'], status['comparator_1'], status['comparator_2'])):
            row = covariates[rsid]
            if row['reference_bases_verified'] != 'True' or row['genotype_status'] != 'exact_alleles_present':
                raise ValueError('Input reference or exact genotype missing')
            inputs.append({'archive_region': label, 'matched_anchor': status['anchor'], 'role': 'anchor' if order == 0 else 'comparator',
                           'rsid': rsid, 'chromosome': 'chr' + row['chromosome'],
                           'position_grch37_1based': int(row['position_grch37_1based']), 'position_grch38_1based': int(row['position_grch38_1based']),
                           'reference_bases': row['reference_bases'], 'alternate_bases': row['alternate_bases'], 'archived_pip': row['pip'],
                           'interval_start_0based': int(windows[label]['interval_start_0based']), 'interval_end_0based_exclusive': int(windows[label]['interval_end_0based_exclusive']),
                           'nearest_tss_gene_id': row['nearest_tss_gene_id'], 'nearest_tss_gene_name': row['nearest_tss_gene_name'],
                           'nearest_tss_distance_bp': int(row['nearest_tss_distance_bp'])})
    if len({r['rsid'] for r in inputs}) != len(inputs) or len(inputs) > protocol['model']['maximum_new_variant_requests']:
        raise ValueError('Repeated variant or request budget exceeded')
    for row in inputs:
        if row['interval_end_0based_exclusive'] - row['interval_start_0based'] != protocol['model']['interval_width'] or not row['interval_start_0based'] <= row['position_grch38_1based'] - 1 < row['interval_end_0based_exclusive']:
            raise ValueError('Invalid model interval')
    files = ['config/psc-controlled-comparison.json', 'config/psc-controlled-comparison-lock.json', 'config/finemap-sources.json',
             'config/psc-comparator-sources.json', 'config/psc-comparator-model-metadata.json', 'config/psc-comparator-genotype-queries.json',
             'data/derived/psc-comparison-gene-universe.csv', 'data/derived/psc-comparator-covariates.csv',
             'data/derived/psc-comparator-anchor-status.csv', 'data/derived/psc-comparator-pair-audit.csv', 'data/derived/psc-comparator-ld-audit.csv',
             'scripts/prepare_matching_covariates.py', 'scripts/match_comparison_variants.py']
    if inputs:
        files.append('data/derived/psc-comparator-matches.csv')
    payload = {'protocol_version': protocol['protocol_version'], 'ready_for_inference': bool(inputs),
               'inputs': inputs, 'gene_universe_by_region': regional_genes, 'primary_track_metadata': tracks,
               'expected_merged_primary_track_names': sorted(grouped), 'model_metadata_sha256': metadata_record['sha256'],
               'frozen_file_sha256': {name: digest(ROOT / name) for name in files}, 'all_anchor_statuses': statuses,
               'limitations': ['Archived fine-mapping run provenance unresolved', 'Unresolved high-PIP identifiers/alleles outside LD screen',
                               'ALL_FOLDS; genomic and donor independence not established', 'Three planned regions; descriptive statistical-priority comparison only']}
    path = ROOT / 'config/psc-matched-comparison-inputs.json'
    if path.exists():
        previous = json.loads(path.read_text())
        timestamp = previous.pop('frozen_at_utc')
        if previous != payload:
            raise ValueError('Existing final comparison lock differs; preserve it and review explicitly')
        print(json.dumps({'status': 'verified-existing-lock', 'requests': len(inputs), 'frozen_at_utc': timestamp}))
    else:
        payload['frozen_at_utc'] = datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(payload, indent=2) + '\n')
        print(json.dumps({'status': 'locked-before-predictions', 'requests': len(inputs), 'frozen_at_utc': payload['frozen_at_utc']}))


if __name__ == '__main__':
    main()
