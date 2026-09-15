"""Select measured regulatory experiments using metadata, before target overlaps."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from prkd2_mechanism_common import ROOT, RAW, PLAN, freeze, now, sha, save_new, encode_request, write_csv

OUT = ROOT / 'data/derived/prkd2-mechanism-encode'
PEAK_PRIORITIES = {
    'IDR thresholded peaks': 0,
    'optimal IDR thresholded peaks': 1,
    'replicated peaks': 2,
    'pseudoreplicated IDR thresholded peaks': 3,
    'peaks': 4,
}
MARKS = ['H3K27ac', 'H3K4me1', 'H3K4me3']


def audit_levels(obj):
    return sorted(k for k, v in obj.get('audit', {}).items() if v)


def assay_group(experiment):
    assay = experiment.get('assay_title')
    if assay in ['DNase-seq', 'ATAC-seq']:
        return assay
    if assay == 'Histone ChIP-seq' and experiment.get('target', {}).get('label') in MARKS:
        return experiment['target']['label']
    return None


def context_failures(experiment):
    failures = []
    ontology = experiment.get('biosample_ontology', {})
    if experiment.get('status') != 'released': failures.append('experiment not released')
    if ontology.get('term_name') != 'CD14-positive monocyte': failures.append('not exact CD14-monocyte context')
    if ontology.get('classification') != 'primary cell': failures.append('not primary cell')
    organisms = {r.get('library', {}).get('biosample', {}).get('organism', {}).get('scientific_name') for r in experiment.get('replicates', [])}
    if organisms != {'Homo sapiens'}: failures.append('human organism not verified for all replicates')
    if 'ERROR' in audit_levels(experiment): failures.append('experiment audit ERROR')
    return failures


def file_failures(file, max_bytes):
    failures = []
    if file.get('status') != 'released': failures.append('file not released')
    if file.get('assembly') != 'GRCh38': failures.append('file not GRCh38')
    if file.get('file_format') != 'bed': failures.append('not BED')
    if file.get('output_type') not in PEAK_PRIORITIES: failures.append('not an eligible peak output')
    if 'ERROR' in audit_levels(file): failures.append('file audit ERROR')
    if not isinstance(file.get('file_size'), int) or not 0 < file['file_size'] <= max_bytes:
        failures.append('unknown file size or transfer cap')
    if not file.get('md5sum'): failures.append('missing publisher MD5')
    return failures


def identity_rows(experiment):
    rows = []
    for replicate in experiment.get('replicates', []):
        library = replicate.get('library', {})
        sample = library.get('biosample', {})
        donor = sample.get('donor', {})
        if isinstance(donor, str): donor = {'@id': donor}
        rows.append({'experiment': experiment['accession'], 'biosample': sample.get('accession'),
                     'library': library.get('accession'), 'donor': donor.get('accession') or donor.get('@id'),
                     'biological_replicate': replicate.get('biological_replicate_number'),
                     'technical_replicate': replicate.get('technical_replicate_number'),
                     'organism': sample.get('organism', {}).get('scientific_name'),
                     'sample_summary': sample.get('summary'), 'treatment': json.dumps(sample.get('treatments', []), sort_keys=True)})
    return rows


def main():
    frozen = freeze()
    lock_path = ROOT / 'config/prkd2-mechanism-encode-selection.json'
    if lock_path.exists():
        lock = json.loads(lock_path.read_text())
        for path, expected in lock['metadata_sha256'].items():
            if sha(ROOT / path) != expected:
                raise ValueError('ENCODE metadata lock changed')
        print(json.dumps({'status': 'verified ENCODE metadata selection', 'experiments': len(lock['selected'])}))
        return
    plan = json.loads(PLAN.read_text())['measured_regulation']
    # An all-assay, still cell-specific query distinguishes empty ATAC search
    # from a provider/network error. The earlier 404 receipt stays preserved.
    all_cd14 = encode_request('all-released-cd14-experiments', 'search/',
                             {'type': 'Experiment', 'status': 'released',
                              'biosample_ontology.term_name': 'CD14-positive monocyte', 'limit': 100})
    if all_cd14 is None or not isinstance(all_cd14.get('total'), int) or all_cd14['total'] > 1000:
        raise ValueError('Bounded CD14 metadata inventory not recovered')
    complete = list(all_cd14['@graph'])
    if len(complete) != all_cd14['total']:
        # The portal rejected its from=100 pagination parameter (receipt retained).
        # Reissue the same filtered metadata query with a sufficient bounded limit.
        full = encode_request('all-released-cd14-experiments-limit200', 'search/',
                              {'type': 'Experiment', 'status': 'released',
                               'biosample_ontology.term_name': 'CD14-positive monocyte', 'limit': 200})
        if full is None or full.get('total') != all_cd14['total']:
            raise ValueError('CD14 metadata total changed during inventory retrieval')
        complete = full['@graph']
    if len(complete) != all_cd14['total'] or len({e['accession'] for e in complete}) != len(complete):
        raise ValueError('Complete unique CD14 inventory not recovered')
    experiments = sorted(complete, key=lambda r: r['accession'])
    counts = Counter(e.get('assay_title') for e in experiments)
    selected, statuses, file_audit, identities = [], [], [], []
    group_counts = Counter()
    total_bytes = 0
    OUT.mkdir(parents=True, exist_ok=True)
    for summary in experiments:
        accession = summary['accession']
        group = assay_group(summary)
        basic = {'experiment': accession, 'assay_title': summary.get('assay_title'),
                 'histone_mark': summary.get('target', {}).get('label'), 'selected_group': group,
                 'search_audit_levels': ';'.join(audit_levels(summary))}
        if group is None:
            statuses.append({**basic, 'status': 'outside prespecified assay/mark scope'})
            continue
        if group_counts[group] >= 2:
            statuses.append({**basic, 'status': 'later accession after two qualifying experiments'})
            continue
        reasons = context_failures(summary)
        if reasons:
            statuses.append({**basic, 'status': '; '.join(reasons)})
            continue
        experiment = encode_request(accession, f'experiments/{accession}/', {'frame': 'embedded'})
        if experiment is None:
            statuses.append({**basic, 'status': 'full metadata request unavailable'})
            continue
        reasons = context_failures(experiment)
        if reasons:
            statuses.append({**basic, 'status': '; '.join(reasons)})
            continue
        candidates = []
        for f in experiment.get('files', []):
            reasons = file_failures(f, plan['limits']['maximum_compressed_bytes_per_peak_file'])
            file_audit.append({'experiment': accession, 'file_accession': f.get('accession'),
                               'assembly': f.get('assembly'), 'file_status': f.get('status'),
                               'file_format': f.get('file_format'), 'output_type': f.get('output_type'),
                               'file_size': f.get('file_size'), 'audit_levels': ';'.join(audit_levels(f)),
                               'eligible': not reasons, 'reason': '; '.join(reasons),
                               'biological_replicates': ';'.join(map(str, f.get('biological_replicates', [])))})
            if not reasons:
                candidates.append(f)
        if not candidates:
            statuses.append({**basic, 'status': 'no eligible released GRCh38 peak file'})
            continue
        file = min(candidates, key=lambda f: (PEAK_PRIORITIES[f['output_type']], -len(f.get('biological_replicates', [])), f['accession']))
        # Resolve the selected file itself for full audit, checksum and download metadata.
        full_file = encode_request(file['accession'], f"files/{file['accession']}/", {'frame': 'embedded'})
        if full_file is None or file_failures(full_file, plan['limits']['maximum_compressed_bytes_per_peak_file']):
            raise ValueError('Selected file fails its direct object check; preserve selection for review')
        if any(full_file.get(k) != file.get(k) for k in ['md5sum', 'assembly', 'file_size', 'output_type', 'status']):
            raise ValueError('Selected file metadata changed between experiment and file object')
        total_bytes += file['file_size']
        if total_bytes > plan['limits']['maximum_total_compressed_peak_bytes']:
            raise ValueError('Total peak transfer cap exceeded before download')
        group_counts[group] += 1
        rep_rows = identity_rows(experiment)
        identities.extend(rep_rows)
        selected.append({'experiment': accession, 'group': group, 'assay_title': experiment['assay_title'],
                         'file_accession': file['accession'], 'file_format': file['file_format'],
                         'file_format_type': file.get('file_format_type'), 'output_type': file['output_type'],
                         'assembly': file['assembly'], 'bytes': file['file_size'], 'publisher_md5': file['md5sum'],
                         'download_url': 'https://www.encodeproject.org' + full_file['href'],
                         'experiment_audit_levels': audit_levels(experiment), 'file_audit_levels': audit_levels(full_file),
                         'biological_replicates': file.get('biological_replicates', []),
                         'biosample_ids': sorted({r['biosample'] for r in rep_rows if r['biosample']}),
                         'donor_ids': sorted({r['donor'] for r in rep_rows if r['donor']}),
                         'experiment_metadata_file': f'data/raw/prkd2-mechanism/encode/{accession}.json',
                         'file_metadata_file': f"data/raw/prkd2-mechanism/encode/{file['accession']}.json"})
        statuses.append({**basic, 'status': 'selected by prespecified accession/metadata order'})
    if len(selected) > plan['limits']['maximum_experiments']:
        raise ValueError('Experiment count cap exceeded')
    write_csv(OUT / 'experiment-screen.csv', statuses)
    if file_audit: write_csv(OUT / 'file-screen.csv', file_audit)
    if identities: write_csv(OUT / 'source-replicate-identities.csv', identities)
    inputs = [p for p in (RAW/'encode').glob('*.json')] + [p for p in (ROOT/'config/prkd2-mechanism-encode-requests').glob('*.json')] + sorted(OUT.glob('*.csv'))
    lock = {'locked_at_utc': now(), 'plan_sha256': frozen['plan_sha256'], 'peak_downloads_or_regional_inspection': False,
            'all_released_cd14_experiments': len(experiments), 'assay_counts': dict(sorted(counts.items())),
            'selected': selected, 'selected_counts': dict(sorted(group_counts.items())), 'total_selected_bytes': total_bytes,
            'empty_ATAC_interpretation': 'Direct ATAC search returned 404; complete successful all-assay CD14 inventory independently counts ATAC experiments.',
            'metadata_source_format': 'Full API JSON parsed and reserialized by the installed ENCODE skill; not byte-identical HTTP wire bodies.',
            'metadata_sha256': {str(p.relative_to(ROOT)): sha(p) for p in inputs}, 'selection_script_sha256': sha(Path(__file__))}
    save_new(lock_path, lock)
    print(json.dumps({'status': 'ENCODE selection locked before peak inspection', 'experiments': len(selected),
                      'groups': lock['selected_counts'], 'bytes': total_bytes, 'assay_counts': lock['assay_counts'],
                      'locked_at_utc': lock['locked_at_utc']}))


if __name__ == '__main__':
    main()
