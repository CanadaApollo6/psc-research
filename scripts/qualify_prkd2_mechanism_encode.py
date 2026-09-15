"""Recover computed ENCODE audit fields omitted by the embedded object frame."""
from __future__ import annotations

import json
from pathlib import Path

from prkd2_mechanism_common import ROOT, RAW, PLAN, now, sha, save_new, encode_request, write_csv


def main():
    selection_path = ROOT / 'config/prkd2-mechanism-encode-selection.json'
    selection = json.loads(selection_path.read_text())
    path = ROOT / 'config/prkd2-mechanism-encode-audit-resolution.json'
    if path.exists():
        record = json.loads(path.read_text())
        for name, expected in record['source_sha256'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Audit source changed')
        print(json.dumps({'status': 'verified computed audit resolution', 'qualified': sum(r['qualified'] for r in record['selected'])}))
        return
    all_summary = json.loads((RAW/'encode/all-released-cd14-experiments-limit200.json').read_text())
    experiments = {r['accession']: r for r in all_summary['@graph']}
    rows, details, sources = [], [], [selection_path, RAW/'encode/all-released-cd14-experiments-limit200.json']
    for selected in selection['selected']:
        accession = selected['file_accession']
        label = 'selected-file-audit-' + accession
        result = encode_request(label, 'search/', {'type': 'File', 'accession': accession, 'limit': 10})
        if result is None or result.get('total') != 1 or len(result.get('@graph', [])) != 1:
            raise ValueError('Computed selected-file audit unavailable')
        file = result['@graph'][0]
        if file['accession'] != accession or file.get('dataset') not in [f"/experiments/{selected['experiment']}/", selected['experiment']]:
            # Search may embed the experiment rather than return its path.
            dataset = file.get('dataset')
            if not isinstance(dataset, dict) or dataset.get('accession') != selected['experiment']:
                raise ValueError('Computed audit belongs to a different experiment')
        for field, expected in [('assembly', 'GRCh38'), ('status', 'released'), ('file_size', selected['bytes']), ('output_type', selected['output_type'])]:
            if file.get(field) != expected:
                raise ValueError('Selected-file metadata differs in audit response')
        experiment = experiments[selected['experiment']]
        exp_levels = sorted(k for k, v in experiment.get('audit', {}).items() if v)
        file_levels = sorted(k for k, v in file.get('audit', {}).items() if v)
        qualified = 'ERROR' not in exp_levels and 'ERROR' not in file_levels
        rows.append({**selected, 'experiment_audit_levels': exp_levels, 'file_audit_levels': file_levels,
                     'qualified': qualified, 'qualification_note': 'Computed audit search checked; embedded object frame omitted audit fields.'})
        for kind, obj in [('experiment', experiment), ('selected_file', file)]:
            for level, entries in obj.get('audit', {}).items():
                for entry in entries:
                    details.append({'experiment': selected['experiment'], 'file_accession': accession,
                                    'object_kind': kind, 'level': level, 'category': entry.get('category'),
                                    'path': entry.get('path'), 'detail': entry.get('detail')})
        sources.extend([RAW/f'encode/{label}.json', RAW/f'encode/{label}-receipt.json',
                        ROOT/f'config/prkd2-mechanism-encode-requests/{label}.json'])
    out = ROOT / 'data/derived/prkd2-mechanism-encode'
    if details: write_csv(out/'computed-audit-details.csv', details)
    summary = [{'experiment': r['experiment'], 'file_accession': r['file_accession'], 'group': r['group'],
                'experiment_audit_levels': ';'.join(r['experiment_audit_levels']),
                'file_audit_levels': ';'.join(r['file_audit_levels']), 'qualified': r['qualified'],
                'donor_ids': ';'.join(r['donor_ids'])} for r in rows]
    write_csv(out/'computed-audit-summary.csv', summary)
    record = {'completed_at_utc': now(), 'plan_sha256': sha(PLAN), 'selection_sha256': sha(selection_path),
              'reason': 'Embedded ENCODE object responses do not include computed audit fields. Recover them through exact File search and the full Experiment search, before peak download or regional peak inspection.',
              'selection_changed': False, 'new_prediction_scores_inspected_by_this_step': False,
              'selected': rows, 'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sources},
              'script_sha256': sha(Path(__file__))}
    save_new(path, record)
    print(json.dumps({'status': 'computed ENCODE audits checked before peaks', 'qualified': sum(r['qualified'] for r in rows),
                      'selected': len(rows), 'warning_or_noncompliant_files': sum(bool(set(r['file_audit_levels']) & {'WARNING','NOT_COMPLIANT'}) for r in rows)}))
    if not all(r['qualified'] for r in rows):
        raise SystemExit('A selected file is unqualified; preserve failure and apply the original selection rule explicitly before analysis')


if __name__ == '__main__':
    main()
