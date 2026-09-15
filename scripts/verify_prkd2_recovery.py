"""Verify the recovery record without fitting a new shared-signal model.

Default checks versioned outputs and independent R coverage calculations.
--with-cache additionally checks source and transferred-byte cache receipts.
"""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

from coloc_common import ROOT, save_json, sha


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(ROOT / path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t' if '.tsv' in str(path) else ','))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-cache', action='store_true')
    args = parser.parse_args()
    primary = json.loads((ROOT / 'reports/prkd2-recovery-analysis.json').read_text())
    secondary = json.loads((ROOT / 'reports/prkd2-recovery-secondary-audit.json').read_text())
    author = json.loads((ROOT / 'reports/prkd2-recovery-ibdverse-analysis.json').read_text())
    distribution = json.loads((ROOT / 'reports/prkd2-recovery-ibdverse-format-audit.json').read_text())
    reports = [primary, secondary, author]
    if any(item['new_H4_posteriors'] != 0 for item in reports):
        raise ValueError('Unexpected posterior in a recovery-only stage')
    outputs = [item for report in reports for item in report['outputs']] + [distribution['output']]
    for script, report in [('analyze_prkd2_recovery.py', primary), ('audit_prkd2_recovery_secondary.py', secondary),
                           ('analyze_prkd2_recovery_ibdverse.py', author), ('audit_prkd2_recovery_ibdverse_format.py', distribution)]:
        if sha((ROOT / 'scripts' / script).read_bytes()) != report['script_sha256']:
            raise ValueError('Changed analysis script: ' + script)
    for item in outputs:
        if sha((ROOT / item['file']).read_bytes()) != item['sha256']:
            raise ValueError('Changed analysis output: ' + item['file'])
    plan = json.loads((ROOT / 'config/prkd2-recovery-plan.json').read_text())
    cohorts = table('data/derived/prkd2-recovery-cohort-screen.csv')
    if {x['dataset_id'] for x in cohorts} != {x['dataset_id'] for x in plan['cohorts']} or len(cohorts) != len(plan['cohorts']):
        raise ValueError('Incomplete or duplicated cohort registry')
    targets = table('data/derived/prkd2-recovery-target-screen.csv')
    if len(targets) != 6 * len(cohorts):
        raise ValueError('Missing target outcomes')
    primary_rows = [x for x in targets if x['target_rsid'] == 'rs112445263']
    measured = [x for x in primary_rows if x['matching_rows'] != '']
    unmeasured = [x for x in primary_rows if x['matching_rows'] == '']
    if len(measured) != 33 or len(unmeasured) != 9 or any(float(x['matching_rows']) != 0 for x in measured):
        raise ValueError('Catalogue availability/absence summary differs')
    queries = {}
    for release in ('r7', 'r8_beta'):
        ledger = json.loads((ROOT / f'config/prkd2-recovery-{release}-queries.json').read_text())
        for q in ledger['queries']:
            queries[q['id']] = q
            if q.get('complete') and q['exact_primary_target_gene_rows'] != 0:
                raise ValueError('Catalogue primary-target count differs')
            if release == 'r7' and q.get('complete') and q['all_gene_point_counts']['46700099'] != 0:
                raise ValueError('Primary variant is present across genes')
    controls = json.loads((ROOT / 'config/prkd2-recovery-controls-audit.json').read_text())
    if len(controls['queries']) != 26 or any(set(q['all_gene_or_site_counts']) != {'46718300'} for q in controls['queries']):
        raise ValueError('Corrected control queries are incomplete')
    original = table('data/derived/prkd2-recovery-original-target-associations.csv')
    d = [x for x in original if x['dataset_id'] == 'DICE-original' and x['target_rsid'] == 'rs112445263']
    if len(d) != 1 or (d[0]['position'], d[0]['ref'], d[0]['alt'], d[0]['source_gene_id']) != ('47203356', 'C', 'A', 'ENSG00000105287.8'):
        raise ValueError('Recovered DICE identity differs')
    if (d[0]['source_pvalue'], d[0]['source_beta'], d[0]['source_se']) != ('0.224', '0.22', ''):
        raise ValueError('DICE source values or missing-SE treatment differs')
    if any(x['dataset_id'] == 'BLUEPRINT-original' and x['target_rsid'] == 'rs112445263' for x in original):
        raise ValueError('Unexpected BLUEPRINT primary association')
    hm_targets = table('data/derived/prkd2-recovery-harmonised-targets.csv')
    if [int(x['harmonised_exact_rows']) for x in hm_targets] != [1, 0, 0, 0, 0, 1]:
        raise ValueError('Harmonised target results differ')
    if any(int(x['unharmonised_source_candidate_rows']) != 0 for x in hm_targets):
        raise ValueError('Unharmonised source target result differs')
    if not all(x['OR_reciprocal_check'] for x in secondary['harmonised_controls']):
        raise ValueError('Harmonised risk-allele comparison failed')
    independent = table('data/derived/prkd2-recovery-independent-coverage.csv')
    expected = {**primary['original_study_overlap'], **{x['dataset_id']: x for x in author['contexts']}}
    if {x['dataset_id'] for x in independent} != set(expected):
        raise ValueError('Independent coverage check omits a source')
    errors = []
    for row in independent:
        result = expected[row['dataset_id']]
        if int(row['shared_edits']) != result['shared_distinct_edits'] or int(row['source_rows']) != result['full_source_gene_rows'] or int(row['distinct_source_rows']) != result['distinct_source_association_rows']:
            raise ValueError('Independent original source universe differs')
        for field, key in [('full_PSC_weight', 'full_PSC_marginal_weight_available'), ('full_RNA_weight', 'RNA_marginal_weight_available')]:
            expected_value = result[key]
            if expected_value is None:
                if row[field] != '':
                    raise ValueError('Independent check invented unavailable evidence')
            else:
                observed = float(row[field])
                errors.append(abs(observed - expected_value))
                if not math.isclose(observed, expected_value, rel_tol=2e-12, abs_tol=2e-12):
                    raise ValueError('Independent R evidence coverage differs')
        if result['multi_signal_model_qualified']:
            raise ValueError('Unexpected model qualification')
    format_checks = table('data/derived/prkd2-recovery-independent-format.csv')
    formats = {x['dataset_id']: x for x in distribution['contexts']}
    if len(format_checks) != 4 or distribution['new_H4_posteriors'] != 0:
        raise ValueError('Incomplete distribution-convention check')
    format_errors = []
    for row in format_checks:
        expected_format = formats[row['dataset_id']]
        if int(row['best_matching_integer_residual_df']) != expected_format['best_matching_integer_residual_df']:
            raise ValueError('Independent inferred residual df differs')
        for key in ('median_abs_log10_P_error', 'p95_abs_log10_P_error', 'max_abs_log10_P_error', 'normal_median_abs_log10_P_error', 'normal_max_abs_log10_P_error'):
            error = abs(float(row[key]) - expected_format[key])
            format_errors.append(error)
            if error > 2e-10:
                raise ValueError('Independent distribution-convention calculation differs')
    if [x['best_matching_integer_residual_df'] for x in distribution['contexts']] != [69, 66, 73, 2]:
        raise ValueError('Unexpected source distribution convention')
    checked, skipped = {}, set()

    def visit(item):
        if isinstance(item, dict):
            if 'file' in item and 'sha256' in item and isinstance(item['file'], str):
                name = item['file']
                if name.startswith('data/raw/') and not args.with_cache:
                    skipped.add(name)
                elif name not in checked:
                    file = ROOT / name
                    if not file.is_file() or sha(file.read_bytes()) != item['sha256']:
                        raise ValueError('Changed or missing receipt: ' + name)
                    checked[name] = item['sha256']
            for value in item.values():
                visit(value)
        elif isinstance(item, list):
            for value in item:
                visit(value)

    for path in sorted((ROOT / 'config').glob('prkd2-recovery-*.json')):
        visit(json.loads(path.read_text()))
    for report in reports + [distribution]:
        visit(report)
    record = {'all_checks_pass': True, 'verifier_sha256': sha(Path(__file__).read_bytes()),
              'cache_receipts_checked': args.with_cache, 'receipt_files_checked': len(checked),
              'raw_receipts_not_checked_in_portable_mode': len(skipped),
              'versioned_analysis_outputs': len(outputs), 'catalogue_contexts': len(cohorts),
              'measured_catalogue_contexts': len(measured), 'unavailable_catalogue_contexts': len(unmeasured),
              'catalogue_target_outcomes': len(targets), 'corrected_control_queries': len(controls['queries']),
              'original_contexts_with_independent_R_coverage': len(independent),
              'maximum_R_python_coverage_difference': max(errors), 'new_H4_posteriors': 0,
              'distribution_contexts_checked': len(format_checks),
              'maximum_R_python_distribution_difference': max(format_errors),
              'checked_file_sha256': checked,
              'limits': ['Receipt checks are not biological replication.', 'Whole original archive CRC checks have a separate record.',
                         'No new multi-signal model is qualified by this recovery audit.']}
    suffix = 'cache-verification' if args.with_cache else 'verification'
    save_json(ROOT / f'reports/prkd2-recovery-{suffix}.json', record)
    print(json.dumps({k: v for k, v in record.items() if k != 'checked_file_sha256'}, indent=2))


if __name__ == '__main__':
    main()
