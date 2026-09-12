"""Reproduce the public GEUVADIS two-junction comparison using donors as units.

No network or model requests. Full public sample-level inputs and the joined
donor table remain in ignored data/raw; aggregate results are versioned.
"""

import csv
import gzip
import hashlib
import io
import json
import re
import platform
from importlib.metadata import version
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from fetch_sources import validate

ROOT = Path(__file__).resolve().parents[1]
EUROPEAN = {'CEU', 'FIN', 'GBR', 'TSI'}
csv.field_size_limit(2_000_000)


class UnestimableError(ValueError):
    """A mathematically unsupported model must remain a reported non-result."""


def read_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))


def write_csv(path, rows):
    if not rows:
        raise ValueError('Cannot publish an empty table without an explicit schema')
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def unique_index(rows, column):
    result = {}
    for row in rows:
        key = row[column]
        if not key or key in result:
            raise ValueError(f'Duplicate or missing {column}')
        result[key] = row
    return result


def parse_junction(text, target):
    rows = read_tsv(text)
    if len(rows) != 1:
        raise ValueError('Expected exactly one complete junction response; missing is not zero')
    row = rows[0]
    actual = (row['DataSource:Type'], row['chromosome'], int(row['start']) - 1,
              int(row['end']), row['strand'])
    expected = ('srav3h:I', target['chromosome'], target['start_0based'],
                target['end_0based_exclusive'], target['biological_strand'])
    if actual != expected or int(row['length']) != actual[3] - actual[2]:
        raise ValueError('Junction build/coordinates/strand/length do not match the fixed endpoint')
    if (row['left_motif'], row['right_motif']) != ('GT', 'AG'):
        raise ValueError('Unexpected motif at the predefined positive-strand endpoints')
    values = row['samples']
    if not values.startswith(','):
        raise ValueError('Malformed sparse sample-count list')
    counts = {}
    for pair in values[1:].split(','):
        sample, count_text = pair.split(':')
        count = int(count_text)
        if not sample.isdigit() or sample in counts or count <= 0:
            raise ValueError('Duplicate sample or nonpositive sparse read count')
        counts[sample] = count
    if len(counts) != int(row['samples_count']) or sum(counts.values()) != int(row['coverage_sum']):
        raise ValueError('Sparse count list is incomplete or disagrees with source totals')
    return counts, {k: row[k] for k in ('snaptron_id', 'start', 'end', 'strand',
                    'left_motif', 'right_motif', 'samples_count', 'coverage_sum')}


def parse_genotypes(text):
    lines = [l for l in text.splitlines() if not l.startswith('##')]
    if len(lines) != 2:
        raise ValueError('Expected one VCF header and one variant')
    header, row = [l.split('\t') for l in lines]
    if header[:2] != ['#CHROM', 'POS'] or len(header) != len(row):
        raise ValueError('Malformed VCF')
    if (row[0], int(row[1]), row[3], row[4]) != ('21', 43855067, 'A', 'C'):
        raise ValueError('VCF does not contain the verified GRCh37 A>C variant')
    if len(set(header[9:])) != len(header[9:]):
        raise ValueError('Duplicate VCF donor IDs')
    gt_index = row[8].split(':').index('GT')
    result = {}
    for donor, value in zip(header[9:], row[9:]):
        gt = value.split(':')[gt_index]
        if '.' in gt:
            result[donor] = None
            continue
        alleles = re.split(r'[|/]', gt)
        if len(alleles) != 2 or not set(alleles) <= {'0', '1'}:
            raise ValueError('Expected a diploid biallelic genotype')
        result[donor] = sum(int(a) for a in alleles)
    return result


def hc3_linear(x, y):
    """OLS coefficients and HC3 covariance; reject unidentified designs."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Non-finite regression data')
    if len(y) <= x.shape[1] or np.linalg.matrix_rank(x) < x.shape[1]:
        raise UnestimableError('Insufficient independent observations for regression')
    inv = np.linalg.inv(x.T @ x)
    beta = inv @ x.T @ y
    resid = y - x @ beta
    leverage = np.einsum('ij,jk,ik->i', x, inv, x)
    if np.any(leverage >= 1 - 1e-10):
        raise UnestimableError('HC3 is undefined for unit leverage (a population group has only one retained donor)')
    omega = (resid / (1 - leverage)) ** 2
    cov = inv @ (x.T @ (x * omega[:, None])) @ inv
    return beta, cov


def fit_ratio(donors, threshold):
    selected = [r for r in donors if r['population'] in EUROPEAN and
                r['dosage_c'] is not None and r['pair_total'] >= threshold]
    populations = sorted({r['population'] for r in selected})
    x = [[1, r['dosage_c']] + [int(r['population'] == p) for p in populations[1:]] for r in selected]
    y = [r['plus29'] / r['pair_total'] for r in selected]
    result = {'minimum_pair_reads': threshold, 'cohort': 'CEU/FIN/GBR/TSI',
            'donors': len(selected), 'unit': 'donor', 'outcome': 'plus29/(canonical+plus29)',
            'status': 'not_estimable', 'reason': '',
            'coefficient_per_C': '', 'hc3_se': '', 'wald_95_lower': '', 'wald_95_upper': '',
            'covariates': 'intercept; C dosage; population indicators',
            'inference': 'exploratory observational association; normal-approximation HC3 interval'}
    try:
        beta, cov = hc3_linear(x, y)
    except UnestimableError as exc:
        result['reason'] = str(exc)
        result['inference'] = 'No allele-effect estimate or confidence interval reported'
        return result
    se = float(np.sqrt(cov[1, 1]))
    result.update(status='estimated', coefficient_per_C=float(beta[1]), hc3_se=se,
                  wald_95_lower=float(beta[1] - 1.96 * se), wald_95_upper=float(beta[1] + 1.96 * se))
    return result


def summarize_groups(donors, threshold):
    rows = []
    for cohort, pops in [('European', EUROPEAN), ('YRI', {'YRI'})]:
        for dosage in (0, 1, 2):
            selected = [r for r in donors if r['population'] in pops and r['dosage_c'] == dosage
                        and r['pair_total'] >= threshold]
            ratios = [r['plus29'] / r['pair_total'] for r in selected]
            rows.append({'cohort': cohort, 'minimum_pair_reads': threshold,
                         'genotype': ('AA', 'AC', 'CC')[dosage], 'dosage_c': dosage,
                         'donors': len(selected), 'canonical_reads': sum(r['canonical'] for r in selected),
                         'plus29_reads': sum(r['plus29'] for r in selected),
                         'donors_with_plus29': sum(r['plus29'] > 0 for r in selected),
                         'mean_donor_ratio': float(np.mean(ratios)) if ratios else '',
                         'median_donor_ratio': float(np.median(ratios)) if ratios else '',
                         'q25_donor_ratio': float(np.quantile(ratios, .25)) if ratios else '',
                         'q75_donor_ratio': float(np.quantile(ratios, .75)) if ratios else ''})
    return rows


def aggregate_donors(by_donor, panel, genotypes):
    donors = []
    for donor, libs in sorted(by_donor.items()):
        a, b = sum(l['canonical'] for l in libs), sum(l['plus29'] for l in libs)
        donors.append({'donor': donor, 'population': panel.get(donor, {}).get('pop', ''),
                       'dosage_c': genotypes.get(donor), 'rna_runs': len(libs),
                       'canonical': a, 'plus29': b, 'pair_total': a + b,
                       'run_accessions': ';'.join(sorted(l['run'] for l in libs)),
                       'sequencing_centers': ';'.join(sorted({l['center'] for l in libs}))})
    return donors


def main():
    run = json.loads((ROOT / 'config/ubash3a-junction-run.json').read_text())
    plan_path = ROOT / run['plan_before_counts']['file']
    if hashlib.sha256(plan_path.read_bytes()).hexdigest() != run['plan_before_counts']['sha256']:
        raise ValueError('The plan recorded before the counts has changed')
    plan = json.loads(plan_path.read_text())
    source_list = json.loads((ROOT / 'config/ubash3a-junction-sources.json').read_text())['sources']
    sources = {s['id']: s for s in source_list}
    for s in source_list:
        validate((ROOT / 'data/raw' / s['file']).read_bytes(), s)

    def source(name, compressed=False):
        s = sources[name]
        content = (ROOT / 'data/raw' / s['file']).read_bytes()
        validate(content, s)
        return (gzip.decompress(content) if compressed else content).decode()

    counts, source_stats = {}, {}
    for target in plan['junctions']:
        name = target['name']
        counts[name], source_stats[name] = parse_junction(source(f'snaptron-{name}-unfiltered'), target)
    metadata = unique_index(read_tsv(source('recount-geuvadis-sra', True)), 'rail_id')
    projects = unique_index(read_tsv(source('recount-geuvadis-project', True)), 'rail_id')
    qc = unique_index(read_tsv(source('recount-geuvadis-qc', True)), 'rail_id')
    ena = unique_index(read_tsv(source('ena-geuvadis-runs')), 'run_accession')
    if set(metadata) != set(projects) or set(metadata) != set(qc):
        raise ValueError('Processed sample inventories disagree')
    direct_rows = []
    for batch in range(1, 5):
        direct_rows.extend(read_tsv(source(f'snaptron-geuvadis-samples-batch{batch}')))
    direct_samples = unique_index(direct_rows, 'rail_id')
    if set(metadata) != set(direct_samples):
        raise ValueError('Snaptron and recount3 processed sample inventories disagree')
    for rail_id, row in metadata.items():
        if (direct_samples[rail_id]['run_acc'], direct_samples[rail_id]['study_acc']) != (row['external_id'], row['study']):
            raise ValueError('Snaptron and recount3 sample identifiers disagree')
    panel = unique_index(read_tsv(source('kgp-phase3-panel')), 'sample')
    g = run['genotype_query']
    genotype_content = (ROOT / g['file']).read_bytes()
    if hashlib.sha256(genotype_content).hexdigest() != g['sha256']:
        raise ValueError('Genotype subset changed')
    genotypes = parse_genotypes(genotype_content.decode())
    by_donor = defaultdict(list)
    flow = Counter({'processed_rna_runs': len(metadata)})
    exclusions = []
    for rail_id, row in metadata.items():
        accession = row['external_id']
        if row['study'] != 'ERP001942' or row['library_strategy'] != 'RNA-Seq':
            raise ValueError('Unexpected study or assay')
        for other in (projects[rail_id], qc[rail_id]):
            if (other['external_id'], other['study']) != (accession, 'ERP001942'):
                raise ValueError('Mismatched processed sample identifiers')
        erow = ena.get(accession)
        match = re.fullmatch(r'GEUV:((?:NA|HG)\d+)', erow['sample_alias']) if erow else None
        if not match or not row['library_name'].startswith(match[1] + '.'):
            flow['excluded_runs_ambiguous_donor'] += 1
            exclusions.append({'run': accession, 'reason': 'ambiguous donor mapping'})
            continue
        if erow['study_accession'] != 'PRJEB3366':
            raise ValueError('ENA study accession changed')
        if float(qc[rail_id]['star.all_mapped_reads_both']) <= 0 or int(qc[rail_id]['junction_count']) <= 0:
            flow['excluded_runs_no_confirmed_processing'] += 1
            exclusions.append({'run': accession, 'reason': 'no confirmed mapped reads/junction processing'})
            continue
        # A verified, processed library absent from a complete sparse junction
        # response has zero observed supporting reads. This is not a missing QTL.
        by_donor[match[1]].append({'run': accession, 'rail_id': rail_id, 'center': row['run_center'],
                                  'canonical': counts['canonical'].get(rail_id, 0),
                                  'plus29': counts['plus29'].get(rail_id, 0)})
    donors = aggregate_donors(by_donor, panel, genotypes)
    flow.update({'unique_rna_donors': len(donors),
                 'donors_with_multiple_runs': sum(r['rna_runs'] > 1 for r in donors),
                 'donors_with_public_genotype': sum(r['dosage_c'] is not None for r in donors),
                 'donors_without_public_genotype': sum(r['dosage_c'] is None for r in donors)})
    for threshold in (10, 20):
        for cohort, pops in [('European', EUROPEAN), ('YRI', {'YRI'})]:
            group = [r for r in donors if r['population'] in pops and r['dosage_c'] is not None]
            flow[f'{cohort}_genotyped_donors'] = len(group)
            flow[f'{cohort}_donors_at_{threshold}_pair_reads'] = sum(r['pair_total'] >= threshold for r in group)
    write_csv(ROOT / 'data/raw/junction-geuvadis-donors.csv', donors)
    write_csv(ROOT / 'reports/ubash3a-junction-genotype-groups.csv',
              summarize_groups(donors, 10) + summarize_groups(donors, 20))
    fits = [fit_ratio(donors, threshold) for threshold in (10, 20)]
    write_csv(ROOT / 'reports/ubash3a-junction-associations.csv', fits)
    result = {'analysis': 'exploratory reanalysis of a previously studied public cohort',
              'software': {'python': platform.python_version(), 'numpy': np.__version__,
                           'matplotlib': version('matplotlib'), 'pysam': version('pysam')},
              'genome_build_junctions': 'GRCh38', 'genotype_build_source': 'GRCh37',
              'flow': dict(flow), 'excluded_runs': exclusions, 'source_junction_totals': source_stats,
              'plan_before_counts': run['plan_before_counts'],
              'genotype_subset_sha256': g['sha256'],
              'donor_table_sha256': hashlib.sha256((ROOT / 'data/raw/junction-geuvadis-donors.csv').read_bytes()).hexdigest(),
              'source_hashes': {s['id']: s['sha256'] for s in source_list},
              'fits': fits,
              'limitations': ['ALL_FOLDS model; not a held-out benchmark', 'GEUVADIS was used previously by Ji et al. 2017',
                              'Lymphoblastoid cells, not CD4 or PSC tissue', 'All-read junction counts; mapping/overhang details not resolved per supporting read',
                              'Two-junction composition, not all transcripts or full intron retention',
                              'Population-adjusted exploratory model; library/batch effects and relatedness not fully modeled']}
    (ROOT / 'reports/ubash3a-junction-audit.json').write_text(json.dumps(result, indent=2) + '\n')
    plot_coverage(donors)
    print(json.dumps({'flow': dict(flow), 'fits': fits}, indent=2))


def plot_coverage(donors):
    import os
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['svg.hashsalt'] = 'psc-ubash3a-junction-coverage'
    import matplotlib.pyplot as plt
    european = [r for r in donors if r['population'] in EUROPEAN and r['dosage_c'] is not None]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={'width_ratios': [1.25, 1]})
    totals = [sum(r['dosage_c'] == d for r in european) for d in (0, 1, 2)]
    retained = [sum(r['dosage_c'] == d and r['pair_total'] >= 10 for r in european) for d in (0, 1, 2)]
    axes[0].bar(np.arange(3), totals, color='#ccd5dc', label='All genotyped European donors')
    axes[0].bar(np.arange(3), retained, color='#287e8a', label='At least 10 junction reads')
    for i, (n, kept) in enumerate(zip(totals, retained)):
        axes[0].text(i, n + 3, f'{kept}/{n} retained', ha='center', fontsize=10)
    axes[0].set(xticks=np.arange(3), xticklabels=['AA', 'AC', 'CC'], ylabel='Donors',
                ylim=(0, max(totals) * 1.25), title='Very few donors meet the fixed coverage rule')
    axes[0].legend(frameon=False, loc='upper right', fontsize=8)
    x = [r['pair_total'] for r in european]
    bins = [0, 1, 2, 3, 5, 10, 20, np.inf]
    depth_counts, _ = np.histogram(x, bins=bins)
    axes[1].bar(np.arange(7), depth_counts, color=['#ccd5dc'] * 5 + ['#287e8a'] * 2)
    axes[1].axvline(4.5, color='#b05935', linestyle='--')
    axes[1].text(4.6, axes[1].get_ylim()[1] * .8, 'Minimum: 10', color='#954724', fontsize=10)
    axes[1].set(xlabel='Combined reads at the two junctions per donor', ylabel='Donors',
                xticks=np.arange(7), xticklabels=['0', '1', '2', '3–4', '5–9', '10–19', '20+'],
                title='Coverage limits the comparison')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.15)
        ax.set_axisbelow(True)
    fig.suptitle('UBASH3A in public GEUVADIS lymphoblastoid-cell RNA-seq', fontsize=14, x=.03, ha='left')
    fig.text(.03, .035, 'Repeated sequencing runs are combined by donor. Sparse coverage is an inconclusive result, not evidence of no effect.', fontsize=9)
    fig.tight_layout(rect=[.01, .09, .99, .92])
    for suffix in ('png', 'svg'):
        fig.savefig(ROOT / f'reports/ubash3a-junction-coverage.{suffix}', dpi=160, facecolor='white',
                    metadata={'Date': None} if suffix == 'svg' else {})
    plt.close(fig)


if __name__ == '__main__':
    main()
