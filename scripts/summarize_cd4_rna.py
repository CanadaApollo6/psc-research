"""Account for every frozen CD4 variant/gene test, including missing evidence."""

import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from fetch_cd4_rna import verify_inputs

ROOT = Path(__file__).resolve().parents[1]


def gene_id(value):
    if not re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', value):
        raise ValueError('Not an unambiguous human Ensembl gene ID')
    return value.split('.')[0]


def exact_allele_rows(rows, chromosome, position, ref, alt):
    expected = (chromosome.removeprefix('chr'), int(position), ref, alt)
    matches = []
    for row in rows:
        key = (row['chromosome'].removeprefix('chr'), int(row['position']), row['ref'], row['alt'])
        declared = 'chr' + key[0] + '_' + str(key[1]) + '_' + key[2] + '_' + key[3]
        if row['variant'] != declared:
            raise ValueError('Source variant ID disagrees with explicit coordinate/allele fields')
        if key == expected:
            matches.append(row)
    return matches


def unique_gene_measurement(rows, target):
    selected = [row for row in rows if gene_id(row['gene_id']) == target and gene_id(row['molecular_trait_id']) == target]
    if not selected:
        return None, 0, ''
    # The archive explicitly duplicates measurements over rsID aliases.
    unique = {tuple(sorted((key, value) for key, value in row.items() if key != 'rsid')) for row in selected}
    if len(unique) != 1:
        raise ValueError('Conflicting measurements for the same exact variant/gene')
    return selected[0], len(selected), ';'.join(sorted({row.get('rsid', '') for row in selected}))


def optional_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def measured_numbers(row, metadata_n):
    beta, se, p, an, ac = [float(row[key]) for key in ('beta', 'se', 'pvalue', 'an', 'ac')]
    if not all(math.isfinite(x) for x in (beta, se, p, an, ac)) or se <= 0 or not 0 <= p <= 1:
        raise ValueError('Nonfinite or invalid molecular-QTL coefficient/uncertainty')
    if an <= 0 or an != int(an) or int(an) % 2 or an > 2 * int(metadata_n) or not 0 <= ac <= an:
        raise ValueError('Allele counts incompatible with autosomal sample metadata')
    r2 = optional_number(row.get('r2'))
    return {'beta_per_alt': beta, 'standard_error': se, 'nominal_pvalue': p, 'wald_95_lower': beta - 1.96 * se,
            'wald_95_upper': beta + 1.96 * se, 'sample_size_from_an': int(an) // 2,
            'imputation_r2': r2, 'imputation_below_0_8': r2 < .8 if r2 is not None else None,
            'median_tpm': optional_number(row.get('median_tpm')), 'pvalue_reported_as_zero': p == 0}


def direction_agreement(direction, beta):
    if direction not in ('increase', 'decrease') or beta == 0:
        return None
    return (beta > 0) == (direction == 'increase')


def missing_pair_reason(exact_rows, annotation, within_cis):
    if not exact_rows:
        return 'exact_variant_absent'
    if annotation is None:
        return 'gene_annotation_absent'
    if not within_cis:
        return 'outside_source_cis_window'
    return 'gene_variant_pair_absent'


def parse_tsv(data):
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))


def main():
    plan, sources = verify_inputs()
    query_path = ROOT / 'config/cd4-rna-queries.json'
    run = json.loads(query_path.read_text())
    lock_path = ROOT / 'config/cd4-rna-input-lock.json'
    if hashlib.sha256(lock_path.read_bytes()).hexdigest() != run['input_lock_sha256']:
        raise ValueError('QTL run input lock changed')
    lock = json.loads(lock_path.read_text())
    for name, digest in lock['frozen_file_sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen QTL inputs changed')
    pairs = list(csv.DictReader((ROOT / 'data/derived/cd4-rna-planned-pairs.csv').open()))
    annotation_rows = parse_tsv(gzip.decompress((ROOT / 'data/raw' / sources['cd4-gene-phenotypes']['file']).read_bytes()))
    annotation = {}
    for row in annotation_rows:
        identifier = gene_id(row['gene_id'])
        if gene_id(row['phenotype_id']) != identifier or identifier in annotation:
            raise ValueError('Ambiguous gene annotation')
        annotation[identifier] = row
    queries = {(q['dataset_id'], q['rsid']): q for q in run['queries']}
    if len(queries) != len(run['queries']):
        raise ValueError('Duplicate query records')
    rows_by_query = {}
    for key, query in queries.items():
        if query['status'] != 'complete':
            continue
        data = (ROOT / query['file']).read_bytes()
        plain = gzip.decompress(data)
        if hashlib.sha256(data).hexdigest() != query['sha256'] or hashlib.sha256(plain).hexdigest() != query['uncompressed_sha256']:
            raise ValueError('Nominal association checksum changed')
        rows_by_query[key] = parse_tsv(plain)
        if len(rows_by_query[key]) != query['rows']:
            raise ValueError('Nominal association row count changed')
    results = []
    threshold = plan['multiplicity']['p_threshold']
    for dataset in plan['datasets']:
        identifier = dataset['dataset_id']
        for pair in pairs:
            source_gene = annotation.get(pair['gene_id'])
            within_cis = (source_gene['chromosome'] == pair['chromosome'].removeprefix('chr') and
                          abs(int(source_gene['phenotype_pos']) - int(pair['position_grch38_1based'])) <= 1_000_000) if source_gene else None
            assay_model = float(pair['model_total_rna_score'] if dataset['cohort_label'] == 'BLUEPRINT' else pair['model_polya_rna_score'])
            result = {**pair, 'dataset_id': identifier, 'cohort': dataset['cohort_label'], 'sample_group': dataset['sample_group'],
                      'source_tissue_label': dataset['tissue_label'], 'context_note': dataset['context_note'],
                      'metadata_sample_size': int(dataset['sample_size']), 'reporting_context': dataset['reporting_context'],
                      'source_annotation_present': source_gene is not None,
                      'source_gene_name': source_gene['gene_name'] if source_gene else '',
                      'source_gene_type': source_gene['gene_type'] if source_gene else '',
                      'source_gene_version': source_gene['gene_version'] if source_gene else '', 'within_source_cis_window': within_cis,
                      'assay_matched_model_raw_score': assay_model, 'status': '', 'absence_note': '',
                      'source_rows_at_position': None, 'exact_variant_source_rows': None, 'alias_row_count': None, 'source_rsid_aliases': '',
                      'source_file': '', 'beta_per_alt': None, 'standard_error': None, 'nominal_pvalue': None,
                      'wald_95_lower': None, 'wald_95_upper': None, 'sample_size_from_an': None,
                      'imputation_r2': None, 'imputation_below_0_8': None, 'median_tpm': None, 'pvalue_reported_as_zero': None,
                      'bonferroni_2280_pass': None, 'primary_direction_agrees': None, 'assay_secondary_direction_agrees': None,
                      'source_beta': '', 'source_se': '', 'source_pvalue': '', 'source_ac': '', 'source_an': '', 'source_maf': '', 'source_ma_samples': '',
                      'source_gene_id': '', 'source_molecular_trait_id': ''}
            key = identifier, pair['rsid']
            query = queries.get(key)
            if not query or query['status'] != 'complete':
                result.update(status='query_unavailable', absence_note='Source query incomplete; no biological inference')
                results.append(result); continue
            raw = rows_by_query[key]
            result.update(source_rows_at_position=len(raw), source_file=query['file'])
            try:
                exact = exact_allele_rows(raw, pair['chromosome'], pair['position_grch38_1based'], pair['reference_bases'], pair['alternate_bases'])
                result['exact_variant_source_rows'] = len(exact)
                measurement, aliases, alias_ids = unique_gene_measurement(exact, pair['gene_id'])
            except ValueError as exc:
                result.update(status='conflicting_or_incompatible_source_rows', absence_note=str(exc))
                results.append(result); continue
            if measurement is None:
                result.update(status=missing_pair_reason(exact, source_gene, within_cis),
                              absence_note='No usable nominal association for this exact planned pair; not a measured null')
                results.append(result); continue
            result.update(alias_row_count=aliases, source_rsid_aliases=alias_ids)
            for field in ['beta', 'se', 'pvalue', 'ac', 'an', 'maf', 'ma_samples', 'gene_id', 'molecular_trait_id']:
                result['source_' + field] = measurement.get(field, '')
            try:
                numbers = measured_numbers(measurement, dataset['sample_size'])
                result.update(numbers, status='measured', bonferroni_2280_pass=numbers['nominal_pvalue'] <= threshold,
                              primary_direction_agrees=direction_agreement(pair['model_consensus_direction'], numbers['beta_per_alt']),
                              assay_secondary_direction_agrees=direction_agreement('increase' if assay_model > 0 else 'decrease' if assay_model < 0 else 'mixed_or_zero', numbers['beta_per_alt']))
            except (ValueError, TypeError) as exc:
                result.update(status='invalid_numerical_measurement', absence_note=str(exc))
            results.append(result)
    frame = pd.DataFrame(results)
    if len(frame) != plan['planned_dataset_variant_gene_tests'] or frame.duplicated(['dataset_id', 'rsid', 'gene_id']).any():
        raise ValueError('Complete planned coverage was not retained')
    reports = ROOT / 'reports'
    frame.to_csv(reports / 'cd4-rna-evidence.csv', index=False)
    hits = frame[frame.bonferroni_2280_pass == True].sort_values(['nominal_pvalue','dataset_id','rsid','gene_id'])
    hits.to_csv(reports / 'cd4-rna-screened-associations.csv', index=False)
    focus = frame[(frame.model_gene_rank.astype(int) == 1) | (frame.is_region_named_gene == 'True')].copy()
    focus['panel'] = np.where(focus.model_gene_rank.astype(int) == 1, 'model_maximum_gene', 'region_named_gene')
    focus.to_csv(reports / 'cd4-rna-focus-genes.csv', index=False)
    coverage = []
    for keys, group in frame.groupby(['dataset_id','cohort','sample_group','archive_region','role'], sort=False):
        record = dict(zip(['dataset_id','cohort','sample_group','archive_region','role'], keys))
        record.update(planned_pairs=len(group), measured_pairs=int((group.status == 'measured').sum()),
                      threshold_passes=int((group.bonferroni_2280_pass == True).sum()),
                      primary_sign_agrees=int((group.primary_direction_agrees == True).sum()),
                      primary_sign_opposes=int((group.primary_direction_agrees == False).sum()),
                      measured_with_mixed_model_tracks=int(((group.status == 'measured') & (group.model_consensus_direction == 'mixed_or_zero')).sum()))
        coverage.append(record)
    pd.DataFrame(coverage).to_csv(reports / 'cd4-rna-coverage.csv', index=False)
    donor_rows = []
    for accession in ['ENCSR033XWU', 'ENCSR411MUF', 'ENCSR463JBR', 'ENCSR545MEZ']:
        data = json.loads((ROOT / 'data/raw' / sources[accession]['file']).read_text())
        for replicate in data.get('replicates', []):
            biosample = replicate.get('library', {}).get('biosample', {})
            donor = biosample.get('donor', {})
            donor_rows.append({'experiment': accession, 'assay': data['assay_title'], 'lab': data['lab']['name'],
                               'biosample': biosample.get('accession', ''), 'donor': donor.get('accession', '') if isinstance(donor, dict) else donor,
                               'public_donor_crossrefs': ';'.join(donor.get('dbxrefs', [])) if isinstance(donor, dict) else '',
                               'qtl_donor_overlap': 'unresolved; accessions alone do not establish non-overlap'})
    pd.DataFrame(donor_rows).to_csv(reports / 'cd4-rna-training-provenance.csv', index=False)
    summary = {'plan_sha256': hashlib.sha256((ROOT / 'config/cd4-rna-plan.json').read_bytes()).hexdigest(),
               'input_lock_sha256': run['input_lock_sha256'], 'query_manifest_sha256': hashlib.sha256(query_path.read_bytes()).hexdigest(),
               'summary_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'datasets': len(plan['datasets']), 'variants': len(plan['variants']), 'planned_variant_gene_pairs': len(pairs),
               'planned_dataset_variant_gene_tests': len(frame), 'complete_queries': sum(q['status'] == 'complete' for q in queries.values()),
               'preserved_nominal_rows': sum(q.get('rows', 0) for q in queries.values()), 'status_counts': dict(Counter(frame.status)),
               'bonferroni_threshold': threshold, 'bonferroni_passes': len(hits),
               'passing_unique_variant_gene_pairs': len(hits.drop_duplicates(['rsid','gene_id'])),
               'passing_genes': sorted(set(hits.gene_name)),
               'passing_primary_sign_agrees': int((hits.primary_direction_agrees == True).sum()),
               'passing_primary_sign_opposes': int((hits.primary_direction_agrees == False).sum()),
               'passing_mixed_model_direction': int((hits.model_consensus_direction == 'mixed_or_zero').sum()),
               'model_maximum_gene_statuses': dict(Counter(focus[focus.panel == 'model_maximum_gene'].status)),
               'no_new_alphagenome_requests': True,
               'limits': ['Both DICE and BLUEPRINT appear in the original PSC study; this is not wholly independent evidence',
                          'DICE cell subsets and conditions share donors; region/cohort are reporting units',
                          'ALL_FOLDS and unresolved exact training-donor overlap preclude a held-out accuracy claim',
                          'Missing annotation or filtered associations are not measured nulls',
                          'Normalized QTL coefficients are not model log changes, disease-risk directions or treatment effects']}
    (reports / 'cd4-rna-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    make_figures(frame, focus, plan)
    make_screened_pair_figure(frame, hits, plan)
    print(json.dumps(summary, indent=2))


def save_figure(fig, name):
    for suffix in ['png', 'svg']:
        path = ROOT / 'reports' / (name + '.' + suffix)
        fig.savefig(path, dpi=160, metadata={'Date': None} if suffix == 'svg' else {})
        if suffix == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')


def make_figures(frame, focus, plan):
    cache = ROOT / 'work/matplotlib'; cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(cache))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['svg.hashsalt'] = 'psc-cd4-rna-followup-v01'
    datasets = plan['datasets']
    labels = [d['cohort_label'] + ' ' + d['sample_group'].replace('CD4_T-cell_', '').replace('_', ' ') + (' *' if d['context_note'] else '') for d in datasets]
    fig, ax = plt.subplots(figsize=(11, 5.7))
    measured = [int(((frame.dataset_id == d['dataset_id']) & (frame.status == 'measured')).sum()) for d in datasets]
    ypos = np.arange(len(datasets))
    ax.barh(ypos, measured, color='#166b82', label='Measured exact variant–gene association')
    ax.barh(ypos, [228-n for n in measured], left=measured, color='#dce0e3', label='Unavailable')
    for y, count in zip(ypos, measured):
        ax.text(count + 3, y, str(count) + ' / 228', va='center', fontsize=9)
    ax.set_yticks(ypos, labels, fontsize=9); ax.invert_yaxis()
    ax.set_xlim(0, 228); ax.set_xlabel('Planned variant–gene pairs per dataset')
    ax.set_title('CD4 RNA evidence: complete accounting of all planned pairs', loc='left', fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False); ax.legend(loc='lower right', fontsize=8)
    fig.text(.01, .015, '* Treg naive/memory metadata labels conflict. DICE datasets share donors. Unavailable does not mean no effect.', fontsize=8)
    fig.tight_layout(rect=(0,.05,1,1)); save_figure(fig, 'cd4-rna-coverage'); plt.close(fig)
    variants = [v['rsid'] for v in plan['variants']]
    fig, axes = plt.subplots(2, 1, figsize=(13, 10.5), sharex=True)
    finite = focus.loc[focus.status == 'measured', 'beta_per_alt'].astype(float)
    limit = float(finite.abs().max()) if len(finite) else 1.0
    if not limit: limit = 1.0
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#e4e7e9')
    for ax, panel, title in zip(axes, ['model_maximum_gene','region_named_gene'], ['Highest model-ranked gene for every variant', 'Region-named gene for every variant']):
        rows = focus[focus.panel == panel]
        grid = np.full((len(variants), len(datasets)), np.nan)
        ylabels = []
        for y, variant in enumerate(variants):
            subset = rows[rows.rsid == variant]
            if subset.empty: raise ValueError('Focus panel lost a planned variant')
            first = subset.iloc[0]
            arrow = {'increase':'↑', 'decrease':'↓', 'mixed_or_zero':'mixed'}[first.model_consensus_direction]
            ylabels.append(variant + (' [P]' if first.role == 'anchor' else ' [C]') + ' · ' + first.gene_name + ' · model ' + arrow)
            for x, dataset in enumerate(datasets):
                selected = subset[subset.dataset_id == dataset['dataset_id']]
                if len(selected) != 1: raise ValueError('Missing/duplicate focus pair')
                row = selected.iloc[0]
                if row.status == 'measured':
                    grid[y, x] = float(row.beta_per_alt)
                    if row.bonferroni_2280_pass == True: ax.text(x,y,'*',ha='center',va='center',fontweight='bold')
        im = ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=-limit, vmax=limit, aspect='auto')
        ax.set_yticks(range(len(variants)), ylabels, fontsize=8)
        ax.set_title(title, loc='left', fontsize=12, fontweight='bold')
        ax.set_xticks(np.arange(-.5, len(datasets)), minor=True); ax.set_yticks(np.arange(-.5, len(variants)), minor=True)
        ax.grid(which='minor', color='white', linewidth=1); ax.tick_params(which='minor', bottom=False, left=False)
    axes[1].set_xticks(range(len(datasets)), labels, rotation=40, ha='right', fontsize=8)
    fig.subplots_adjust(left=.32, right=.90, bottom=.21, top=.93, hspace=.24)
    colorbar = fig.colorbar(im, cax=fig.add_axes([.92,.30,.018,.5])); colorbar.set_label('Measured beta per ALT allele (normalized RNA)', fontsize=9)
    fig.suptitle('Measured RNA directions for the complete variant set', x=.02, ha='left', fontsize=15, fontweight='bold')
    fig.text(.02,.04,'[P] prioritized; [C] comparison. Gray: unavailable. Cell *: passes the fixed 2,280-test correction.\nModel arrows require agreement of both CD4 tracks. DICE columns share donors; label * marks conflicting Treg metadata.', fontsize=9)
    save_figure(fig, 'cd4-rna-focus-genes'); plt.close(fig)


def make_screened_pair_figure(frame, hits, plan):
    """Show every dataset for every pair crossing the fixed threshold anywhere."""
    selected_pairs = hits[['rsid', 'gene_id', 'gene_name']].drop_duplicates().sort_values(['rsid','gene_id'])
    if selected_pairs.empty:
        return
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(selected_pairs), figsize=(6 * len(selected_pairs), 6.3), sharey=True, squeeze=False)
    datasets = plan['datasets']
    labels = [d['cohort_label'] + ' ' + d['sample_group'].replace('CD4_T-cell_', '').replace('_', ' ') + (' *' if d['context_note'] else '') for d in datasets]
    for ax, selected in zip(axes[0], selected_pairs.itertuples()):
        subset = frame[(frame.rsid == selected.rsid) & (frame.gene_id == selected.gene_id)]
        for y, dataset in enumerate(datasets):
            row = subset[subset.dataset_id == dataset['dataset_id']].iloc[0]
            if row.status != 'measured':
                ax.text(0,y,'unavailable',ha='left',va='center',fontsize=8,color='#888888'); continue
            color = '#166b82'
            ax.errorbar(row.beta_per_alt, y, xerr=1.96 * row.standard_error, fmt='o', color=color,
                        markerfacecolor=color if row.bonferroni_2280_pass == True else 'white', capsize=3, markersize=6)
        ax.axvline(0,color='#777777',linewidth=1)
        ax.set_title(selected.gene_name + '\n' + selected.rsid + ' · model tracks disagree', loc='left', fontsize=12, fontweight='bold')
        ax.set_yticks(range(len(datasets)), labels, fontsize=9); ax.invert_yaxis()
        ax.set_xlabel('Measured beta per ALT allele\n(normalized RNA; descriptive 95% Wald interval)', fontsize=9)
        ax.spines[['top','right']].set_visible(False); ax.grid(axis='x',alpha=.15)
    # Shared axes need one direction assignment, irrespective of panel count.
    axes[0,0].set_ylim(len(datasets)-.5, -.5)
    fig.suptitle('Two variant–gene pairs pass the fixed screen in at least one CD4 dataset', x=.02, ha='left', fontsize=13, fontweight='bold')
    fig.text(.02,.025,'Filled points: pass the 2,280-test correction. All datasets shown; intervals are not multiplicity-adjusted.\nDICE subsets share donors. * Treg cell-state labels conflict. These are molecular associations, not treatment effects.', fontsize=9)
    fig.tight_layout(rect=(0,.12,1,.92),w_pad=2)
    save_figure(fig,'cd4-rna-screened-pairs'); plt.close(fig)


if __name__ == '__main__':
    main()
