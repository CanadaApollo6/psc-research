"""Availability and full-evidence overlap audit; never emits an H4 posterior."""

import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from coloc_common import ROOT, ChainMap, gzip_bytes, save_json, sha
from prkd2_expanded_common import load_references, map_edit, qtl_identity
from verify_prkd2_statistics import binary_bf, quantitative_bf, weights

QUERY = ROOT / 'data/derived/prkd2-recovery-queries'
OUT = ROOT / 'data/derived'


def frame_file(name):
    return pd.read_csv(QUERY / name, sep='\t', dtype=str, keep_default_na=False)


def parse_dice(lines):
    rows = []
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 8:
            raise ValueError('Malformed original DICE VCF row')
        info = dict(x.split('=', 1) for x in fields[7].split(';'))
        if info['Gene'].split('.')[0] != 'ENSG00000105287':
            raise ValueError('Unexpected gene in retained DICE input')
        rows.append({'source_row': len(rows), 'chromosome': fields[0], 'position': int(fields[1]),
                     'ref': fields[3], 'alt': fields[4], 'rsid': fields[2], 'source_gene_id': info['Gene'],
                     'source_pvalue': info['Pvalue'], 'source_beta': info['Beta'], 'source_se': '',
                     'source_statistic': info['Statistic'], 'source_fdr': info['FDR'], 'source_af': '',
                     'se_status': 'not supplied; rounded Statistic is not substituted for SE'})
    return rows


def parse_blueprint(lines):
    rows = []
    for line in lines:
        fields = line.split()
        if len(fields) != 9:
            raise ValueError('Malformed original BLUEPRINT row')
        locus, ref, alt = fields[0].split('_')
        chrom, position = locus.split(':')
        if fields[2].split('.')[0] != 'ENSG00000105287':
            raise ValueError('Unexpected gene in retained BLUEPRINT input')
        rows.append({'source_row': len(rows), 'chromosome': chrom, 'position': int(position),
                     'ref': ref, 'alt': alt, 'rsid': fields[1], 'source_gene_id': fields[2],
                     'source_pvalue': fields[3], 'source_beta': fields[4], 'source_se': fields[8],
                     'source_statistic': '', 'source_fdr': fields[6], 'source_af': fields[7],
                     'source_bonferroni_pvalue': fields[5], 'se_status': 'explicit source SE'})
    return rows


def save_csv(frame, path, compressed=False):
    plain = frame.to_csv(index=False, lineterminator='\n').encode()
    data = gzip_bytes(plain) if compressed else plain
    path.write_bytes(data)
    return {'file': str(path.relative_to(ROOT)), 'rows': len(frame), 'sha256': sha(data)}


def original_identities(rows, refs, forward, reverse, plan):
    result = []
    for row in rows:
        item = dict(row)
        hits = forward.map_base('chr19', row['position'] - 1)
        item['mapped_anchor_grch38'] = hits[0][1] + 1 if len(hits) == 1 and hits[0][0] == 'chr19' else ''
        item['in_fixed_region'] = (isinstance(item['mapped_anchor_grch38'], int) and
                                   plan['start_grch38_0based'] < item['mapped_anchor_grch38'] <= plan['end_grch38_exclusive'])
        item['snp'] = ''
        item['canonical_beta_sign'] = ''
        try:
            mapping = map_edit(refs['GRCh37'], refs['GRCh38'], forward, reverse, row['position'], row['ref'], row['alt'])
            item.update({'snp': mapping['snp'], 'canonical_beta_sign': mapping['dosage_sign'],
                         'normalized_grch38': ':'.join(map(str, mapping['normalized_grch38'])), 'identity_status': 'verified-full-edit'})
        except ValueError as exc:
            item['identity_status'] = str(exc)
        result.append(item)
    frame = pd.DataFrame(result)
    for key, group in frame[frame.snp.ne('')].groupby('snp'):
        if len(group) < 2:
            continue
        values = group[['source_beta', 'source_se', 'source_pvalue', 'canonical_beta_sign']].drop_duplicates()
        if len(values) > 1:
            frame.loc[group.index, 'identity_status'] = 'conflicting equivalent edits'
            frame.loc[group.index, 'snp'] = ''
    return frame


def overlap(frame, gwas, full_gwas):
    distinct = frame.drop_duplicates(subset=[c for c in frame.columns if c != 'source_row'])
    eligible = frame[frame.snp.ne('') & frame.in_fixed_region]
    common = set(gwas.snp) & set(eligible.snp)
    source_rows = gwas[gwas.snp.isin(common)].source_row
    gm = float(weights(binary_bf(full_gwas))[full_gwas.source_row.isin(source_rows)].sum())
    result = {'full_source_gene_rows': len(frame), 'source_rows_in_fixed_region': int(frame.in_fixed_region.sum()),
              'distinct_source_association_rows': len(distinct), 'exact_duplicate_source_rows': len(frame) - len(distinct),
              'verified_edit_rows': int(frame.snp.ne('').sum()), 'shared_distinct_edits': len(common),
              'full_PSC_marginal_weight_available': gm, 'RNA_marginal_weight_available': None,
              'RNA_weight_denominator': 'Every distinct original PRKD2 association, including unmapped edits and rows outside the fixed window; collapse only exact duplicates in all fields except source row index',
              'quantitative_prior_SD_in_source_units': .15,
              'coverage_gate': None, 'multi_signal_model_qualified': False}
    try:
        beta = pd.to_numeric(distinct.source_beta, errors='raise').to_numpy(float)
        se = pd.to_numeric(distinct.source_se, errors='raise').to_numpy(float)
        if not np.isfinite(beta).all() or not np.isfinite(se).all() or (se <= 0).any():
            raise ValueError('Missing or invalid original beta/SE in full denominator')
        bf = quantitative_bf(pd.DataFrame({'beta': beta, 'varbeta': se ** 2}), 1)
        qm = float(weights(bf)[distinct.snp.isin(common) & distinct.in_fixed_region].sum())
        result['RNA_marginal_weight_available'] = qm
        result['coverage_gate'] = len(common) >= 100 and gm >= .9 and qm >= .9
        result['RNA_weight_status'] = 'diagnostic marginal overlap only; source units and original mixed-model analysis retained'
    except (ValueError, TypeError):
        result['RNA_weight_status'] = 'unavailable: no complete explicit finite beta/SE vector; no reconstructed or imputed SE'
    result['model_gate_reason'] = 'Original and reprocessed effects cannot be spliced together. No full original-source component vectors or matching study covariance have been recovered.'
    return result


def target_status(frame, target, gene='ENSG00000105287'):
    if frame is None:
        return {'status': 'not-accessible', 'matching_rows': None}
    if not len(frame):
        return {'status': 'no-exact-target-gene-row', 'matching_rows': 0}
    keep = (frame.position.astype(int).eq(target['position']) & frame.ref.eq(target['ref']) & frame.alt.eq(target['alt']))
    if 'chromosome' in frame:
        keep &= frame.chromosome.astype(str).str.removeprefix('chr').eq(target['chromosome'])
    if 'gene_id' in frame:
        keep &= frame.gene_id.str.split('.').str[0].eq(gene)
    selected = frame[keep]
    if selected.empty:
        return {'status': 'no-exact-target-gene-row', 'matching_rows': 0}
    measurements = selected[['beta', 'se', 'pvalue']].drop_duplicates()
    return {'status': 'measured' if len(measurements) == 1 else 'multiple-distinct-measurements', 'matching_rows': len(selected),
            'beta': ';'.join(measurements.beta), 'se': ';'.join(measurements.se), 'pvalue': ';'.join(measurements.pvalue)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--available-only', action='store_true', help='Exploratory preview; final report requires every planned query outcome')
    args = parser.parse_args()
    plan_path = ROOT / 'config/prkd2-recovery-plan.json'
    followup_path = ROOT / 'config/prkd2-recovery-followup-plan.json'
    plan = json.loads(plan_path.read_text())
    followup = json.loads(followup_path.read_text())
    refs = load_references()
    forward = ChainMap.from_file(ROOT / 'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    targets = followup['targets']
    registry = []
    for target in targets:
        mapped = map_edit(refs['GRCh38'], refs['GRCh37'], reverse, forward, target['position'], target['ref'], target['alt'])
        identity = qtl_identity(f"chr19_{target['position']}_{target['ref']}_{target['alt']}", refs['GRCh38'])
        registry.append({**target, 'canonical_edit_grch38': identity['snp'],
                         'position_grch37': mapped['normalized_grch38'][0], 'ref_grch37': mapped['normalized_grch38'][1],
                         'alt_grch37': mapped['normalized_grch38'][2]})
    outputs = [save_csv(pd.DataFrame(registry), OUT / 'prkd2-recovery-targets.csv')]
    gwas = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-PRKD2.csv.gz')
    full_gwas = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz')
    primary_mask = full_gwas.SNP.eq(plan['primary_missing_rsid'])
    if primary_mask.sum() != 1:
        raise ValueError('Primary missing GWAS variant is not unique')
    primary_weight = float(weights(binary_bf(full_gwas))[primary_mask].sum())
    prior_weight = next(float(r['full_trait_bf_weight']) for r in plan['prior_missing_evidence'] if r['trait'] == 'GWAS')
    if not np.isclose(primary_weight, prior_weight, rtol=1e-12, atol=1e-14):
        raise ValueError('Completed primary missing-variant weight changed')
    originals, original_targets = {}, []
    for name, filename, parser_fn in [('DICE-original', 'original-DICE-PRKD2.vcf.gz', parse_dice),
                                      ('BLUEPRINT-original', 'original-BLUEPRINT-PRKD2.txt.gz', parse_blueprint)]:
        if not (QUERY / filename).exists() and args.available_only:
            continue
        with gzip.open(QUERY / filename, 'rt', newline='') as stream:
            rows = parser_fn(stream)
        frame = original_identities(rows, refs, forward, reverse, plan)
        outputs.append(save_csv(frame, OUT / ('prkd2-recovery-' + name + '-alleles.csv.gz'), True))
        originals[name] = overlap(frame, gwas, full_gwas)
        originals[name]['source_coordinate_min'] = int(frame.position.min())
        originals[name]['source_coordinate_max'] = int(frame.position.max())
        for target in registry:
            matched = frame[frame.snp.eq(target['canonical_edit_grch38'])]
            for _, row in matched.iterrows():
                original_targets.append({'dataset_id': name, 'target_rsid': target['rsid'], 'target_edit_grch38': target['canonical_edit_grch38'],
                                         'role': target['role'], **{k: row[k] for k in ['source_row', 'rsid', 'position', 'ref', 'alt', 'source_gene_id', 'source_beta', 'source_se', 'source_pvalue', 'source_statistic', 'source_af', 'canonical_beta_sign', 'identity_status']}})
    outputs.append(save_csv(pd.DataFrame(original_targets), OUT / 'prkd2-recovery-original-target-associations.csv'))
    all_queries = {}
    for release in ['r7', 'r8_beta']:
        p = ROOT / ('config/prkd2-recovery-' + release + '-queries.json')
        if p.exists():
            all_queries.update({q['id']: q for q in json.loads(p.read_text())['queries']})
    cohort_rows, target_rows = [], []
    for ds in plan['cohorts']:
        dsid = ds['dataset_id']
        q = all_queries.get(dsid)
        frame = None
        if ds['prior_analysis']:
            frame = pd.read_csv(ROOT / ('data/derived/prkd2-query-rows/' + dsid + '-PRKD2.tsv.gz'), sep='\t', dtype=str, keep_default_na=False)
            status = 'prior-complete-region; missing-required-PSC-variant'
            scope = 'complete prior PRKD2 regional rows'
        elif q and q.get('complete'):
            if ds['release'] == 'r7':
                frame = frame_file(dsid + '-points-all-genes.tsv.gz')
                control = QUERY / (dsid + '-controls.tsv.gz')
                if control.exists():
                    frame = pd.concat([frame, frame_file(control.name)], ignore_index=True)
                elif not args.available_only:
                    raise ValueError('Required true rs313839 control query missing: ' + dsid)
                scope = 'five exact missing-evidence points plus separately verified rs313839 control; indel aliases outside queried positions not excluded'
            else:
                frame = frame_file(dsid + '-PRKD2-region.tsv.gz')
                scope = 'complete PRKD2 gene rows in the fixed region'
            status = q['status']
        elif q and q.get('status') == 'access-or-format-failure':
            status = q['status']
            scope = q['error']
        elif args.available_only:
            continue
        else:
            raise ValueError('Planned dataset has no completed outcome: ' + dsid)
        primary = target_status(frame, targets[0])
        cohort_rows.append({**ds, 'query_status': status, 'query_scope': scope,
                            'primary_variant_status': primary['status'], 'primary_variant_gene_rows': primary['matching_rows'],
                            'maximum_PSC_marginal_weight_if_primary_absent': (1 - primary_weight) if primary['matching_rows'] == 0 else None})
        for target in registry:
            target_rows.append({'dataset_id': dsid, 'study_label': ds['study_label'], 'sample_group': ds['sample_group'],
                                'target_edit_grch38': target['canonical_edit_grch38'], 'target_rsid': target['rsid'], 'role': target['role'],
                                **target_status(frame, target)})
    outputs.append(save_csv(pd.DataFrame(cohort_rows), OUT / 'prkd2-recovery-cohort-screen.csv'))
    outputs.append(save_csv(pd.DataFrame(target_rows), OUT / 'prkd2-recovery-target-screen.csv'))
    reference = frame_file('1000G-30x-20201028-reference.tsv.gz')
    refrows = []
    for target in registry:
        matches = reference[(reference.POS.astype(int) == target['position']) & reference.REF.eq(target['ref']) & reference.ALT.eq(target['alt'])]
        refrows.append({**target, 'reference_exact_site_rows': len(matches),
                        'reference_status': 'present' if len(matches) else 'absent-at-exact-site',
                        'filter_values': ';'.join(matches.FILTER), 'site_ids': ';'.join(matches.ID)})
    outputs.append(save_csv(pd.DataFrame(refrows), OUT / 'prkd2-recovery-reference-panel-audit.csv'))
    result = {'analysis_status': 'preview' if args.available_only else 'complete', 'new_H4_posteriors': 0,
              'plan_sha256': sha(plan_path.read_bytes()), 'followup_plan_sha256': sha(followup_path.read_bytes()),
              'script_sha256': sha(Path(__file__).read_bytes()), 'cohorts_planned': len(plan['cohorts']),
              'cohort_outcomes': len(cohort_rows), 'outcome_counts': dict(Counter(r['query_status'] for r in cohort_rows)),
              'original_study_overlap': originals, 'outputs': outputs,
              'limitations': ['Metadata-listed but inaccessible datasets are not biological nulls.',
                              'Reannotations share source cohorts; study labels are not independent replications.',
                              'Original study estimates cannot repair missing entries in reprocessed component vectors.',
                              'The high-coverage reference check does not prove which internal Catalogue filter removed a variant.']}
    save_json(ROOT / 'reports/prkd2-recovery-analysis.json', result)
    print(json.dumps({k: result[k] for k in ['analysis_status', 'cohort_outcomes', 'outcome_counts', 'original_study_overlap']}, indent=2))


if __name__ == '__main__':
    main()
