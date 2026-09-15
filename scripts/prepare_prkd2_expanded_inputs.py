"""Recover exact full alleles without altering the completed SNP-only analysis."""

import csv
import gzip
import json
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from prkd2_common import ROOT, ChainMap, complement, reciprocal_base, save_json, sha, verify_plan
from prkd2_expanded_common import load_references, map_edit, normalize_edit, qtl_identity
from prepare_coloc_inputs import save_frame


OUT = ROOT / 'data/derived/prkd2-expanded-inputs'
LD = ROOT / 'data/raw/prkd2-expanded-ld'


def literal_key(position, a, b):
    return int(position), tuple(sorted([a, b]))


def main():
    plan = verify_plan()
    extension = json.loads((ROOT / 'config/prkd2-expanded-plan.json').read_text())
    refs = load_references()
    ref37, ref38 = refs['GRCh37'], refs['GRCh38']
    forward = ChainMap.from_file(ROOT / 'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    OUT.mkdir(exist_ok=True)
    LD.mkdir(exist_ok=True)
    original_audit = pd.read_csv(ROOT / 'data/derived/prkd2-gwas-alignment.csv.gz', dtype=str, keep_default_na=False)
    original = pd.read_csv(ROOT / 'data/derived/prkd2-inputs/GWAS-PRKD2.csv.gz')
    original_by_row = {int(r['source_row']): r for r in original.to_dict('records')}
    source_keys = {literal_key(r['pos'], r['allele_0'], r['allele_1']) for r in original_audit.to_dict('records')}
    reference_matches = defaultdict(list)
    with gzip.open(ROOT / 'data/raw/prkd2-expanded-EUR-dosages.tsv.gz', 'rt') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        if len(header) != 508:
            raise ValueError('503 reference donors required')
        for row in reader:
            key = literal_key(row[1], row[3], row[4])
            if key not in source_keys:
                continue
            values = np.array(row[5:], dtype=np.int8)
            if not np.isin(values, [-1, 0, 1, 2]).all():
                raise ValueError('Invalid public dosage')
            called = values >= 0
            frequency = float(values[called].sum() / (2 * called.sum())) if called.any() else float('nan')
            reference_matches[key].append({'position': int(row[1]), 'source_id': row[2], 'ref': row[3], 'alt': row[4],
                                           'values': values, 'called': int(called.sum()), 'alt_frequency': frequency})
    gwas_rows, gwas_audit, normalization_checks = [], [], []
    vectors = {}
    for raw in original_audit.to_dict('records'):
        source_row = int(raw['source_row'])
        audit = {'source_row': source_row, 'source_id': raw['SNP'], 'position_grch37': int(raw['pos']),
                 'allele_0': raw['allele_0'], 'allele_1': raw['allele_1'], 'original_status': raw['status'],
                 'status': '', 'snp': ''}
        if raw['status'] not in ['eligible', 'not_biallelic_snv']:
            audit['status'] = 'original_exclusion_preserved'
            gwas_audit.append(audit)
            continue
        matches = reference_matches[literal_key(raw['pos'], raw['allele_0'], raw['allele_1'])]
        is_snv = len(raw['allele_0']) == len(raw['allele_1']) == 1
        try:
            if is_snv:
                reference_allele = ref37.get(int(raw['pos']))
                alleles = {raw['allele_0'], raw['allele_1']}
                if reference_allele not in alleles:
                    raise ValueError('GWAS SNV alleles do not include the GRCh37 reference')
                alternative = next(a for a in alleles if a != reference_allele)
                mapped = map_edit(ref37, ref38, forward, reverse, int(raw['pos']), reference_allele, alternative)
                row = original_by_row[source_row].copy()
                if row['snp'] != mapped['snp']:
                    raise ValueError('Original SNP identity changed')
                row['variant_class'] = 'SNV'
            else:
                if len(matches) != 1:
                    raise ValueError('Exact original-position full-allele VCF identity is absent or ambiguous')
                reference = matches[0]
                risk_frequency = reference['alt_frequency'] if raw['allele_1'] == reference['alt'] else 1 - reference['alt_frequency']
                difference = abs(risk_frequency - float(raw['freq_1_controls']))
                audit.update(reference_risk_frequency=risk_frequency, source_frequency_difference=difference,
                             reference_source_id=reference['source_id'])
                if not math.isfinite(difference) or difference > .15:
                    raise ValueError('Recovered-allele frequency disagreement exceeds 0.15')
                mapped = map_edit(ref37, ref38, forward, reverse, reference['position'], reference['ref'], reference['alt'])
                scope = plan['platform_counts'][raw['platform']]
                p, maf, odds = float(raw['p']), min(float(raw['freq_1_controls']), 1 - float(raw['freq_1_controls'])), float(raw['or'])
                if not (0 < p <= 1 and .01 - 1e-12 <= maf <= .5 and odds >= 1):
                    raise ValueError('Recovered allele fails original numerical or MAF rules')
                row = {'snp': mapped['snp'], 'position': mapped['normalized_grch38'][0], 'pvalue': p, 'MAF': maf,
                       'N': scope['N'], 'cases': scope['cases'], 'controls': scope['controls'],
                       'case_fraction': scope['cases'] / scope['N'],
                       'canonical_sign': mapped['dosage_sign'] * (1 if raw['allele_1'] == reference['alt'] else -1),
                       'source_id': raw['SNP'], 'source_or': raw['or'], 'source_se': raw['se'],
                       'platform': raw['platform'], 'source_row': source_row, 'variant_class': 'non_SNV'}
                for build, edit in [('GRCh37', (reference['position'], reference['ref'], reference['alt'])),
                                    ('GRCh38', mapped['normalized_grch38'])]:
                    normalization_checks.append({'id': f'gwas-{source_row}-{build}', 'build': build, 'position': edit[0], 'ref': edit[1], 'alt': edit[2],
                                                 'expected': normalize_edit(refs[build], *edit)})
            pos = row['position']
            region = plan['regions'][0]
            if not region['start_grch38_0based'] < pos <= region['end_grch38_exclusive']:
                raise ValueError('Normalized edit lies outside the fixed analysis window')
            audit.update(status='eligible', snp=row['snp'], reciprocal_bases=mapped['reciprocal_bases'],
                         strand=mapped['strand'], normalized_grch37=json.dumps(mapped['normalized_grch37']),
                         normalized_grch38=json.dumps(mapped['normalized_grch38']))
            if len(matches) == 1 and matches[0]['called'] == 503 and 0 < matches[0]['alt_frequency'] < 1:
                reference = matches[0]
                genotype_mapping = map_edit(ref37, ref38, forward, reverse, reference['position'], reference['ref'], reference['alt'])
                if genotype_mapping['snp'] != row['snp']:
                    raise ValueError('Reference dosage edit differs from association edit')
                values = reference['values'].astype(float)
                if genotype_mapping['dosage_sign'] == -1:
                    values = 2 - values
                vectors[source_row] = values
            gwas_rows.append(row)
        except (ValueError, KeyError) as exc:
            audit.update(status='unavailable', reason=str(exc))
        gwas_audit.append(audit)
    duplicated = {key for key, count in Counter(r['snp'] for r in gwas_rows).items() if count > 1}
    for audit in gwas_audit:
        if audit['status'] == 'eligible' and audit['snp'] in duplicated:
            audit.update(status='conflicting_or_duplicate_normalized_GWAS_edit')
    gwas_rows = [r for r in gwas_rows if r['snp'] not in duplicated]
    gwas = pd.DataFrame(gwas_rows).sort_values(['position', 'snp'])
    save_frame(gwas, OUT / 'GWAS-PRKD2.csv.gz')
    # Every original numeric row remains in the full evidence denominator.
    full_gwas = pd.read_csv(ROOT / 'data/derived/prkd2-inputs/GWAS-full-evidence.csv.gz')
    save_frame(full_gwas, OUT / 'GWAS-full-evidence.csv.gz')
    save_frame(pd.DataFrame(gwas_audit), ROOT / 'data/derived/prkd2-expanded-gwas-alignment.csv.gz')

    maps, dataset_audits = [], []
    for ds in ['QTD000021', 'QTD000504']:
        full = pd.read_csv(ROOT / ('data/derived/prkd2-inputs/' + ds + '-full-evidence.csv.gz'))
        groups = defaultdict(list)
        for row in full.to_dict('records'):
            item = {'dataset_id': ds, 'source_variant': row['source_variant'], 'status': '', 'snp': ''}
            try:
                identity = qtl_identity(row['source_variant'], ref38)
                item.update(identity, status='eligible')
                numeric = {**row, **identity, 'beta': row['beta'] * identity['canonical_beta_sign']}
                groups[identity['snp']].append(numeric)
                fields = row['source_variant'].split('_')
                if len(fields[2]) != 1 or len(fields[3]) != 1:
                    normalization_checks.append({'id': ds + '-' + row['source_variant'], 'build': 'GRCh38',
                          'position': int(fields[1]), 'ref': fields[2], 'alt': fields[3],
                          'expected': (identity['position'], identity['normalized_ref'], identity['normalized_alt'])})
            except ValueError as exc:
                item.update(status='unavailable', reason=str(exc))
            maps.append(item)
        kept = []
        equivalent = 0
        conflicting = set()
        for key, rows in groups.items():
            # Representation-only duplicates must agree on every available numerical field.
            signatures = {tuple(r[k] for k in ['beta', 'varbeta', 'pvalue', 'MAF', 'N']) for r in rows}
            if len(signatures) != 1:
                conflicting.add(key)
                continue
            kept.append({**rows[0], 'source_variant': rows[0]['source_variant'],
                         'source_variant_aliases': ';'.join(sorted(r['source_variant'] for r in rows))})
            equivalent += len(rows) - 1
        for item in maps:
            if item['dataset_id'] == ds and item['snp'] in conflicting:
                item.update(status='conflicting_equivalent_edit')
        qtl = pd.DataFrame(kept).sort_values(['position', 'snp'])
        save_frame(qtl, OUT / (ds + '-PRKD2.csv.gz'))
        save_frame(full, OUT / (ds + '-full-evidence.csv.gz'))
        dataset_audits.append({'dataset_id': ds, 'numeric_full_variants': len(full), 'eligible_distinct_edits': len(qtl),
            'equivalent_alias_rows_collapsed': equivalent, 'conflicting_edit_groups': len(conflicting),
            'gwas_overlap': int(qtl.snp.isin(gwas.snp).sum())})
    save_frame(pd.DataFrame(maps), ROOT / 'data/derived/prkd2-expanded-qtl-alignment.csv.gz')
    selected = gwas[gwas.source_row.isin(vectors)].copy()
    x = np.array([vectors[int(i)] for i in selected.source_row], dtype=float)
    x -= x.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(x, axis=1)
    if (norms <= 0).any():
        raise ValueError('Degenerate dosage vector')
    x /= norms[:, None]
    matrix = x @ x.T
    if not np.isfinite(matrix).all() or not np.allclose(np.diag(matrix), 1, atol=1e-12):
        raise ValueError('Invalid expanded signed LD')
    matrix_path = LD / 'GWAS-PRKD2.f64'
    matrix.astype('<f8').tofile(matrix_path)
    save_frame(selected, OUT / 'GWAS-PRKD2-ld.csv.gz')
    old_ld = pd.read_csv(ROOT / 'data/derived/prkd2-inputs/GWAS-PRKD2-ld.csv.gz')
    selected_index = {s: i for i, s in enumerate(selected.snp)}
    if not set(old_ld.snp) <= set(selected.snp):
        raise ValueError('A previously qualified SNP was lost from reference LD')
    indices = [selected_index[s] for s in old_ld.snp]
    previous = np.fromfile(ROOT / 'data/raw/prkd2-ld/GWAS-PRKD2.f64', dtype='<f8').reshape(len(old_ld), len(old_ld))
    ld_error = float(np.max(np.abs(previous - matrix[np.ix_(indices, indices)])))
    if ld_error > 1e-12:
        raise ValueError('Expanded source does not reproduce every original SNP LD entry')
    save_json(LD / 'GWAS-PRKD2.json', {'matrix_sha256': sha(matrix_path.read_bytes()), 'n_variants': len(selected),
        'n_reference_donors': 503, 'coding': 'Lexical-maximum allele for SNVs; normalized ALT for non-SNV edits',
        'in_sample': False, 'original_SNP_LD_entries_compared': len(old_ld) ** 2, 'original_SNP_LD_max_abs_difference': ld_error})
    save_json(ROOT / 'work/prkd2/expanded-normalization-checks.json', normalization_checks)
    summary = {'plan_sha256': sha((ROOT / 'config/prkd2-expanded-plan.json').read_bytes()), 'datasets': dataset_audits,
        'GWAS_eligible_variants': len(gwas), 'GWAS_variant_classes': gwas.variant_class.value_counts().to_dict(),
        'GWAS_status_counts': dict(Counter(r['status'] for r in gwas_audit)), 'GWAS_LD_variants': len(selected),
        'qtl_status_counts': dict(Counter(r['status'] for r in maps)), 'normalization_checks_queued': len(normalization_checks),
        'original_SNP_LD_max_abs_difference': ld_error,
        'inputs': {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in sorted(OUT.iterdir())}}
    save_json(ROOT / 'reports/prkd2-expanded-input-audit.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k != 'inputs'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
