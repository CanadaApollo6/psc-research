"""Align complete PRKD2 source regions and retain every exclusion explicitly."""

import csv
import gzip
import json
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from prkd2_common import ROOT, ChainMap, complement, reciprocal_base, save_json, sha, verify_plan
from prepare_coloc_inputs import canonical_key, save_frame, unique_alias_rows
from prepare_coloc_platform_inputs import platform_gwas
from summarize_cd4_rna import gene_id, measured_numbers


OUT = ROOT / 'data/derived/prkd2-inputs'
LD = ROOT / 'data/raw/prkd2-ld'


def queries(kind):
    record = json.loads((ROOT / ('config/prkd2-' + kind + '-queries.json')).read_text())
    assert record['plan_sha256'] == sha((ROOT / 'config/prkd2-plan.json').read_bytes())
    assert record['source_manifest_sha256'] == sha((ROOT / 'config/prkd2-sources.json').read_bytes())
    for row in record['queries']:
        data = (ROOT / row['file']).read_bytes()
        assert sha(data) == row['sha256'] and sha(gzip.decompress(data)) == row['uncompressed_sha256']
    return {r['id']: r for r in record['queries']}


def valid_gwas_evidence(row, counts):
    """Compute full-source evidence even for an unmatchable indel or identity."""
    try:
        p, f, odds = float(row['p']), float(row['freq_1_controls']), float(row['or'])
        if not all(math.isfinite(x) for x in [p, f, odds]) or not (0 < p <= 1 and 0 < f < 1 and odds > 0):
            raise ValueError('Invalid source statistics')
        scope = counts[row['platform']]
        return {'evidence_valid': True, 'pvalue': p, 'MAF': min(f, 1 - f),
                'N': scope['N'], 'case_fraction': scope['cases'] / scope['N']}
    except (ValueError, KeyError, TypeError) as exc:
        return {'evidence_valid': False, 'evidence_error': str(exc)}


def main():
    plan = verify_plan()
    region = plan['regions'][0]
    runs = {kind: queries(kind) for kind in ['gwas', 'qtl', 'reference']}
    if [len(runs[k]) for k in ['gwas', 'qtl', 'reference']] != [1, 2, 1]:
        raise ValueError('All planned regional sources are required')
    OUT.mkdir(exist_ok=True)
    LD.mkdir(exist_ok=True)
    forward = ChainMap.from_file(ROOT / 'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    genotypes = pd.read_csv(ROOT / runs['reference']['EUR-PRKD2']['file'], sep='\t', dtype={'chromosome': str})
    if genotypes.shape[1] != 508:
        raise ValueError('Exactly 503 selected EUR donor columns required')
    groups = defaultdict(list)
    ref_audit = []
    for i, row in genotypes.iloc[:, :5].iterrows():
        item = dict(row, row_index=i)
        hit = reciprocal_base(forward, reverse, 'chr' + row.chromosome, int(row.position_grch37) - 1)
        if hit is None:
            item['status'] = 'nonunique_or_nonreciprocal_mapping'
            ref_audit.append(item)
            continue
        chrom, pos0, strand = hit
        if not region['start_grch38_0based'] <= pos0 < region['end_grch38_exclusive']:
            item['status'] = 'outside_fixed_region'
            ref_audit.append(item)
            continue
        a, b = (complement(row.ref), complement(row.alt)) if strand == '-' else (row.ref, row.alt)
        key, _, effect = canonical_key(chrom, pos0 + 1, a, b)
        values = genotypes.iloc[i, 5:].to_numpy(dtype=np.int8)
        if not np.isin(values, [-1, 0, 1, 2]).all():
            raise ValueError('Invalid public reference dosage')
        called = values >= 0
        frequency = float(values[called].sum() / (2 * called.sum())) if called.any() else float('nan')
        if b != effect:
            frequency = 1 - frequency
        item.update(status='mapped', snp=key, position=pos0 + 1, effect_allele=effect,
                    effect_frequency=frequency, called_donors=int(called.sum()), flip_reference_dosage=b != effect)
        groups[key].append(item)
    reference = {}
    for key, rows in groups.items():
        if len(rows) == 1:
            reference[key] = rows[0]
        else:
            for row in rows:
                row['status'] = 'duplicate_reference_mapping'
        ref_audit.extend(rows)

    with gzip.open(ROOT / runs['gwas']['GWAS-PRKD2']['file'], 'rt') as stream:
        lines = stream.read().splitlines()
    fields = lines[0].lstrip('#').split()
    raw_gwas = [dict(zip(fields, line.split(), strict=True)) for line in lines[1:]]
    audited = []
    for index, row in enumerate(raw_gwas):
        item = platform_gwas(row, region, forward, reverse, reference, plan['platform_counts'])
        item.update(source_row=index, **valid_gwas_evidence(row, plan['platform_counts']))
        audited.append(item)
    counts = Counter(r['snp'] for r in audited if r['status'] == 'eligible')
    for row in audited:
        if row['status'] == 'eligible' and counts[row['snp']] > 1:
            row['status'] = 'duplicate_gwas_mapping'
    gwas = pd.DataFrame([{'snp': r['snp'], 'position': r['position'], 'pvalue': float(r['p']),
                         'MAF': r['control_maf'], 'N': r['N'], 'cases': r['cases'], 'controls': r['controls'],
                         'case_fraction': r['cases'] / r['N'], 'canonical_sign': r['canonical_sign'],
                         'source_id': r['SNP'], 'source_or': r['or'], 'source_se': r['se'],
                         'platform': r['platform'], 'source_row': r['source_row']}
                        for r in audited if r['status'] == 'eligible']).sort_values(['position', 'snp'])
    save_frame(gwas, OUT / 'GWAS-PRKD2.csv.gz')
    save_frame(pd.DataFrame(audited), ROOT / 'data/derived/prkd2-gwas-alignment.csv.gz')
    save_frame(pd.DataFrame(ref_audit), ROOT / 'data/derived/prkd2-reference-alignment.csv.gz')
    save_frame(pd.DataFrame([r for r in audited if r['evidence_valid']]), OUT / 'GWAS-full-evidence.csv.gz')

    qtl_audit = []
    dataset_audits = []
    for dataset in plan['datasets']:
        ds = dataset['dataset_id']
        with gzip.open(ROOT / runs['qtl'][ds + '-PRKD2']['file'], 'rt') as stream:
            source_rows = list(csv.DictReader(stream, delimiter='\t'))
        groups = defaultdict(list)
        for row in source_rows:
            groups[row['variant']].append(row)
        kept, full = [], []
        for variant, rows in groups.items():
            item = {'dataset_id': ds, 'source_variant': variant, 'source_alias_rows': len(rows), 'status': '', 'snp': ''}
            try:
                row, aliases, rsids = unique_alias_rows(rows)
                if gene_id(row['gene_id']) != region['gene_id'] or gene_id(row['molecular_trait_id']) != region['gene_id']:
                    raise ValueError('Molecular trait is not exact PRKD2')
                if variant != 'chr' + row['chromosome'] + '_' + row['position'] + '_' + row['ref'] + '_' + row['alt']:
                    raise ValueError('Variant columns disagree')
                if row['chromosome'] != '19' or not region['start_grch38_0based'] < int(row['position']) <= region['end_grch38_exclusive']:
                    raise ValueError('Variant is outside the frozen region')
                numbers = measured_numbers(row, dataset['sample_size'])
                maf = float(row['maf'])
                if not math.isfinite(maf) or not 0 < maf <= .5:
                    raise ValueError('Invalid source MAF')
                numeric = {'source_variant': variant, 'position': int(row['position']), 'source_rsids': rsids,
                           'beta': numbers['beta_per_alt'], 'varbeta': numbers['standard_error'] ** 2,
                           'pvalue': numbers['nominal_pvalue'], 'MAF': maf, 'N': numbers['sample_size_from_an'],
                           'source_beta': row['beta'], 'source_se': row['se']}
                # Full BF denominator includes well-formed non-SNV records too.
                full.append(numeric.copy())
                key, _, effect = canonical_key('chr19', row['position'], row['ref'], row['alt'])
                numeric.update(snp=key, beta=numeric['beta'] * (1 if row['alt'] == effect else -1))
                item.update(status='eligible', snp=key, source_rsids=rsids,
                            source_beta=row['beta'], source_se=row['se'], source_pvalue=row['pvalue'])
                kept.append(numeric)
            except (ValueError, TypeError) as exc:
                item.update(status='unavailable', reason=str(exc))
            qtl_audit.append(item)
        qtl = pd.DataFrame(kept).sort_values(['position', 'snp'])
        if qtl.snp.duplicated().any():
            raise ValueError('Duplicate canonical QTL identity')
        save_frame(qtl, OUT / (ds + '-PRKD2.csv.gz'))
        save_frame(pd.DataFrame(full).sort_values(['position', 'source_variant']), OUT / (ds + '-full-evidence.csv.gz'))
        dataset_audits.append({'dataset_id': ds, 'source_rows': len(source_rows), 'source_unique_variants': len(groups),
                               'numeric_full_variants': len(full), 'eligible_snvs': len(qtl),
                               'gwas_overlap': int(qtl.snp.isin(gwas.snp).sum()), 'sample_sizes': sorted(set(map(int, qtl.N)))})
        old_file = ROOT / ('data/derived/qtl-query-rows/' + ds + '-chr19_46718300_C_G.tsv')
        old = [r for r in csv.DictReader(old_file.open(), delimiter='\t') if gene_id(r['gene_id']) == region['gene_id']]
        new = [r for r in source_rows if r['variant'] == 'chr19_46718300_C_G']
        if sorted(json.dumps(r, sort_keys=True) for r in old) != sorted(json.dumps(r, sort_keys=True) for r in new):
            raise ValueError('Earlier rs313839 QTL row is not reproduced exactly')
    save_frame(pd.DataFrame(qtl_audit), ROOT / 'data/derived/prkd2-qtl-alignment.csv.gz')

    selected, vectors = [], []
    for row in gwas.to_dict('records'):
        entry = reference.get(row['snp'])
        if not entry or entry['called_donors'] != 503 or not 0 < entry['effect_frequency'] < 1:
            continue
        values = genotypes.iloc[entry['row_index'], 5:].to_numpy(dtype=np.float64)
        if entry['flip_reference_dosage']:
            values = 2 - values
        selected.append(row)
        vectors.append(values)
    selected = pd.DataFrame(selected)
    x = np.asarray(vectors, dtype=np.float64)
    x -= x.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(x, axis=1)
    if (norm == 0).any():
        raise ValueError('Constant reference genotype')
    x /= norm[:, None]
    ld = x @ x.T
    if not np.isfinite(ld).all() or not np.allclose(ld, ld.T, atol=1e-12) or not np.allclose(np.diag(ld), 1, atol=1e-12):
        raise ValueError('Invalid signed LD')
    matrix = LD / 'GWAS-PRKD2.f64'
    ld.astype('<f8').tofile(matrix)
    save_frame(selected, OUT / 'GWAS-PRKD2-ld.csv.gz')
    save_json(LD / 'GWAS-PRKD2.json', {'matrix_sha256': sha(matrix.read_bytes()),
              'snp_file_sha256': sha((OUT / 'GWAS-PRKD2-ld.csv.gz').read_bytes()), 'n_variants': len(selected),
              'n_reference_donors': 503, 'coding': 'Signed r for lexical-maximum allele dosage',
              'matrix_type': 'Gram matrix of centered unit-norm dosage vectors; PSD by construction',
              'reference_build': 'GRCh37 with unique reciprocal GRCh38 mapping', 'in_sample': False})

    sources = {s['id']: s for s in json.loads((ROOT / 'config/prkd2-sources.json').read_text())['sources']}
    cs = []
    for ds in ['QTD000021', 'QTD000504']:
        with gzip.open(ROOT / 'data/raw' / sources[ds + '-credible-sets']['file'], 'rt') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if gene_id(row['molecular_trait_id']) == region['gene_id']:
                    cs.append({'dataset_id': ds, **row})
    save_frame(pd.DataFrame(cs), ROOT / 'data/derived/prkd2-published-credible-sets.csv.gz')
    audit = {'datasets': dataset_audits, 'gwas_source_rows': len(raw_gwas), 'gwas_eligible_snvs': len(gwas),
             'gwas_ld_variants': len(selected), 'gwas_status_counts': dict(Counter(r['status'] for r in audited)),
             'qtl_status_counts': dict(Counter(r['status'] for r in qtl_audit)),
             'reference_status_counts': dict(Counter(r['status'] for r in ref_audit)),
             'published_credible_set_rows': len(cs), 'prior_rs313839_qtl_rows_exactly_reproduced': True,
             'plan_sha256': sha((ROOT / 'config/prkd2-plan.json').read_bytes()),
             'preparation_script_sha256': sha(__import__('pathlib').Path(__file__).read_bytes()),
             'inputs': {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in sorted(OUT.iterdir())}}
    save_json(ROOT / 'reports/prkd2-input-audit.json', audit)
    print(json.dumps({k: v for k, v in audit.items() if k not in ['inputs']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
