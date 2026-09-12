"""Coverage and exact BCL2L11 RNA effects in the separately frozen OneK1K cohort.

No H0-H4 colocalisation posterior is calculated. RNA effects remain
observational associations in correlated donor/cell-group measurements.
"""

import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from analyze_coverage_repair import read_json, source_rows
from coloc_common import ROOT, save_json, sha
from coverage_repair_common import component_weights, exact_qtl_status
from prepare_coloc_inputs import canonical_key, save_frame, unique_alias_rows
from summarize_coloc import gwas_logbf


def eligible_qtl(row, gwas_keys, plan):
    """Apply the alternative plan before constructing the common SNP set."""
    if row['gene_id'].split('.')[0] != plan['gene_id'] or row['molecular_trait_id'].split('.')[0] != plan['gene_id']:
        raise ValueError('wrong_gene_or_trait')
    position = int(row['position'])
    if row['chromosome'] != '2' or not plan['start_grch38_0based'] < position <= plan['end_grch38_exclusive']:
        raise ValueError('outside_fixed_region')
    ref, alt = row['ref'], row['alt']
    if row['variant'] != f'chr2_{position}_{ref}_{alt}': raise ValueError('variant_fields_disagree')
    try: key, _, effect = canonical_key('chr2', position, ref, alt)
    except ValueError as error: raise ValueError('not_biallelic_SNV') from error
    if {ref, alt} in [{'A', 'T'}, {'C', 'G'}] and key not in gwas_keys:
        raise ValueError('nonshared_palindromic_frequency_unverified')
    beta, se, p, maf = (float(row[c]) for c in ['beta', 'se', 'pvalue', 'maf'])
    if not all(math.isfinite(x) for x in [beta, se, p, maf]) or se <= 0 or not 0 < p <= 1 or not .01 <= maf <= .5:
        raise ValueError('invalid_statistics_or_MAF_below_0.01')
    r2 = row['r2']
    if r2 not in ['', 'NA', 'nan', 'None', '.']:
        value = float(r2)
        if not math.isfinite(value) or not 0 <= value <= 1 or value < .8: raise ValueError('reported_R2_below_0.8_or_invalid')
    shrinkage = .15 ** 2 / (.15 ** 2 + se ** 2)
    logbf = .5 * (np.log1p(-shrinkage) + shrinkage * (beta / se) ** 2)
    return dict(snp=key, position=position, source_variant=row['variant'], source_beta=beta,
                beta_canonical=beta * (1 if alt == effect else -1), se=se, pvalue=p, maf=maf,
                logbf=float(logbf), source_an=row['an'], source_ma_samples=row['ma_samples'])


def main():
    plan = read_json('config/coverage-repair-alternative-plan.json')
    for name, digest in plan['frozen_file_sha256'].items():
        if sha((ROOT / name).read_bytes()) != digest: raise ValueError('Alternative input changed')
    queries = read_json('config/coverage-repair-parquet-queries.json')['queries']
    for query in queries:
        if 'result' not in query or sha((ROOT / query['result']['file']).read_bytes()) != query['result']['sha256']:
            raise ValueError('Incomplete or changed alternative query')
    if {q['id'] for q in queries} != {d['dataset_id'] for d in plan['datasets']}: raise ValueError('Incomplete cohort contexts')
    targets = read_json('data/derived/coverage-repair-exact-identities.json')['bcl2l11']
    gwas = pd.read_csv(ROOT / 'data/derived/coloc-platform-inputs/GWAS-BCL2L11.csv.gz')
    keys = set(gwas.snp)
    gwas['weight'] = component_weights(gwas_logbf(gwas))
    # The previous report independently checked these same logBFs against R coloc.
    if not np.isclose(gwas.weight.sum(), 1): raise ValueError('Incomplete full GWAS evidence')
    output, effects, audit, cs_tables = [], [], [], []
    for ds in plan['datasets']:
        identifier = ds['dataset_id']
        rows = source_rows(f'data/derived/coverage-repair-queries/{identifier}-alternative.tsv.gz')
        grouped = defaultdict(list)
        for row in rows: grouped[row['variant']].append(row)
        kept, exclusions = [], Counter()
        for variant, aliases in grouped.items():
            try:
                row, count, rsids = unique_alias_rows(aliases)
                result = eligible_qtl(row, keys, plan)
                kept.append(result); status = 'eligible'
            except ValueError as exc:
                status = str(exc); exclusions[status] += 1
            audit.append(dict(dataset_id=identifier, source_variant=variant, source_alias_rows=len(aliases), status=status))
        frame = pd.DataFrame(kept)
        if frame.snp.duplicated().any(): raise ValueError('Conflicting canonical QTL alleles')
        if frame.empty: raise ValueError('No eligible QTL SNPs')
        frame['weight'] = component_weights(frame.logbf)
        shared_qtl = frame.snp.isin(keys)
        shared_gwas = gwas.snp.isin(frame.snp)
        gwas_mass = float(gwas.loc[shared_gwas, 'weight'].sum())
        qtl_mass = float(frame.loc[shared_qtl, 'weight'].sum())
        full_cs = pd.read_parquet(ROOT / f'data/raw/repair-r8-{identifier}-cs.parquet')
        cs = full_cs[full_cs.gene_id.str.split('.').str[0].eq(plan['gene_id'])].copy()
        cs.insert(0, 'dataset_id', identifier); cs_tables.append(cs)
        missing = gwas.loc[~shared_gwas].nlargest(10, 'weight')
        source_frame = pd.DataFrame(rows)
        output.append(dict(dataset_id=identifier, sample_group=ds['sample_group'], source_metadata_sample_size=int(ds['sample_size']),
                           source_row_an=sorted(map(int, set(source_frame.an))), source_nominal_rows=len(rows),
                           source_unique_variants=len(grouped), eligible_qtl_snps=len(frame), excluded_unique_variants=dict(exclusions),
                           full_eligible_gwas_snps=len(gwas), shared_snps=int(shared_qtl.sum()),
                           gwas_bf_coverage=gwas_mass, qtl_bf_coverage=qtl_mass,
                           joint_coverage_gate=bool(shared_qtl.sum() >= 100 and min(gwas_mass, qtl_mass) >= .9),
                           published_target_gene_cs_rows=len(cs), minimum_nominal_p=float(frame.pvalue.min()),
                           missing_gwas_top10=missing[['snp', 'source_id', 'weight']].to_dict('records'),
                           status='coverage_and_direct_effects_only_no_colocalisation_posterior'))
        for target in targets:
            status, at_position, exact, measured = exact_qtl_status(rows, target['position_grch38'], target['ref'], target['alt'], plan['gene_id'])
            record = dict(dataset_id=identifier, sample_group=ds['sample_group'], source_id=target['source_id'],
                          position_grch38=target['position_grch38'], ref=target['ref'], alt=target['alt'], status=status)
            if measured:
                row, aliases, rsids = unique_alias_rows(measured)
                risk = target['gwas']['allele_1']
                record.update(beta_per_alt=float(row['beta']), se=float(row['se']), nominal_pvalue=float(row['pvalue']),
                              maf=float(row['maf']), minor_allele_carriers=int(row['ma_samples']), source_an=int(row['an']),
                              r2=row['r2'], source_metadata_sample_size=int(ds['sample_size']), risk_allele=risk,
                              beta_per_source_risk_allele=float(row['beta']) * (1 if risk == row['alt'] else -1),
                              eligible_for_coverage=row['variant'] in set(frame.source_variant), alias_rows=aliases)
            effects.append(record)
        save_frame(frame.sort_values(['position', 'snp']), ROOT / f'data/derived/coverage-repair-queries/{identifier}-eligible.csv.gz')
    save_json(ROOT / 'data/derived/coverage-repair-alternative-coverage.json', output)
    save_frame(pd.DataFrame(effects), ROOT / 'data/derived/coverage-repair-alternative-target-effects.csv')
    save_frame(pd.DataFrame(audit), ROOT / 'data/derived/coverage-repair-alternative-variant-audit.csv.gz')
    save_frame(pd.concat(cs_tables, ignore_index=True), ROOT / 'data/derived/coverage-repair-alternative-credible-sets.csv.gz')
    print(pd.DataFrame(output)[['dataset_id', 'shared_snps', 'gwas_bf_coverage', 'qtl_bf_coverage', 'joint_coverage_gate', 'published_target_gene_cs_rows']].to_string(index=False))
    print(pd.DataFrame(effects).to_string(index=False))


if __name__ == '__main__': main()
