"""Summarize public model evidence without fitting a new association or posterior."""
import csv
import gzip
import json
from pathlib import Path

import pandas as pd

from fetch_ibdverse_model_audit import ROOT, RAW, freeze, save, sha

OUT = ROOT / 'data/derived/ibdverse-model-audit'


def write_csv(name, rows):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows: raise ValueError('Unexpected empty audit table')
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return {'file': str(path.relative_to(ROOT)), 'rows': len(rows), 'sha256': sha(path)}


def rows(path, sep='\t'):
    with Path(path).open() as f: return list(csv.DictReader(f, delimiter=sep))


def sample_counts():
    samples = pd.read_csv(RAW / 'E-MTAB-16986-restored-sample_metadata.tsv.data', sep='\t', dtype=str)
    mapping = pd.read_csv(RAW / 'E-MTAB-16986-restored-sample_individual_mapping.tsv.data', sep='\t', dtype=str)
    individuals = pd.read_csv(RAW / 'E-MTAB-16986-restored-individual_metadata.tsv.data', sep='\t', dtype=str)
    if mapping.sanger_sample_id.duplicated().any() or individuals.individual_id.duplicated().any():
        raise ValueError('Ambiguous published sample/donor mapping')
    joined = samples.merge(mapping, on='sanger_sample_id', how='outer', validate='one_to_one', indicator=True)
    if not joined['_merge'].eq('both').all(): raise ValueError('Sample metadata join incomplete')
    joined = joined.drop(columns='_merge').merge(individuals, on='individual_id', how='left', validate='many_to_one', indicator=True)
    if not joined['_merge'].eq('both').all(): raise ValueError('Individual metadata join incomplete')
    return [{'biopsy_type': tissue, 'samples': len(part), 'distinct_public_individual_labels': part.individual_id.nunique(),
             'CD_individuals': part.loc[part.disease_status.eq('CD'), 'individual_id'].nunique(),
             'healthy_individuals': part.loc[part.disease_status.eq('Healthy'), 'individual_id'].nunique(),
             'qualification': 'Collection metadata; not final eQTL model N'}
            for tissue, part in joined.groupby('biopsy_type', sort=True)]


def main():
    frozen = freeze()
    plan = json.loads((ROOT / 'config/ibdverse-model-audit-plan.json').read_text())
    sources = json.loads((ROOT / 'config/ibdverse-model-audit-sources.json').read_text())['sources']
    for source in sources:
        for file_key, hash_key in [('file', 'sha256'), ('receipt', 'receipt_sha256')]:
            if source.get(file_key) and sha(ROOT / source[file_key]) != source[hash_key]:
                raise ValueError('Source receipt or bytes changed: ' + source['id'])
    outputs = []
    file_inventory = []
    for accession in ['S-BSST2922', 'E-MTAB-16986']:
        table = json.loads((RAW / ('biostudies-' + accession + '-file-table.json')).read_text())
        if len(table['data']) != table['recordsTotal'] or table['recordsFiltered'] != table['recordsTotal']:
            raise ValueError('Incomplete public file inventory')
        file_inventory += [{'accession': accession, 'filename': row['path'], 'bytes': row['Size'],
                            'source_description': row.get('Description', '')} for row in table['data']]
    outputs.append(write_csv('public-file-inventory.csv', file_inventory))
    outputs.append(write_csv('collection-sample-counts.csv', sample_counts()))
    archive = rows(ROOT / 'data/derived/prkd2-recovery-ibdverse-archive-inventory.csv', ',')
    if sha(ROOT / 'data/derived/prkd2-recovery-ibdverse-archive-inventory.csv') != frozen['prior_tracked_sha256']['data/derived/prkd2-recovery-ibdverse-archive-inventory.csv']:
        raise ValueError('Prior ZIP inventory changed')
    basenames = pd.Series([Path(r['filename']).name for r in archive]).value_counts().sort_index()
    outputs.append(write_csv('base-eqtl-member-types.csv', [{'basename': name, 'members': int(n)} for name, n in basenames.items()]))

    ld_source = RAW / 'IBDverse-LD-clumped-eQTLs.data'
    ld = pd.read_csv(ld_source, sep='\t', dtype=str, keep_default_na=False)
    repo_ld = gzip.decompress((RAW / 'ibdverse-code-data__clumped_all.txt.gz.data').read_bytes())
    if repo_ld != ld_source.read_bytes(): raise ValueError('Released and repository clumping exports differ')
    gene_ld = ld[ld.phenotype_id.eq(plan['target_gene'])]
    outputs.append(write_csv('PRKD2-published-clump-membership.csv', gene_ld.to_dict('records')))
    clumps = gene_ld.set_index('variant_id').qtl_clump_index.to_dict()
    h5ad = json.loads((ROOT / 'reports/ibdverse-model-audit-h5ad-metadata.json').read_text())
    metadata_counts = {r['context']: r for r in h5ad['counts']}
    summaries, comparisons = [], []
    for context in plan['contexts']:
        name = context['source_context']; label = context['label']
        gene_rows = rows(OUT / (name + '-Cis_eqtls_qval-PRKD2.tsv'))
        independent = rows(OUT / (name + '-Cis_eqtls_independent-PRKD2.tsv'))
        if len(gene_rows) != 1: raise ValueError('Expected one complete source gene-summary row')
        gene = gene_rows[0]
        if any(r['phenotype_id'] != plan['target_gene'] for r in [gene, *independent]):
            raise ValueError('Off-target gene in selected rows')
        number = name.split('Myeloid_')[1].split('_')[0]
        original = pd.read_csv(ROOT / f'data/derived/prkd2-recovery-queries/IBDverse-original-Myeloid-{number}-PRKD2.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
        match = original[original.variant_id.eq(gene['variant_id'])]
        if len(match) != 1: raise ValueError('Gene-summary lead not uniquely present in prior nominal vector')
        old = match.iloc[0]
        n_rna = metadata_counts[name.removeprefix('dMean__').removesuffix('_all')]['donor_labels_with_at_least_five_cells']
        inferred_pcs = n_rna - 2 - 5 - context['previous_inferred_residual_df']
        summaries.append({'context': name, 'cell_label': label, 'gene_id': gene['phenotype_id'],
                          'source_tested_variants': gene['num_var'], 'lead_variant': gene['variant_id'],
                          'beta_ALT': gene['slope'], 'SE': gene['slope_se'], 'nominal_P_gene_summary': gene['pval_nominal'],
                          'permutation_P': gene['pval_perm'], 'beta_approximation_P': gene['pval_beta'],
                          'gene_FDR_qval': gene['qval'], 'gene_nominal_threshold': gene['pval_nominal_threshold'],
                          'selected_independent_rows': len(independent), 'selected_ranks': '|'.join(r['rank'] for r in independent),
                          'source_beta_fit_true_df': gene['true_df'],
                          'previous_inferred_regression_residual_df': context['previous_inferred_residual_df'],
                          'RNA_eligible_public_donor_labels': n_rna,
                          'expression_PCs_if_RNA_count_equals_model_N_and_five_genotype_PCs_only': inferred_pcs,
                          'actual_executed_model_N_observed': False, 'actual_selected_PCs_observed': False,
                          'published_global_clump': clumps.get(gene['variant_id'], ''),
                          'qualification': 'One selected hit is not proof of one biological causal effect; PC number is conditional inference only'})
        for field in ['slope', 'slope_se', 'af', 'ma_count', 'ma_samples', 'pval_nominal']:
            comparisons.append({'context': name, 'variant_id': gene['variant_id'], 'field': field,
                                'gene_summary_original': gene[field], 'nominal_vector_original': str(old[field]),
                                'absolute_difference': abs(float(gene[field]) - float(old[field])),
                                'relative_difference': abs(float(gene[field]) - float(old[field])) / max(abs(float(old[field])), 1e-300)})
    outputs.append(write_csv('PRKD2-source-model-summary.csv', summaries))
    outputs.append(write_csv('lead-source-field-comparison.csv', comparisons))

    browser = rows(RAW / 'IBDverse-browser-filtered-export.data', ',')
    outputs.append(write_csv('browser-filtered-export-as-published.csv', browser))
    if len(browser) != 1 or browser[0]['annotation'] != 'Myeloid_0_blood':
        raise ValueError('Browser snapshot selection differs')
    classical = summaries[0]
    portal_comparison = {'portal_filters': {'gene': 'PRKD2', 'tissue': 'Blood', 'major_population': 'Myeloid', 'annotation_level': 'Cell type', 'maximum_displayed_qvalue': 0.05},
        'portal_selected_rows': len(browser), 'portal_variant': browser[0]['variant_id'],
        'portal_qval_as_published': browser[0]['qval'], 'archive_gene_FDR_qval': classical['gene_FDR_qval'],
        'archive_permutation_P': classical['permutation_P'],
        'matches_permutation_P_at_three_significant_digits': float(browser[0]['qval']) == float(format(float(classical['permutation_P']), '.3g')),
        'matches_gene_FDR_at_three_significant_digits': float(browser[0]['qval']) == float(format(float(classical['gene_FDR_qval']), '.3g')),
        'interpretation': 'Preserve the public browser label, but use the complete archived qval field as the gene FDR. Omission from a globally clumped browser is not an absence of an original intermediate-monocyte eQTL.'}

    repositories = [
        {'repository': 'andersonlab/IBDVerse-sc-eQTL-code', 'ref_interpreted': 'current main; publication points to repository', 'commit': json.loads((RAW / 'ibdverse-code-current-retry1.data').read_text())['sha'], 'executed_run_commit_confirmed': False},
        {'repository': 'wtsi-hgi/QTLight', 'ref_interpreted': 'v1.4 branch; README reports v1.4.0 but no exact v1.4.0 release/tag found', 'commit': json.loads((RAW / 'QTLight-v1.4-commit.data').read_text())['sha'], 'executed_run_commit_confirmed': False},
        {'repository': 'broadinstitute/tensorqtl', 'ref_interpreted': 'v1.0.9 tag; version explicitly reported in Reporting Summary', 'commit': json.loads((RAW / 'tensorqtl-v1.0.9-tag.data').read_text())['object']['sha'], 'executed_run_commit_confirmed': False}]
    for name in ['snakemake_colocalisation_sceqtl', 'ld_clump_tqtl']:
        repositories.append({'repository': 'andersonlab/' + name, 'ref_interpreted': 'current HEAD of explicitly linked repository',
                             'commit': json.loads((RAW / (name + '-current.data')).read_text())['sha'], 'executed_run_commit_confirmed': False})
    outputs.append(write_csv('code-repository-versions.csv', repositories))
    missing = [
        {'input': 'Actual model sample/covariate records for each fixed context', 'specific_files_or_fields': 'Covariates.tsv header plus row labels; final genotype_phenotype_mapping.tsv; final sample-QC inclusion counts; duplicate-ID handling', 'public_status': 'RNA eligibility 94/93 observed, genotype model N not directly documented', 'needed_for': 'Confirm participant N, covariate count, row/rank conventions and allele-dosage sample matching'},
        {'input': 'Selected expression PC and execution version records', 'specific_files_or_fields': 'optim_pcs.txt; optimise_nPCs-FDR0pt05.txt; Nextflow trace/config/command; container digest; exact QTLight and TensorQTL revisions', 'public_status': '18/20 PCs are consistent conditional inferences; actual selected values and run manifests missing', 'needed_for': 'Verify actual model instead of transferring defaults or an unlabeled branch snapshot'},
        {'input': 'Original RNA covariance or complete per-signal vectors', 'specific_files_or_fields': 'Signed dosage LD/covariance with exact ordered GRCh38 variant/REF/ALT manifest and donor/covariate provenance; or complete per-signal posterior/BF vectors, conditioning variants and model parameters', 'public_status': 'Nominal rows and one selected independent hit per context recovered; global clump membership has no signed correlations or complete component weights', 'needed_for': 'Separate overlapping RNA signals and test shared PSC cause under qualified assumptions'},
        {'input': 'PSC model/covariance provenance', 'specific_files_or_fields': 'Per-variant platform sample counts and cross-platform overlap/covariance or adequate study-matched LD; documented coefficient/SE scale; original exact allele/QC metadata', 'public_status': 'Previously documented platform Ns 11386, 3504 and 14890 do not reconstruct per-variant covariance; RNA-side recovery does not resolve this', 'needed_for': 'Qualify a new shared-signal model across the full fixed locus'}]
    outputs.append(write_csv('remaining-required-inputs.csv', missing))
    auxiliary_bytes = sum(r['bytes'] for filename in ['ibdverse-model-audit-member-queries.json', 'ibdverse-model-audit-coloc-queries.json']
                          for q in json.loads((ROOT / 'config' / filename).read_text())['queries'] for r in q['ranges'])
    total = sum(s.get('bytes', 0) or 0 for s in sources) + auxiliary_bytes
    if total > plan['retrieval']['initial_total_new_download_budget_bytes']: raise ValueError('Stage byte budget exceeded')
    record = {'audit_plan_sha256': sha(ROOT / 'config/ibdverse-model-audit-plan.json'), 'baseline_commit': frozen['baseline_commit'],
              'source_receipts': len(sources), 'successful_source_receipts': sum(bool(s.get('ok')) for s in sources),
              'retained_successful_source_and_range_bytes': total, 'public_file_inventory_rows': len(file_inventory),
              'base_eQTL_archive_members': len(archive), 'selected_source_summary': summaries,
              'RNA_observation_metadata': {k: v for k, v in h5ad.items() if k != 'schema'},
              'global_clump_rows': len(ld), 'PRKD2_global_clump_rows': len(gene_ld),
              'published_clump_file_equals_repository_decoded_file': True,
              'browser_comparison': portal_comparison, 'outputs': outputs,
              'ready_for_new_shared_cause_posterior': False, 'new_H4_posteriors': 0,
              'decision': 'Public audit complete. Exact executed model records and adequate RNA/PSC covariance or complete RNA components remain unqualified.',
              'limitations': ['Source branch snapshots and declared versions are not executed-run manifests.', 'RNA-eligible donor counts and conditional PC arithmetic are not direct genotype-model metadata.', 'Conditional lead tables and global clumps are not complete causal components.', 'The base-coloc ZIP has three disease result streams; only headers were read, not every association. It is an IBD/CD/UC release, not PSC.', 'This bounded audit does not prove that no suitable data exist anywhere.']}
    save(ROOT / 'reports/ibdverse-model-audit.json', record)
    print(json.dumps({'ready': False, 'new_H4_posteriors': 0, 'outputs': len(outputs), 'retained_source_bytes': total,
                      'RNA_donor_counts': [r['RNA_eligible_public_donor_labels'] for r in summaries],
                      'conditionally_implied_expression_PCs': [r['expression_PCs_if_RNA_count_equals_model_N_and_five_genotype_PCs_only'] for r in summaries]}, indent=2))


if __name__ == '__main__': main()
