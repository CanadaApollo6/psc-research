"""Freeze the two-gene, ten-context regional colocalisation follow-up."""

import csv
import gzip
import json
from datetime import datetime, timezone

from coloc_common import ROOT, ChainMap, save_json, sha, verify_plan


def main():
    path = ROOT / 'config/coloc-plan.json'
    if path.exists():
        plan = verify_plan()
        print(json.dumps({'status': 'verified-existing-plan', 'frozen_at_utc': plan['frozen_at_utc']}))
        return
    source_plan = json.loads((ROOT / 'config/cd4-rna-plan.json').read_text())
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    genes = {}
    with gzip.open(ROOT / 'data/raw/cd4-gene-phenotypes.tsv.gz', 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['gene_id'] in ['ENSG00000170525', 'ENSG00000153094']:
                genes[row['gene_id']] = row
    regions = []
    for gene, anchor in [('ENSG00000170525', 'rs7923054'), ('ENSG00000153094', 'rs72837826')]:
        row = genes[gene]
        chrom, tss = 'chr' + row['chromosome'], int(row['phenotype_pos'])
        start0, end0 = tss - 1 - 1_000_000, tss + 1_000_000
        old_start, old_end = reverse.enclosing_destination(chrom, start0, end0)
        regions.append({'gene_id': gene, 'gene_name': row['gene_name'], 'anchor': anchor,
                        'chromosome': chrom, 'tss_grch38_1based': tss, 'start_grch38_0based': start0,
                        'end_grch38_exclusive': end0, 'start_grch37_0based': max(0, old_start-1000),
                        'end_grch37_exclusive': old_end+1000})
    frozen = ['config/cd4-rna-plan.json', 'config/cd4-rna-sources.json', 'reports/cd4-rna-summary.json',
              'reports/cd4-rna-evidence.csv', 'data/raw/cd4-gene-phenotypes.tsv.gz',
              'data/raw/coloc-hg19-to-hg38-chain.gz', 'data/raw/coloc-hg38-to-hg19-chain.gz',
              'data/raw/GCST004030-original-header.txt', 'data/raw/GCST004030.json', 'data/raw/Ji2017.html']
    plan = {
        'analysis': 'PSC and CD4 RNA regional colocalisation follow-up v0.1',
        'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
        'selection_status': 'Two genes selected after the prior CD4 RNA screen; exploratory follow-up, not an independent benchmark or registered protocol',
        'frozen_file_sha256': {p: sha((ROOT / p).read_bytes()) for p in frozen},
        'regions': regions, 'datasets': source_plan['datasets'], 'planned_gene_dataset_pairs': 20,
        'regional_rule': 'Full source gene TSS +/- 1 Mb cis window; no result-dependent lead-SNP window',
        'gwas': {
            'accession': 'GCST004030', 'source_build': 'GRCh37',
            'url': 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004030/ipscsg2016.result.combined.full.with_header.txt',
            'source_bytes': 662022204, 'cases': 2871, 'controls': 12019, 'n': 14890,
            'scope': 'Discovery GWAS, not the 4796-case combined replication sample or later fine-mapping run',
            'primary_platform': 'both',
            'other_platforms': 'Preserve rows; exclude from calculations because per-platform sample sizes are not fixed in the source metadata',
            'effect_input': 'Primary marginal ABF uses reported P, control MAF and known case fraction; avoids assuming the undocumented scale of source SE accompanying OR',
            'direction': 'For signed-z reference-LD sensitivity, use inverse-normal two-sided P and source risk-allele orientation; preserve original OR/SE separately',
            'minimum_maf': 0.01,
        },
        'alignment': {
            'variants': 'Biallelic single-nucleotide substitutions only; retain all excluded indels, duplicate conflicts and missing mappings in the audit',
            'build': 'Unique reciprocal UCSC chain mapping, exact chromosome/position and allele set, explicit strand complement only from chain orientation',
            'palindromic': 'Require exact coordinate/alleles plus GWAS risk/control-frequency versus reference EUR frequency difference <= 0.15; otherwise exclude',
            'qtl': 'Exact GRCh38 REF/ALT; numeric Ensembl gene-version removal only; collapse rsID-only aliases, exclude conflicting rows',
            'missing': 'Unavailable is never a null result; do not substitute LD proxies',
        },
        'methods': {
            'preferred': 'Signal-wise coloc using complete published SuSiE log Bayes factors or coherent signed LD with documented sample/allele provenance',
            'baseline': 'Official coloc.abf 5.2.3; at most one causal variant per trait; all 20 pairs retained regardless of result',
            'qtl_baseline': 'Reported beta and SE, quantitative trait, sdY=1 for inverse-normalized expression; sensitivity estimates sdY from MAF/SE/N',
            'priors': {'p1': 0.0001, 'p2': 0.0001, 'p12_primary': 0.000001, 'p12_sensitivity': [0.0000001, 0.00001]},
            'reference_ld_sensitivity': 'If complete published BFs or in-sample LD are unavailable, fit GWAS and BLUEPRINT only with signed LD from 503 EUR phase-3 donors. This remains approximate external-LD sensitivity. DICE is not assumed ancestry-matched to EUR.',
            'susie': {'L': 10, 'coverage': 0.95, 'min_abs_corr': 0.5, 'estimate_residual_variance': False, 'max_iter': 200, 'seed': 20260912},
            'ld_qc': 'Exact SNP/allele order; symmetric finite signed r, unit diagonal and PSD; preserve diagnostics, convergence and missing credible sets. No arbitrary ridge chosen for favorable results.',
            'coverage_gate': 'At least 100 shared variants; report full pre-intersection counts and univariate-BF mass retained for both traits. Interpret high posterior only when each retained mass >=0.90; never hide missing strong variants.',
            'interpretation': 'Report H0-H4 for all priors. Shared-signal support requires H4>=0.9 at primary prior and >=0.8 at both sensitivities, adequate coverage and a usable multi-signal assessment; otherwise label assumption-dependent, weak or unavailable.',
            'not_permitted': 'No multiplication by AlphaGenome scores, no favorable-track selection, no back-calculation of incomplete BF vectors from credible-set-only PIPs, no clinical or mediation claim',
        },
        'transfer_bounds': {'minimum_seconds_between_indexed_queries': 2, 'maximum_gwas_range_requests': 80,
                            'maximum_gwas_range_bytes': 10000000, 'maximum_nominal_rows_per_region': 2000000,
                            'maximum_nominal_uncompressed_bytes_per_region': 750000000,
                            'maximum_reference_uncompressed_bytes_per_region': 2000000000,
                            'maximum_reference_rows_per_region': 200000, 'maximum_published_bf_probe_bytes_per_file': 1000000},
        'limits': ['Both QTL cohorts and PSC data are reused from prior publications', 'DICE contexts share donors',
                   'Treg QTD000464/469 naive/memory metadata conflict remains unresolved',
                   'Original fine-mapping association inputs and in-sample LD are not in the public archive',
                   'Susceptibility colocalisation does not establish mediation, progression or treatment benefit'],
    }
    save_json(path, plan)
    print(json.dumps({'status': 'frozen-before-regional-effects', 'frozen_at_utc': plan['frozen_at_utc'], 'pairs': 20, 'regions': regions}, indent=2))


if __name__ == '__main__':
    main()
