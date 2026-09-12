"""Freeze the complete CD4 follow-up design before querying new QTL effects."""

import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path, delimiter=','):
    return list(csv.DictReader(path.open(), delimiter=delimiter))


def main():
    original_sources = json.loads((ROOT / 'config/qtl-followup-sources.json').read_text())['sources']
    by_source = {s['id']: s for s in original_sources}
    source_ids = ['eqtl-r7-metadata', 'eqtl-tabix-paths', 'eqtl-columns', 'eqtl-data-access', 'eqtl-methods',
                  'eqtl-study-design', 'DICE-paper-html', 'blueprint-paper', 'alphagenome-supplement-tables']
    for name in source_ids:
        source = by_source[name]
        if digest(ROOT / 'data/raw' / source['file']) != source['sha256']:
            raise ValueError('Original source changed: ' + name)
    metadata = rows(ROOT / 'data/raw' / by_source['eqtl-r7-metadata']['file'], '\t')
    paths = {r['dataset_id']: r for r in rows(ROOT / 'data/raw' / by_source['eqtl-tabix-paths']['file'], '\t')}
    groups = {'Tfh_memory', 'Th17_memory', 'Th1_memory', 'Th2_memory', 'Th1-17_memory',
              'Treg_memory', 'Treg_naive', 'CD4_T-cell_naive', 'CD4_T-cell_anti-CD3-CD28'}
    selected = [r for r in metadata if r['quant_method'] == 'ge' and
                (r['study_id'] == 'QTS000002' and r['sample_group'] == 'T-cell' or
                 r['study_id'] == 'QTS000026' and r['sample_group'] in groups)]
    if len(selected) != 10 or len({r['dataset_id'] for r in selected}) != 10:
        raise ValueError('Unexpected CD4 dataset inventory')
    sources, datasets = [by_source[name] for name in source_ids], []
    for row in selected:
        data = paths[row['dataset_id']]
        if any(data[k] != row[k] for k in data if k in row):
            raise ValueError('Metadata tables disagree')
        identifier = row['dataset_id']
        url = data['ftp_path'].replace('ftp://', 'https://', 1)
        if not url.endswith('.all.tsv.gz'):
            raise ValueError('Expected unfiltered nominal gene-level archive')
        primary = row['sample_group'] in ('T-cell', 'CD4_T-cell_naive')
        datasets.append({**row, 'url': url, 'cohort_label': 'BLUEPRINT' if row['study_id'] == 'QTS000002' else 'DICE',
                         'reporting_context': 'closest_unstimulated' if primary else 'activation' if 'anti-CD3' in row['sample_group'] else 'specialized_subset',
                         'assay_matched_model_track': 'CL:0000624 total RNA-seq' if row['study_id'] == 'QTS000002' else 'CL:0000624 polyA plus RNA-seq',
                         'context_note': 'Sample-group naive/memory label conflicts with tissue_label; preserve both and do not resolve cell state by inference' if identifier in ('QTD000464', 'QTD000469') else '',
                         'independence_note': 'Source cohort reused in PSC fine-mapping; not independent replication of that evidence' if row['study_id'] == 'QTS000002' else 'Separate named cohort from BLUEPRINT; DICE conditions and subsets share donors; exact model-training donor overlap unresolved'})
        for suffix, extra in [('index', {'url': url + '.tbi', 'max_bytes': 12_000_000}),
                              ('header', {'url': url, 'max_bytes': 65536, 'request_headers': {'Range': 'bytes=0-65535'}})]:
            sid = identifier + '-' + suffix
            sources.append(by_source[sid] if sid in by_source else {'id': sid, 'file': 'cd4-' + sid + ('.tbi' if suffix == 'index' else '.bgz'), **extra})
    sources.append({'id': 'cd4-phenotype-record', 'file': 'cd4-phenotype-record.json',
                    'url': 'https://zenodo.org/api/records/7808390', 'max_bytes': 5_000_000})
    for accession in ['ENCSR033XWU', 'ENCSR411MUF', 'ENCSR463JBR', 'ENCSR545MEZ']:
        sources.append({'id': accession, 'file': 'cd4-' + accession + '.json',
                        'url': f'https://www.encodeproject.org/experiments/{accession}/?format=json', 'max_bytes': 5_000_000})
    locked = json.loads((ROOT / 'config/psc-matched-comparison-inputs.json').read_text())
    ranks = rows(ROOT / 'reports/matched-comparison-gene-rankings.csv')
    scores = rows(ROOT / 'reports/matched-comparison-primary-track-scores.csv')
    score_by_key = {(r['rsid'], r['gene_id'], r['track_name']): float(r['raw_score']) for r in scores}
    if len(score_by_key) != len(scores):
        raise ValueError('Repeated model track score')
    rank_by_key = {(r['rsid'], r['gene_id']): r for r in ranks}
    pairs = []
    for variant in locked['inputs']:
        for gene in locked['gene_universe_by_region'][variant['archive_region']]:
            rank = rank_by_key[variant['rsid'], gene['gene_id']]
            pair = {key: variant[key] for key in ['archive_region', 'matched_anchor', 'role', 'rsid', 'chromosome', 'position_grch38_1based', 'reference_bases', 'alternate_bases']}
            values = [score_by_key[variant['rsid'], gene['gene_id'], track] for track in locked['expected_merged_primary_track_names']]
            direction = 'increase' if all(v > 0 for v in values) else 'decrease' if all(v < 0 for v in values) else 'mixed_or_zero'
            pair.update(gene_id=gene['gene_id'], gene_name=gene['gene_name'], gene_biotype=gene['biotype'], model_gene_rank=rank['model_rank'],
                        model_gene_effect=rank['gene_effect'], model_consensus_direction=direction,
                        model_total_rna_score=score_by_key[variant['rsid'], gene['gene_id'], 'CL:0000624 total RNA-seq'],
                        model_polya_rna_score=score_by_key[variant['rsid'], gene['gene_id'], 'CL:0000624 polyA plus RNA-seq'],
                        is_region_named_gene=gene['gene_name'] == variant['archive_region'])
            pairs.append(pair)
    if len(pairs) != 228 or len({(p['rsid'], p['gene_id']) for p in pairs}) != 228:
        raise ValueError('The complete 228-pair universe is required')
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(pairs[0]), lineterminator='\n')
    writer.writeheader(); writer.writerows(pairs)
    pair_bytes = buffer.getvalue().encode()
    pair_path = ROOT / 'data/derived/cd4-rna-planned-pairs.csv'
    immutable = ['config/psc-matched-comparison-inputs.json', 'reports/matched-comparison-gene-rankings.csv',
                 'reports/matched-comparison-primary-track-scores.csv', 'config/qtl-followup-sources.json']
    payload = {
        'protocol_version': '0.1', 'purpose': 'Measured CD4 gene-expression evidence for every variant/gene in the completed matched model comparison',
        'context': 'Follow-up selected after model results; new target QTL coefficients have not been queried or inspected. This is a local analysis freeze, not a registry preregistration.',
        'datasets': datasets, 'variants': locked['inputs'], 'variant_gene_pairs': 228, 'planned_dataset_variant_gene_tests': 2280,
        'model_requests_planned': 0, 'planned_pairs_sha256': hashlib.sha256(pair_bytes).hexdigest(),
        'frozen_file_sha256': {name: digest(ROOT / name) for name in immutable},
        'selection': 'All nine DICE CD4-related bulk gene-expression subsets in the pinned release-7 metadata, plus BLUEPRINT CD4; no effect-based dataset selection. Other cohorts and modalities are outside this bounded follow-up.',
        'alignment': 'Exact GRCh38 chromosome, 1-based position, REF and ALT. ALT is the source effect allele. Match stable Ensembl gene ID with only numeric version suffix removed; require gene-level molecular_trait_id. No rsID-only matching, complement guessing, LD proxies or gene-symbol substitution.',
        'duplicates': 'Collapse rows differing only in rsID aliases; retain every original row and all aliases. Conflicting measurements for one exact variant/gene key are unavailable and explicitly reported.',
        'missingness': 'Account for all 2280 planned cells. Separate failed query, absent exact-allele variant, absent gene annotation, absent pair, conflicting duplicate, and invalid numerical measurement. Missing or filtered rows never imply a null effect.',
        'numerical_checks': 'Require finite beta, finite positive SE, P in [0,1], positive even AN at most twice metadata sample size, and AC in [0,AN]. Preserve optional imputation quality and median TPM; missing optional fields stay unknown. Flag R2 below 0.8 separately without deleting rows.',
        'primary_direction': 'Compare ALT coefficient sign with model direction only where both fixed CD4 tracks have the same nonzero sign. Keep mixed/zero model tracks explicitly unavailable for the primary direction comparison. Signs alone are descriptive, including non-significant estimates.',
        'assay_secondary': 'Also report the fixed assay-matched score: total RNA for BLUEPRINT, polyA RNA for DICE, based on original library methods. This secondary comparison does not replace mixed primary directions.',
        'uncertainty': 'Preserve source beta, SE and nominal P. Report descriptive beta +/- 1.96 SE intervals; QTL beta is on a normalized-expression scale, not model log change or measured percent change.',
        'multiplicity': {'method': 'Bonferroni across all 2280 planned dataset-variant-gene tests', 'family_alpha': 0.05, 'p_threshold': 0.05 / 2280,
                         'interpretation': 'Conservative family-wise screening; missing tests remain in the planned denominator. This flags molecular associations, not causal variants or model accuracy.'},
        'reporting': 'Complete coverage/effect table; all corrected-threshold associations; coverage by dataset/region/role; model-maximum gene and region-named gene panels for every variant. No strongest-only reporting, meta-analysis, accuracy percentage or treatment inference.',
        'independence': 'Region and cohort are reporting units. DICE datasets share donors. BLUEPRINT reuses evidence underlying PSC molecular fine-mapping. ALL_FOLDS does not hold out these loci. Exact donor overlap with AlphaGenome training remains unestablished; ENCODE metadata will be checked without private donor matching.',
        'transfer_bounds': {'queries': 90, 'minimum_seconds_between_queries': 2.0, 'maximum_rows_per_position': 5000, 'maximum_uncompressed_bytes_per_position': 5_000_000},
    }
    plan_path = ROOT / 'config/cd4-rna-plan.json'
    if plan_path.exists():
        previous = json.loads(plan_path.read_text()); timestamp = previous.pop('frozen_at_utc')
        if previous != payload or pair_path.read_bytes() != pair_bytes:
            raise ValueError('Existing CD4 plan or planned pairs changed')
        print(json.dumps({'status': 'verified-existing-plan', 'frozen_at_utc': timestamp, 'tests': 2280}))
        return
    payload['frozen_at_utc'] = datetime.now(timezone.utc).isoformat()
    pair_path.write_bytes(pair_bytes)
    plan_path.write_text(json.dumps(payload, indent=2) + '\n')
    (ROOT / 'config/cd4-rna-sources.json').write_text(json.dumps({'sources': sources}, indent=2) + '\n')
    print(json.dumps({'status': 'frozen-before-new-QTL-effects', 'frozen_at_utc': payload['frozen_at_utc'], 'datasets': len(datasets), 'tests': 2280}))


if __name__ == '__main__':
    main()
