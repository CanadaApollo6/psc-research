"""Full-allele and complete-evidence audit of four original monocyte exports."""

import json
from pathlib import Path

import pandas as pd

from analyze_prkd2_recovery import overlap, save_csv
from coloc_common import ROOT, save_json, sha
from prkd2_expanded_common import load_references, qtl_identity


def main():
    plan = json.loads((ROOT / 'config/prkd2-recovery-plan.json').read_text())
    author_path = ROOT / 'config/prkd2-recovery-author-plan.json'
    authors = json.loads(author_path.read_text())
    queries = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-queries.json').read_text())
    reference = load_references()['GRCh38']
    targets = pd.read_csv(ROOT / 'data/derived/prkd2-recovery-targets.csv', keep_default_na=False)
    gwas = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-PRKD2.csv.gz')
    full_gwas = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz')
    summaries, matches, outputs = [], [], []
    for selected in authors['ibdverse']:
        q = next(q for q in queries['queries'] if q['id'] == selected['id'])
        if not q.get('result', {}).get('deflate_eof_verified'):
            raise ValueError('Incomplete original IBDverse member: ' + q['id'])
        receipt = q['result']['outputs'][0]
        source = ROOT / receipt['file']
        if sha(source.read_bytes()) != receipt['sha256']:
            raise ValueError('Changed original IBDverse extract')
        frame = pd.read_csv(source, sep='\t', dtype=str, keep_default_na=False)
        if len(frame) != receipt['rows'] or not frame.phenotype_id.str.split('.').str[0].eq(plan['gene_id']).all():
            raise ValueError('Unexpected source gene or row count')
        rows = []
        for index, row in frame.iterrows():
            fields = row.variant_id.split(':')
            if len(fields) != 4 or fields[0] != 'chr19' or not fields[1].isdigit():
                raise ValueError('Non-explicit original IBDverse allele identity')
            position = int(fields[1])
            item = {'source_row': index, **row.to_dict(), 'position': position,
                    'source_beta': row.slope, 'source_se': row.slope_se, 'source_pvalue': row.pval_nominal,
                    'source_af': row.af, 'snp': '', 'canonical_beta_sign': '', 'identity_status': '',
                    'in_fixed_region': plan['start_grch38_0based'] < position <= plan['end_grch38_exclusive']}
            try:
                mapped = qtl_identity('_'.join(fields), reference)
                item.update({'snp': mapped['snp'], 'canonical_beta_sign': mapped['canonical_beta_sign'], 'identity_status': 'verified-full-edit'})
            except ValueError as exc:
                item['identity_status'] = str(exc)
            rows.append(item)
        mapped = pd.DataFrame(rows)
        if len(mapped) == 0:
            raise ValueError('PRKD2 has no source rows; do not calculate coverage')
        for key, group in mapped[mapped.snp.ne('')].groupby('snp'):
            if len(group[['source_beta', 'source_se', 'source_pvalue', 'canonical_beta_sign']].drop_duplicates()) > 1:
                mapped.loc[group.index, 'snp'] = ''
                mapped.loc[group.index, 'identity_status'] = 'conflicting equivalent edits'
        path = ROOT / f"data/derived/prkd2-recovery-{q['id']}-alleles.csv.gz"
        outputs.append(save_csv(mapped, path, True))
        summary = {'dataset_id': q['id'], 'context': selected['label'], 'source_chromosome_rows': q['result']['source_rows'],
                   'source_build': 'GRCh38', 'effect_allele': 'ALT dosage, explicit in primary methods',
                   'cohort_context': 'Original IBDverse blood; Crohn disease and healthy participants, not a PSC case-control expression comparison',
                   **overlap(mapped, gwas, full_gwas)}
        summary['model_gate_reason'] = 'No complete original-source component vectors or adequate study-matched LD/covariance recovered. Catalogue sample counts are not silently assigned to this different original export.'
        summary['RNA_weight_status'] = 'Diagnostic Gaussian-ABF overlap in original inverse-normal expression units; finite-sample distribution check is separate.'
        summary['primary_variant_present'] = bool(mapped.snp.eq(targets.iloc[0].canonical_edit_grch38).any())
        summary['source_coordinate_min'] = int(mapped.position.min())
        summary['source_coordinate_max'] = int(mapped.position.max())
        summaries.append(summary)
        for _, target in targets.iterrows():
            found = mapped[mapped.snp.eq(target.canonical_edit_grch38)]
            if len(found):
                for _, row in found.iterrows():
                    matches.append({'dataset_id': q['id'], 'context': selected['label'], 'target_rsid': target.rsid,
                                    'target_edit_grch38': target.canonical_edit_grch38, 'status': 'present',
                                    'variant_id': row.variant_id, 'source_beta': row.source_beta, 'source_se': row.source_se,
                                    'source_pvalue': row.source_pvalue, 'source_af': row.source_af, 'canonical_beta_sign': row.canonical_beta_sign})
            else:
                matches.append({'dataset_id': q['id'], 'context': selected['label'], 'target_rsid': target.rsid,
                                'target_edit_grch38': target.canonical_edit_grch38, 'status': 'absent-from-complete-source-PRKD2-vector'})
    outputs.append(save_csv(pd.DataFrame(summaries), ROOT / 'data/derived/prkd2-recovery-ibdverse-summary.csv'))
    outputs.append(save_csv(pd.DataFrame(matches), ROOT / 'data/derived/prkd2-recovery-ibdverse-targets.csv'))
    result = {'analysis_status': 'complete', 'author_plan_sha256': sha(author_path.read_bytes()),
              'script_sha256': sha(Path(__file__).read_bytes()), 'new_H4_posteriors': 0,
              'contexts': summaries, 'outputs': outputs,
              'limitations': ['Four contexts share one study family.', 'Source nominal vectors are distinct from Catalogue reanalysis.',
                              'Marginal overlap is not multi-signal colocalization.', 'No original-source study-matched covariance is supplied by these nominal files.']}
    save_json(ROOT / 'reports/prkd2-recovery-ibdverse-analysis.json', result)
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
