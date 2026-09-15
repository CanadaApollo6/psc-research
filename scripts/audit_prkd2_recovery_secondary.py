"""Audit the harmonised PSC copy and bounded author-source availability."""

import json
import math
from pathlib import Path

import pandas as pd

from analyze_prkd2_recovery import save_csv
from coloc_common import ROOT, save_json, sha
from prkd2_expanded_common import load_references, qtl_identity


def main():
    plan = json.loads((ROOT / 'config/prkd2-recovery-plan.json').read_text())
    targets = pd.read_csv(ROOT / 'data/derived/prkd2-recovery-targets.csv', keep_default_na=False)
    refs = load_references()
    source = ROOT / 'data/derived/prkd2-recovery-queries/PSC-harmonised-union-region.tsv.gz'
    query = json.loads((ROOT / 'config/prkd2-recovery-harmonised-query.json').read_text())
    if not query['complete'] or sha(source.read_bytes()) != query['output']['sha256']:
        raise ValueError('Incomplete or changed harmonised source')
    frame = pd.read_csv(source, sep='\t', dtype=str, keep_default_na=False)
    controls = []
    for rsid in ('rs112445263', 'rs313839'):
        target = targets[targets.rsid.eq(rsid)].iloc[0]
        match = frame[frame.hm_rsid.eq(rsid)]
        if len(match) != 1:
            raise ValueError('Nonunique harmonised build control')
        row = match.iloc[0]
        expected = f"19_{target.position}_{target.ref}_{target.alt}"
        if row.hm_variant_id != expected or int(row.hm_pos) != target.position:
            raise ValueError('Harmonised control disagrees with GRCh38 reference')
        if int(row.hm_pos) == target.position_grch37:
            raise ValueError('Harmonised control did not discriminate the builds')
        controls.append({'rsid': rsid, 'verified_build': 'GRCh38', 'hm_variant_id': row.hm_variant_id,
                         'source_effect_allele': row.effect_allele, 'source_odds_ratio': row.odds_ratio,
                         'harmonised_effect_allele': row.hm_effect_allele, 'harmonised_odds_ratio': row.hm_odds_ratio,
                         'OR_reciprocal_check': math.isclose(float(row.odds_ratio) * float(row.hm_odds_ratio), 1, rel_tol=1e-12)})
    rows = []
    for index, row in frame.iterrows():
        item = {'retained_union_row': index, **row.to_dict(), 'canonical_edit_grch38': '', 'identity_status': '',
                'unharmonised_source_edit_candidates': '', 'unharmonised_source_identity_status': 'not-needed'}
        try:
            variant = 'chr' + row.hm_variant_id
            identity = qtl_identity(variant, refs['GRCh38'])
            item['canonical_edit_grch38'] = identity['snp']
            item['identity_status'] = 'verified-full-edit'
        except (ValueError, IndexError) as exc:
            item['identity_status'] = str(exc)
        if row.hm_variant_id == 'NA':
            candidates = set()
            failures = []
            for ref, alt in [(row.other_allele, row.effect_allele), (row.effect_allele, row.other_allele)]:
                try:
                    variant = f'chr{row.chromosome}_{row.base_pair_location}_{ref}_{alt}'
                    candidates.add(qtl_identity(variant, refs['GRCh38'])['snp'])
                except (ValueError, IndexError) as exc:
                    failures.append(str(exc))
            item['unharmonised_source_edit_candidates'] = '|'.join(sorted(candidates))
            item['unharmonised_source_identity_status'] = 'verified-full-source-edit' if len(candidates) == 1 else 'ambiguous' if candidates else '; '.join(sorted(set(failures)))
        rows.append(item)
    identities = pd.DataFrame(rows)
    output = ROOT / 'data/derived/prkd2-recovery-harmonised-identity-audit.csv.gz'
    outputs = [save_csv(identities, output, True)]
    exact = []
    for _, target in targets.iterrows():
        matches = identities[identities.canonical_edit_grch38.eq(target.canonical_edit_grch38)]
        fallback = identities[identities.unharmonised_source_edit_candidates.str.split('|').map(lambda x: target.canonical_edit_grch38 in x)]
        exact.append({'target_rsid': target.rsid, 'role': target.role, 'canonical_edit_grch38': target.canonical_edit_grch38,
                      'harmonised_exact_rows': len(matches), 'source_rsids': ';'.join(matches.variant_id),
                      'unharmonised_source_candidate_rows': len(fallback),
                      'status': 'present' if len(matches) else 'source-candidate-recovered' if len(fallback) else 'not-recovered-from-complete-harmonised-regional-copy'})
    outputs.append(save_csv(pd.DataFrame(exact), ROOT / 'data/derived/prkd2-recovery-harmonised-targets.csv'))
    numeric_pos = pd.to_numeric(identities.hm_pos, errors='coerce')
    regional = identities[numeric_pos.gt(plan['start_grch38_0based']) & numeric_pos.le(plan['end_grch38_exclusive'])]
    prior = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz', keep_default_na=False)
    prior_ids = set(prior.SNP.str.strip())
    current_ids = set(regional.variant_id.str.strip())
    sources = json.loads((ROOT / 'config/prkd2-recovery-extra-sources.json').read_text())['sources']
    conditional, inventory = [], []
    for context in ('UT', 'IFN', 'LPS24'):
        receipt = next(x for x in sources if x['id'] == f'nassiri-{context}-conditional-profile')
        raw = ROOT / receipt['file']
        if sha(raw.read_bytes()) != receipt['sha256']:
            raise ValueError('Changed conditional author export')
        data = pd.read_csv(raw, sep='\t', dtype=str, keep_default_na=False)
        if not {'fwd_pval', 'bwd_pval', 'rank', 'phe_id', 'var_id'} <= set(data):
            raise ValueError('Unexpected author conditional-export schema')
        gene = data[data.phe_id.eq('PRKD2')]
        source_file = ROOT / f'data/derived/prkd2-recovery-queries/Nassiri-original-{context}-PRKD2-conditional.tsv.gz'
        from coloc_common import gzip_bytes
        encoded = gzip_bytes(gene.to_csv(sep='\t', index=False, lineterminator='\n').encode())
        source_file.write_bytes(encoded)
        outputs.append({'file': str(source_file.relative_to(ROOT)), 'rows': len(gene), 'sha256': sha(encoded)})
        for _, row in gene.iterrows():
            conditional.append({'context': context, **row.to_dict()})
        inventory.append({'family': 'Nassiri_2025', 'context': context, 'full_export_rows': len(data),
                          'PRKD2_rows': len(gene), 'primary_rsid_PRKD2_rows': int(gene.var_id.eq('rs112445263').sum()),
                          'full_nominal_vector': False, 'reason': 'Forward/backward conditional output; no explicit allele pair; absence is not a null association.'})
    outputs.append(save_csv(pd.DataFrame(conditional), ROOT / 'data/derived/prkd2-recovery-nassiri-conditional-rows.csv'))
    natri = next(x for x in sources if x['id'] == 'natri-limix-immune-head')
    inventory.append({'family': 'Natri_2024', 'context': 'Monocyte and Inflammatory monocyte; secondary lung context',
                      'full_nominal_vector': None, 'reason': 'Public LIMIX immune TAR.GZ is unindexed and exceeds the frozen transfer cap; no variant absence claim.',
                      'archive_bytes': int(natri['headers']['Content-Length']), 'url': natri['url']})
    outputs.append(save_csv(pd.DataFrame(inventory), ROOT / 'data/derived/prkd2-recovery-secondary-author-availability.csv'))
    result = {'script_sha256': sha(Path(__file__).read_bytes()), 'new_H4_posteriors': 0,
              'harmonised_source_rows': query['source_rows'], 'harmonised_retained_union_rows': len(frame),
              'harmonised_fixed_region_rows': len(regional), 'original_fixed_region_rows': len(prior),
              'unharmonised_retained_union_rows': int(identities.hm_variant_id.eq('NA').sum()),
              'unharmonised_verified_source_edits': int(identities.unharmonised_source_identity_status.eq('verified-full-source-edit').sum()),
              'harmonised_controls': controls, 'target_results': exact,
              'original_rsids_missing_from_harmonised_fixed_window': sorted(prior_ids - current_ids),
              'harmonised_rsids_new_to_original_fixed_window': sorted(current_ids - prior_ids),
              'secondary_author_availability': inventory, 'outputs': outputs,
              'interpretation': 'The harmonised file is a representation of the same PSC study, not an independent cohort. Conditional author exports cannot supply a complete marginal association vector.'}
    save_json(ROOT / 'reports/prkd2-recovery-secondary-audit.json', result)
    print(json.dumps({k: result[k] for k in ['harmonised_fixed_region_rows', 'original_fixed_region_rows', 'harmonised_controls', 'target_results', 'secondary_author_availability']}, indent=2))


if __name__ == '__main__':
    main()
