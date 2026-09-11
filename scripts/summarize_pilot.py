"""Apply the prespecified pilot ranking and report every selected variant."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LFC = 'GeneMaskLFCScorer(requested_output=RNA_SEQ)'
ACTIVE = 'GeneMaskActiveScorer(requested_output=RNA_SEQ)'
SPLICE = 'GeneMaskSplicingScorer(requested_output=SPLICE_SITES, width=None)'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_genes(scores, candidates, ontology):
    primary = scores[(scores['variant_scorer'] == LFC) & (scores['ontology_curie'] == ontology)].copy()
    primary['gene_id'] = primary['gene_id'].str.split('.').str[0]
    primary = primary[primary['gene_id'].isin(candidates['gene_id'])]
    if primary.empty or not np.isfinite(primary.raw_score).all():
        raise ValueError('Missing or non-finite primary scores')
    primary['absolute_score'] = primary['raw_score'].abs()
    ranked = primary.groupby('gene_id', as_index=False).agg(
        median_absolute_lfc=('absolute_score', 'median'),
        median_signed_lfc=('raw_score', 'median'),
        primary_track_count=('raw_score', 'size'),
        minimum_signed_lfc=('raw_score', 'min'), maximum_signed_lfc=('raw_score', 'max'),
    )
    ranked = ranked.merge(candidates[['gene_id', 'gene_name', 'biotype', 'distance_to_variant']], on='gene_id', validate='one_to_one')
    ranked = ranked.sort_values(['median_absolute_lfc', 'gene_id'], ascending=[False, True]).reset_index(drop=True)
    ranked['model_rank'] = ranked.index + 1
    baseline = ranked.sort_values(['distance_to_variant', 'gene_id']).iloc[0]
    return ranked, baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', help='Local run directory containing run.json and score files')
    args = parser.parse_args()
    folder = Path(args.run_directory)
    manifest = json.loads((folder / 'run.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Report requires all four prespecified requests; failed requests must remain visible')
    if digest(folder / 'inputs.json') != manifest['input_sha256'] or digest(folder / 'protocol.json') != manifest['protocol_sha256']:
        raise ValueError('Run inputs or protocol changed')
    protocol = json.loads((folder / 'protocol.json').read_text())
    inputs = json.loads((folder / 'inputs.json').read_text())
    candidates = pd.read_csv(ROOT / 'data/derived/benchmark-candidate-genes.csv')
    summary, rankings, tracks, splice_rows = [], [], [], []
    for row in inputs:
        request = next(entry for entry in manifest['requests'] if entry['rsid'] == row['rsid'])
        path = folder / request['file']
        if digest(path) != request['sha256']:
            raise ValueError('Prediction file checksum changed')
        scores = pd.read_csv(path)
        genes = candidates[candidates.rsid == row['rsid']]
        primary = protocol['primary_ontology_by_gene'][row['target_gene']]
        ranking, baseline = rank_genes(scores, genes, primary)
        target_rows = ranking[ranking.gene_id == row['target_gene_id']]
        if len(target_rows) != 1:
            raise ValueError('Expected target missing or duplicated in model gene scores')
        target = target_rows.iloc[0]
        expected = row['expected_alt_minus_ref_expression_sign']
        actual = int(np.sign(target['median_signed_lfc']))
        agreement = 'excluded: unresolved source direction' if expected is None else 'inconclusive: zero score' if actual == 0 else 'matches source sign' if actual == expected else 'opposes source sign'
        active = scores[(scores.variant_scorer == ACTIVE) & (scores.ontology_curie == primary) & (scores.gene_id.str.split('.').str[0] == row['target_gene_id'])]
        summary.append({
            'rsid': row['rsid'], 'target_gene': row['target_gene'], 'allele_change': row['reference_bases'] + '>' + row['alternate_bases'],
            'primary_ontology': primary, 'target_rank': int(target.model_rank), 'evaluable_genes': len(ranking),
            'annotated_candidates': len(genes), 'candidates_not_scored': len(genes) - len(ranking),
            'median_absolute_lfc': target.median_absolute_lfc, 'median_signed_lfc': target.median_signed_lfc,
            'primary_track_count': int(target.primary_track_count), 'tracks_disagree_in_sign': bool(target.minimum_signed_lfc < 0 < target.maximum_signed_lfc),
            'expected_alt_minus_ref_sign': expected, 'direction_comparison': agreement,
            'nearest_tss_evaluable_gene': baseline.gene_name, 'nearest_tss_full_reference_gene': row['nearest_tss_gene'],
            'target_median_active_score': float(active.raw_score.median()) if len(active) else None,
        })
        rankings.append(ranking.assign(rsid=row['rsid'], target_gene=row['target_gene']))
        selected = [primary] + protocol['supporting_ontologies_by_gene'][row['target_gene']]
        target_mask = scores.gene_id.str.split('.').str[0] == row['target_gene_id']
        tracks.append(scores[target_mask & (scores.variant_scorer == LFC) & scores.ontology_curie.isin(selected)].assign(rsid=row['rsid'], primary=lambda frame: frame.ontology_curie == primary))
        splice_rows.append(scores[target_mask & (scores.variant_scorer == SPLICE)].assign(rsid=row['rsid']))
    reports = ROOT / 'reports'
    pd.DataFrame(summary).to_csv(reports / 'pilot-target-summary.csv', index=False)
    pd.concat(rankings).to_csv(reports / 'pilot-gene-rankings.csv', index=False)
    pd.concat(tracks).to_csv(reports / 'pilot-target-track-scores.csv', index=False)
    pd.concat(splice_rows).to_csv(reports / 'pilot-target-splice-scores.csv', index=False)
    (reports / 'pilot-run-provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == '__main__':
    main()
