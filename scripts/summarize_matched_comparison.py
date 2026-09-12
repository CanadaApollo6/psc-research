"""Report every fixed matched comparison, including unavailable model outputs."""

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LFC = 'GeneMaskLFCScorer(requested_output=RNA_SEQ)'


def score_variant_table(table, genes, track_names, ontology):
    subset = table[(table['variant_scorer'] == LFC) & (table['ontology_curie'] == ontology)].copy()
    subset['gene_id'] = subset['gene_id'].str.split('.').str[0]
    ids = {g['gene_id'] for g in genes}
    subset = subset[subset['gene_id'].isin(ids)]
    expected = {(identifier, track) for identifier in ids for track in track_names}
    observed = list(zip(subset['gene_id'], subset['track_name']))
    if len(observed) != len(set(observed)) or set(observed) != expected:
        raise ValueError('Missing, duplicate or unexpected fixed gene/track scores')
    if not np.isfinite(subset['raw_score']).all():
        raise ValueError('Nonfinite fixed gene/track score')
    subset['absolute_score'] = subset['raw_score'].abs()
    ranked = subset.groupby('gene_id', as_index=False).agg(gene_effect=('absolute_score', 'median'),
                                                         signed_median=('raw_score', 'median'), tracks=('raw_score', 'size'))
    names = {g['gene_id']: g['gene_name'] for g in genes}
    ranked['gene_name'] = ranked['gene_id'].map(names)
    ranked = ranked.sort_values(['gene_effect', 'gene_id'], ascending=[False, True]).reset_index(drop=True)
    ranked['model_rank'] = ranked.index + 1
    return ranked, subset


def matched_contrast(anchor, comparator1, comparator2):
    values = [anchor, comparator1, comparator2]
    if any(value is None or not np.isfinite(value) for value in values):
        return None
    return float(anchor - np.median([comparator1, comparator2]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    folder = args.run_directory
    run = json.loads((folder / 'run.json').read_text())
    locked = json.loads((folder / 'inputs.json').read_text())
    protocol = json.loads((folder / 'protocol.json').read_text())
    for name, expected in [('inputs.json', run['input_sha256']), ('protocol.json', run['protocol_sha256'])]:
        if hashlib.sha256((folder / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Saved input or protocol changed')
    variants, ranks, track_rows = [], [], []
    for row in locked['inputs']:
        output = {**row, 'status': '', 'variant_effect': None, 'max_gene_id': '', 'max_gene_name': '',
                  'max_gene_signed_median': None, 'nearest_tss_model_rank': None,
                  'fixed_gene_count': len(locked['gene_universe_by_region'][row['archive_region']])}
        request = next((r for r in run['requests'] if r['rsid'] == row['rsid']), None)
        if not request or request['status'] != 'complete':
            output['status'] = 'unavailable: model request incomplete'
            variants.append(output)
            continue
        path = folder / request['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != request['sha256']:
            raise ValueError('Prediction checksum changed')
        table = pd.read_csv(path)
        try:
            ranked, scores = score_variant_table(table, locked['gene_universe_by_region'][row['archive_region']],
                                                 locked['expected_merged_primary_track_names'], protocol['model']['primary_ontology'])
            top = ranked.iloc[0]
            nearest = ranked[ranked.gene_id == row['nearest_tss_gene_id']]
            if len(nearest) != 1:
                raise ValueError('Nearest-TSS gene not in the fixed scored universe')
            output.update(status='complete', variant_effect=float(top.gene_effect), max_gene_id=top.gene_id, max_gene_name=top.gene_name,
                          max_gene_signed_median=float(top.signed_median), nearest_tss_model_rank=int(nearest.iloc[0].model_rank))
            ranks.append(ranked.assign(rsid=row['rsid'], archive_region=row['archive_region'], role=row['role'], matched_anchor=row['matched_anchor']))
            track_rows.append(scores[['gene_id', 'track_name', 'track_strand', 'raw_score']].assign(rsid=row['rsid'], archive_region=row['archive_region']))
        except ValueError as exc:
            output['status'] = 'unavailable: ' + str(exc)
        variants.append(output)
    by_id = {r['rsid']: r for r in variants}
    contrasts = []
    for status in locked['all_anchor_statuses']:
        result = {'archive_region': status['archive_region'], 'anchor': status['anchor'], 'status': status['status'],
                  'anchor_effect': None, 'comparator_1': status['comparator_1'], 'comparator_1_effect': None,
                  'comparator_2': status['comparator_2'], 'comparator_2_effect': None, 'comparator_median_effect': None, 'anchor_contrast': None}
        if status['status'] == 'matched':
            anchor, first, second = [by_id[status[k]] for k in ['anchor', 'comparator_1', 'comparator_2']]
            contrast = matched_contrast(anchor['variant_effect'], first['variant_effect'], second['variant_effect'])
            result.update(status='complete' if contrast is not None else 'unavailable: incomplete fixed model outputs',
                          anchor_effect=anchor['variant_effect'], comparator_1_effect=first['variant_effect'], comparator_2_effect=second['variant_effect'],
                          comparator_median_effect=float(np.median([first['variant_effect'], second['variant_effect']])) if contrast is not None else None,
                          anchor_contrast=contrast)
        contrasts.append(result)
    regions = []
    for label in protocol['locus_selection']['selected_archive_labels']:
        regional = [r for r in contrasts if r['archive_region'] == label]
        values = [r['anchor_contrast'] for r in regional if r['anchor_contrast'] is not None]
        regions.append({'archive_region': label, 'planned_anchors': len(regional), 'evaluable_anchors': len(values),
                        'region_median_contrast': float(np.median(values)) if values else None,
                        'status': 'complete' if len(values) == len(regional) else 'partial' if values else 'unavailable'})
    reports = ROOT / 'reports'
    pd.DataFrame(variants).to_csv(reports / 'matched-comparison-variant-effects.csv', index=False)
    pd.DataFrame(contrasts).to_csv(reports / 'matched-comparison-anchor-contrasts.csv', index=False)
    pd.DataFrame(regions).to_csv(reports / 'matched-comparison-region-contrasts.csv', index=False)
    if ranks:
        pd.concat(ranks).to_csv(reports / 'matched-comparison-gene-rankings.csv', index=False)
        pd.concat(track_rows).to_csv(reports / 'matched-comparison-primary-track-scores.csv', index=False)
    available_regions = [r['region_median_contrast'] for r in regions if r['region_median_contrast'] is not None]
    summary = {'run_directory': str(folder), 'model_requests': len(run['requests']), 'complete_requests': sum(r['status'] == 'complete' for r in run['requests']),
               'evaluable_variants': sum(v['status'] == 'complete' for v in variants), 'evaluable_regions': len(available_regions),
               'equal_weight_cross_region_median_contrast': float(np.median(available_regions)) if available_regions else None,
               'input_sha256': run['input_sha256'], 'primary_gene_effect': protocol['analysis']['gene_effect'],
               'variant_effect': protocol['analysis']['variant_effect'], 'regions': regions,
               'summary_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'evaluated_variant_gene_pairs': sum(len(rank) for rank in ranks),
               'evaluated_primary_track_scores': sum(len(rows) for rows in track_rows),
               'interpretation_limit': 'Descriptive model-score contrasts against uncertain statistical-priority labels. No independent accuracy, disease-risk direction or treatment-effect inference.'}
    (reports / 'matched-comparison-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (reports / 'matched-comparison-run-provenance.json').write_text(json.dumps(run, indent=2) + '\n')
    if available_regions:
        cache = ROOT / 'work/matplotlib'
        cache.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault('MPLCONFIGDIR', str(cache))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        matplotlib.rcParams['svg.hashsalt'] = 'psc-matched-comparison-v0.1'
        fig, ax = plt.subplots(figsize=(10, 4.3))
        for index, row in enumerate(contrasts):
            if row['anchor_contrast'] is None:
                ax.text(index, .04, 'No valid\ncomparison pair', ha='center', va='bottom',
                        transform=ax.get_xaxis_transform(), color='#666666', fontsize=8)
                continue
            values = [row['anchor_effect'], row['comparator_1_effect'], row['comparator_2_effect']]
            ax.plot([index] * 3, values, color='#b5bdc5', linewidth=1, zorder=1)
            ax.scatter([index], values[:1], color='#1565a8', s=65, label='Prioritized variant' if index == 0 else None, zorder=3)
            ax.scatter([index-.10, index+.10], values[1:], color='#737d85', marker='s', s=38, label='Matched comparison variants' if index == 0 else None, zorder=3)
        ax.set_xticks(range(len(contrasts)), [r['archive_region']+'\n'+r['anchor'] for r in contrasts], fontsize=9)
        ax.set_ylabel('Maximum gene effect\n(median absolute RNA score across CD4 tracks)')
        ax.set_title('Fixed matched comparison: all five prioritized variants', loc='left', fontweight='bold')
        ax.set_ylim(bottom=0)
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y', alpha=.18)
        ax.legend(frameon=False, fontsize=9)
        fig.text(.01, .015, 'Model predictions; low-PIP comparators are not proven negatives. Regions are the reporting unit.', fontsize=8, color='#555555')
        fig.tight_layout(rect=(0,.05,1,1))
        for extension in ['png','svg']:
            output_path = reports / ('matched-comparison-effects.'+extension)
            fig.savefig(output_path, dpi=170, metadata={'Date': None} if extension == 'svg' else {})
            if extension == 'svg':
                output_path.write_text('\n'.join(line.rstrip() for line in output_path.read_text().splitlines()) + '\n')
        plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
