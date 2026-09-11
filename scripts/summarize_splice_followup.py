"""Summarize saved splice predictions, preserving the post-request coordinate correction."""

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fetch_sources import validate
from splice_context import donor_contexts, junction_in_scope, position_index

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aligned_delta(reference, alternate, keys):
    """Align biological features rather than assuming returned rows share order."""
    if reference.duplicated(keys).any() or alternate.duplicated(keys).any():
        raise ValueError('Duplicate feature keys would double count predictions')
    merged = reference.merge(alternate, on=keys, how='outer', validate='one_to_one', indicator=True)
    merged['delta_alt_minus_ref'] = merged['alternate'] - merged['reference']
    return merged.rename(columns={'_merge': 'prediction_availability'})


def read_values(folder, record):
    values = np.load(folder / record['values_file'], allow_pickle=False)['values']
    metadata = pd.read_csv(folder / record['metadata_file']).fillna('')
    if list(values.shape) != record['shape'] or values.shape[1] != len(metadata):
        raise ValueError('Saved data shape or metadata mismatch')
    if not np.isfinite(values).all():
        raise ValueError('Non-finite prediction')
    return values, metadata


def point_records(folder, record, positions):
    values, metadata = read_values(folder, record)
    result = []
    for label, position in positions.items():
        index = position_index(position, record['start_0based'], record['end_0based_exclusive'], record['resolution'])
        for track, meta in metadata.iterrows():
            result.append({
                'site': label, 'position_0based': position, 'position_1based': position + 1,
                'modality': record['modality'], 'track_name': meta['name'], 'strand': meta['strand'],
                'assay': meta.get('Assay title', ''), 'ontology_curie': meta.get('ontology_curie', ''),
                'reference' if record['allele'] == 'REF' else 'alternate': float(values[index, track]),
            })
    return pd.DataFrame(result)


def junction_records(folder, record, donor, acceptor):
    values, metadata = read_values(folder, record)
    coordinates = pd.read_csv(folder / record['coordinates_file'])
    if len(coordinates) != len(values):
        raise ValueError('Saved junction count mismatch')
    rows = []
    for index, junction in coordinates.iterrows():
        start, end = int(junction.start_0based), int(junction.end_0based_exclusive)
        if not junction_in_scope(start, end, junction.strand, donor, acceptor):
            continue
        for track, meta in metadata.iterrows():
            rows.append({
                'chromosome': junction.chromosome, 'start_0based': start, 'end_0based_exclusive': end,
                'strand': junction.strand, 'track_name': meta['name'], 'assay': meta['Assay title'],
                'ontology_curie': meta['ontology_curie'],
                'reference' if record['allele'] == 'REF' else 'alternate': float(values[index, track]),
            })
    return pd.DataFrame(rows)


def make_figure(points, destination):
    # Keep matplotlib's cache inside the ignored research workspace.
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'work/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['svg.hashsalt'] = 'psc-ubash3a-splice-followup'
    import matplotlib.pyplot as plt
    wanted = ['annotated_donor_corrected', 'largest_donor_gain_in_fixed_window']
    sites = points[(points.modality == 'splice_sites') & (points.track_name == 'donor') & (points.strand == '+')].set_index('site').loc[wanted]
    usage = points[(points.modality == 'splice_site_usage') & (points.strand == '+') & points.site.isin(wanted)].copy()
    usage['order'] = usage.site.map({wanted[0]: 0, wanted[1]: 1})
    usage = usage.sort_values(['order', 'assay'])
    colors = ['#274C77', '#C76732']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.7), gridspec_kw={'width_ratios': [1, 1.55]})
    for ax, table in [(axes[0], sites), (axes[1], usage)]:
        x = np.arange(len(table))
        ax.bar(x - .18, table.reference, width=.36, color=colors[0], label='Reference A')
        ax.bar(x + .18, table.alternate, width=.36, color=colors[1], label='Alternate C')
        ax.set_xticks(x)
        ax.set_ylim(0, 1.06)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.17)
        ax.set_axisbelow(True)
        if ax is axes[1]:
            for index, value in enumerate(table.alternate):
                if value < .02:
                    ax.text(index + .18, value + .022, f'{value:.4f}', ha='center', fontsize=8, color=colors[1])
    axes[0].set_xticklabels(['Annotated\ndonor', 'Alternative donor\n29 nt downstream'])
    axes[0].set_ylabel('Predicted splice-site probability')
    axes[0].set_title('A  |  Donor recognition', loc='left', fontsize=12)
    labels = [('Annotated' if row.site == wanted[0] else '+29 nt') + '\n' + ('PolyA RNA' if row.assay.startswith('polyA') else 'Total RNA') for row in usage.itertuples()]
    axes[1].set_xticklabels(labels, fontsize=9)
    axes[1].set_ylabel('Predicted splice-site usage')
    axes[1].set_title('B  |  CD4 T-cell assay tracks', loc='left', fontsize=12)
    fig.suptitle('UBASH3A rs1893592: a predicted shift in splice-donor use', fontsize=14, x=.07, ha='left')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(.98,.91), ncol=2, frameon=False, fontsize=9)
    fig.text(.07, .015, 'Exploratory model outputs, not measured RNA fractions. Coordinate correction is documented in the audit.\nAlphaGenome ALL_FOLDS; one A>C request; assay tracks are not independent donors.', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.075, right=.985, top=.75, bottom=.22, wspace=.35)
    fig.savefig(destination.with_suffix('.png'), dpi=180)
    svg = destination.with_suffix('.svg')
    fig.savefig(svg, metadata={'Date': None})
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    folder = args.run_directory.resolve()
    run = json.loads((folder / 'run.json').read_text())
    if run['status'] != 'complete' or sha(folder / 'protocol.json') != run['protocol_sha256']:
        raise ValueError('Incomplete run or changed request protocol')
    for name, expected in run['output_hashes'].items():
        if sha(folder / name) != expected:
            raise ValueError(f'Output checksum mismatch: {name}')
    original = json.loads((folder / 'protocol.json').read_text())
    corrected = json.loads((ROOT / 'config/ubash3a-splice-followup.json').read_text())
    if original['protocol_version'] not in {'0.1', '0.2'}:
        raise ValueError('Review the analysis for this protocol version')
    if sha(ROOT / 'config/mechanism-sources.json') != run['mechanism_sources_sha256_at_request']:
        raise ValueError('Request-time source manifest has changed')
    for key in ['variant', 'interval', 'requested_outputs', 'ontology_terms', 'client_commit', 'model_version']:
        if original[key] != corrected[key]:
            raise ValueError('Correction must not change the prediction request')
    for manifest_name in ['mechanism-sources.json', 'splice-audit-sources.json']:
        for source in json.loads((ROOT / 'config' / manifest_name).read_text())['sources']:
            validate((ROOT / 'data/raw' / source['file']).read_bytes(), source)
    gene = json.loads((ROOT / 'data/raw/UBASH3A-transcripts-GRCh38.json').read_text())
    contexts = donor_contexts(gene, original['variant']['position_1based'])
    pd.DataFrame(contexts).to_csv(ROOT / 'data/derived/ubash3a-transcript-boundaries.csv', index=False)
    records = {(r['allele'], r['modality']): r for r in run['outputs']}
    ref, meta = read_values(folder, records['REF', 'splice_sites'])
    alt, alt_meta = read_values(folder, records['ALT', 'splice_sites'])
    if not meta.equals(alt_meta):
        raise ValueError('Cannot scan donor window without aligned track metadata')
    donor_track = meta.index[(meta['name'] == 'donor') & (meta.strand == '+')].tolist()
    if len(donor_track) != 1:
        raise ValueError('Missing or duplicated positive donor track')
    track = donor_track[0]
    start = records['REF', 'splice_sites']['start_0based']
    window = original['diagnostic_window']
    lo, hi = window['start_0based'] - start, window['end_0based_exclusive'] - start
    gain_position = int(np.argmax(alt[lo:hi, track] - ref[lo:hi, track]) + lo + start)
    positions = {
        ('original_v0.1_index_error' if original['protocol_version'] == '0.1' else 'request_endpoint_v0.2'): original['primary_endpoint']['position_0based'],
        'annotated_donor_corrected': corrected['primary_endpoint']['position_0based'],
        'largest_donor_gain_in_fixed_window': gain_position,
    }
    keys = ['site', 'position_0based', 'position_1based', 'modality', 'track_name', 'strand', 'assay', 'ontology_curie']
    points = []
    windows = []
    for modality in ['splice_sites', 'splice_site_usage']:
        tables = [point_records(folder, records[label, modality], positions) for label in ['REF', 'ALT']]
        points.append(aligned_delta(*tables, keys))
        all_positions = {f'offset_{i-corrected["primary_endpoint"]["position_0based"]:+d}': i for i in range(window['start_0based'], window['end_0based_exclusive'])}
        tables = [point_records(folder, records[label, modality], all_positions) for label in ['REF', 'ALT']]
        windows.append(aligned_delta(*tables, keys))
    points = pd.concat(points, ignore_index=True)
    points.to_csv(ROOT / 'reports/ubash3a-splice-site-comparison.csv', index=False)
    pd.concat(windows, ignore_index=True).to_csv(ROOT / 'reports/ubash3a-splice-local-window.csv', index=False)
    donor = original['annotation']['upstream_exon_end_1based']
    acceptor = original['annotation']['downstream_exon_start_1based'] - 1
    tables = [junction_records(folder, records[label, 'splice_junctions'], donor, acceptor) for label in ['REF', 'ALT']]
    junctions = aligned_delta(*tables, ['chromosome', 'start_0based', 'end_0based_exclusive', 'strand', 'track_name', 'assay', 'ontology_curie'])
    junctions.to_csv(ROOT / 'reports/ubash3a-splice-junction-comparison.csv', index=False)
    summary = {
        'run_directory': str(folder.relative_to(ROOT)), 'run_manifest_sha256': sha(folder / 'run.json'),
        'original_protocol_sha256': run['protocol_sha256'], 'corrected_protocol_sha256': sha(ROOT / 'config/ubash3a-splice-followup.json'),
        'correction': corrected['amendment'] if original['protocol_version'] == '0.1' else 'This run used the corrected v0.2 coordinate from the outset; the earlier project amendment is recorded in its protocol.',
        'coordinate_source': corrected['coordinate_source'],
        'positions_0based': positions, 'alternative_donor_distance_nt': gain_position - positions['annotated_donor_corrected'],
        'selection_of_alternative_donor': 'Largest positive donor-probability difference in the original fixed 101-base diagnostic window; exploratory endpoint',
        'transcripts_at_boundary': len(contexts), 'junction_track_comparisons': len(junctions),
        'source_manifests': {name: sha(ROOT / 'config' / name) for name in ['mechanism-sources.json', 'splice-audit-sources.json', 'association-sources.json']},
        'request_provenance': run,
    }
    (ROOT / 'reports/splice-followup-provenance.json').write_text(json.dumps(summary, indent=2) + '\n')
    make_figure(points, ROOT / 'reports/ubash3a-splice-followup')
    print(json.dumps({'site_comparisons': len(points), 'junction_comparisons': len(junctions), 'alternative_donor_distance_nt': summary['alternative_donor_distance_nt'], 'status': 'summarized'}))


if __name__ == '__main__':
    main()
