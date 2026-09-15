"""Bounded published-training-source crosswalk; absence never proves independence."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from prkd2_mechanism_common import ROOT, RAW, encode_request, sha, write_csv

SOURCE = ROOT / 'data/raw/qtl-alphagenome-supplement-tables.data'
SOURCE_SHA = '833cb78b6ae6fe39415cfff296ac00c48d800326139f13eb307531a1cc133154'
DEFAULT = ROOT / 'data/derived/prkd2-mechanism-training'


def accessions(value, prefix):
    return sorted(set(re.findall(prefix + r'[A-Z0-9]{6}', str(value))))


def specimen_ids(experiment):
    biosamples, donors = set(), set()
    unresolved = 0
    for rep in experiment.get('replicates', []):
        bs = rep.get('library', {}).get('biosample', {})
        if not isinstance(bs, dict):
            unresolved += 1
            continue
        biosamples.update(accessions(bs.get('accession', bs.get('@id', '')), 'ENCBS'))
        ds = accessions(bs.get('donor', {}), 'ENCDO')
        donors.update(ds)
        unresolved += not bool(ds)
    return sorted(biosamples), sorted(donors), unresolved


def run(output, offline=False):
    if sha(SOURCE) != SOURCE_SHA:
        raise ValueError('Published training workbook changed')
    with SOURCE.open('rb') as stream:
        tracks = pd.read_excel(stream, sheet_name='Suppl Table 2 Track metadata (f').fillna('')
    tracks.insert(0, 'workbook_row_1based', range(2, len(tracks) + 2))
    human = tracks[tracks.organism == 'human']
    contexts = human[human.ontology_curie.isin(['CL:0001054', 'CL:0000576']) &
                     human.output_type.isin(['RNA_SEQ', 'DNASE', 'ATAC', 'CAGE', 'CHIP_HISTONE'])]
    # All monocyte histone marks are included in this provenance audit, even
    # though only the three prospectively chosen marks are interpreted.
    exps = sorted({a for s in contexts['Experiment accession'] for a in accessions(s, 'ENCSR')})
    all_exps = {a for s in human['Experiment accession'] for a in accessions(s, 'ENCSR')}
    all_files = {a for s in human['File accession'] for a in accessions(s, 'ENCFF')}
    metadata, rows = {}, []
    for exp in exps:
        existing = RAW / 'encode' / (exp + '.json')
        label = exp if existing.exists() else 'training-' + exp
        p = RAW / 'encode' / (label + '.json')
        if offline:
            receipt_path = p.with_name(label + '-receipt.json')
            if not p.exists() or not receipt_path.exists():
                raise ValueError('Missing cached training metadata or receipt: ' + exp)
            receipt = json.loads(receipt_path.read_text())
            if not receipt['ok'] or sha(p) != receipt['sha256']:
                raise ValueError('Training source receipt does not verify')
            data = json.loads(p.read_text())
        else:
            data = encode_request(label, 'experiments/' + exp + '/', {'frame': 'embedded'})
        if not data:
            rows.append({'experiment': exp, 'biosample_ids': '', 'donor_ids': '', 'unresolved_replicates': '', 'status': 'metadata unavailable'})
            continue
        if data['accession'] != exp:
            raise ValueError('Experiment identity mismatch')
        bs, ds, unresolved = specimen_ids(data)
        rows.append({'experiment': exp, 'biosample_ids': ';'.join(bs), 'donor_ids': ';'.join(ds),
                     'unresolved_replicates': unresolved, 'status': 'resolved' if ds and not unresolved else 'partial'})
        metadata[str(p.relative_to(ROOT))] = sha(p)
    by_exp = {r['experiment']: r for r in rows}
    selection = json.loads((ROOT / 'config/prkd2-mechanism-encode-selection.json').read_text())
    overlap = []
    for s in selection['selected']:
        bs_hits, ds_hits = [], []
        for exp, row in by_exp.items():
            if set(s['biosample_ids']) & set(row['biosample_ids'].split(';')):
                bs_hits.append(exp)
            if set(s['donor_ids']) & set(row['donor_ids'].split(';')):
                ds_hits.append(exp)
        direct = s['experiment'] in all_exps
        overlap.append({'experiment': s['experiment'], 'file_accession': s['file_accession'], 'group': s['group'],
                        'source_donor_ids': ';'.join(s['donor_ids']), 'exact_experiment_in_all_human_training_rows': direct,
                        'exact_peak_file_in_all_human_training_rows': s['file_accession'] in all_files,
                        'same_biosample_training_experiments': ';'.join(bs_hits),
                        'same_donor_training_experiments': ';'.join(ds_hits),
                        'independence_status': ('published training experiment reused' if direct else
                                                'donor reused in published monocyte training sources' if ds_hits else
                                                'independence unresolved; no match in bounded donor audit')})
    output.mkdir(parents=True, exist_ok=True)
    contexts.to_csv(output / 'published-monocyte-training-tracks.csv', index=False, lineterminator='\n')
    write_csv(output / 'training-experiment-specimens.csv', rows)
    write_csv(output / 'selected-peak-training-overlap.csv', overlap)
    result = {'source_file': str(SOURCE.relative_to(ROOT)), 'source_sha256': SOURCE_SHA,
              'source_url': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-10014-0/MediaObjects/41586_2025_10014_MOESM3_ESM.xlsx',
              'published_rows': len(tracks), 'human_rows': len(human), 'selected_context_rows': len(contexts),
              'audited_training_experiments': exps, 'exact_training_experiment_matches': sum(r['exact_experiment_in_all_human_training_rows'] for r in overlap),
              'same_donor_matches': sum(bool(r['same_donor_training_experiments']) for r in overlap),
              'source_metadata_sha256': metadata,
              'cage_context_limitation': 'The generic monocyte CAGE track lists both CD14 monocytes and monocyte-derived endothelial progenitor cells, donors 1-3. Donor equivalence to ENCODE is unresolved.',
              'limits': ['Exact file/experiment membership screened against every published human training row.',
                         'Biosample/donor crosswalk bounded to the published primary/generic monocyte RNA, DNase, ATAC, CAGE and all histone rows.',
                         'All replicates of each training experiment screened conservatively; per-file contributor assignment may be narrower.',
                         'Published training metadata is not a receipt for the unexposed serving-model build. ALL_FOLDS is not held-out-locus validation.',
                         'No match is not evidence of independent donors or experimental validation of an allele effect.']}
    (output / 'training-audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['human_rows', 'exact_training_experiment_matches', 'same_donor_matches']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-directory', type=Path, default=DEFAULT)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    run(args.output_directory, args.offline)
