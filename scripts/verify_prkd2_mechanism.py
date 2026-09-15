"""Independent reconstruction from original matrices, source DNA and whole BEDs."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.ipc as ipc

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'data/predictions/20260915T185853Z-prkd2-mechanism'
RESULTS = ROOT / 'data/derived/prkd2-mechanism-results'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads((ROOT / path).read_text())


def check_hash(path, expected):
    if digest(path) != expected:
        raise AssertionError('Checksum mismatch: ' + str(path))


def main():
    inputs = load('config/prkd2-mechanism-inputs.json')
    execution = load('config/prkd2-mechanism-execution.json')
    freeze = load('config/prkd2-mechanism-freeze.json')
    run = json.loads((RUN / 'run.json').read_text())
    assert run['status'] == 'complete' and len(run['scores']) == 3 and len(run['anchor_tracks']) == 2
    for path, key in [('config/prkd2-mechanism-inputs.json', 'input_sha256'), ('config/prkd2-mechanism-plan.json', 'plan_sha256'),
                      ('config/prkd2-mechanism-encode-selection.json', 'encode_selection_sha256'),
                      ('scripts/run_prkd2_mechanism.py', 'runner_sha256')]:
        check_hash(ROOT/path, execution[key])
    for path, value in inputs['frozen_file_sha256'].items(): check_hash(ROOT/path, value)
    assert freeze['frozen_at_utc'] < inputs['locked_at_utc'] < execution['created_before_predictions_at_utc'] < run['started_at_utc']
    # Reconstruct the exact decoded model reference and compare the separate source.
    refs = load('config/prkd2-mechanism-reference-sources.json')
    for r in refs['sources']: check_hash(ROOT/r['file'], r['sha256'])
    model_ref = next(r for r in refs['sources'] if r['id'] == 'model-GRCh38.p13-common-window.range')
    model_seq = b''.join((ROOT/model_ref['file']).read_bytes().split()).upper()
    audit = load('data/derived/prkd2-mechanism/reference-audit.json')
    old = audit['ensembl_reference_source']
    check_hash(ROOT/old['file'], old['sha256'])
    ensembl_seq = (ROOT/old['file']).read_text().strip().upper()
    offset = audit['start0'] - (old['start_1based'] - 1)
    assert model_seq.decode() == ensembl_seq[offset:offset+audit['width']]
    assert hashlib.sha256(model_seq).hexdigest() == audit['sequence_sha256']
    for v in inputs['variants']:
        i = v['position_grch38_1based'] - 1 - audit['start0']
        assert model_seq[i:i+1].decode() == v['reference']
        assert model_seq[i-1:i+2].decode() == v['triplet']
        assert v['alternate'] != v['reference'] and len(v['alternate']) == 1
    # Rebuild gene membership from the source annotations, not the derived list.
    gene_json = next(r for r in refs['sources'] if r['id'] == 'ensembl-common-window-genes.json')
    ensembl_ids = set()
    for gene in json.loads((ROOT/gene_json['file']).read_text()):
        tss1 = gene['start'] if gene['strand'] == 1 else gene['end']
        if gene['seq_region_name'] == '19' and audit['start0'] < tss1 <= audit['end0']:
            ensembl_ids.add(gene['gene_id'].split('.')[0])
    model_ids = set()
    gencode = ROOT/'data/raw'/audit['model_gencode_source']['file']
    check_hash(gencode, audit['model_gencode_source']['sha256'])
    with pa.memory_map(str(gencode),'r') as stream:
        reader = ipc.open_file(stream)
        for i in range(reader.num_record_batches):
            batch = reader.get_batch(i).select(['Chromosome','Feature','Start','End','Strand','gene_id']).to_pandas()
            tss0 = np.where(batch.Strand == '+', batch.Start, batch.End)
            use = (batch.Chromosome == 'chr19') & (batch.Feature == 'transcript') & (tss0 >= audit['start0']) & (tss0 < audit['end0'])
            model_ids.update(g.split('.')[0] for g in batch.loc[use,'gene_id'])
    assert ensembl_ids & model_ids == {g['gene_id'] for g in inputs['fixed_gene_universe']}
    # Independent join from each original gene x track matrix, with float32
    # round-trip equality. Do not invoke tidy_scores or the summary script.
    numerical_cells, quantile_cells, raw_arrays = 0, 0, 0
    primary_values = {}
    for entry in run['scores']:
        check_hash(RUN/entry['file'], entry['sha256'])
        tidy = pd.read_csv(RUN/entry['file'], low_memory=False).fillna('')
        expected_rows = 0
        for component in entry['components']:
            check_hash(RUN/component['values_file'], component['values_sha256'])
            for axis in ['genes', 'tracks']:
                check_hash(RUN/component[axis]['file'], component[axis]['sha256'])
            z = np.load(RUN/component['values_file'])
            values = z['scores']
            assert list(values.shape) == component['shape'] and np.isfinite(values).all()
            track = pd.read_csv(RUN/component['tracks']['file']).fillna('')
            assert not track.duplicated(['name', 'strand']).any()
            if values.shape[0] == 0:
                assert not (tidy.variant_scorer == component['scorer']).any()
                continue
            try: genes = pd.read_csv(RUN/component['genes']['file']).fillna('')
            except pd.errors.EmptyDataError: genes = pd.DataFrame({'gene_id': ['']})
            gids = [str(g).split('.')[0] for g in genes.gene_id] if 'gene_id' in genes else ['']
            assert len(gids) == values.shape[0]
            keys = [(gid, t['name'], t['strand']) for gid in gids for t in track.to_dict('records')]
            subset = tidy[tidy.variant_scorer == component['scorer']].set_index(['gene_id','track_name','track_strand'])
            assert subset.index.is_unique and len(subset) == len(keys)
            matched = subset.loc[keys]
            assert np.array_equal(matched.raw_score.to_numpy(dtype=np.float32), values.ravel().astype(np.float32))
            numerical_cells += values.size
            expected_rows += values.size
            if 'quantiles' in z:
                assert np.array_equal(matched.quantile_score.to_numpy(dtype=np.float32), z['quantiles'].ravel().astype(np.float32), equal_nan=True)
                quantile_cells += values.size
            if component['scorer'] == 'GeneMaskLFCScorer(requested_output=RNA_SEQ)':
                ti = track.index[track.name == 'CL:0001054 polyA plus RNA-seq'].tolist()
                assert len(ti) == 1
                primary_values[entry['rsid']] = dict(zip(gids, values[:, ti[0]].astype(float)))
        assert expected_rows == len(tidy) == entry['rows']
    for entry in run['anchor_tracks']:
        for output in entry['outputs']:
            check_hash(RUN/output['values_file'], output['values_sha256'])
            check_hash(RUN/output['metadata_file'], output['metadata_sha256'])
            values = np.load(RUN/output['values_file'])['values']
            assert list(values.shape) == output['shape'] and np.isfinite(values).all()
            assert values.shape[0] * output['resolution'] == audit['width']
            raw_arrays += 1
    ranking = pd.read_csv(RESULTS/'primary-gene-ranks.csv')
    universe = {g['gene_id'] for g in inputs['fixed_gene_universe']}
    for rsid, source in primary_values.items():
        target = ranking[ranking.rsid == rsid]
        assert set(target.gene_id) == universe and len(target) == 50
        ranked = sorted(universe, key=lambda gid: (-abs(source[gid]), gid))
        for row in target.itertuples():
            assert np.isclose(row.signed_rna_median, source[row.gene_id], rtol=1e-7, atol=1e-12)
            assert row.model_rank == ranked.index(row.gene_id) + 1
    # Parse every original peak row independently. Bounds convert directly
    # to one-based inclusive coordinates, distinct from the production parser.
    peaks = load('data/derived/prkd2-mechanism-encode/peak-analysis.json')
    summary = pd.read_csv(ROOT/'data/derived/prkd2-mechanism-encode/peak-overlap-summary.csv')
    total_bed_rows, comparisons = 0, 0
    for receipt in peaks['files']:
        path = ROOT/receipt['file']
        check_hash(path, receipt['sha256'])
        with path.open('rb') as stream: md5 = hashlib.file_digest(stream,'md5').hexdigest()
        assert md5 == receipt['publisher_md5']
        counters = {p['label']: [0,0,None] for p in peaks['points']}
        count = 0
        with gzip.open(path,'rt') as stream:
            for line in stream:
                if not line.strip() or line.startswith(('#','track ','browser ')): continue
                cells = line.split('\t')
                left1, right1 = int(cells[1])+1, int(cells[2])
                assert 1 <= left1 <= right1
                count += 1
                if cells[0] != 'chr19' or right1 <= audit['start0'] or left1 > audit['end0']: continue
                for p in peaks['points']:
                    q = p['point1']
                    distance = max(left1-q, q-right1, 0)
                    c = counters[p['label']]
                    c[0] += left1 <= q <= right1
                    c[1] += left1 <= q+250 and right1 >= q-250
                    c[2] = distance if c[2] is None else min(c[2], distance)
        assert count == receipt['total_source_rows']
        total_bed_rows += count
        for label, (exact, nearby, nearest) in counters.items():
            row = summary[(summary.file_accession == receipt['file_accession']) & (summary.label == label)].iloc[0]
            assert exact == row.exact_base_peak_count and nearby == row.within_250bp_peak_count
            assert nearest == row.nearest_regional_peak_distance_bp
            comparisons += 1
    # External LD from explicit centered sums, independent of np.corrcoef.
    source = load('data/derived/prkd2-mechanism-results/reference-LD-provenance.json')['source']
    wanted = {str(v['position_grch37_1based']): v for v in inputs['variants']}
    dosages = {}
    with gzip.open(ROOT/source['file'],'rt') as stream:
        reader = csv.reader(stream,delimiter='\t'); next(reader)
        for row in reader:
            if row[1] in wanted:
                v = wanted[row[1]]
                assert row[0] == '19' and row[3:5] == [v['reference'],v['alternate']]
                dosages[v['rsid']] = np.asarray(row[5:],dtype=float)
    ld = pd.read_csv(RESULTS/'reference-LD.csv')
    for row in ld.itertuples():
        a,b = dosages[row.variant_a],dosages[row.variant_b]
        assert len(a) == len(b) == row.joint_donors == 503
        a,b = a-a.mean(),b-b.mean()
        r = (a@b)/np.sqrt((a@a)*(b@b))
        assert np.isclose(r,row.r_alt_dosage,rtol=1e-12,atol=1e-12)
        assert np.isclose(r*r,row.r_squared,rtol=1e-12,atol=1e-12)
    allowed = {'README.md','docs/evidence-log.md','docs/research-plan.md'}
    changed = [p for p,h in freeze['prior_tracked_sha256'].items() if digest(ROOT/p) != h]
    assert set(changed) <= allowed, 'An earlier scientific output changed: ' + str(changed)
    result = {'status': 'PASS', 'model_numerical_cells_exact_float32': numerical_cells, 'quantile_cells_exact_float32': quantile_cells,
              'raw_prediction_arrays_verified': raw_arrays, 'fixed_gene_rank_rows': len(ranking),
              'whole_BED_rows_independently_parsed': total_bed_rows, 'independent_point_file_checks': comparisons,
              'reference_LD_cells_independently_verified': len(ld), 'reference_sequence_bases': len(model_seq),
              'source_edit_identities_checked': len(inputs['variants']), 'source_gene_universe_reconstructed': len(ensembl_ids & model_ids),
              'prior_tracked_files': len(freeze['prior_tracked_sha256']), 'intentional_navigation_updates': changed,
              'verifier_sha256': digest(Path(__file__)), 'model_run_sha256': digest(RUN/'run.json'),
              'limits': 'Checks numerical transport, joins, coordinates and deterministic summaries, not biological causality or model accuracy.'}
    (ROOT/'reports/prkd2-mechanism-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
