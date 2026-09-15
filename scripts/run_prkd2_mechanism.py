"""Execute the locked PRKD2 mechanism predictions; all outputs are cached once."""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.metadata
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

from alphagenome_access import key_is_configured, load_key
from prkd2_mechanism_common import ROOT, PLAN, now, sha, save_new


def validate_inputs():
    plan = json.loads(PLAN.read_text())
    inputs_path = ROOT / 'config/prkd2-mechanism-inputs.json'
    inputs = json.loads(inputs_path.read_text())
    for file, expected in inputs['frozen_file_sha256'].items():
        if sha(ROOT / file) != expected:
            raise ValueError('Frozen input changed: ' + file)
    if not inputs['ready_for_inference'] or inputs['plan_sha256'] != sha(PLAN):
        raise ValueError('Model inputs not ready')
    encode_path = ROOT / 'config/prkd2-mechanism-encode-selection.json'
    encode = json.loads(encode_path.read_text())
    if encode['plan_sha256'] != sha(PLAN):
        raise ValueError('ENCODE selection used a different plan')
    for file, expected in encode['metadata_sha256'].items():
        if sha(ROOT / file) != expected:
            raise ValueError('ENCODE metadata changed before predictions')
    origin = json.loads(importlib.metadata.distribution('alphagenome').read_text('direct_url.json'))
    model = plan['model']
    if alphagenome.__version__ != model['client_version'] or origin.get('vcs_info', {}).get('commit_id') != model['client_commit']:
        raise ValueError('Pinned AlphaGenome client required')
    return plan, inputs, inputs_path, encode_path


def redacted_error(exc):
    return exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__


def verify_model_metadata(client):
    lock = json.loads((ROOT / 'config/prkd2-mechanism-model-metadata.json').read_text())
    live = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    for file in lock['files']:
        frame = getattr(live, file['modality'])
        if frame is None:
            raise ValueError('Frozen model output disappeared')
        buffer = io.StringIO()
        frame.to_csv(buffer, index=False)
        if hashlib.sha256(buffer.getvalue().encode()).hexdigest() != file['sha256']:
            raise ValueError('Live model metadata differs from frozen metadata')


def write_score_components(folder, rsid, scores):
    """Retain original numerical matrices and axes for independent CSV checks."""
    records = []
    for i, adata in enumerate(scores):
        prefix = f'{rsid}-component-{i:02d}'
        axes = {}
        for label, frame in [('genes', adata.obs), ('tracks', adata.var)]:
            path = folder / f'{prefix}-{label}.csv'
            frame.to_csv(path, index=False)
            axes[label] = {'file': path.name, 'sha256': sha(path), 'rows': len(frame)}
        arrays = {'scores': np.asarray(adata.X)}
        if 'quantiles' in adata.layers:
            arrays['quantiles'] = np.asarray(adata.layers['quantiles'])
        path = folder / (prefix + '.npz')
        np.savez_compressed(path, **arrays)
        records.append({'scorer': str(adata.uns['variant_scorer']), 'modality': adata.uns['variant_scorer'].requested_output.name,
                        'variant': str(adata.uns.get('variant')), 'interval': str(adata.uns.get('interval')),
                        'shape': list(arrays['scores'].shape), 'values_file': path.name, 'values_sha256': sha(path), **axes})
    return records


def save_tracks(folder, rsid, allele, output, modalities):
    records = []
    for modality in modalities:
        data = getattr(output, modality.lower())
        if data is None or data.interval is None or data.values.ndim != 2:
            raise ValueError('Requested track output missing or malformed')
        if not np.isfinite(data.values).all():
            raise ValueError('Nonfinite predicted track values')
        prefix = f'{rsid}-{allele}-{modality.lower()}'
        meta_path, values_path = folder / (prefix + '-metadata.csv'), folder / (prefix + '.npz')
        data.metadata.to_csv(meta_path, index=False)
        np.savez_compressed(values_path, values=data.values)
        records.append({'rsid': rsid, 'allele': allele, 'modality': modality.lower(),
                        'metadata_file': meta_path.name, 'metadata_sha256': sha(meta_path),
                        'values_file': values_path.name, 'values_sha256': sha(values_path),
                        'shape': list(data.values.shape), 'dtype': str(data.values.dtype),
                        'chromosome': data.interval.chromosome, 'start0': data.interval.start, 'end0': data.interval.end,
                        'resolution': data.resolution})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    plan, inputs, input_path, encode_path = validate_inputs()
    model = plan['model']
    anchor_rows = [r for r in inputs['variants'] if r['role'] == 'anchor']
    anchor_ontologies = [model['primary_ontology'], 'CL:0000576']
    metadata = inputs['selected_track_metadata']
    modalities = [m for m in ['RNA_SEQ', 'DNASE', 'ATAC', 'CAGE', 'CHIP_HISTONE']
                  if any(t['modality'] == m.lower() and t['ontology_curie'] in anchor_ontologies for t in metadata)]
    if not args.run:
        print(json.dumps({'status': 'dry-run; no prediction requests', 'key_configured': key_is_configured(),
                          'score_requests': len(inputs['variants']), 'anchor_track_requests': len(anchor_rows),
                          'anchor_track_modalities': modalities, 'scorers': model['scorers']}))
        return
    execution_path = ROOT / 'config/prkd2-mechanism-execution.json'
    if execution_path.exists():
        execution = json.loads(execution_path.read_text())
        if execution['input_sha256'] != sha(input_path) or execution['runner_sha256'] != sha(Path(__file__)):
            raise ValueError('Prediction execution inputs/code changed')
        folder = ROOT / execution['run_directory']
    else:
        folder = ROOT / 'data/predictions' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-prkd2-mechanism')
        folder.mkdir(parents=True, exist_ok=False)
        execution = {'created_before_predictions_at_utc': now(), 'plan_sha256': sha(PLAN),
                     'input_sha256': sha(input_path), 'encode_selection_sha256': sha(encode_path),
                     'runner_sha256': sha(Path(__file__)), 'run_directory': str(folder.relative_to(ROOT)),
                     'client_version': alphagenome.__version__, 'client_commit': model['client_commit'],
                     'model_selection': model['selection'], 'server_build_identifier': 'Not exposed by this client response',
                     'score_requests': len(inputs['variants']), 'anchor_track_requests': len(anchor_rows),
                     'anchor_track_modalities': modalities, 'anchor_track_ontologies': anchor_ontologies,
                     'raw_RNA_track_note': 'Raw predictions retain all returned strands. Gene-score API merges stranded names and drops duplicate unstranded names.',
                     'model_metadata_sha256': sha(ROOT/'config/prkd2-mechanism-model-metadata.json')}
        save_new(execution_path, execution)
        (folder/'inputs.json').write_bytes(input_path.read_bytes())
        (folder/'plan.json').write_bytes(PLAN.read_bytes())
    manifest_path = folder / 'run.json'
    run = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        'status': 'running', 'started_at_utc': now(), 'execution_sha256': sha(execution_path), 'scores': [], 'anchor_tracks': []}
    def checkpoint():
        manifest_path.write_text(json.dumps(run, indent=2, allow_nan=False) + '\n')
    checkpoint()
    try:
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[model['selection']], timeout=60)
        verify_model_metadata(client)
    except Exception as exc:
        run.update(status='metadata-preflight-failed', error_code=redacted_error(exc))
        checkpoint()
        raise SystemExit('Model metadata preflight failed: ' + redacted_error(exc)) from None
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in model['scorers']]
    for row in inputs['variants']:
        previous = next((r for r in run['scores'] if r['rsid'] == row['rsid']), None)
        if previous:
            if previous['status'] == 'complete' and sha(folder/previous['file']) == previous['sha256']:
                continue
            raise ValueError('Prior failed/incomplete score request requires an explicit execution repair; no silent retry')
        entry = {'rsid': row['rsid'], 'role': row['role'], 'status': 'running', 'started_at_utc': now()}
        run['scores'].append(entry)
        checkpoint()
        interval = genome.Interval('chr19', row['interval_start0'], row['interval_end0'])
        variant = genome.Variant('chr19', row['position_grch38_1based'], row['reference'], row['alternate'], name=row['rsid'])
        try:
            scores = client.score_variant(interval=interval, variant=variant, variant_scorers=scorers,
                                          organism=dna_client.Organism.HOMO_SAPIENS, merge_stranded_gene_tracks=True)
            components = write_score_components(folder, row['rsid'], scores)
            table = variant_scorers.tidy_scores(scores, include_extended_metadata=True)
            if table is None or table.empty or not np.isfinite(table.raw_score).all():
                raise ValueError('Missing or nonfinite prediction scores')
            path = folder / (row['rsid'] + '-scores.csv.gz')
            table.to_csv(path, index=False, compression='gzip')
            entry.update(status='complete', rows=len(table), file=path.name, sha256=sha(path), components=components)
        except Exception as exc:
            entry.update(status='failed', error_code=redacted_error(exc))
        entry['finished_at_utc'] = now()
        checkpoint()
        print(json.dumps({k: v for k, v in entry.items() if k != 'components'}), flush=True)
    for row in anchor_rows:
        previous = next((r for r in run['anchor_tracks'] if r['rsid'] == row['rsid']), None)
        if previous:
            if previous['status'] == 'complete':
                for output in previous['outputs']:
                    if sha(folder/output['values_file']) != output['values_sha256'] or sha(folder/output['metadata_file']) != output['metadata_sha256']:
                        raise ValueError('Cached anchor output changed')
                continue
            raise ValueError('Prior failed/incomplete track request requires explicit execution repair')
        entry = {'rsid': row['rsid'], 'status': 'running', 'started_at_utc': now()}
        run['anchor_tracks'].append(entry)
        checkpoint()
        try:
            result = client.predict_variant(
                interval=genome.Interval('chr19', row['interval_start0'], row['interval_end0']),
                variant=genome.Variant('chr19', row['position_grch38_1based'], row['reference'], row['alternate'], name=row['rsid']),
                organism=dna_client.Organism.HOMO_SAPIENS, ontology_terms=anchor_ontologies,
                requested_outputs=[dna_client.OutputType[name] for name in modalities])
            outputs = save_tracks(folder, row['rsid'], 'REF', result.reference, modalities)
            outputs += save_tracks(folder, row['rsid'], 'ALT', result.alternate, modalities)
            entry.update(status='complete', outputs=outputs)
        except Exception as exc:
            entry.update(status='failed', error_code=redacted_error(exc))
        entry['finished_at_utc'] = now()
        checkpoint()
        print(json.dumps({k: v for k, v in entry.items() if k != 'outputs'}), flush=True)
    complete = all(r['status'] == 'complete' for r in run['scores'] + run['anchor_tracks'])
    run.update(status='complete' if complete else 'incomplete', finished_at_utc=now())
    checkpoint()
    print(json.dumps({'status': run['status'], 'directory': str(folder.relative_to(ROOT))}))
    if not complete:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
