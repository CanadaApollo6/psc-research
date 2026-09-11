"""Make one prespecified UBASH3A variant request; keep full outputs out of Git."""

import argparse
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client
import numpy as np
import pandas as pd

from alphagenome_access import key_is_configured, load_key
from fetch_sources import validate
from run_alphagenome_pilot import sha
from splice_context import donor_contexts, position_index

ROOT = Path(__file__).resolve().parents[1]


def verify_inputs(protocol):
    for filename in ['reference-sources.json', 'mechanism-sources.json']:
        manifest = json.loads((ROOT / 'config' / filename).read_text())
        for source in manifest['sources']:
            validate((ROOT / 'data/raw' / source['file']).read_bytes(), source)
    gene = json.loads((ROOT / 'data/raw' / protocol['annotation']['source']).read_text())
    contexts = donor_contexts(gene, protocol['variant']['position_1based'])
    donor = protocol['primary_endpoint']['position_0based']
    expected = [row for row in contexts if row['model_donor_position_0based'] == donor and row['intron_offset_1based'] == 3]
    if not expected or not any(row['is_canonical'] for row in expected):
        raise ValueError('Prespecified donor is not supported by the canonical annotation')
    if any(row['downstream_exon_start_1based'] != protocol['annotation']['downstream_exon_start_1based'] for row in expected):
        raise ValueError('Downstream exon boundary changed')
    if (ROOT / 'data/raw/model-reference-rs1893592.txt').read_text().strip() != protocol['variant']['reference']:
        raise ValueError('Reference allele disagrees with the model reference')
    origin = json.loads(importlib.metadata.distribution('alphagenome').read_text('direct_url.json') or '{}')
    if alphagenome.__version__ != protocol['client_version'] or origin.get('vcs_info', {}).get('commit_id') != protocol['client_commit']:
        raise ValueError('Install the recorded AlphaGenome client commit')
    return contexts


def save_output(folder, label, output):
    records = []
    for modality in ['splice_sites', 'splice_site_usage', 'splice_junctions']:
        data = getattr(output, modality)
        if data is None or data.interval is None:
            raise ValueError(f'Missing {modality} data or coordinates')
        prefix = f'{label}-{modality}'
        data.metadata.to_csv(folder / f'{prefix}-metadata.csv', index=False)
        np.savez_compressed(folder / f'{prefix}.npz', values=data.values)
        row = {'allele': label, 'modality': modality, 'values_file': f'{prefix}.npz', 'metadata_file': f'{prefix}-metadata.csv',
               'shape': list(data.values.shape), 'chromosome': data.interval.chromosome,
               'start_0based': data.interval.start, 'end_0based_exclusive': data.interval.end}
        if modality == 'splice_junctions':
            coordinates = [{'chromosome': j.chromosome, 'start_0based': j.start, 'end_0based_exclusive': j.end, 'strand': j.strand} for j in data.junctions]
            name = f'{prefix}-coordinates.csv.gz'
            pd.DataFrame(coordinates).to_csv(folder / name, index=False, compression='gzip')
            row['coordinates_file'] = name
        else:
            row['resolution'] = data.resolution
            if data.resolution != 1:
                raise ValueError('Expected nucleotide-resolution splicing output')
        records.append(row)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    path = ROOT / 'config/ubash3a-splice-followup.json'
    protocol = json.loads(path.read_text())
    contexts = verify_inputs(protocol)
    v = protocol['variant']
    interval = genome.Interval(v['chromosome'], protocol['interval']['start_0based'], protocol['interval']['end_0based_exclusive'])
    position_index(protocol['primary_endpoint']['position_0based'], interval.start, interval.end)
    if not args.run:
        print(json.dumps({'status': 'dry-run; no predictions requested', 'key_configured': key_is_configured(), 'requests': 1,
                          'protocol_sha256': sha(path), 'annotated_transcripts': len(contexts)}, indent=2))
        return
    folder = ROOT / 'data/predictions' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-ubash3a-splice')
    folder.mkdir(parents=True, exist_ok=False)
    (folder / 'protocol.json').write_bytes(path.read_bytes())
    manifest = {'status': 'running', 'started_at_utc': datetime.now(timezone.utc).isoformat(), 'request_count': 1,
                'client_version': alphagenome.__version__, 'client_commit': protocol['client_commit'],
                'model_version_selection': protocol['model_version'], 'server_model_build_identifier': 'Not exposed by client response',
                'protocol_sha256': sha(path), 'runner_sha256': sha(Path(__file__)),
                'mechanism_sources_sha256_at_request': sha(ROOT / 'config/mechanism-sources.json'),
                'reference_sources_sha256': sha(ROOT / 'config/reference-sources.json'), 'outputs': []}
    manifest_path = folder / 'run.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    try:
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[protocol['model_version']], timeout=60)
        result = client.predict_variant(
            interval=interval,
            variant=genome.Variant(v['chromosome'], v['position_1based'], v['reference'], v['alternate'], name=v['rsid']),
            organism=dna_client.Organism.HOMO_SAPIENS,
            requested_outputs=[dna_client.OutputType[name] for name in protocol['requested_outputs']],
            ontology_terms=protocol['ontology_terms'],
        )
        manifest['outputs'] = save_output(folder, 'REF', result.reference) + save_output(folder, 'ALT', result.alternate)
        manifest['output_hashes'] = {p.name: sha(p) for p in sorted(folder.iterdir()) if p.name != 'run.json'}
        manifest['status'] = 'complete'
    except Exception as exc:
        # Do not serialize arbitrary RPC errors: provider metadata might include credentials.
        manifest['status'] = 'failed'
        manifest['error_code'] = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
    manifest['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'status': manifest['status'], 'run_directory': str(folder.relative_to(ROOT)), 'error_code': manifest.get('error_code')}), flush=True)
    if manifest['status'] != 'complete':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
