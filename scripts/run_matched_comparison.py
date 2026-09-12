"""Run only the locked complete matched groups; default is a credential-free dry run."""

import argparse
import hashlib
import importlib.metadata
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

from alphagenome_access import load_key

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    input_path = ROOT / 'config/psc-matched-comparison-inputs.json'
    locked = json.loads(input_path.read_text())
    protocol = json.loads((ROOT / 'config/psc-controlled-comparison.json').read_text())
    for name, sha in locked['frozen_file_sha256'].items():
        if digest(ROOT / name) != sha:
            raise ValueError('Frozen comparison file changed: ' + name)
    if not locked['ready_for_inference'] or not locked['inputs']:
        print(json.dumps({'status': 'no eligible matched groups', 'requests': 0}))
        return
    model = protocol['model']
    origin = json.loads(importlib.metadata.distribution('alphagenome').read_text('direct_url.json') or '{}')
    if alphagenome.__version__ != model['client_version'] or origin.get('vcs_info', {}).get('commit_id') != model['client_commit']:
        raise ValueError('Pinned AlphaGenome client is required')
    if not args.run:
        print(json.dumps({'status': 'dry-run', 'requests': len(locked['inputs']), 'input_sha256': digest(input_path), 'model': model['model_version']}))
        return
    try:
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[model['model_version']], timeout=45)
        frame = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS).rna_seq
        buffer = io.StringIO()
        frame.to_csv(buffer, index=False)
        if hashlib.sha256(buffer.getvalue().encode()).hexdigest() != locked['model_metadata_sha256']:
            raise ValueError('Live RNA track metadata differs from the frozen snapshot')
    except Exception as exc:
        code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
        raise SystemExit('Pre-inference metadata check failed: ' + code) from None
    folder = ROOT / 'data/predictions' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-matched-comparison')
    folder.mkdir(parents=True, exist_ok=False)
    (folder / 'inputs.json').write_bytes(input_path.read_bytes())
    (folder / 'protocol.json').write_bytes((ROOT / 'config/psc-controlled-comparison.json').read_bytes())
    run = {'started_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'running', 'input_sha256': digest(input_path),
           'protocol_sha256': digest(folder / 'protocol.json'), 'client_version': alphagenome.__version__,
           'client_commit': model['client_commit'], 'model_version_selection': model['model_version'],
           'server_model_build_identifier': 'Not exposed by this client; selection and request times recorded',
           'runner_sha256': digest(Path(__file__)), 'requests': []}
    for row in locked['inputs']:
        record = {'rsid': row['rsid'], 'archive_region': row['archive_region'], 'role': row['role'], 'status': 'running'}
        run['requests'].append(record)
        (folder / 'run.json').write_text(json.dumps(run, indent=2) + '\n')
        try:
            interval = genome.Interval(row['chromosome'], row['interval_start_0based'], row['interval_end_0based_exclusive'])
            variant = genome.Variant(row['chromosome'], row['position_grch38_1based'], row['reference_bases'], row['alternate_bases'], name=row['rsid'])
            scores = client.score_variant(interval=interval, variant=variant,
                                           variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
                                           organism=dna_client.Organism.HOMO_SAPIENS, merge_stranded_gene_tracks=True)
            table = variant_scorers.tidy_scores([scores], include_extended_metadata=True)
            if table is None or table.empty:
                raise ValueError('Missing scores')
            path = folder / (row['rsid'] + '-scores.csv.gz')
            table.to_csv(path, index=False, compression='gzip')
            record.update(status='complete', rows=len(table), file=path.name, sha256=digest(path))
        except Exception as exc:
            code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
            record.update(status='failed', error_code=code)
        (folder / 'run.json').write_text(json.dumps(run, indent=2) + '\n')
        print(json.dumps(record), flush=True)
    run.update(status='complete' if all(r['status'] == 'complete' for r in run['requests']) else 'incomplete', finished_at_utc=datetime.now(timezone.utc).isoformat())
    (folder / 'run.json').write_text(json.dumps(run, indent=2) + '\n')
    print(json.dumps({'run_directory': str(folder.relative_to(ROOT)), 'status': run['status']}), flush=True)


if __name__ == '__main__':
    main()
