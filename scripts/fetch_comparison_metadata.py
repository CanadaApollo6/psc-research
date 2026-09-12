"""Verify or retrieve the pinned RNA track metadata; never request predictions."""

import argparse
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    record = json.loads((ROOT / 'config/psc-comparator-model-metadata.json').read_text())
    protocol = json.loads((ROOT / 'config/psc-controlled-comparison.json').read_text())
    path = ROOT / 'data/raw' / record['file']
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
            raise SystemExit('Cached RNA metadata differs from the frozen snapshot')
        print(json.dumps({'status': 'verified-cached-metadata', 'prediction_requests': 0}))
        return
    if args.offline:
        raise SystemExit('Pinned RNA metadata is absent; an API metadata request is needed')
    import alphagenome
    from alphagenome.models import dna_client
    from alphagenome_access import load_key
    model = protocol['model']
    origin = json.loads(importlib.metadata.distribution('alphagenome').read_text('direct_url.json') or '{}')
    if alphagenome.__version__ != model['client_version'] or origin.get('vcs_info', {}).get('commit_id') != model['client_commit']:
        raise SystemExit('Pinned AlphaGenome client is required')
    try:
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[model['model_version']], timeout=45)
        frame = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS).rna_seq
        buffer = io.StringIO()
        frame.to_csv(buffer, index=False)
        data = buffer.getvalue().encode()
        if hashlib.sha256(data).hexdigest() != record['sha256']:
            raise ValueError('Metadata differs')
    except Exception as exc:
        code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
        raise SystemExit('RNA metadata retrieval or frozen-hash check failed: ' + code) from None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    print(json.dumps({'status': 'retrieved-pinned-metadata', 'prediction_requests': 0}))


if __name__ == '__main__':
    main()
