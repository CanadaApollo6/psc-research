"""Run the fixed small pilot; store predictions and run provenance locally."""

import argparse
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

from alphagenome_access import key_is_configured, load_key
from normalize_benchmark import verify_sources

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def request_for(row):
    interval = genome.Interval(row['chromosome'], row['interval_start_0based'], row['interval_end_0based_exclusive'])
    variant = genome.Variant(chromosome=row['chromosome'], position=row['position_grch38_1based'], reference_bases=row['reference_bases'], alternate_bases=row['alternate_bases'], name=row['rsid'])
    if interval.width != 1048576 or not interval.start <= variant.start < interval.end:
        raise ValueError('Unexpected input interval or coordinate conversion')
    return interval, variant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true', help='Make the four authenticated prediction requests')
    args = parser.parse_args()
    protocol_path = ROOT / 'config/alphagenome-pilot.json'
    input_path = ROOT / 'data/derived/benchmark-variants.json'
    protocol = json.loads(protocol_path.read_text())
    variants = json.loads(input_path.read_text())
    verify_sources()
    if alphagenome.__version__ != protocol['client_version']:
        raise ValueError('AlphaGenome client version differs from the fixed protocol')
    origin = json.loads(importlib.metadata.distribution('alphagenome').read_text('direct_url.json') or '{}')
    installed_commit = origin.get('vcs_info', {}).get('commit_id')
    if installed_commit != protocol['client_commit']:
        raise ValueError('Install the pinned client commit from requirements-alphagenome.lock')
    if len(variants) != 4:
        raise ValueError('This pilot is limited to the four selected variants')
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in protocol['scorers']]
    requests = [request_for(row) for row in variants]
    if not args.run:
        print(json.dumps({'status': 'dry-run; no predictions requested', 'key_configured': key_is_configured(), 'variants': [row['rsid'] for row in variants], 'scorers': [scorer.name for scorer in scorers], 'protocol_sha256': sha(protocol_path), 'input_sha256': sha(input_path)}, indent=2))
        return
    metadata_path = ROOT / 'data/predictions/metadata/inventory.json'
    if not metadata_path.is_file():
        raise ValueError('Fetch the model metadata before running predictions')
    started = datetime.now(timezone.utc).isoformat()
    identifier = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    folder = ROOT / 'data/predictions' / identifier
    folder.mkdir(parents=True, exist_ok=False)
    manifest = {
        'started_at_utc': started, 'client_version': alphagenome.__version__,
        'client_commit': installed_commit, 'model_version_selection': protocol['model_version'],
        'server_model_build_identifier': 'Not exposed by this client response; selection and retrieval time recorded',
        'protocol_sha256': sha(protocol_path), 'input_sha256': sha(input_path),
        'reference_manifest_sha256': sha(ROOT / 'config/reference-sources.json'),
        'metadata_inventory_sha256': sha(metadata_path), 'requests': [], 'status': 'running',
    }
    (folder / 'protocol.json').write_bytes(protocol_path.read_bytes())
    (folder / 'inputs.json').write_bytes(input_path.read_bytes())
    manifest_path = folder / 'run.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[protocol['model_version']], timeout=30)
    for row, (interval, variant) in zip(variants, requests):
        entry = {'rsid': row['rsid'], 'status': 'running'}
        manifest['requests'].append(entry)
        try:
            scores = client.score_variant(interval=interval, variant=variant, variant_scorers=scorers, organism=dna_client.Organism.HOMO_SAPIENS, merge_stranded_gene_tracks=True)
            table = variant_scorers.tidy_scores([scores], include_extended_metadata=True)
            if table is None or table.empty:
                raise ValueError('No scores returned')
            output = folder / f"{row['rsid']}-scores.csv.gz"
            table.to_csv(output, index=False, compression='gzip')
            entry.update(status='complete', rows=len(table), file=output.name, sha256=sha(output))
        except Exception as exc:
            code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
            entry.update(status='failed', error_code=code)
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        print(json.dumps(entry), flush=True)
    manifest['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
    manifest['status'] = 'complete' if all(item['status'] == 'complete' for item in manifest['requests']) else 'incomplete'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'run_directory': str(folder.relative_to(ROOT)), 'status': manifest['status']}), flush=True)
    if manifest['status'] != 'complete':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
