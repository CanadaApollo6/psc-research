"""Freeze this follow-up and fetch metadata without requesting predictions."""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.metadata
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/prkd2-mechanism'
PLAN = ROOT / 'config/prkd2-mechanism-plan.json'
SKILL = Path('/home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/encode-skill/scripts/rest_request.py')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def save_new(path, value):
    path = Path(path)
    data = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError('Refusing to replace an existing source/freeze: ' + str(path))
    else:
        path.write_bytes(data)


def freeze():
    path = ROOT / 'config/prkd2-mechanism-freeze.json'
    if path.exists():
        record = json.loads(path.read_text())
        if record['plan_sha256'] != sha(PLAN):
            raise ValueError('Frozen plan changed')
        return record
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    baseline = {name: sha(ROOT / name) for name in names if name}
    record = {'frozen_at_utc': now(), 'plan_sha256': sha(PLAN),
              'baseline_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(),
              'prior_tracked_sha256': baseline,
              'new_prediction_or_peak_overlap_inspected': False,
              'known_prior_results': 'See plan prior_knowledge; this is exploratory follow-up, not preregistration.'}
    save_new(path, record)
    return record


def metadata():
    freeze()
    path = ROOT / 'config/prkd2-mechanism-model-metadata.json'
    if path.exists():
        record = json.loads(path.read_text())
        for row in record['files']:
            if sha(ROOT / row['file']) != row['sha256']:
                raise ValueError('Metadata checksum mismatch')
        print(json.dumps({'status': 'verified cached metadata', 'files': len(record['files'])}))
        return
    import alphagenome
    from alphagenome.models import dna_client
    from alphagenome_access import load_key
    plan = json.loads(PLAN.read_text())['model']
    dist = importlib.metadata.distribution('alphagenome')
    origin = json.loads(dist.read_text('direct_url.json'))
    if alphagenome.__version__ != plan['client_version'] or origin.get('vcs_info', {}).get('commit_id') != plan['client_commit']:
        raise ValueError('Pinned AlphaGenome client required')
    try:
        client = dna_client.create(load_key(), model_version=dna_client.ModelVersion[plan['selection']], timeout=45)
        result = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    except Exception as exc:
        code = exc.code().name if callable(getattr(exc, 'code', None)) else type(exc).__name__
        raise SystemExit('Metadata request failed: ' + code) from None
    out = RAW / 'model-metadata'
    out.mkdir(parents=True, exist_ok=True)
    files = []
    for field in dataclasses.fields(result):
        frame = getattr(result, field.name)
        if frame is None:
            continue
        dest = out / (field.name + '.csv')
        if dest.exists():
            raise ValueError('Unpinned metadata output already exists')
        frame.to_csv(dest, index=False)
        files.append({'modality': field.name, 'file': str(dest.relative_to(ROOT)),
                      'sha256': sha(dest), 'rows': len(frame), 'columns': list(frame.columns)})
    record = {'retrieved_at_utc': now(), 'model': plan['selection'],
              'client_version': alphagenome.__version__, 'client_commit': plan['client_commit'],
              'server_build': 'Not exposed by this response', 'prediction_requests': 0, 'files': files}
    save_new(path, record)
    print(json.dumps(record))


def encode_search():
    freeze()
    folder = RAW / 'encode'
    folder.mkdir(parents=True, exist_ok=True)
    queries = [('dnase', {'assay_title': 'DNase-seq'}), ('atac', {'assay_title': 'ATAC-seq'}),
               ('histone', {'assay_title': 'Histone ChIP-seq'})]
    for label, extra in queries:
        req = {'base_url': 'https://www.encodeproject.org', 'path': 'search/',
               'params': {'type': 'Experiment', 'status': 'released',
                          'biosample_ontology.term_name': 'CD14-positive monocyte',
                          'limit': 100, 'format': 'json', **extra},
               'headers': {'Accept': 'application/json'}, 'record_path': '@graph',
               'max_items': 3, 'max_depth': 2, 'timeout_sec': 45, 'save_raw': True,
               'raw_output_path': str(folder / f'search-{label}.json')}
        request_path = ROOT / f'config/prkd2-mechanism-encode-{label}-request.json'
        save_new(request_path, req)
        receipt_path = folder / f'search-{label}-receipt.json'
        if receipt_path.exists():
            prior = json.loads(receipt_path.read_text())
            if prior['result'].get('ok'):
                print(json.dumps({'query': label, 'status': 'successful receipt already exists'}))
                continue
            attempt = 1
            while receipt_path.exists():
                if json.loads(receipt_path.read_text())['result'].get('ok'):
                    break
                receipt_path = folder / f'search-{label}-receipt-retry{attempt}.json'
                attempt += 1
            if receipt_path.exists():
                print(json.dumps({'query': label, 'status': 'successful retry already exists'}))
                continue
        start = now()
        p = subprocess.run([sys.executable, str(SKILL)], input=json.dumps(req), text=True, capture_output=True)
        try:
            result = json.loads(p.stdout)
        except json.JSONDecodeError:
            result = {'ok': False, 'error': {'code': 'non_json_script_output'}, 'exit_code': p.returncode}
        save_new(receipt_path, {'started_at_utc': start, 'finished_at_utc': now(),
                               'request_sha256': sha(request_path), 'skill_script_sha256': sha(SKILL), 'result': result})
        print(json.dumps({'query': label, 'result': result}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['freeze', 'model-metadata', 'encode-search'])
    args = parser.parse_args()
    if args.stage == 'freeze':
        record = freeze()
        print(json.dumps({'frozen_at_utc': record['frozen_at_utc'], 'plan_sha256': record['plan_sha256'],
                          'baseline_files': len(record['prior_tracked_sha256'])}))
    elif args.stage == 'model-metadata':
        metadata()
    else:
        encode_search()


if __name__ == '__main__':
    main()
