"""Retrieve bounded public matching inputs; pin GET/POST requests and responses."""

import argparse
import concurrent.futures
import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'config/psc-comparator-sources.json'


def retrieve(spec, offline=False, refresh=False):
    path = ROOT / 'data/raw' / spec['file']
    if Path(spec['file']).name != spec['file']:
        raise ValueError('Source filename must be a basename')
    if path.exists() and 'sha256' in spec:
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != spec['sha256']:
            raise ValueError('Pinned matching input changed: ' + spec['id'])
        if not refresh:
            return spec
    if offline:
        raise FileNotFoundError(spec['file'])
    body = json.dumps(spec['request_json'], separators=(',', ':')).encode() if 'request_json' in spec else None
    headers = {'User-Agent': 'PSC-public-research/0.5', **spec.get('request_headers', {})}
    if body:
        headers.update({'Content-Type': 'application/json', 'Accept': 'application/json'})
    request = urllib.request.Request(spec['url'], data=body, headers=headers, method=spec.get('method', 'GET'))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                if spec.get('expected_content_range') and (response.status != 206 or not response.headers.get('Content-Range', '').startswith(spec['expected_content_range'])):
                    raise ValueError('Reference response did not honor the exact byte range')
                data = response.read(spec['max_bytes'] + 1)
                if len(data) > spec['max_bytes']:
                    raise ValueError('Public transfer cap exceeded: ' + spec['id'])
                break
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt == 2:
                raise
            time.sleep(min(10, max(1, int(exc.headers.get('Retry-After', 2)))))
    digest = hashlib.sha256(data).hexdigest()
    if 'sha256' in spec and digest != spec['sha256']:
        raise ValueError('Remote matching input changed; preserve pin: ' + spec['id'])
    if not path.exists():
        temporary = path.with_suffix(path.suffix + '.part')
        temporary.write_bytes(data)
        temporary.replace(path)
    elif hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError('Unpinned existing file differs: ' + spec['id'])
    return {**spec, 'sha256': digest, 'bytes': len(data),
            'retrieved_at_utc': spec.get('retrieved_at_utc', datetime.now(timezone.utc).isoformat())}


def acquire(specs, offline=False, refresh=False):
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {'snapshot_date': '2026-09-12', 'sources': []}
    known = {s['id']: s for s in manifest['sources']}
    for spec in specs:
        if spec['id'] in known:
            for key, value in spec.items():
                if key not in ('sha256', 'bytes', 'retrieved_at_utc') and known[spec['id']].get(key) != value:
                    raise ValueError('Request changed for an existing source ID')
    requested = [known.get(s['id'], s) for s in specs]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(retrieve, spec, offline, refresh): spec for spec in requested}
        errors = []
        for future in concurrent.futures.as_completed(futures):
            spec = futures[future]
            try:
                result = future.result()
                known[result['id']] = result
                manifest['sources'] = sorted(known.values(), key=lambda s: s['id'])
                temporary = MANIFEST.with_suffix('.json.part')
                temporary.write_text(json.dumps(manifest, indent=2) + '\n')
                temporary.replace(MANIFEST)
            except Exception as exc:
                errors.append({'id': spec['id'], 'error': type(exc).__name__, 'detail': str(exc)})
        if errors:
            raise ValueError(json.dumps(errors))
    print(json.dumps({'verified_sources': len(requested), 'pinned_total': len(known)}), flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--requests', type=Path, help='JSON list of new public request specifications')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    if args.offline and args.refresh:
        parser.error('Cannot combine offline and refresh')
    specs = json.loads(args.requests.read_text()) if args.requests else json.loads(MANIFEST.read_text())['sources']
    acquire(specs, args.offline, args.refresh)


if __name__ == '__main__':
    main()
