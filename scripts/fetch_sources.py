"""Retrieve small public research inputs and verify their pinned output hashes."""

import argparse
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(content):
    return hashlib.sha256(content).hexdigest()


def validate(content, source):
    if len(content) > source['max_bytes']:
        raise ValueError(f"Size limit exceeded for {source['id']}")
    if sha256(content) != source['sha256']:
        raise ValueError(f"Checksum changed for {source['id']}; review the source before updating its pin")


def fetch(source, directory, offline=False):
    name = source['file']
    if Path(name).name != name:
        raise ValueError('Source filenames must not include a directory')
    target = directory / name
    if target.exists():
        content = target.read_bytes()
        validate(content, source)
        return {'id': source['id'], 'status': 'verified-cache', 'bytes': len(content)}
    if offline:
        raise FileNotFoundError(f'Missing cached source: {name}')
    headers = {'User-Agent': 'PSC-public-research/0.2', **source.get('request_headers', {})}
    request = urllib.request.Request(source['url'], headers=headers)
    with urllib.request.urlopen(request, timeout=45) as response:
        if 'Range' in headers and response.status != 206:
            raise ValueError('Reference server did not honor the byte-range request')
        content = response.read(source['max_bytes'] + 1)
    if len(content) > source['max_bytes']:
        raise ValueError(f"Transfer limit exceeded for {source['id']}")
    if source.get('compression') == 'gzip':
        with gzip.GzipFile(fileobj=io.BytesIO(content)) as stream:
            content = stream.read(source['max_bytes'] + 1)
    validate(content, source)
    temporary = target.with_suffix(target.suffix + '.part')
    temporary.write_bytes(content)
    temporary.replace(target)
    return {'id': source['id'], 'status': 'downloaded', 'bytes': len(content)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true', help='Only verify cached files')
    parser.add_argument('--manifest', default='config/sources.json', help='Source manifest relative to repository root')
    args = parser.parse_args()
    sources = json.loads((ROOT / args.manifest).read_text())['sources']
    directory = ROOT / 'data/raw'
    directory.mkdir(parents=True, exist_ok=True)
    results = []
    failed = False
    for source in sources:
        try:
            result = fetch(source, directory, offline=args.offline)
        except Exception as exc:
            failed = True
            result = {'id': source['id'], 'status': 'failed', 'error': str(exc)}
        results.append(result)
        print(json.dumps(result))
    log = {'checked_at_utc': datetime.now(timezone.utc).isoformat(), 'results': results}
    (directory / 'retrieval-log.json').write_text(json.dumps(log, indent=2) + '\n')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
