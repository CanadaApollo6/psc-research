"""Source receipts and safe, immutable inputs for the PRKD2 mechanism follow-up."""
from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests

from prepare_prkd2_mechanism_sources import ROOT, RAW, PLAN, SKILL, freeze, now, sha, save_new


def write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def update_ledger(path, entry):
    path = Path(path)
    ledger = json.loads(path.read_text()) if path.exists() else {'sources': []}
    previous = next((r for r in ledger['sources'] if r['id'] == entry['id']), None)
    if previous is not None:
        if previous != entry:
            raise ValueError('Source ID changed: ' + entry['id'])
    else:
        ledger['sources'].append(entry)
        ledger['sources'].sort(key=lambda r: r['id'])
        path.write_text(json.dumps(ledger, indent=2, allow_nan=False) + '\n')


def fetch_public(label, url, max_bytes=5000000, headers=None, body=None, content_range=None):
    """Bounded non-ENCODE metadata/reference transfer; API secrets are never used."""
    if 'encodeproject.org' in urlparse(url).netloc:
        raise ValueError('Use the ENCODE skill wrapper for ENCODE API requests')
    ledger_path = ROOT / 'config/prkd2-mechanism-reference-sources.json'
    path = RAW / 'references' / label
    known = json.loads(ledger_path.read_text())['sources'] if ledger_path.exists() else []
    match = next((r for r in known if r['id'] == label), None)
    if match:
        if match['url'] != url or match.get('request_headers', {}) != (headers or {}) or match.get('request_body') != body:
            raise ValueError('Reference request changed')
        if sha(path) != match['sha256']:
            raise ValueError('Pinned reference changed')
        return path
    if path.exists():
        raise ValueError('Unpinned reference already exists: ' + str(path))
    requested_at = now()
    with requests.request('POST' if body is not None else 'GET', url,
                          headers={'User-Agent': 'PSC-public-research/1.0', **(headers or {})},
                          json=body, timeout=45, stream=True) as response:
        response.raise_for_status()
        if content_range and (response.status_code != 206 or response.headers.get('Content-Range') != content_range):
            raise ValueError('Exact reference byte range not honored')
        data = bytearray()
        for block in response.iter_content(65536):
            data.extend(block)
            if len(data) > max_bytes:
                raise ValueError('Bounded reference transfer exceeded')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        safe_headers = {k: response.headers[k] for k in ['Content-Type', 'Content-Length', 'Content-Range', 'ETag', 'Last-Modified'] if k in response.headers}
        entry = {'id': label, 'url': url, 'request_headers': headers or {}, 'request_body': body,
                 'requested_at_utc': requested_at, 'retrieved_at_utc': now(), 'status_code': response.status_code,
                 'response_headers': safe_headers, 'max_bytes': max_bytes, 'bytes': len(data),
                 'file': str(path.relative_to(ROOT)), 'sha256': sha(path)}
    update_ledger(ledger_path, entry)
    return path


def encode_request(label, path, params=None):
    """Run all ENCODE metadata requests through the required installed skill."""
    folder = RAW / 'encode'
    destination = folder / (label + '.json')
    receipt = folder / (label + '-receipt.json')
    request = {'base_url': 'https://www.encodeproject.org', 'path': path,
               'params': {'format': 'json', **(params or {})}, 'headers': {'Accept': 'application/json'},
               'timeout_sec': 45, 'save_raw': True, 'raw_output_path': str(destination), 'max_items': 2, 'max_depth': 1}
    request_path = ROOT / 'config/prkd2-mechanism-encode-requests' / (label + '.json')
    save_new(request_path, request)
    if receipt.exists():
        record = json.loads(receipt.read_text())
        if record['request_sha256'] != sha(request_path):
            raise ValueError('ENCODE request changed')
        if not record['ok']:
            return None
        if sha(destination) != record['sha256']:
            raise ValueError('ENCODE source changed')
        return json.loads(destination.read_text())
    start = now()
    result = subprocess.run([sys.executable, str(SKILL)], input=json.dumps(request), capture_output=True, text=True)
    try:
        summary = json.loads(result.stdout)
    except json.JSONDecodeError:
        summary = {'ok': False, 'error': {'code': 'non_json_skill_output'}, 'exit_code': result.returncode}
    record = {'requested_at_utc': start, 'retrieved_at_utc': now(), 'ok': summary.get('ok', False),
              'request_sha256': sha(request_path), 'skill_script_sha256': sha(SKILL), 'summary': summary}
    if record['ok']:
        record.update(file=str(destination.relative_to(ROOT)), sha256=sha(destination), bytes=destination.stat().st_size)
    save_new(receipt, record)
    print(json.dumps({'encode_object': label, 'ok': record['ok'], 'bytes': record.get('bytes'),
                      'error_code': summary.get('error', {}).get('code')}), flush=True)
    return json.loads(destination.read_text()) if record['ok'] else None


def verify_receipt(path, expected):
    if sha(path) != expected:
        raise ValueError('Input checksum mismatch: ' + str(path))


def clean_records(frame):
    return json.loads(frame.to_json(orient='records'))
