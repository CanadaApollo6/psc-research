"""Resume the slow original BLUEPRINT acquisition with verified bounded ranges.

At most two public HTTP requests are in flight, with starts separated by two
seconds. Identical URL, ETag, size and byte positions are required throughout.
"""

import concurrent.futures
import hashlib
import json
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, save_json, sha


def main():
    sources = json.loads((ROOT / 'config/prkd2-recovery-sources.json').read_text())['sources']
    head = next(s for s in sources if s['id'] == 'blueprint-original-monocytes-head')
    size = int(head['headers']['Content-Length'])
    archive = ROOT / 'data/raw/prkd2-recovery-original-BLUEPRINT-monocytes.txt.gz'
    path = ROOT / 'config/prkd2-recovery-blueprint-range-acquisition.json'
    directory = ROOT / 'data/raw/prkd2-recovery-blueprint-ranges'
    directory.mkdir(exist_ok=True)
    pins = {'source_head': head, 'reader_sha256': sha(Path(__file__).read_bytes()),
            'chunk_bytes': 8 * 1024 * 1024, 'maximum_simultaneous_requests': 2,
            'minimum_seconds_between_request_starts': 2}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'ranges': [], 'created_at_utc': datetime.now(timezone.utc).isoformat()}
    if any(record[k] != v for k, v in pins.items()):
        raise ValueError('Changed range acquisition plan')
    if archive.exists():
        if archive.stat().st_size != size or sha(archive.read_bytes()) != record['archive_sha256']:
            raise ValueError('Changed completed archive')
        print(json.dumps({'status': 'verified-cache'}))
        return
    abandoned = archive.with_suffix('.part')
    if abandoned.exists():
        kept = abandoned.with_name(abandoned.name + '.initial-stream')
        if kept.exists():
            raise ValueError('Initial stream already retained')
        abandoned.rename(kept)
        record['initial_stream_attempt'] = {'file': str(kept.relative_to(ROOT)), 'bytes': kept.stat().st_size,
                                          'sha256': sha(kept.read_bytes()), 'status': 'Interrupted for a faster bounded acquisition; no association rows parsed'}
    if size + record.get('initial_stream_attempt', {}).get('bytes', 0) > 2000000000:
        raise ValueError('Total acquisition cap exceeded')
    save_json(path, record)
    lock = threading.Lock()
    last_start = [0.0]

    def get_chunk(start):
        end = min(size - 1, start + pins['chunk_bytes'] - 1)
        prior = next((r for r in record['ranges'] if r['start'] == start), None)
        if prior:
            if sha((ROOT / prior['file']).read_bytes()) != prior['sha256']:
                raise ValueError('Changed cached range')
            return prior
        with lock:
            delay = max(0, 2 - (time.monotonic() - last_start[0]))
            time.sleep(delay)
            last_start[0] = time.monotonic()
        request = urllib.request.Request(head['url'], headers={'Range': f'bytes={start}-{end}', 'User-Agent': 'PSC-public-research/PRKD2-recovery-v1'})
        with urllib.request.urlopen(request, timeout=60) as response:
            if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{size}' or response.headers.get('ETag') != head['headers']['ETag']:
                raise ValueError('Range or archive identity mismatch')
            data = response.read(end - start + 2)
        if len(data) != end - start + 1:
            raise ValueError('Incomplete or oversized range')
        target = directory / f'{start}-{end}.data'
        target.write_bytes(data)
        return {'start': start, 'end': end, 'bytes': len(data), 'sha256': sha(data),
                'file': str(target.relative_to(ROOT)), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}

    starts = list(range(0, size, pins['chunk_bytes']))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(get_chunk, starts):
            if not any(r['start'] == result['start'] for r in record['ranges']):
                record['ranges'].append(result)
                save_json(path, record)
            if len(record['ranges']) % 8 == 0:
                print(json.dumps({'verified_bytes': sum(r['bytes'] for r in record['ranges']), 'total': size}), flush=True)
    record['ranges'].sort(key=lambda r: r['start'])
    if [r['start'] for r in record['ranges']] != starts:
        raise ValueError('Incomplete range inventory')
    digest = hashlib.sha256()
    part = archive.with_suffix('.assembling')
    with part.open('wb') as out:
        for item in record['ranges']:
            data = (ROOT / item['file']).read_bytes()
            if sha(data) != item['sha256']:
                raise ValueError('Range changed during assembly')
            out.write(data)
            digest.update(data)
    if part.stat().st_size != size:
        raise ValueError('Assembled archive size mismatch')
    prefix = next(s for s in sources if s['id'] == 'blueprint-original-monocytes-prefix')
    with part.open('rb') as stream:
        if sha(stream.read(65536)) != prefix['sha256']:
            raise ValueError('Independent source prefix mismatch')
    initial = record.get('initial_stream_attempt')
    if initial:
        with part.open('rb') as stream:
            if sha(stream.read(initial['bytes'])) != initial['sha256']:
                raise ValueError('Range data disagrees with original stream')
    part.replace(archive)
    record.update({'archive_sha256': digest.hexdigest(), 'archive_bytes': size,
                   'completed_at_utc': datetime.now(timezone.utc).isoformat()})
    save_json(path, record)
    print(json.dumps({'status': 'complete', 'bytes': size, 'sha256': digest.hexdigest()}), flush=True)


if __name__ == '__main__':
    main()
