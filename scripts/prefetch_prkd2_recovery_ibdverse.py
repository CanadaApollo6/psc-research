"""Fetch the frozen ZIP ranges with two workers and one manifest writer.

This changes transfer scheduling only. The original extractor then verifies
each cached block, the complete DEFLATE stream, byte count and ZIP CRC.
"""

import concurrent.futures
import json
import struct
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, save_json, sha
from fetch_prkd2_recovery_ibdverse import AuthorRangeReader


def main():
    plan_path = ROOT / 'config/prkd2-recovery-author-plan.json'
    plan = json.loads(plan_path.read_text())
    path = ROOT / 'config/prkd2-recovery-ibdverse-queries.json'
    record = json.loads(path.read_text())
    if record['plan_sha256'] != sha(plan_path.read_bytes()):
        raise ValueError('Changed frozen author plan')
    if record['reader_sha256'] != sha((ROOT / 'scripts/fetch_prkd2_recovery_ibdverse.py').read_bytes()):
        raise ValueError('Changed original archive reader')
    metadata = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-zip-metadata.json').read_text())['queries'][0]
    tasks = []
    for selected in plan['ibdverse']:
        query = next((q for q in record['queries'] if q['id'] == selected['id']), None)
        if query is None:
            query = {**selected, **{k: metadata[k] for k in ('url', 'remote_size', 'etag', 'last_modified')},
                     'ranges': [], 'started_at_utc': datetime.now(timezone.utc).isoformat()}
            record['queries'].append(query)
        reader = AuthorRangeReader(query, record, path, False)
        member = selected['member']
        offset = int(member['local_header_offset'])
        fields = struct.unpack('<4s5H3L2H', reader.fetch(offset, 30))
        if fields[0] != b'PK\x03\x04' or fields[3] != 8 or fields[2] & 1:
            raise ValueError('Unexpected ZIP header')
        n, extra = fields[-2:]
        if reader.fetch(offset + 30, n + extra)[:n].decode() != member['filename']:
            raise ValueError('Wrong archive member')
        start, length = offset + 30 + n + extra, int(member['compressed_bytes'])
        if length + 30 + n + extra > plan['transfer_cap_per_context_bytes']:
            raise ValueError('Source exceeds frozen transfer cap')
        for relative in range(0, length, 8388608):
            size = min(8388608, length - relative)
            identifier = f'{start + relative}-{start + relative + size - 1}'
            previous = next((r for r in query['ranges'] if r['range'] == identifier), None)
            if previous:
                if sha((ROOT / previous['file']).read_bytes()) != previous['sha256']:
                    raise ValueError('Changed cached range')
            else:
                tasks.append((query['id'], start + relative, size))
    lock, last = threading.Lock(), [0.0]

    def transfer(task):
        qid, start, length = task
        query = next(q for q in record['queries'] if q['id'] == qid)
        identifier = f'{start}-{start + length - 1}'
        with lock:
            time.sleep(max(0, 2 - (time.monotonic() - last[0])))
            last[0] = time.monotonic()
        req = urllib.request.Request(query['url'], headers={'Range': 'bytes=' + identifier, 'User-Agent': 'PSC-public-research/PRKD2-recovery-v1'})
        with urllib.request.urlopen(req, timeout=60) as response:
            if response.status != 206 or response.headers.get('Content-Range') != f"bytes {identifier}/{query['remote_size']}" or response.headers.get('ETag') != query['etag'] or response.headers.get('Last-Modified') != query['last_modified']:
                raise ValueError('Changed remote range identity')
            data = response.read(length + 1)
        if len(data) != length:
            raise ValueError('Truncated byte range')
        output = ROOT / f'data/raw/prkd2-recovery-zip-ranges/{qid}-{identifier}.data'
        output.write_bytes(data)
        return qid, {'range': identifier, 'file': str(output.relative_to(ROOT)), 'bytes': length,
                     'sha256': sha(data), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}

    record['range_prefetch'] = {'script_sha256': sha(Path(__file__).read_bytes()), 'workers': 2,
                                'minimum_request_start_spacing_seconds': 2,
                                'reason': 'Preserve the original frozen extractor; change transfer scheduling only.'}
    save_json(path, record)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(transfer, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            qid, result = future.result()
            query = next(q for q in record['queries'] if q['id'] == qid)
            query['ranges'].append(result)
            save_json(path, record)
            print(json.dumps({'id': qid, 'verified_bytes': sum(r['bytes'] for r in query['ranges'])}), flush=True)


if __name__ == '__main__':
    main()
