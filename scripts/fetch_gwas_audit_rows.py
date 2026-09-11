"""Retrieve two exact PSC GWAS rows using bounded public byte ranges.

The original uncompressed GWAS file is coordinate sorted. Binary search avoids
downloading its 662 MB. All fetched windows, including discarded fragments, are
hashed and retained; selected rows are copied verbatim into the derived table.
This is a targeted lookup, not a complete import or a fine-mapping analysis.
"""

import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004030/ipscsg2016.result.combined.full.with_header.txt'
TOTAL = 662022204
WINDOW = 32768
TARGETS = [('rs313839', 19, 47221557), ('rs1893592', 21, 43855067)]


def complete_rows(content):
    # Every range may start/end inside a line. Keep only complete interior lines.
    rows = []
    for line in content.decode('ascii').split('\n')[1:-1]:
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 13 or not fields[0].isdigit():
            raise ValueError('Unexpected GWAS row schema')
        rows.append(((int(fields[0]), int(fields[2])), fields[1], line))
    if any(a[0] > b[0] for a, b in zip(rows, rows[1:])):
        raise ValueError('GWAS window is not coordinate sorted')
    return rows


def main():
    manifest_path = ROOT / 'config/association-sources.json'
    manifest = json.loads(manifest_path.read_text())
    cache = {row['file']: row for row in manifest['sources']}
    records = []

    def fetch(start, end):
        start, end = max(0, int(start)), min(TOTAL - 1, int(end))
        name = f'GCST004030-bytes-{start}-{end}.txt'
        path = ROOT / 'data/raw' / name
        if name in cache:
            content = path.read_bytes()
            if hashlib.sha256(content).hexdigest() != cache[name]['sha256']:
                raise ValueError('Cached GWAS range checksum mismatch')
        else:
            if len(records) >= 40:
                raise ValueError('Targeted request budget exceeded')
            headers = {'Range': f'bytes={start}-{end}'}
            req = urllib.request.Request(URL, headers={'User-Agent': 'PSC-public-research/0.3', **headers})
            with urllib.request.urlopen(req, timeout=30) as response:
                if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{TOTAL}':
                    raise ValueError('Unexpected public-file byte range or source size')
                content = response.read(end - start + 2)
            if len(content) != end - start + 1:
                raise ValueError('Truncated or oversized response')
            path.write_bytes(content)
            source = {'id': name[:-4], 'file': name, 'url': URL, 'request_headers': headers,
                      'sha256': hashlib.sha256(content).hexdigest(), 'max_bytes': 200000,
                      'bytes_after_decompression': len(content), 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}
            records.append(source)
            cache[name] = source
        return complete_rows(content)

    found = []
    try:
        for rsid, chromosome, position in TARGETS:
            target = chromosome, position
            low, high = 0, TOTAL - WINDOW
            for step in range(26):
                if high - low < WINDOW:
                    rows = fetch(low - WINDOW, high + 2 * WINDOW)
                else:
                    mid = (low + high) // 2
                    rows = fetch(mid, mid + WINDOW - 1)
                matches = [r for r in rows if r[0] == target and r[1] == rsid]
                if len(matches) == 1:
                    found.append(matches[0][2])
                    print(json.dumps({'rsid': rsid, 'status': 'exact original row found'}), flush=True)
                    break
                if not rows or high - low < WINDOW:
                    raise ValueError(f'Could not verify exact row for {rsid}')
                if target < rows[0][0]:
                    high = mid
                elif target > rows[-1][0]:
                    low = mid
                else:
                    raise ValueError(f'Coordinate present but requested rsID missing: {rsid}')
            else:
                raise ValueError('Search budget exhausted')
    finally:
        manifest['sources'].extend(records)
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    header = next(line for line in (ROOT / 'data/raw/GCST004030-original-header.txt').read_text().splitlines() if line.startswith('#chr '))
    output = ROOT / 'data/derived/gwas-audit-original-rows.txt'
    output.write_text(header + '\n' + '\n'.join(found) + '\n')
    print(json.dumps({'selected_rows': len(found), 'new_requests': len(records), 'new_bytes': sum(r['bytes_after_decompression'] for r in records)}))


if __name__ == '__main__':
    main()
