"""Verify recorded QTL subsets, or repeat the eight indexed public-data queries.

Refresh compares with the frozen result and never overwrites it. Indexes and
headers are pinned separately; no checksum of the entire remote file is claimed.
"""

import argparse
import hashlib
import json
import time
import zlib
from pathlib import Path

from fetch_sources import fetch

ROOT = Path(__file__).resolve().parents[1]


def canonical_subset(header, rows):
    return (header + '\n' + '\n'.join(rows) + ('\n' if rows else '')).encode()


def source_map():
    return {s['id']: s for s in json.loads(
        (ROOT / 'config/qtl-followup-sources.json').read_text())['sources']}


def verified_queries():
    manifest = json.loads((ROOT / 'config/qtl-followup-queries.json').read_text())
    for query in manifest['queries']:
        content = (ROOT / query['file']).read_bytes()
        if hashlib.sha256(content).hexdigest() != query['sha256']:
            raise ValueError(f"Recorded subset changed: {query['dataset']['dataset_id']}")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--offline', action='store_true', help='Verify recorded subsets only (the default)')
    mode.add_argument('--refresh', action='store_true', help='Repeat remote queries and require exact agreement')
    args = parser.parse_args()
    manifest = verified_queries()
    if not args.refresh:
        print(json.dumps({'status': 'verified-recorded-subsets', 'queries': len(manifest['queries']),
                          'rows': sum(q['rows'] for q in manifest['queries'])}))
        return

    import pysam
    if pysam.__version__ != manifest['pysam_version']:
        raise ValueError('Install the recorded pysam version before refreshing')
    sources = source_map()
    for query in manifest['queries']:
        index = sources[query['index_source_id']]
        header_source = sources[query['header_source_id']]
        for source in (index, header_source):
            fetch(source, ROOT / 'data/raw')
        header_bytes = (ROOT / 'data/raw' / header_source['file']).read_bytes()
        header = zlib.decompressobj(31).decompress(header_bytes).decode().splitlines()[0]
        chromosome = query['query']['variant'].split('_')[0].removeprefix('chr')
        with pysam.TabixFile(query['url'], index=str(ROOT / 'data/raw' / index['file'])) as table:
            available = chromosome in table.contigs
            rows = list(table.fetch(chromosome, query['query_start0'], query['query_end0'])) if available else []
        content = canonical_subset(header, rows)
        if available != query['contig_present'] or hashlib.sha256(content).hexdigest() != query['sha256']:
            raise ValueError(f"Remote result changed for {query['dataset']['dataset_id']}; preserve the recorded result and review")
        print(json.dumps({'dataset_id': query['dataset']['dataset_id'], 'status': 'matched-remote', 'rows': len(rows)}), flush=True)
        time.sleep(2)  # The archive requests a couple of seconds between tabix queries.


if __name__ == '__main__':
    main()
