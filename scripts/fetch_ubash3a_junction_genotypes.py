"""Fetch or verify the single public variant required for the junction analysis."""

import argparse
import hashlib
import json
from pathlib import Path

from fetch_sources import fetch

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    if args.offline and args.refresh:
        parser.error('--offline and --refresh cannot be combined')
    run = json.loads((ROOT / 'config/ubash3a-junction-run.json').read_text())
    query = run['genotype_query']
    target = ROOT / query['file']
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest() != query['sha256']:
            raise ValueError('Recorded genotype subset changed')
        if not args.refresh:
            print('Verified single-variant genotype subset')
            return
    elif args.offline:
        raise FileNotFoundError(target)
    import pysam
    if pysam.__version__ != query['pysam_version']:
        raise ValueError('Use the recorded pysam version')
    sources = json.loads((ROOT / 'config/ubash3a-junction-sources.json').read_text())['sources']
    index = next(s for s in sources if s['id'] == query['index_source_id'])
    fetch(index, ROOT / 'data/raw')
    with pysam.TabixFile(query['url'], index=str(ROOT / 'data/raw' / index['file'])) as table:
        header = list(table.header)
        rows = list(table.fetch(query['chromosome'], *query['query_interval0']))
    content = ('\n'.join(header + rows) + '\n').encode()
    if hashlib.sha256(content).hexdigest() != query['sha256']:
        raise ValueError('Remote genotype subset differs; preserve the recorded result and review')
    if not target.exists():
        target.write_bytes(content)
    print('Verified remote single-variant query')


if __name__ == '__main__':
    main()
