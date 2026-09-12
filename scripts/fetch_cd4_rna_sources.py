"""Retrieve and pin the public CD4 metadata, indexes and bounded header ranges."""

import argparse
import json
from pathlib import Path

from fetch_comparison_sources import retrieve

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'config/cd4-rna-sources.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    (ROOT / 'data/raw').mkdir(exist_ok=True)
    for index, source in enumerate(manifest['sources']):
        result = retrieve(source, offline=args.offline)
        manifest['sources'][index] = result
        encoded = json.dumps(manifest, indent=2) + '\n'
        if encoded != MANIFEST.read_text():
            temporary = MANIFEST.with_suffix('.json.part')
            temporary.write_text(encoded); temporary.replace(MANIFEST)
        print(json.dumps({'source': result['id'], 'bytes': result.get('bytes'), 'status': 'verified'}), flush=True)


if __name__ == '__main__':
    main()
