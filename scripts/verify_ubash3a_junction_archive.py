"""Cross-check two recorded web-service rows against their indexed raw archive."""

import csv
import hashlib
import io
import json
from pathlib import Path

import pysam

from fetch_sources import fetch, validate

ROOT = Path(__file__).resolve().parents[1]
csv.field_size_limit(2_000_000)


def main():
    manifest = json.loads((ROOT / 'config/ubash3a-junction-sources.json').read_text())
    sources = {s['id']: s for s in manifest['sources']}
    index = sources['snaptron-junction-index']
    fetch(index, ROOT / 'data/raw')
    expected = {}
    for name in ('canonical', 'plus29'):
        s = sources[f'snaptron-{name}-unfiltered']
        content = (ROOT / 'data/raw' / s['file']).read_bytes()
        validate(content, s)
        rows = list(csv.DictReader(io.StringIO(content.decode()), delimiter='\t'))
        if len(rows) != 1:
            raise ValueError('Expected one complete saved web-service row')
        row = rows[0]
        expected[(row['chromosome'], row['start'], row['end'])] = ('\t'.join(list(row.values())[1:]), name)
    url = 'https://snaptron.cs.jhu.edu/data/srav3h/junctions.bgz'
    with pysam.TabixFile(url, index=str(ROOT / 'data/raw' / index['file'])) as table:
        rows = list(table.fetch('chr21', 42434954, 42434984))
    observed = {}
    for row in rows:
        fields = row.split('\t')
        key = tuple(fields[1:4])
        if key in expected:
            if key in observed or row != expected[key][0]:
                raise ValueError('Raw archive and recorded web-service row disagree')
            observed[key] = row
    if set(observed) != set(expected):
        raise ValueError('A target junction is missing from the raw archive')
    result = {'status': 'exact-match', 'url': url, 'genome_build': 'GRCh38',
              'query_chromosome': 'chr21', 'query_interval0': [42434954, 42434984],
              'index_sha256': index['sha256'], 'pysam_version': pysam.__version__,
              'regional_rows': len(rows), 'matched_target_rows': len(observed),
              'target_row_sha256': {expected[key][1]: hashlib.sha256(row.encode()).hexdigest()
                                    for key, row in sorted(observed.items())},
              'scope': 'Exact target rows and local index verified; no full remote archive hash claimed'}
    (ROOT / 'reports/ubash3a-junction-archive-check.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
