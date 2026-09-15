"""Small, cached public-source receipts for the PRKD2 recovery audit.

This module never treats an access failure as evidence that a variant is absent.
Association archives require separate bounded readers and a frozen query plan.
"""

import argparse
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, save_json, sha

MANIFEST = ROOT / 'config/prkd2-recovery-sources.json'
DIRECTORY = ROOT / 'data/raw/prkd2-recovery-sources'


def fetch(identifier, url, *, method='GET', byte_range=None, limit=3000000):
    DIRECTORY.mkdir(exist_ok=True)
    record = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {'sources': []}
    spec = {'id': identifier, 'url': url, 'method': method, 'byte_range': byte_range, 'byte_limit': limit}
    prior = next((s for s in record['sources'] if s['id'] == identifier), None)
    if prior:
        if any(prior[k] != v for k, v in spec.items()):
            raise ValueError('Source identity reused for another request')
        if prior.get('file') and sha((ROOT / prior['file']).read_bytes()) != prior['sha256']:
            raise ValueError('Cached source changed')
        return prior
    receipt = {**spec, 'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}
    headers = {'User-Agent': 'PSC-public-research/PRKD2-recovery-v1'}
    if byte_range:
        headers['Range'] = byte_range
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method=method, headers=headers), timeout=35) as response:
            excluded = {'set-cookie', 'cookie', 'authorization', 'proxy-authorization', 'x-api-key'}
            receipt.update({'status': response.status, 'resolved_url': response.url,
                            'headers': {k: v for k, v in response.headers.items() if k.lower() not in excluded}})
            omitted = [k for k in response.headers if k.lower() in excluded]
            if omitted:
                receipt['headers_omitted'] = omitted
            if byte_range and response.status != 206:
                raise ValueError('Server did not honor range; archive not downloaded')
            data = response.read(limit + 1)
            if len(data) > limit:
                raise ValueError('Source exceeds byte cap')
            target = DIRECTORY / (identifier + '.data')
            target.write_bytes(data)
            receipt.update({'file': str(target.relative_to(ROOT)), 'bytes': len(data), 'sha256': sha(data)})
    except (OSError, ValueError) as exc:
        receipt.update({'error': str(exc), 'status': getattr(exc, 'code', receipt.get('status'))})
    record['sources'].append(receipt)
    save_json(MANIFEST, record)
    print(json.dumps({k: receipt.get(k) for k in ['id', 'status', 'bytes', 'error']}), flush=True)
    time.sleep(2)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--id')
    parser.add_argument('--url')
    parser.add_argument('--method', default='GET')
    parser.add_argument('--range')
    args = parser.parse_args()
    if args.id:
        fetch(args.id, args.url, method=args.method, byte_range=args.range)
        return
    requests = [
        ('catalogue-current-commit', 'https://api.github.com/repos/eQTL-Catalogue/eQTL-Catalogue-resources/commits/master'),
        ('catalogue-data-access', 'https://www.ebi.ac.uk/eqtl/Data_access/'),
        ('catalogue-studies', 'https://www.ebi.ac.uk/eqtl/Studies/'),
        ('dice-downloads', 'https://dice-database.org/downloads'),
        ('blueprint-original-qtls', 'https://blueprint-dev.bioinfo.cnio.es/WP10/qtls'),
        ('blueprint-epivar-index', 'https://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint_Epivar/'),
        ('psc-gcst004030-index', 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004030/'),
    ]
    for identifier, url in requests:
        fetch(identifier, url)


if __name__ == '__main__':
    main()
