"""Pin source documents, indexes and BF prefixes for the frozen PRKD2 scope."""

import json
import urllib.request
from datetime import datetime, timezone

from prkd2_common import ROOT, save_json, sha, verify_plan


def main():
    plan = verify_plan()
    path = ROOT / 'config/prkd2-sources.json'
    if path.exists():
        record = json.loads(path.read_text())
        for item in record['sources']:
            assert sha((ROOT / 'data/raw' / item['file']).read_bytes()) == item['sha256']
        print(json.dumps({'status': 'verified-cache', 'sources': len(record['sources'])}))
        return
    initial = json.loads((ROOT / 'work/prkd2/source-preparation.json').read_text())
    requests = [
        {'id': 'match-kgp-chr19-index', 'url': 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'max_bytes': 6000000},
        {'id': 'prkd2-eqtl-data-access', 'url': 'https://www.ebi.ac.uk/eqtl/Data_access/', 'max_bytes': 2000000},
        {'id': 'prkd2-coloc-susie-documentation', 'url': 'https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html', 'max_bytes': 2000000},
    ]
    for dataset in plan['datasets']:
        ds = dataset['dataset_id']
        requests.extend([
            {'id': ds + '-credible-sets', 'url': dataset['ftp_cs_path'].replace('ftp://', 'https://'), 'max_bytes': 50000000},
            {'id': ds + '-lbf-header', 'url': dataset['ftp_lbf_path'].replace('ftp://', 'https://'), 'max_bytes': 65536, 'range': 'bytes=0-65535'},
        ])
    progress = ROOT / 'work/prkd2/source-download-progress.json'
    record = json.loads(progress.read_text()) if progress.exists() else {
        'plan_sha256': sha((ROOT / 'config/prkd2-plan.json').read_bytes()),
        'sources': initial['sources'], 'failed_attempts': [],
    }
    for request in requests:
        known = next((s for s in record['sources'] if s['id'] == request['id']), None)
        if known:
            assert sha((ROOT / 'data/raw' / known['file']).read_bytes()) == known['sha256']
            continue
        headers = {'User-Agent': 'PSC-public-research/PRKD2-v1'}
        if 'range' in request:
            headers['Range'] = request['range']
        started = datetime.now(timezone.utc).isoformat()
        try:
            with urllib.request.urlopen(urllib.request.Request(request['url'], headers=headers), timeout=60) as response:
                data = response.read(request['max_bytes'] + 1)
                if len(data) > request['max_bytes']:
                    raise ValueError('Source exceeds fixed byte limit')
                if 'range' in request and (response.status != 206 or not response.headers.get('Content-Range', '').startswith('bytes 0-65535/')):
                    raise ValueError('Server did not honor exact compressed prefix')
                receipt = {'http_status': response.status, 'content_range': response.headers.get('Content-Range'),
                           'content_length': response.headers.get('Content-Length'), 'etag': response.headers.get('ETag'),
                           'last_modified': response.headers.get('Last-Modified')}
            filename = 'prkd2-source-' + request['id'] + '.data'
            (ROOT / 'data/raw' / filename).write_bytes(data)
            record['sources'].append({**request, **receipt, 'file': filename, 'bytes': len(data), 'sha256': sha(data), 'retrieved_at_utc': started})
            save_json(progress, record)
            print(json.dumps({'source': request['id'], 'bytes': len(data), 'status': 'pinned'}), flush=True)
        except Exception as exc:
            record['failed_attempts'].append({'request': request, 'at_utc': started, 'error': str(exc)})
            save_json(progress, record)
            raise
    record['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    record['source_fetcher_sha256'] = sha(__import__('pathlib').Path(__file__).read_bytes())
    record['archive_access'] = 'Official documented indexed archive downloads; no deprecated eQTL REST API used.'
    save_json(path, record)
    print(json.dumps({'status': 'complete', 'sources': len(record['sources'])}), flush=True)


if __name__ == '__main__':
    main()
