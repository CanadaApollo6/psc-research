"""Read each complete public BF archive; retain all PRKD2 component vectors."""

import argparse
import gzip
import io
import json
import urllib.request
from datetime import datetime, timezone

from fetch_coloc_signal_bfs import CountedReader
from prkd2_common import ROOT, gzip_bytes, save_json, sha, verify_plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', choices=['QTD000021', 'QTD000504'])
    args = parser.parse_args()
    plan = verify_plan()
    ds = args.dataset
    path = ROOT / ('config/prkd2-bf-query-' + ds + '.json')
    if path.exists():
        old = json.loads(path.read_text())
        data = (ROOT / old['file']).read_bytes()
        assert old['archive_eof_verified'] and sha(data) == old['sha256']
        assert sha(gzip.decompress(data)) == old['uncompressed_sha256']
        print(json.dumps({'dataset': ds, 'status': 'verified-cache'}))
        return
    sources = {s['id']: s for s in json.loads((ROOT / 'config/prkd2-sources.json').read_text())['sources']}
    source = sources[ds + '-lbf-header']
    cap = plan['transfer_bounds']['maximum_published_bf_compressed_bytes_per_file']
    record = {'dataset_id': ds, 'gene_id': 'ENSG00000105287', 'url': source['url'],
              'started_at_utc': datetime.now(timezone.utc).isoformat(),
              'plan_sha256': sha((ROOT / 'config/prkd2-plan.json').read_bytes()),
              'reader_sha256': sha(__import__('pathlib').Path(__file__).read_bytes()),
              'stream_reader_sha256': sha((ROOT / 'scripts/fetch_coloc_signal_bfs.py').read_bytes())}
    save_json(ROOT / ('work/prkd2/bf-start-' + ds + '.json'), record)
    try:
        req = urllib.request.Request(source['url'], headers={'User-Agent': 'PSC-public-research/PRKD2-v1'})
        with urllib.request.urlopen(req, timeout=60) as response:
            size = int(response.headers['Content-Length'])
            if response.status != 200 or size > cap or size != int(source['content_range'].split('/')[-1]):
                raise ValueError('BF size/status differs from frozen source prefix')
            reader = CountedReader(response, cap, ds)
            selected = []
            scanned = 0
            with gzip.GzipFile(fileobj=io.BufferedReader(reader, buffer_size=1024 * 1024), mode='rb') as stream:
                header = stream.readline()
                if not header.startswith(b'molecular_trait_id\tregion\tvariant\t'):
                    raise ValueError('Unexpected published BF header')
                for line in stream:
                    scanned += 1
                    if line.startswith(b'ENSG00000105287\t'):
                        selected.append(line)
            if reader.count != size or sha(reader.prefix) != source['sha256']:
                raise ValueError('Incomplete archive or altered pinned prefix')
        if not selected:
            raise ValueError('No full PRKD2 BF vectors in published archive')
        plain = header + b''.join(selected)
        data = gzip_bytes(plain)
        name = 'data/derived/prkd2-query-rows/' + ds + '-PRKD2-published-bfs.tsv.gz'
        (ROOT / name).write_bytes(data)
        record.update(finished_at_utc=datetime.now(timezone.utc).isoformat(), compressed_archive_bytes=size,
                      compressed_archive_sha256=reader.digest.hexdigest(), archive_eof_verified=True,
                      all_archive_rows_scanned=scanned, retained_gene_rows=len(selected), file=name,
                      sha256=sha(data), uncompressed_sha256=sha(plain), bytes=len(data))
        save_json(path, record)
        print(json.dumps({'dataset': ds, 'status': 'complete', 'all_gene_rows_retained': len(selected), 'archive_bytes': size}), flush=True)
    except Exception as exc:
        record.update(failed_at_utc=datetime.now(timezone.utc).isoformat(), error=str(exc))
        save_json(ROOT / ('work/prkd2/bf-failure-' + ds + '.json'), record)
        raise


if __name__ == '__main__':
    main()
