"""Reproduce the IBDverse archive inventory from its cached ZIP directory."""

import csv
import io
import json
import zipfile

from coloc_common import ROOT, sha
from fetch_prkd2_recovery_ibdverse import AuthorRangeReader


def main():
    path = ROOT / 'config/prkd2-recovery-ibdverse-zip-metadata.json'
    manifest = json.loads(path.read_text())
    reader = AuthorRangeReader(manifest['queries'][0], manifest, path, True)
    with zipfile.ZipFile(reader) as archive:
        rows = [{'filename': item.filename, 'uncompressed_bytes': item.file_size,
                 'compressed_bytes': item.compress_size, 'compression_method': item.compress_type,
                 'local_header_offset': item.header_offset, 'crc32': item.CRC}
                for item in archive.infolist()]
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    data = stream.getvalue().encode()
    expected = manifest['result']
    if len(rows) != expected['members'] or sha(data) != expected['inventory_sha256']:
        # Accept a newline representation only when it matches the frozen hash.
        data = data.replace(b'\r\n', b'\n')
    if sha(data) != expected['inventory_sha256'] or data != (ROOT / expected['inventory_file']).read_bytes():
        raise ValueError('ZIP central-directory inventory replay differs')
    print(json.dumps({'status': 'verified', 'members': len(rows), 'sha256': sha(data)}))


if __name__ == '__main__':
    main()
