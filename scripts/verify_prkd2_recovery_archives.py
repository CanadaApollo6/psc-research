"""Independent gzip/ripgrep extraction of the two entire original-study files."""

import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from coloc_common import ROOT, save_json, sha


def main():
    records = []
    for study, suffix in [('dice', 'vcf'), ('blueprint', 'txt')]:
        manifest = json.loads((ROOT / ('config/prkd2-recovery-' + study + '-stream.json')).read_text())
        archive = ROOT / manifest['archive']
        digest = hashlib.sha256()
        with archive.open('rb') as stream:
            while block := stream.read(16 * 1024 * 1024):
                digest.update(block)
        if digest.hexdigest() != manifest['compressed_sha256']:
            raise ValueError('Original archive checksum changed')
        unzip = subprocess.Popen(['gzip', '-dc', str(archive)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        selector = subprocess.Popen(['rg', '--text', '--fixed-strings', '-e', 'ENSG00000105287.', '-e', 'rs112445263', '-'],
                                    stdin=unzip.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        unzip.stdout.close()
        selected, selector_error = selector.communicate()
        _, unzip_error = unzip.communicate()
        if unzip.returncode != 0 or selector.returncode != 0:
            raise ValueError('Independent whole-archive reader failed: ' + (unzip_error + selector_error).decode())
        gene, variant = [], []
        for line in selected.splitlines(keepends=True):
            fields = line.decode().split()
            if fields[0].startswith('#'):
                continue
            if study == 'dice':
                info = dict(x.split('=', 1) for x in fields[7].split(';'))
                gene_id, rsid = info['Gene'].split('.')[0], fields[2]
            else:
                gene_id, rsid = fields[2].split('.')[0], fields[1]
            if gene_id == 'ENSG00000105287':
                gene.append(line)
            if rsid == 'rs112445263':
                variant.append(line)
        for index, rows in enumerate([gene, variant]):
            expected = manifest['outputs'][index]
            body = b''.join(line for line in gzip.decompress((ROOT / expected['file']).read_bytes()).splitlines(keepends=True) if not line.startswith(b'#'))
            if b''.join(rows) != body or len(rows) != expected['rows']:
                raise ValueError('Independent original-source extraction disagrees: ' + study)
        records.append({'study': study, 'archive_sha256': digest.hexdigest(), 'gene_rows': len(gene), 'variant_rows_all_genes': len(variant),
                        'all_source_rows_reread': True, 'gzip_eof_and_crc_verified': True, 'exact_extraction_match': True})
        print(json.dumps(records[-1]), flush=True)
    save_json(ROOT / 'reports/prkd2-recovery-archive-verification.json', {
        'all_checks_pass': True, 'method': 'Independent external gzip decompression and ripgrep selection, followed by exact identifier parsing; every compressed archive reread through EOF',
        'completed_at_utc': datetime.now(timezone.utc).isoformat(), 'verifier_sha256': sha(Path(__file__).read_bytes()), 'studies': records})


if __name__ == '__main__':
    main()
