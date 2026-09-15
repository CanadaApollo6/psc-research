"""Recover only the four frozen gene/conditional-lead ZIP members, with CRC checks."""
import argparse
import csv
import json
import struct
import zlib
from pathlib import Path

from fetch_ibdverse_model_audit import ROOT, RAW, freeze, save, sha
from fetch_coverage_repair_parquet import RangeReader


def decode_member(data, member):
    decoder = zlib.decompressobj(-15)
    plain = decoder.decompress(data) + decoder.flush()
    if not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
        raise ValueError('ZIP member did not finish at exact DEFLATE boundary')
    if len(plain) != int(member['uncompressed_bytes']):
        raise ValueError('ZIP uncompressed size differs')
    if zlib.crc32(plain) != int(member['crc32']):
        raise ValueError('ZIP CRC differs')
    return plain


def extract_gene(plain, gene):
    lines = plain.decode().splitlines()
    names = lines[0].split('\t')
    if len(set(names)) != len(names) or 'phenotype_id' not in names or 'variant_id' not in names:
        raise ValueError('Unexpected source columns')
    keep = []
    for line in lines[1:]:
        values = line.split('\t')
        if len(values) != len(names):
            raise ValueError('Malformed source row')
        if values[names.index('phenotype_id')].split('.')[0] == gene:
            keep.append(line)
    return ('\n'.join([lines[0], *keep]) + '\n').encode(), len(lines) - 1, len(keep)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    frozen = freeze()
    plan = json.loads((ROOT / 'config/ibdverse-model-audit-plan.json').read_text())
    inventory = ROOT / 'data/derived/prkd2-recovery-ibdverse-archive-inventory.csv'
    if sha(inventory) != frozen['prior_tracked_sha256'][str(inventory.relative_to(ROOT))]:
        raise ValueError('Previous complete ZIP inventory changed')
    contexts = [r['source_context'] for r in plan['contexts']]
    names = ['Cis_eqtls_qval.tsv', 'Cis_eqtls_independent.tsv']
    targets = [r for r in csv.DictReader(inventory.open())
               if r['filename'].split('/')[0] in contexts and Path(r['filename']).name in names]
    if len(targets) != 4:
        raise ValueError('Frozen four-member selection changed')
    parent = json.loads((ROOT / 'config/prkd2-recovery-ibdverse-zip-metadata.json').read_text())['queries'][0]
    path = ROOT / 'config/ibdverse-model-audit-member-queries.json'
    pins = {'plan_sha256': sha(ROOT / 'config/ibdverse-model-audit-plan.json'), 'inventory_sha256': sha(inventory)}
    record = json.loads(path.read_text()) if path.exists() else {**pins, 'queries': []}
    if any(record[k] != v for k, v in pins.items()):
        raise ValueError('Frozen ZIP selection changed')
    for member in targets:
        context, name = member['filename'].split('/')[0], Path(member['filename']).stem
        identifier = context + '-' + name
        query = next((r for r in record['queries'] if r['id'] == identifier), None)
        if query is None:
            query = {'id': identifier, 'member': member, 'ranges': [],
                     **{k: parent[k] for k in ['url', 'remote_size', 'etag', 'last_modified']}}
            record['queries'].append(query)
        if query['member'] != member:
            raise ValueError('Member definition changed')
        reader = RangeReader(query, record, path, args.offline)
        reader.directory = RAW / 'zip-ranges'
        reader.directory.mkdir(parents=True, exist_ok=True)
        offset, length = int(member['local_header_offset']), int(member['compressed_bytes'])
        if length > plan['retrieval']['maximum_new_selected_archive_member_compressed_bytes']:
            raise ValueError('Member exceeds frozen transfer cap')
        fields = struct.unpack('<4s5H3L2H', reader.fetch(offset, 30))
        if fields[0] != b'PK\x03\x04' or fields[3] != 8 or fields[2] & 1:
            raise ValueError('Invalid, unsupported or encrypted ZIP header')
        name_length, extra_length = fields[-2:]
        actual_name = reader.fetch(offset + 30, name_length + extra_length)[:name_length].decode()
        if actual_name != member['filename']:
            raise ValueError('Local and central ZIP filenames differ')
        start = offset + 30 + name_length + extra_length
        plain = decode_member(reader.fetch(start, length), member)
        whole_path = RAW / (identifier + '.tsv')
        whole_path.write_bytes(plain)
        selected, count, n_selected = extract_gene(plain, plan['target_gene'])
        output = ROOT / 'data/derived/ibdverse-model-audit' / (identifier + '-PRKD2.tsv')
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() and output.read_bytes() != selected:
            raise ValueError('Source-row replay differs')
        output.write_bytes(selected)
        result = {'source_rows': count, 'source_bytes': len(plain), 'source_sha256': sha(whole_path),
                  'source_file': str(whole_path.relative_to(ROOT)), 'crc32_verified': True,
                  'deflate_eof_verified': True, 'selected_rows': n_selected,
                  'file': str(output.relative_to(ROOT)), 'sha256': sha(output),
                  'transfer_bytes': sum(r['bytes'] for r in query['ranges'])}
        if 'result' in query and query['result'] != result:
            raise ValueError('Member result replay differs')
        query['result'] = result
        save(path, record)
        print(json.dumps({'id': identifier, **result}), flush=True)


if __name__ == '__main__':
    main()
