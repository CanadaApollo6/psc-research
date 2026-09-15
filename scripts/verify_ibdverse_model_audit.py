"""Independent, offline verification of the bounded IBDverse public-source audit.

Uses standard ZIP handling and a separate categorical-code Counter aggregation.
It does not import the primary member, metadata or analysis implementation.
"""
import csv
import difflib
import gzip
import hashlib
import io
import json
import subprocess
import zipfile
import zlib
from collections import Counter
from pathlib import Path
from urllib.parse import unquote, urlsplit

import h5py
import numpy as np

from verify_prkd2_recovery_ibdverse import CachedArchive

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/ibdverse-model-audit'
OUT = ROOT / 'data/derived/ibdverse-model-audit'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads((ROOT / path).read_text())


def source_rows(path, separator='\t'):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=separator))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def verify_sources():
    sources = read_json('config/ibdverse-model-audit-sources.json')['sources']
    checked = set()
    for source in sources:
        for key, hash_key in [('file', 'sha256'), ('receipt', 'receipt_sha256')]:
            if source.get(key):
                path = ROOT / source[key]
                check(digest(path.read_bytes()) == source[hash_key], 'Source hash differs: ' + source['id'])
                checked.add(source[key])
    trees = {}
    for source in sources:
        if source.get('ok') and source['id'].endswith('-tree'):
            tree = read_json(source['file'])
            check(not tree['truncated'], 'Incomplete code tree')
            parts = urlsplit(source['url']).path.split('/')
            repository, commit = '/'.join(parts[2:4]), parts[6]
            trees[(repository, commit)] = {r['path']: r for r in tree['tree']}
    blobs = 0
    for source in sources:
        if not source.get('ok') or urlsplit(source.get('url', '')).netloc != 'raw.githubusercontent.com':
            continue
        parts = unquote(urlsplit(source['url']).path).split('/')
        repository, commit, filename = '/'.join(parts[1:3]), parts[3], '/'.join(parts[4:])
        if (repository, commit) not in trees:
            continue
        blob = (ROOT / source['file']).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(blob)).encode() + b'\0' + blob).hexdigest()
        check(actual == trees[(repository, commit)][filename]['sha'], 'Git blob differs: ' + filename)
        blobs += 1
    return {'unique_source_or_receipt_files_checked': len(checked), 'Git_blobs_matching_pinned_trees': blobs,
            'complete_pinned_repository_trees': len(trees)}


def verify_members():
    parent = read_json('config/prkd2-recovery-ibdverse-zip-metadata.json')['queries'][0]
    queries = read_json('config/ibdverse-model-audit-member-queries.json')['queries']
    ranges = parent['ranges'] + [r for q in queries for r in q['ranges']]
    outcomes = []
    with zipfile.ZipFile(CachedArchive(ranges, parent['remote_size'])) as archive:
        inventory = source_rows(ROOT / 'data/derived/prkd2-recovery-ibdverse-archive-inventory.csv', ',')
        check(len(archive.infolist()) == len(inventory), 'Full base-eQTL inventory count differs')
        for item, row in zip(archive.infolist(), inventory):
            expected = (row['filename'], int(row['uncompressed_bytes']), int(row['compressed_bytes']),
                        int(row['compression_method']), int(row['local_header_offset']), int(row['crc32']))
            check((item.filename, item.file_size, item.compress_size, item.compress_type, item.header_offset, item.CRC) == expected,
                  'Full base-eQTL inventory metadata differs')
        for query in queries:
            info = archive.getinfo(query['member']['filename'])
            # Reading the entire member also triggers standard-library ZIP CRC verification.
            plain = archive.read(info)
            expected = query['result']
            check(digest(plain) == expected['source_sha256'], 'Independent full-member digest differs')
            lines = plain.splitlines()
            check(len(lines) - 1 == expected['source_rows'], 'Independent full-member row count differs')
            names = lines[0].split(b'\t')
            gene_index = names.index(b'phenotype_id')
            selected = []
            for line in lines[1:]:
                fields = line.split(b'\t')
                check(len(fields) == len(names), 'Malformed source row during independent scan')
                if fields[gene_index].partition(b'.')[0] == b'ENSG00000105287':
                    selected.append(line)
            extracted = b'\n'.join([lines[0], *selected]) + b'\n'
            check(extracted == (ROOT / expected['file']).read_bytes(), 'Independent PRKD2 extraction differs')
            check(len(selected) == expected['selected_rows'], 'Independent target count differs')
            outcomes.append({'id': query['id'], 'source_rows': len(lines) - 1, 'PRKD2_rows': len(selected),
                             'full_source_sha256': digest(plain), 'standard_ZIP_CRC_checked': True})
    return {'complete_archive_inventory_members': len(inventory), 'members': outcomes}


def verify_coloc_headers():
    ledger = read_json('config/ibdverse-model-audit-coloc-queries.json')
    query = ledger['queries'][0]
    ranges = query['ranges']
    results = []
    with zipfile.ZipFile(CachedArchive(ranges, query['remote_size'])) as archive:
        inventory = [{'filename': item.filename, 'uncompressed_bytes': item.file_size,
                      'compressed_bytes': item.compress_size, 'compression_method': item.compress_type,
                      'local_header_offset': item.header_offset, 'crc32': item.CRC} for item in archive.infolist()]
        actual = source_rows(ROOT / ledger['inventory']['file'], ',')
        check([{k: str(v) for k, v in row.items()} for row in inventory] == actual, 'Coloc inventory differs')
        rebuilt = io.StringIO(newline='')
        writer = csv.DictWriter(rebuilt, fieldnames=list(inventory[0]))
        writer.writeheader()
        writer.writerows(inventory)
        data = rebuilt.getvalue().encode()
        if digest(data) != ledger['inventory']['sha256']:
            data = data.replace(b'\r\n', b'\n')
        check(digest(data) == ledger['inventory']['sha256'], 'Coloc inventory byte reconstruction differs')
        for expected in ledger['member_headers']:
            # Standard ZIP inflates only the nested gzip prefix. This is not a full-file CRC check.
            with archive.open(expected['member']) as stream:
                prefix = stream.read(16384)
            text = zlib.decompressobj(31).decompress(prefix, 2_000_000)
            check(b'\n' in text, 'Incomplete nested gzip header')
            columns = next(csv.reader(io.StringIO(text.split(b'\n', 1)[0].decode()), delimiter='\t'))
            check(columns == expected['columns'], 'Coloc header differs')
            check(expected['complete_member_scanned'] is False, 'Partial stream mislabeled as complete')
            results.append({'member': expected['member'], 'columns': columns, 'complete_member_scanned': False})
    return results


def verify_rna_counts():
    plan = read_json('config/ibdverse-model-audit-h5ad-plan.json')
    ledger = read_json('config/ibdverse-model-audit-h5ad-queries.json')
    amendment = read_json('config/ibdverse-model-audit-h5ad-mapping-amendment.json')
    check(digest((ROOT / 'config/ibdverse-model-audit-h5ad-plan.json').read_bytes()) == amendment['original_plan_sha256'], 'Amended plan differs')
    check(digest((ROOT / 'reports/ibdverse-model-audit-h5ad-schema.json').read_bytes()) == amendment['observed_schema_sha256'], 'Frozen schema differs')
    mapping_rows = source_rows(ROOT / amendment['mapping_file'])
    sample_to_individual = {r['sanger_sample_id']: r['individual_id'] for r in mapping_rows}
    check(len(sample_to_individual) == len(mapping_rows), 'Duplicate sample mapping')
    records = [{**r, 'range': f"{r['start']}-{r['start'] + r['bytes'] - 1}"} for r in ledger['blocks']]
    names = ['sanger_sample_id', 'predicted_labels', 'tissue']
    with h5py.File(CachedArchive(records, plan['source_bytes']), 'r') as source:
        categories, codes = [], []
        for name in names:
            node = source['obs'][name]
            check(set(node.keys()) == {'codes', 'categories'}, 'Expected categorical encoding changed')
            labels = [v.decode() if isinstance(v, bytes) else str(v) for v in node['categories'][:]]
            indices = node['codes'][:]
            check(np.all(indices >= 0) and np.all(indices < len(labels)), 'Missing/invalid category in observed counts')
            categories.append(labels)
            codes.append(indices)
        check(len({len(v) for v in codes}) == 1, 'Metadata row counts differ')
        grouped = Counter(zip(*codes))
    target_counts = {name: Counter() for name in plan['counts_only_for']}
    unmatched = 0
    for (sample, cell_type, tissue), count in grouped.items():
        individual = sample_to_individual.get(categories[0][sample])
        if individual is None:
            unmatched += count
            continue
        context = categories[1][cell_type] + '_' + categories[2][tissue]
        if context in target_counts:
            target_counts[context][individual] += count
    expected = read_json('reports/ibdverse-model-audit-h5ad-metadata.json')
    check(unmatched == expected['unmapped_cells_all_annotations'] == 0, 'Unmatched observation mapping')
    check(sum(grouped.values()) == expected['observation_rows'], 'Total observation count differs')
    outcomes = []
    for row in expected['counts']:
        counts = target_counts[row['context']]
        actual = {'context': row['context'], 'cells': sum(counts.values()), 'donor_labels_with_any_cells': len(counts),
                  'donor_labels_with_at_least_five_cells': sum(n >= plan['minimum_cells_per_donor'] for n in counts.values())}
        check(all(row[k] == v for k, v in actual.items()), 'Independent RNA-eligibility aggregation differs')
        outcomes.append(actual)
    check(expected['actual_executed_model_N_recovered'] is False and expected['actual_selected_PCs_recovered'] is False,
          'RNA counts upgraded into observed model metadata')
    return {'method': 'Categorical-code tuple Counter then public sample map; primary implementation uses expanded pandas columns and groupby.',
            'observation_rows': sum(grouped.values()), 'unmatched_cells': unmatched, 'counts': outcomes,
            'actual_executed_model_N_verified': False}


def verify_source_fields():
    summary = source_rows(OUT / 'PRKD2-source-model-summary.csv', ',')
    differences = source_rows(OUT / 'lead-source-field-comparison.csv', ',')
    common_fields, comparisons, nonidentical = 0, 0, []
    for row in summary:
        context = row['context']
        gene = source_rows(OUT / (context + '-Cis_eqtls_qval-PRKD2.tsv'))[0]
        number = context.split('Myeloid_')[1].split('_')[0]
        with gzip.open(ROOT / f'data/derived/prkd2-recovery-queries/IBDverse-original-Myeloid-{number}-PRKD2.tsv.gz', 'rt') as stream:
            nominal = [r for r in csv.DictReader(stream, delimiter='\t') if r['variant_id'] == gene['variant_id']]
        check(len(nominal) == 1, 'Prior nominal lead identity is not unique')
        for comparison in (r for r in differences if r['context'] == context):
            field = comparison['field']
            a, b = gene[field], nominal[0][field]
            check(a == comparison['gene_summary_original'] and b == comparison['nominal_vector_original'], 'Source numeric string was changed')
            difference = abs(float(a) - float(b))
            relative = difference / abs(float(b))
            check(difference == float(comparison['absolute_difference']) and relative == float(comparison['relative_difference']), 'Source difference was not retained exactly')
            comparisons += 1
            if difference:
                nonidentical.append({'context': context, 'field': field, 'absolute_difference': difference,
                                     'relative_difference': relative, 'harmonized_or_overwritten': False})
        independent = source_rows(OUT / (context + '-Cis_eqtls_independent-PRKD2.tsv'))
        check(len(independent) == int(row['selected_independent_rows']), 'Conditional hit count differs')
        lead = next(r for r in independent if r['variant_id'] == gene['variant_id'])
        for key in set(gene) & set(lead):
            check(gene[key] == lead[key], 'Gene and independent-hit field disagreement: ' + key)
            common_fields += 1
        check(row['actual_executed_model_N_observed'] == row['actual_selected_PCs_observed'] == 'False',
              'Conditional PC inference upgraded into observed model metadata')
    check(comparisons == len(differences), 'Uncompared source fields')
    return {'lead_nominal_fields_compared': comparisons, 'gene_and_conditional_hit_fields_identical': common_fields,
            'nonidentical_gene_summary_vs_nominal_fields_preserved': nonidentical}


def verify_previous_outputs():
    frozen = read_json('config/ibdverse-model-audit-freeze.json')
    unchanged, navigation = 0, []
    for filename, expected in frozen['prior_tracked_sha256'].items():
        actual = (ROOT / filename).read_bytes()
        if digest(actual) == expected:
            unchanged += 1
            continue
        check(filename in {'README.md', 'docs/evidence-log.md'}, 'Prior output changed: ' + filename)
        original = subprocess.check_output(['git', 'show', frozen['baseline_commit'] + ':' + filename], cwd=ROOT)
        check(digest(original) == expected, 'Frozen baseline differs from Git')
        operations = difflib.SequenceMatcher(a=original.splitlines(), b=actual.splitlines(), autojunk=False).get_opcodes()
        check(all(tag in {'equal', 'insert'} for tag, *_ in operations), 'Navigation edit replaced prior content')
        navigation.append(filename)
    return {'baseline_commit': frozen['baseline_commit'], 'prior_files_unchanged': unchanged,
            'navigation_files_with_insertions_only': navigation, 'prior_scientific_content_preserved': True}


def main():
    result = {'all_checks_pass': True, 'verifier_sha256': digest(Path(__file__).read_bytes()),
              'source_hashes_and_code': verify_sources(), 'base_eQTL_archive': verify_members(),
              'coloc_headers_only': verify_coloc_headers(), 'RNA_metadata': verify_rna_counts(),
              'source_field_checks': verify_source_fields(), 'prior_preservation': verify_previous_outputs()}
    summary = read_json('reports/ibdverse-model-audit.json')
    check(summary['ready_for_new_shared_cause_posterior'] is False and summary['new_H4_posteriors'] == 0, 'Audit readiness gate changed')
    for row in summary['outputs']:
        check(digest((ROOT / row['file']).read_bytes()) == row['sha256'], 'Derived output hash differs')
    result['new_H4_posteriors'] = 0
    path = ROOT / 'reports/ibdverse-model-audit-verification.json'
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
