"""Map complete published component vectors by independently checked full edits."""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from prkd2_common import ROOT, save_json, sha
from prkd2_expanded_common import load_references, qtl_identity
from prepare_coloc_inputs import save_frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', choices=['QTD000021', 'QTD000504'])
    args = parser.parse_args()
    ds = args.dataset
    receipt = json.loads((ROOT / ('config/prkd2-bf-query-' + ds + '.json')).read_text())
    assert receipt['archive_eof_verified']
    raw = (ROOT / receipt['file']).read_bytes()
    assert sha(raw) == receipt['sha256'] and sha(gzip.decompress(raw)) == receipt['uncompressed_sha256']
    rows = pd.read_csv(ROOT / receipt['file'], sep='\t', dtype=str, keep_default_na=False)
    if set(rows.molecular_trait_id) != {'ENSG00000105287'} or rows.variant.duplicated().any():
        raise ValueError('Incorrect or duplicated full molecular trait vectors')
    reference = load_references()['GRCh38']
    records = []
    for row in rows.to_dict('records'):
        try:
            identity = qtl_identity(row['variant'], reference)
            records.append({**row, 'snp': identity['snp'], 'normalization_status': 'mapped'})
        except ValueError as exc:
            records.append({**row, 'snp': 'unmatched_source:' + row['variant'], 'normalization_status': str(exc)})
    duplicate = {key for key, count in Counter(r['snp'] for r in records).items() if count > 1}
    for row in records:
        if row['snp'] in duplicate:
            row.update(snp='unmatched_equivalent:' + row['variant'], normalization_status='duplicate_normalized_full_vector_identity')
    mapped = pd.DataFrame(records)
    nominal = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-qtl-alignment.csv.gz', dtype=str, keep_default_na=False)
    nominal = nominal[nominal.dataset_id == ds].set_index('source_variant')
    if set(mapped.variant) != set(nominal.index):
        raise ValueError('Full BF variants differ from the independently normalized nominal variant universe')
    for row in records:
        expected = nominal.loc[row['variant']]
        if expected['status'] != 'eligible' or expected['snp'] != row['snp']:
            raise ValueError('Published and nominal normalization mapping disagrees')
    # Keep published numeric strings unchanged; no renormalization or posterior fitting here.
    for column in rows.columns:
        if not mapped[column].equals(rows[column]):
            raise ValueError('An original published column changed')
    out = ROOT / 'data/derived/prkd2-expanded-published-bfs'
    out.mkdir(exist_ok=True)
    output = out / (ds + '.csv.gz')
    save_frame(mapped, output)
    record = {'dataset_id': ds, 'rows': len(mapped), 'full_vector_source_sha256': receipt['sha256'],
              'all_original_columns_preserved': True, 'same_complete_variant_universe_as_nominal_source': True,
              'all_mappings_agree_with_independently_checked_nominal_edits': True,
              'normalization_status_counts': mapped.normalization_status.value_counts().to_dict(),
              'output_file': str(output.relative_to(ROOT)), 'output_sha256': sha(output.read_bytes()),
              'source_manifest_sha256': sha((ROOT / 'config/prkd2-expanded-sources.json').read_bytes()),
              'normalization_verification_sha256': sha((ROOT / 'reports/prkd2-normalization-verification.json').read_bytes()),
              'mapper_sha256': sha(Path(__file__).read_bytes())}
    save_json(ROOT / ('config/prkd2-expanded-bf-mapping-' + ds + '.json'), record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
