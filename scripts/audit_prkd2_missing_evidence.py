"""Account for omitted association weight and directly query dominant GWAS gaps."""

import csv
import gzip
import json
import time
import zlib
from datetime import datetime, timezone

import pandas as pd
import pysam

from prkd2_common import ROOT, gzip_bytes, save_json, sha
from verify_prkd2_statistics import binary_bf, quantitative_bf, weights
from prkd2_expanded_common import load_references
from prepare_coloc_inputs import save_frame


def main():
    folder = ROOT / 'data/derived/prkd2-expanded-inputs'
    g = pd.read_csv(folder / 'GWAS-PRKD2.csv.gz')
    gf = pd.read_csv(folder / 'GWAS-full-evidence.csv.gz')
    gf['full_trait_bf_weight'] = weights(binary_bf(gf))
    names = dict(zip(g.source_row, g.snp))
    gf['matched_edit'] = gf.source_row.map(names).fillna('')
    rows, query_plan = [], []
    reference = load_references()['GRCh38']
    for ds in ['QTD000021', 'QTD000504']:
        q = pd.read_csv(folder / (ds + '-PRKD2.csv.gz'))
        qf = pd.read_csv(folder / (ds + '-full-evidence.csv.gz'))
        qf['full_trait_bf_weight'] = weights(quantitative_bf(qf, 1))
        qf['matched_edit'] = qf.source_variant.map(dict(zip(q.source_variant, q.snp))).fillna('')
        common = set(g.snp) & set(q.snp)
        for trait, frame in [('GWAS', gf), ('QTL', qf)]:
            missing = frame[~frame.matched_edit.isin(common)].sort_values('full_trait_bf_weight', ascending=False)
            for rank, row in enumerate(missing.head(20).to_dict('records'), 1):
                item = {'dataset_id': ds, 'trait': trait, 'omitted_rank': rank, 'full_trait_bf_weight': row['full_trait_bf_weight'],
                        'matched_edit': row['matched_edit'], 'source_id': row['SNP'] if trait == 'GWAS' else row['source_variant'],
                        'pvalue': row['pvalue']}
                rows.append(item)
            # Explicit rule: all omitted, already identity-qualified GWAS variants with >=5% full-trait weight.
            if trait == 'GWAS':
                for row in missing[missing.full_trait_bf_weight >= .05].to_dict('records'):
                    key = row['matched_edit'].split(':')
                    if len(key) != 4 or len(key[2]) != 1 or len(key[3]) != 1:
                        raise ValueError('Dominant targeted source gap needs a full-allele request definition')
                    ref = reference.get(int(key[1]))
                    if ref not in key[2:]:
                        raise ValueError('Target reference not in allele pair')
                    alt = next(a for a in key[2:] if a != ref)
                    query_plan.append({'dataset_id': ds, 'source_id': row['SNP'], 'chromosome': '19', 'position': int(key[1]),
                                       'ref': ref, 'alt': alt, 'full_gwas_bf_weight': row['full_trait_bf_weight']})
    save_frame(pd.DataFrame(rows), ROOT / 'data/derived/prkd2-missing-evidence.csv')
    plan_path = ROOT / 'config/prkd2-missing-evidence-query-plan.json'
    if not plan_path.exists():
        save_json(plan_path, {'frozen_at_utc': datetime.now(timezone.utc).isoformat(), 'selection': 'All identity-qualified omitted GWAS variants with >=5% full-source marginal BF weight after the frozen expansion', 'queries': query_plan})
    else:
        assert json.loads(plan_path.read_text())['queries'] == query_plan
    sources = {s['id']: s for s in json.loads((ROOT / 'config/prkd2-sources.json').read_text())['sources']}
    record_path = ROOT / 'config/prkd2-missing-evidence-queries.json'
    if record_path.exists():
        record = json.loads(record_path.read_text())
        for query in record['queries']:
            assert sha((ROOT / query['file']).read_bytes()) == query['sha256']
        print(json.dumps({'status': 'verified-cache', 'queries': len(record['queries'])}))
        return
    receipts = []
    for request in query_plan:
        ds = request['dataset_id']
        source = sources[ds + '-index']
        url = source['url'].removesuffix('.tbi')
        header = zlib.decompressobj(31).decompress((ROOT / 'data/raw' / sources[ds + '-header']['file']).read_bytes()).decode().splitlines()[0]
        with pysam.TabixFile(url, index=str(ROOT / 'data/raw' / source['file'])) as table:
            lines = list(table.fetch('19', request['position'] - 1, request['position']))
        parsed = [dict(zip(header.split('\t'), line.split('\t'), strict=True)) for line in lines]
        exact = [r for r in parsed if r['ref'] == request['ref'] and r['alt'] == request['alt']]
        plain = (header + '\n' + '\n'.join(lines) + ('\n' if lines else '')).encode()
        path = ROOT / ('data/derived/prkd2-query-rows/' + ds + '-dominant-missing-variant.tsv.gz')
        path.write_bytes(gzip_bytes(plain))
        receipts.append({**request, 'queried_at_utc': datetime.now(timezone.utc).isoformat(), 'url': url,
                         'rows_at_position_all_genes': len(parsed), 'rows_exact_alleles_all_genes': len(exact),
                         'rows_exact_alleles_PRKD2': sum(r['gene_id'].split('.')[0] == 'ENSG00000105287' for r in exact),
                         'observed_allele_pairs': sorted({r['ref'] + '>' + r['alt'] for r in parsed}),
                         'file': str(path.relative_to(ROOT)), 'sha256': sha(path.read_bytes()), 'uncompressed_sha256': sha(plain)})
        time.sleep(2)
    save_json(record_path, {'query_plan_sha256': sha(plan_path.read_bytes()), 'queries': receipts,
                           'filter_cause': 'Not inferred from absence in aggregate files', 'reader_sha256': sha(__import__('pathlib').Path(__file__).read_bytes())})
    print(json.dumps(receipts, indent=2))


if __name__ == '__main__':
    main()
