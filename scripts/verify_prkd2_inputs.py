"""Independently check raw association fields, risk signs and all signed LD entries."""

import csv
import gzip
import json
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri

from prkd2_common import ROOT, save_json, sha


def near(a, b, name):
    if not np.allclose(np.asarray(a, float), np.asarray(b, float), atol=2e-12, rtol=2e-12):
        raise ValueError('Independent input check failed: ' + name)


def main():
    plan = json.loads((ROOT / 'config/prkd2-plan.json').read_text())
    ns = {'x': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(ROOT / 'data/raw/coloc-Ji2017-MOESM245_ESM.xlsx') as archive:
        tree = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
    cells = {c.get('r'): float(c.find('x:v', ns).text) for c in tree.findall('.//x:c', ns)
             if c.get('r') in ['B13', 'C13', 'D13', 'E13', 'F13', 'G13']}
    assert cells['B13'] + cells['C13'] == plan['platform_counts']['omni']['cases']
    assert cells['D13'] == plan['platform_counts']['omni']['controls']
    assert sum(cells[k] for k in ['E13', 'F13', 'G13']) == plan['platform_counts']['affy']['N']
    assert sum(plan['platform_counts'][p]['cases'] for p in ['omni', 'affy']) == 2871
    assert sum(plan['platform_counts'][p]['controls'] for p in ['omni', 'affy']) == 12019

    source = pd.read_csv(ROOT / 'data/derived/prkd2-query-rows/GWAS-PRKD2.txt.gz', sep=r'\s+').rename(columns={'#chr': 'chr'})
    full = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz').sort_values('source_row')
    assert list(full.source_row) == list(range(len(source))) and len(source) == 4727
    for column in source.columns:
        if pd.api.types.is_numeric_dtype(source[column]):
            near(source[column], full[column], 'full source GWAS field ' + column)
        else:
            assert list(source[column]) == list(full[column]), column
    g = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-PRKD2.csv.gz')
    for row in g.to_dict('records'):
        original = source.iloc[int(row['source_row'])]
        counts = plan['platform_counts'][original.platform]
        near([row['pvalue'], row['MAF'], row['N'], row['cases'], row['controls'], row['case_fraction']],
             [original.p, min(original.freq_1_controls, 1-original.freq_1_controls), counts['N'], counts['cases'],
              counts['controls'], counts['cases']/counts['N']], 'prepared GWAS row')
    qtl_checks = {}
    for ds, n in [('QTD000021', 191), ('QTD000504', 91)]:
        raw = pd.read_csv(ROOT / ('data/derived/prkd2-query-rows/' + ds + '-PRKD2.tsv.gz'), sep='\t')
        numeric = ['position', 'beta', 'se', 'pvalue', 'maf', 'an', 'ac']
        assert (raw.groupby('variant')[numeric].nunique() == 1).all().all()
        unique = raw.drop_duplicates('variant').set_index('variant')
        fullq = pd.read_csv(ROOT / ('data/derived/prkd2-expanded-inputs/' + ds + '-full-evidence.csv.gz'))
        assert set(fullq.source_variant) == set(unique.index)
        original = unique.loc[fullq.source_variant]
        near(fullq.beta, original.beta, 'unoriented full RNA beta')
        near(fullq.varbeta, original.se ** 2, 'RNA variance')
        near(fullq.pvalue, original.pvalue, 'RNA P')
        near(fullq.MAF, original.maf, 'RNA MAF')
        near(fullq.N, original.an / 2, 'RNA sample count')
        assert set(fullq.N) == {n}
        q = pd.read_csv(ROOT / ('data/derived/prkd2-expanded-inputs/' + ds + '-PRKD2.csv.gz'))
        for row in q.to_dict('records'):
            original = unique.loc[row['source_variant']]
            effect_sign = 1
            if len(original.ref) == len(original.alt) == 1:
                effect_sign = 1 if original.alt == max(original.ref, original.alt) else -1
            near(row['beta'], original.beta * effect_sign, 'RNA effect orientation')
        qtl_checks[ds] = {'all_source_rows': len(raw), 'distinct_full_variants': len(unique), 'prepared_variants': len(q)}

    # Re-read donor dosages directly; do not reuse the preparation matrices or
    # its orientation function as input to this independent correlation check.
    lookup = defaultdict(list)
    with gzip.open(ROOT / 'data/raw/prkd2-expanded-EUR-dosages.tsv.gz', 'rt') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        assert len(header[5:]) == len(set(header[5:])) == 503
        for row in reader:
            lookup[(int(row[1]), tuple(sorted(row[3:5])))].append(row)
    ga = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-gwas-alignment.csv.gz').set_index('source_row')
    ldrows = pd.read_csv(ROOT / 'data/derived/prkd2-expanded-inputs/GWAS-PRKD2-ld.csv.gz')
    arrays = []
    for row in ldrows.to_dict('records'):
        original = source.iloc[int(row['source_row'])]
        matches = lookup[(int(original.pos), tuple(sorted([original.allele_0, original.allele_1])))]
        assert len(matches) == 1
        reference = matches[0]
        values = np.array(reference[5:], dtype=float)
        assert np.isin(values, [0, 1, 2]).all()
        audit = ga.loc[int(row['source_row'])]
        _, ref, alt = json.loads(audit.normalized_grch38)
        allele_sign = -1 if len(ref) == len(alt) == 1 and alt < ref else 1
        near(row['canonical_sign'], allele_sign * (1 if original.allele_1 == reference[4] else -1), 'PSC risk orientation')
        arrays.append(values if allele_sign == 1 else 2-values)
    calculated = np.corrcoef(np.asarray(arrays))
    observed = np.fromfile(ROOT / 'data/raw/prkd2-expanded-ld/GWAS-PRKD2.f64', dtype='<f8').reshape(len(arrays), len(arrays))
    max_ld_error = float(np.max(np.abs(calculated-observed)))
    assert max_ld_error < 2e-12
    for n in [11386, 3504, 14890]:
        diagnostics = pd.read_csv(ROOT / f'reports/prkd2-expanded-LD-diagnostics-N{n}.csv').set_index('snp').loc[ldrows.snp]
        raw_z = -ndtri(ldrows.pvalue.to_numpy()/2)*ldrows.canonical_sign.to_numpy()
        adjusted_z = raw_z * np.sqrt((n-1)/(n-2+raw_z**2))
        near(diagnostics.z, adjusted_z, 'R finite-N-adjusted signed Z')
    record = {'all_checks_pass': True, 'publisher_table_cells': cells, 'GWAS_source_rows_all_fields_checked': len(source),
        'GWAS_prepared_rows_checked': len(g), 'QTL_checks': qtl_checks, 'reference_donors': 503,
        'independently_reconstructed_LD_variants': len(arrays), 'LD_entries_checked': len(arrays)**2,
        'maximum_LD_difference': max_ld_error, 'signed_Z_vectors_checked': 3,
        'kriging_Z_convention': 'Finite-N adjustment z * sqrt((n-1)/(n-2+z*z)); the original fit input remains the raw signed Z',
        'method': 'Raw-source pandas/CSV joins, publisher XLSX XML cells, independent direct-dosage correlation and allele orientation',
        'verifier_sha256': sha(Path(__file__).read_bytes())}
    save_json(ROOT / 'reports/prkd2-input-verification.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
