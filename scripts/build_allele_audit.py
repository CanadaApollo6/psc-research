"""Derive a separate allele audit without changing the frozen pilot inputs."""

import csv
import json
from pathlib import Path

from fetch_gwas_audit_rows import complete_rows
from fetch_sources import validate

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifests = [json.loads((ROOT / 'config' / name).read_text()) for name in ['association-sources.json', 'mechanism-sources.json']]
    sources = {row['id']: row for m in manifests for row in m['sources']}
    for source in sources.values():
        validate((ROOT / 'data/raw' / source['file']).read_bytes(), source)
    original = (ROOT / 'data/derived/gwas-audit-original-rows.txt').read_text()
    source_rows = {}
    for source in sources.values():
        if source['id'].startswith('GCST004030-bytes-'):
            for coordinate, rsid, line in complete_rows((ROOT / 'data/raw' / source['file']).read_bytes()):
                if rsid in {'rs1893592', 'rs313839'}:
                    source_rows[rsid] = (line, source['id'])
    input_rows = list(csv.DictReader(original.splitlines(), delimiter=' '))
    output_rows = []
    for row in input_rows:
        line, source_id = source_rows[row['SNP']]
        if line != ' '.join(row.values()):
            raise ValueError('Derived GWAS row is not an exact copy of a pinned source line')
        if float(row['or']) <= 1 or float(row['freq_1_cases']) <= float(row['freq_1_controls']):
            raise ValueError('The risk-allele interpretation requires review')
        output_rows.append({
            'rsid': row['SNP'], 'assembly': 'GRCh37', 'chromosome': row['#chr'], 'position_1based': row['pos'],
            'risk_allele_original_gwas': row['allele_1'], 'other_allele_original_gwas': row['allele_0'],
            'risk_frequency_cases': row['freq_1_cases'], 'risk_frequency_controls': row['freq_1_controls'],
            'odds_ratio_original_gwas': row['or'], 'standard_error_as_published': row['se'], 'p_value_original_gwas': row['p'],
            'source_range_id': source_id,
            'goode2024_supplement_risk_allele_as_published': 'A',
            'goode2020_thesis_risk_allele_as_published': 'G' if row['SNP'] == 'rs313839' else 'A',
            'working_risk_allele': row['allele_1'],
            'status': 'Original GWAS and population frequencies support C; contradictory thesis and molecular-QTL labeling remain unresolved' if row['SNP'] == 'rs313839' else 'A risk / C protective is consistent across these sources',
        })
    output = ROOT / 'data/derived/allele-audit.csv'
    with output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(output_rows)
    populations = json.loads((ROOT / 'data/raw/GRCh38-rs313839-populations.json').read_text())['populations']
    population_rows = [row for row in populations if row['population'] in {'1000GENOMES:phase_3:CEU', '1000GENOMES:phase_3:EUR'}]
    with (ROOT / 'data/derived/prkd2-population-alleles.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['population', 'allele', 'frequency', 'allele_count'], lineterminator='\n')
        writer.writeheader()
        writer.writerows(sorted(population_rows, key=lambda row: (row['population'], row['allele'])))
    print(json.dumps({'audited_variants': len(output_rows), 'population_allele_rows': len(population_rows), 'pilot_inputs_changed': False}))


if __name__ == '__main__':
    main()
