"""Validate published benchmark variants against pinned public references."""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_json(name):
    return json.loads((ROOT / 'data/raw' / name).read_text())


def verify_sources():
    for name in ['sources.json', 'reference-sources.json']:
        manifest = json.loads((ROOT / 'config' / name).read_text())
        for source in manifest['sources']:
            content = (ROOT / 'data/raw' / source['file']).read_bytes()
            if hashlib.sha256(content).hexdigest() != source['sha256']:
                raise ValueError(f"Unverified source: {source['file']}")


def primary_mapping(record, assembly, chromosome):
    matches = [item for item in record['mappings'] if item['assembly_name'] == assembly and item['seq_region_name'] == chromosome and item['coord_system'] == 'chromosome']
    if len(matches) != 1:
        raise ValueError('Expected exactly one primary chromosome mapping')
    match = matches[0]
    if match['strand'] != 1 or match['start'] != match['end']:
        raise ValueError('Pilot only accepts forward-strand single-base mappings')
    return match


def validate_alleles(mapping, alternate, sequence_base, model_base):
    alleles = mapping['allele_string'].split('/')
    if any(base not in 'ACGT' or len(base) != 1 for base in alleles):
        raise ValueError('Unexpected single-nucleotide alleles')
    reference = alleles[0]
    if reference != sequence_base or reference != model_base:
        raise ValueError('Reference-base mismatch')
    if alternate == reference or alternate not in alleles:
        raise ValueError('Selected alternate is not a distinct recorded allele')
    return reference, alleles


def interpret_direction(reference, alternate, evidence):
    if not evidence:
        return None, None, 'No direct risk-allele evidence for this rsID; proxy direction not transferred'
    risks = {row['gwas_risk_allele'] for row in evidence}
    if len(risks) != 1:
        return None, None, 'Conflicting risk-allele evidence'
    risk = risks.pop()
    if risk not in {reference, alternate}:
        return risk, None, 'Published risk allele absent from selected alleles; direction excluded'
    if any(row['functional_risk_allele_as_published'] != risk for row in evidence):
        return risk, None, 'Functional effect allele differs from the risk label'
    signs = {1 if float(row['functional_beta']) > 0 else -1 if float(row['functional_beta']) < 0 else 0 for row in evidence}
    if len(signs) != 1 or 0 in signs:
        return risk, None, 'Expression effect varies across source contexts or includes a zero'
    sign = signs.pop() * (1 if risk == alternate else -1)
    return risk, sign, 'Aligned to the directly tabulated allele; source is a benchmark, not independent validation'


def candidate_genes(records, chromosome, position, width):
    start0 = position - 1 - width // 2
    end0 = start0 + width
    if start0 < 0:
        raise ValueError('Pilot window crosses chromosome start')
    result = []
    for gene in records:
        if gene['seq_region_name'] != chromosome or gene['strand'] not in [-1, 1]:
            continue
        tss = gene['start'] if gene['strand'] == 1 else gene['end']
        if start0 <= tss - 1 < end0:
            result.append({'gene_id': gene['gene_id'], 'gene_version': gene.get('version'), 'gene_name': gene.get('external_name', gene['gene_id']), 'biotype': gene['biotype'], 'strand': gene['strand'], 'tss_1based': tss, 'distance_to_variant': abs(tss - position)})
    result.sort(key=lambda gene: (gene['distance_to_variant'], gene['gene_id']))
    if not result or len({gene['gene_id'] for gene in result}) != len(result):
        raise ValueError('Missing or duplicate candidate genes')
    return start0, end0, result


def write_csv(path, records):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def main():
    verify_sources()
    config = json.loads((ROOT / 'config/benchmark-variants.json').read_text())
    with (ROOT / 'data/derived/psc-signal-summaries.csv').open() as stream:
        signals = list(csv.DictReader(stream))
    with (ROOT / 'data/derived/psc-molecular-colocalisation.csv').open() as stream:
        molecular = list(csv.DictReader(stream))
    normalized, all_genes = [], []
    for selection in config['variants']:
        rsid, target = selection['rsid'], selection['target_gene']
        rows = [row for row in signals if row['candidate_gene_as_published'] == target and rsid in [row['lead_gwas_snp'], row['highest_pp_snp']]]
        if len(rows) != 1:
            raise ValueError(f'Ambiguous published row for {rsid}')
        source = rows[0]
        positions = {int(source[field]) for identifier, field in [('lead_gwas_snp', 'lead_position_b37'), ('highest_pp_snp', 'highest_pp_position_b37')] if source[identifier] == rsid}
        if len(positions) != 1:
            raise ValueError(f'Conflicting source positions for {rsid}')
        chromosome = source['chromosome']
        mapped = {}
        for build in ['GRCh37', 'GRCh38']:
            record = read_json(f'{build}-{rsid}.json')
            if record['name'] != rsid:
                raise ValueError(f'Renamed variant requires review: {rsid}')
            mapped[build] = primary_mapping(record, build, chromosome)
        if mapped['GRCh37']['start'] != positions.pop():
            raise ValueError(f'Published coordinate mismatch: {rsid}')
        if mapped['GRCh37']['allele_string'] != mapped['GRCh38']['allele_string']:
            raise ValueError(f'Alleles differ between builds: {rsid}')
        base37 = read_json(f'GRCh37-{rsid}-sequence.json')['seq'].upper()
        base38 = read_json(f'GRCh38-{rsid}-sequence.json')['seq'].upper()
        model_base = (ROOT / f'data/raw/model-reference-{rsid}.txt').read_text().strip().upper()
        if base37 != base38:
            raise ValueError(f'Reference changes between builds: {rsid}')
        reference, alleles = validate_alleles(mapped['GRCh38'], selection['alternate_bases'], base38, model_base)
        evidence = [row for row in molecular if row['lead_gwas_snp'] == rsid and row['egene_as_published'] == target and row['molecular_qtl_type'] == 'eQTL']
        risk, expected, direction_note = interpret_direction(reference, selection['alternate_bases'], evidence)
        position = mapped['GRCh38']['start']
        start0, end0, genes = candidate_genes(read_json(f'GRCh38-{rsid}-genes.json'), chromosome, position, config['interval_width'])
        target_genes = [gene for gene in genes if gene['gene_name'] == target]
        if len(target_genes) != 1:
            raise ValueError(f'Missing/ambiguous target gene annotation for {rsid}')
        normalized.append({
            'rsid': rsid, 'target_gene': target, 'target_gene_id': target_genes[0]['gene_id'],
            'role': selection['role'], 'chromosome': 'chr' + chromosome,
            'position_grch37_1based': mapped['GRCh37']['start'], 'position_grch38_1based': position,
            'reference_bases': reference, 'alternate_bases': selection['alternate_bases'],
            'all_database_alleles': '/'.join(alleles), 'published_risk_allele': risk,
            'expected_alt_minus_ref_expression_sign': expected, 'direction_note': direction_note,
            'interval_start_0based': start0, 'interval_end_0based_exclusive': end0,
            'nearest_tss_gene': genes[0]['gene_name'], 'nearest_tss_distance_bp': genes[0]['distance_to_variant'],
            'candidate_gene_count': len(genes), 'ensembl_release': read_json('GRCh38-release.json')['releases'][0],
            'reference_status': 'Verified against both builds and AlphaGenome GRCh38.p13 FASTA',
        })
        all_genes.extend({'rsid': rsid, 'target_gene': target, **gene} for gene in genes)
    derived = ROOT / 'data/derived'
    write_csv(derived / 'benchmark-variants.csv', normalized)
    write_csv(derived / 'benchmark-candidate-genes.csv', all_genes)
    (derived / 'benchmark-variants.json').write_text(json.dumps(normalized, indent=2) + '\n')
    vcf = ['##fileformat=VCFv4.3', '##reference=GRCh38.p13', '##INFO=<ID=TARGET,Number=1,Type=String,Description="Published benchmark gene">', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    for row in normalized:
        vcf.append('\t'.join(map(str, [row['chromosome'], row['position_grch38_1based'], row['rsid'], row['reference_bases'], row['alternate_bases'], '.', 'PASS', 'TARGET=' + row['target_gene']])))
    (derived / 'benchmark-variants.vcf').write_text('\n'.join(vcf) + '\n')
    report = {'variants_verified': len(normalized), 'loci': len({row['target_gene'] for row in normalized}), 'direction_tests_supported_by_source': sum(row['expected_alt_minus_ref_expression_sign'] is not None for row in normalized), 'unresolved': [{'rsid': row['rsid'], 'reason': row['direction_note']} for row in normalized if row['expected_alt_minus_ref_expression_sign'] is None], 'model_predictions_run_by_this_script': False}
    (ROOT / 'reports/benchmark-input-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
