"""Build the complete mapping audit, common gene universe and pre-genotype pool."""

import csv
import gzip
import io
import json
import math
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.ipc as ipc

from audit_finemap_archive import csv_bytes
from prepare_controlled_comparison import assess_mapping, mapped_record, noncoding_context

ROOT = Path(__file__).resolve().parents[1]


def read_csv(name):
    return list(csv.DictReader(io.StringIO((ROOT / name).read_text())))


def mutation_class(ref, alt):
    if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
        raise ValueError('Expected distinct single-base alleles')
    return 'transition' if {ref, alt} in ({'A', 'G'}, {'C', 'T'}) else 'transversion'


def model_genes(windows):
    """Match the official exon extractor's transcript-TSS membership rule."""
    genes = {w['archive_region']: set() for w in windows}
    with pa.memory_map(str(ROOT / 'data/raw/match-model-gencode-v46.feather'), 'r') as stream:
        reader = ipc.open_file(stream)
        for index in range(reader.num_record_batches):
            batch = reader.get_batch(index)
            transcripts = batch.filter(pc.equal(batch['Feature'], 'transcript')).select(['Chromosome', 'Start', 'End', 'Strand', 'gene_id_nopatch'])
            for row in transcripts.to_pylist():
                # Official AlphaGenome helper uses End (exclusive), not End-1,
                # for its negative-strand transcript-TSS point query.
                tss0 = row['End'] if row['Strand'] == '-' else row['Start']
                for window in windows:
                    if row['Chromosome'] == window['chromosome'] and int(window['interval_start_0based']) <= tss0 < int(window['interval_end_0based_exclusive']):
                        genes[window['archive_region']].add(row['gene_id_nopatch'])
    return genes


def gene_universe(windows):
    model_ids = model_genes(windows)
    result = []
    for w in windows:
        label = w['archive_region']
        records = json.loads((ROOT / f'data/raw/match-GRCh38-{label}-genes.json').read_text())
        seen = set()
        for gene in records:
            if gene['seq_region_name'] != w['chromosome'].removeprefix('chr') or gene['strand'] not in (-1, 1):
                continue
            tss = gene['start'] if gene['strand'] == 1 else gene['end']
            if not int(w['interval_start_0based']) <= tss - 1 < int(w['interval_end_0based_exclusive']):
                continue
            identifier = gene['gene_id']
            if identifier in seen:
                raise ValueError('Duplicate gene annotation')
            seen.add(identifier)
            result.append({'archive_region': label, 'gene_id': identifier, 'gene_name': gene.get('external_name', identifier),
                           'ensembl_gene_version': gene.get('version', ''), 'biotype': gene['biotype'], 'strand': gene['strand'],
                           'tss_grch38_1based': tss, 'present_in_model_transcript_tss_universe': identifier in model_ids[label]})
    return sorted(result, key=lambda r: (r['archive_region'], r['gene_id']))


def nearest(genes, label, position):
    candidates = [g for g in genes if g['archive_region'] == label and g['present_in_model_transcript_tss_universe']]
    if not candidates:
        raise ValueError('Empty fixed gene universe')
    gene = min(candidates, key=lambda g: (abs(g['tss_grch38_1based'] - position), g['gene_id']))
    return gene['gene_id'], gene['gene_name'], abs(gene['tss_grch38_1based'] - position)


def static_pair(anchor, comparator, limits):
    return (anchor['mutation_class'] == comparator['mutation_class']
            and anchor['noncoding_context'] == comparator['noncoding_context']
            and abs(anchor['position_grch38_1based'] - comparator['position_grch38_1based']) <= limits['maximum_distance_from_anchor_bp']
            and abs(math.log10(anchor['nearest_tss_distance_bp'] + 1) - math.log10(comparator['nearest_tss_distance_bp'] + 1)) <= .5)


def bulk_records(manifest, label, build):
    records, requested = {}, set()
    specs = [s for s in manifest['sources'] if s['id'].startswith(f'match-{build}-{label}-variants-')]
    for spec in specs:
        payload = json.loads((ROOT / 'data/raw' / spec['file']).read_text())
        ids = set(spec['request_json']['ids'])
        if ids & requested or not isinstance(payload, dict) or 'error' in payload:
            raise ValueError('Duplicate batch identifiers or failed response')
        requested |= ids
        # Ensembl may key a merged rsID's response by its replacement ID.
        # Keep raw responses pinned, but do not use an unrequested replacement
        # as the original identifier or overwrite another batch's direct query.
        records.update({key: value for key, value in payload.items() if key in ids})
    return records, requested


def main():
    protocol = json.loads((ROOT / 'config/psc-controlled-comparison.json').read_text())
    manifest = json.loads((ROOT / 'config/psc-comparator-sources.json').read_text())
    windows = read_csv('data/derived/psc-comparison-windows.csv')
    genes = gene_universe(windows)
    anchors = []
    for source in read_csv('data/derived/psc-comparison-candidates.csv'):
        if source['reference_eligible'] != 'True':
            continue
        row = {key: source[key] for key in ['reference_bases', 'alternate_bases', 'noncoding_context']}
        row.update(archive_region=source['dataset_id'].split('/')[1], rsid=source['variant_id_as_published'],
                   chromosome=source['chromosome'], position_grch37_1based=int(source['position_grch37_1based']),
                   position_grch38_1based=int(source['position_grch38_1based']), pip=source['snp_prob_as_published'], role=source['role'],
                   most_severe_consequence=source['most_severe_consequence_grch38'], mutation_class=mutation_class(source['reference_bases'], source['alternate_bases']))
        row['nearest_tss_gene_id'], row['nearest_tss_gene_name'], row['nearest_tss_distance_bp'] = nearest(genes, row['archive_region'], row['position_grch38_1based'])
        anchors.append(row)
    with gzip.open(ROOT / 'data/derived/psc-finemap-variants.csv.gz', 'rt') as stream:
        all_rows = list(csv.DictReader(stream))
    audit, potential, queries, requests = [], [], [], []
    for w in windows:
        label, chrom = w['archive_region'], w['chromosome'].removeprefix('chr')
        low = [r for r in all_rows if r['dataset_id'] == 'gwas/' + label and float(r['snp_prob_as_published']) <= protocol['comparators']['maximum_archived_marginal_pip'] and r['identifier_type'] == 'rsid']
        maps = {b: bulk_records(manifest, label, b) for b in ('GRCh37', 'GRCh38')}
        expected = {r['variant_id_as_published'] for r in low}
        if any(requested != expected for records, requested in maps.values()):
            raise ValueError('Bulk query inventory does not cover the full literal-rsID low-PIP pool')
        for source in low:
            rsid = source['variant_id_as_published']
            row = {'archive_region': label, 'rsid': rsid, 'pip': source['snp_prob_as_published'], 'status': '',
                   'chromosome': chrom, 'position_grch37_1based': '', 'position_grch38_1based': '',
                   'reference_bases': '', 'alternate_bases': '', 'most_severe_consequence': '', 'noncoding_context': '',
                   'mutation_class': '', 'nearest_tss_gene_id': '', 'nearest_tss_gene_name': '', 'nearest_tss_distance_bp': '',
                   'possible_anchor_ids': ''}
            try:
                records = {b: maps[b][0].get(rsid) for b in maps}
                if any(not r or 'error' in r for r in records.values()):
                    raise ValueError('Missing or failed database mapping')
                mapped = {b: mapped_record(records[b], rsid, b, chrom) for b in records}
                reason = assess_mapping(mapped['GRCh37'], mapped['GRCh38'])
                if reason:
                    raise ValueError(reason)
                ref, alt = mapped['GRCh38']['allele_string'].split('/')
                context = noncoding_context(records['GRCh38']['most_severe_consequence'])
                if not context:
                    raise ValueError('Consequence outside fixed noncoding classes')
                position = mapped['GRCh38']['start']
                if not int(w['interval_start_0based']) <= position - 1 < int(w['interval_end_0based_exclusive']):
                    raise ValueError('Outside fixed region window')
                gene_id, gene_name, distance = nearest(genes, label, position)
                row.update(position_grch37_1based=mapped['GRCh37']['start'], position_grch38_1based=position,
                           reference_bases=ref, alternate_bases=alt, most_severe_consequence=records['GRCh38']['most_severe_consequence'],
                           noncoding_context=context, mutation_class=mutation_class(ref, alt),
                           nearest_tss_gene_id=gene_id, nearest_tss_gene_name=gene_name, nearest_tss_distance_bp=distance)
                possible = [a['rsid'] for a in anchors if a['archive_region'] == label and a['role'] == 'anchor_pending_comparators' and static_pair(a, row, protocol['comparators'])]
                if not possible:
                    raise ValueError('No anchor passes fixed mutation/context/TSS/physical-distance criteria')
                row.update(status='potential_comparator_pending_bases_genotypes', possible_anchor_ids=';'.join(possible))
                potential.append({k: row[k] for k in anchors[0] if k != 'role'} | {'role': 'comparator_pool'})
            except ValueError as exc:
                row['status'] = str(exc)
            audit.append(row)
        region = [r for r in anchors + potential if r['archive_region'] == label]
        start37 = min(r['position_grch37_1based'] for r in region) - 1
        end37 = max(r['position_grch37_1based'] for r in region)
        if end37 - start37 > 1_500_000:
            raise ValueError('Genotype query span exceeds bound')
        queries.append({'archive_region': label, 'chromosome': chrom, 'start_grch37_0based': start37, 'end_grch37_0based_exclusive': end37,
                        'variant_ids': [r['rsid'] for r in region], 'index_source_id': f'match-kgp-chr{chrom}-index',
                        'url': f'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'})
        for build, server, start, end in [('GRCh37', 'https://grch37.rest.ensembl.org', start37, end37),
                                         ('GRCh38', 'https://rest.ensembl.org', int(w['interval_start_0based']), int(w['interval_end_0based_exclusive']))]:
            name = f'match-{build}-{label}-reference'
            requests.append({'id': name, 'file': name + '.json', 'url': f'{server}/sequence/region/human/{chrom}:{start+1}..{end}:1?content-type=application/json',
                             'max_bytes': 2000000, 'region': {'chromosome': chrom, 'start0': start, 'end0': end}})
        fai = {line.split('\t')[0]: list(map(int, line.split('\t')[1:])) for line in (ROOT / 'data/raw/GRCh38.p13.genome.fa.fai').read_text().splitlines()}
        _, offset, bases, width = fai['chr' + chrom]
        def byte_at(pos):
            return offset + (pos // bases) * width + pos % bases
        byte_start, byte_end = byte_at(int(w['interval_start_0based'])), byte_at(int(w['interval_end_0based_exclusive']) - 1)
        name = f'match-model-{label}-reference'
        requests.append({'id': name, 'file': name + '.txt', 'url': 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/GRCh38.p13.genome.fa',
                         'request_headers': {'Range': f'bytes={byte_start}-{byte_end}'}, 'expected_content_range': f'bytes {byte_start}-{byte_end}/',
                         'max_bytes': byte_end - byte_start + 1, 'region': {'chromosome': chrom, 'start0': int(w['interval_start_0based']), 'end0': int(w['interval_end_0based_exclusive'])}})
    (ROOT / 'data/derived/psc-comparison-gene-universe.csv').write_bytes(csv_bytes(genes))
    (ROOT / 'data/derived/psc-comparator-mapping-audit.csv').write_bytes(csv_bytes(audit))
    (ROOT / 'data/derived/psc-comparator-genotype-targets.csv').write_bytes(csv_bytes(anchors + potential))
    (ROOT / 'data/derived/psc-comparator-query-plan.json').write_text(json.dumps(queries, indent=2) + '\n')
    (ROOT / 'data/raw/match-reference-requests.json').write_text(json.dumps(requests, indent=2) + '\n')
    print(json.dumps({'low_pip_rsid_records_audited': len(audit), 'potential_comparators': len(potential), 'high_pip_genotype_targets': len(anchors),
                      'gene_universe_by_region': {w['archive_region']: sum(g['archive_region'] == w['archive_region'] and g['present_in_model_transcript_tss_universe'] for g in genes) for w in windows}}))


if __name__ == '__main__':
    main()
