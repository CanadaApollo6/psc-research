"""Check target splice motifs and a conditional explanation of reversed labels.

The truth table reproduces the bit-flag expression in regtools 0.6.0 source;
it does not execute that binary or prove the historical Lepik run parameters.
"""

import json
import gzip
from pathlib import Path

from fetch_sources import validate
from summarize_ubash3a_junctions import read_tsv, write_csv

ROOT = Path(__file__).resolve().parents[1]


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def motifs(sequence, sequence_start0, intron_start0, intron_end0):
    left, right = intron_start0 - sequence_start0, intron_end0 - sequence_start0
    if left < 0 or right > len(sequence) or right - left < 4:
        raise ValueError('Junction outside verified reference interval')
    donor, acceptor = sequence[left:left + 2], sequence[right - 2:right]
    return donor, acceptor, reverse_complement(acceptor), reverse_complement(donor)


def regtools_strand(flag, mode):
    if mode not in (1, 2):
        raise ValueError('Only RF=1 and FR=2 bit-flag modes are modeled')
    reversed_read, mate_reversed = (flag >> 4) % 2, (flag >> 5) % 2
    first, second = (flag >> 6) % 2, (flag >> 7) % 2
    base = int(not (mode - 1))
    first_strand = base ^ first ^ reversed_read
    second_strand = base ^ second ^ mate_reversed
    return ('+' if first_strand else '-') if first_strand == second_strand else '?'


def main():
    source_list = json.loads((ROOT / 'config/ubash3a-junction-sources.json').read_text())['sources']
    sources = {s['id']: s for s in source_list}

    def source(name):
        s = sources[name]
        content = (ROOT / 'data/raw' / s['file']).read_bytes()
        validate(content, s)
        return content

    reference = json.loads(source('ubash-intron-GRCh38'))
    if reference['id'] != 'chromosome:GRCh38:21:42434944:42437504:1' or len(reference['seq']) != 2561:
        raise ValueError('Unexpected reference sequence interval')
    plan = json.loads((ROOT / 'config/ubash3a-public-junction-plan.json').read_text())
    seq = reference['seq']
    if seq[42434957 - 42434944] != 'A':
        raise ValueError('Reference does not match rs1893592-A')
    annotation = read_tsv(gzip.decompress(source('QTD000377-phenotype-metadata')).decode())
    rows = []
    for j in plan['junctions']:
        a, b = j['start_0based'], j['end_0based_exclusive']
        d, ac, rd, ra = motifs(seq, 42434943, a, b)
        phenotype = f'21:{a}:{b+1}:clu_31404_-'
        matches = [r for r in annotation if r['phenotype_id'] == phenotype and r['gene_id'] == 'ENSG00000160185']
        if len(matches) != 1 or matches[0]['strand'] != '1':
            raise ValueError('Published target annotation differs')
        rows.append({'junction': j['name'], 'chromosome': 'chr21', 'build': 'GRCh38',
                     'start_0based': a, 'end_0based_exclusive': b,
                     'positive_strand_donor': d, 'positive_strand_acceptor': ac,
                     'negative_strand_donor': rd, 'negative_strand_acceptor': ra,
                     'Lepik_source_cluster_strand': '-', 'UBASH3A_annotation_strand': '+',
                     'interpretation': 'GT-AG reference motif supports positive-strand boundary; source label preserved'})
    write_csv(ROOT / 'reports/ubash3a-junction-motifs.csv', rows)
    extractor = source('regtools-v06-junction-cc').decode()
    module = source('eqtl-leafcutter-module').decode()
    alignment = source('rnaseq-align-module').decode()
    if ('1 = first-strand/RF, 2, = second-strand/FR' not in extractor or
            'leafcutter_strand = 2' not in module or "'--rna-strandness RF'" not in alignment):
        raise ValueError('Pinned strand conventions differ; review the conditional truth table')
    truth = []
    for strand, flags in [('+', (83, 163)), ('-', (99, 147))]:
        for ordinal, flag in enumerate(flags, 1):
            truth.append({'assumed_library': 'first-strand/RF', 'true_RNA_strand': strand,
                          'read_in_pair': ordinal, 'SAM_flag': flag,
                          'regtools_s1_RF': regtools_strand(flag, 1),
                          'regtools_s2_FR': regtools_strand(flag, 2),
                          'catalogue_reverse_stranded_branch_selects': 2,
                          'scope': 'synthetic source-code truth table; actual historical flags unavailable'})
    write_csv(ROOT / 'reports/ubash3a-strand-flag-audit.csv', truth)
    permuted = read_tsv(gzip.decompress(source('lepik-splice-permuted')).decode())
    target = [r for r in permuted if r['molecular_trait_object_id'] == 'clu_31404_-']
    if len(target) != 1:
        raise ValueError('Expected one permuted lead record for the target cluster')
    write_csv(ROOT / 'data/derived/ubash3a-lepik-cluster-lead.csv', target)
    print(json.dumps({'positive_GT_AG_junctions': sum(r['positive_strand_donor'] == 'GT' and r['positive_strand_acceptor'] == 'AG' for r in rows),
                      'truth_table_rows': len(truth), 'lead_records': len(target)}))


if __name__ == '__main__':
    main()
