"""Coordinate checks for a positive-strand, annotated splice-donor experiment."""

def donor_contexts(gene, position_1based):
    if gene['strand'] != 1 or gene['assembly_name'] != 'GRCh38':
        raise ValueError('This audit requires a positive-strand GRCh38 gene')
    rows = []
    for transcript in gene['Transcript']:
        exons = sorted(transcript['Exon'], key=lambda exon: exon['start'])
        for index, (upstream, downstream) in enumerate(zip(exons, exons[1:]), 1):
            if upstream['end'] < position_1based < downstream['start']:
                rows.append({
                    'transcript_id': transcript['id'],
                    'transcript_version': transcript['version'],
                    'biotype': transcript['biotype'],
                    'is_canonical': transcript.get('is_canonical', 0),
                    'upstream_exon_number': index,
                    'upstream_exon_id': upstream['id'],
                    'upstream_exon_end_1based': upstream['end'],
                    'junction_start_0based': upstream['end'],
                    'model_donor_position_0based': upstream['end'] - 1,
                    'intron_offset_1based': position_1based - upstream['end'],
                    'downstream_exon_start_1based': downstream['start'],
                    'intron_length_bp': downstream['start'] - upstream['end'] - 1,
                })
    return rows


def position_index(position_0based, start_0based, end_0based_exclusive, resolution=1):
    if resolution != 1:
        raise ValueError('A nucleotide-level endpoint requires 1-bp output resolution')
    if not start_0based <= position_0based < end_0based_exclusive:
        raise ValueError('Endpoint is outside the predicted interval')
    return position_0based - start_0based


def junction_in_scope(start, end, strand, donor, downstream_exon_start):
    return strand == '+' and (
        start == donor or end == downstream_exon_start
        or (start < donor and end > downstream_exon_start)
    )
