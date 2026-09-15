"""Full-allele identities for the separately frozen PRKD2 coverage extension."""

import json
from pathlib import Path

from prkd2_common import ROOT, ChainMap, complement, reciprocal_base, sha
from coverage_repair_common import Reference, normalize, reciprocal_edit
from prepare_coloc_inputs import canonical_key


def load_references():
    record = json.loads((ROOT / 'config/prkd2-expanded-sources.json').read_text())
    assert record['plan_sha256'] == sha((ROOT / 'config/prkd2-expanded-plan.json').read_bytes())
    refs = {}
    for item in record['files']:
        if item.get('build') in ['GRCh37', 'GRCh38']:
            data = (ROOT / item['file']).read_bytes()
            assert sha(data) == item['sha256']
            sequence = data.decode().strip().upper()
            assert len(sequence) == item['length'] and sha(sequence.encode()) == item['sequence_sha256']
            refs[item['build']] = Reference('chr19', item['start_1based'], sequence)
    return refs


def normalize_edit(reference, position, ref, alt):
    if len(ref) == len(alt) == 1:
        if ref == alt or set(ref + alt) - set('ACGT') or reference.get(position) != ref:
            raise ValueError('Invalid SNV or reference-base disagreement')
        return position, ref, alt
    start = max(reference.start, position - 10000)
    end = min(reference.start + len(reference.sequence), position + len(ref) + 10000)
    local = Reference(reference.chromosome, start, reference.get(start, end - start))
    return normalize(local, position, ref, alt)


def key_and_sign(position, ref, alt):
    if len(ref) == len(alt) == 1:
        key, _, effect = canonical_key('chr19', position, ref, alt)
        return key, 1 if alt == effect else -1
    return f'chr19:{position}:{ref}:{alt}', 1


def qtl_identity(variant, reference):
    fields = variant.split('_')
    if len(fields) != 4 or fields[0] != 'chr19':
        raise ValueError('Not an explicit chromosome-19 edit')
    position, ref, alt = normalize_edit(reference, int(fields[1]), fields[2], fields[3])
    key, sign = key_and_sign(position, ref, alt)
    return {'snp': key, 'position': position, 'normalized_ref': ref, 'normalized_alt': alt,
            'canonical_beta_sign': sign, 'source_variant': variant}


def map_edit(ref37, ref38, forward, reverse, position, ref, alt):
    norm37 = normalize_edit(ref37, position, ref, alt)
    if len(ref) == len(alt) == 1:
        hit = reciprocal_base(forward, reverse, 'chr19', position - 1)
        if hit is None:
            raise ValueError('Nonunique or nonreciprocal SNV mapping')
        chrom, p0, strand = hit
        if chrom != 'chr19':
            raise ValueError('Wrong target chromosome')
        mapped_ref, mapped_alt = (complement(ref), complement(alt)) if strand == '-' else (ref, alt)
        mapped = (p0 + 1, mapped_ref, mapped_alt)
        bases = 1
    else:
        mapped, bases, strand = reciprocal_edit(ref37, ref38, forward, reverse, *norm37, flank=100)
    norm38 = normalize_edit(ref38, *mapped)
    key, dosage_sign = key_and_sign(*norm38)
    return {'snp': key, 'normalized_grch37': norm37, 'normalized_grch38': norm38,
            'reciprocal_bases': bases, 'strand': strand, 'dosage_sign': dosage_sign}
